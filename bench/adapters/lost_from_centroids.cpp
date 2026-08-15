/**
 * LOST adapter: frozen centroids in → pyramid star-ID + DQM attitude out.
 *
 * Usage:
 *   lost_from_centroids CENTROIDS.tsv DATABASE.dat WIDTH HEIGHT FOV_DEG
 *
 * CENTROIDS.tsv lines: x y [intensity]
 * Coordinates use top-left origin, matching LOST Star.position convention.
 */
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "attitude-estimators.hpp"
#include "attitude-utils.hpp"
#include "camera.hpp"
#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"

using namespace lost;

namespace {

std::unique_ptr<unsigned char[]> LoadDatabase(const std::string &path, long *lengthOut) {
    std::fstream fs;
    fs.open(path, std::fstream::in | std::fstream::binary);
    fs.seekg(0, fs.end);
    long length = fs.tellg();
    fs.seekg(0, fs.beg);
    if (fs.fail() || length <= 0) {
        std::cerr << "ERROR: cannot read database: " << path << std::endl;
        std::exit(2);
    }
    auto buffer = std::unique_ptr<unsigned char[]>(new unsigned char[length]);
    fs.read(reinterpret_cast<char *>(buffer.get()), length);
    *lengthOut = length;
    return buffer;
}

Stars LoadCentroids(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        std::cerr << "ERROR: cannot read centroids: " << path << std::endl;
        std::exit(2);
    }
    Stars stars;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        decimal x = 0;
        decimal y = 0;
        int intensity = 1000;
        if (!(iss >> x >> y)) {
            continue;
        }
        iss >> intensity;
        if (intensity < 0) {
            intensity = 0;
        }
        stars.push_back(Star(x, y, DECIMAL(1.0), DECIMAL(1.0), intensity));
    }
    // Pyramid expects brightest first for its stepped scan; sort by magnitude descending.
    std::sort(stars.begin(), stars.end(), [](const Star &a, const Star &b) {
        return a.magnitude > b.magnitude;
    });
    return stars;
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 6 || argc > 8) {
        std::cerr << "Usage: lost_from_centroids CENTROIDS.tsv DATABASE.dat "
                     "WIDTH HEIGHT FOV_DEG "
                     "[ANGULAR_TOL_DEG=0.05] [MAX_MISMATCH_PROB=0.0001]\n";
        return 2;
    }

    const std::string centroidsPath = argv[1];
    const std::string databasePath = argv[2];
    const int width = std::atoi(argv[3]);
    const int height = std::atoi(argv[4]);
    const decimal fovDeg = STR_TO_DECIMAL(argv[5]);
    const decimal angularTolDeg =
        (argc >= 7) ? STR_TO_DECIMAL(argv[6]) : DECIMAL(0.05);
    const decimal maxMismatchProb =
        (argc >= 8) ? STR_TO_DECIMAL(argv[7]) : DECIMAL(0.0001);

    if (width <= 0 || height <= 0 || fovDeg <= DECIMAL(0.0) ||
        angularTolDeg <= DECIMAL(0.0) || maxMismatchProb <= DECIMAL(0.0)) {
        std::cerr << "ERROR: invalid width/height/fov/tolerance/mismatch\n";
        return 2;
    }

    Stars stars = LoadCentroids(centroidsPath);
    if (stars.size() < 4) {
        std::cout << "attitude_known 0\n";
        std::cout << "num_centroids " << stars.size() << "\n";
        std::cout << "num_star_ids 0\n";
        return 0;
    }

    long dbLength = 0;
    std::unique_ptr<unsigned char[]> database = LoadDatabase(databasePath, &dbLength);

    MultiDatabase multiDatabase(database.get());
    const unsigned char *catalogBuffer =
        multiDatabase.SubDatabasePointer(kCatalogMagicValue);
    if (catalogBuffer == NULL) {
        std::cerr << "ERROR: database has no catalog\n";
        return 3;
    }
    DeserializeContext catalogDes(catalogBuffer);
    Catalog catalog = DeserializeCatalog(&catalogDes, NULL, NULL);

    const decimal focalLength = FovToFocalLength(DegToRad(fovDeg), width);
    Camera camera(focalLength, width, height);

    PyramidStarIdAlgorithm pyramid(
        DegToRad(angularTolDeg),
        1000,
        maxMismatchProb,
        1000);
    DavenportQAlgorithm attitudeAlgo;

    const auto starIdBegin = std::chrono::steady_clock::now();
    StarIdentifiers starIds = pyramid.Go(database.get(), stars, catalog, camera);
    const auto starIdEnd = std::chrono::steady_clock::now();
    const long long starIdNs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(starIdEnd - starIdBegin).count();

    std::cout << "num_centroids " << stars.size() << "\n";
    std::cout << "num_star_ids " << starIds.size() << "\n";
    std::cout << "starid_average_ns " << starIdNs << "\n";

    if (starIds.size() < 2) {
        std::cout << "attitude_known 0\n";
        return 0;
    }

    const auto attitudeBegin = std::chrono::steady_clock::now();
    Attitude attitude = attitudeAlgo.Go(camera, stars, catalog, starIds);
    const auto attitudeEnd = std::chrono::steady_clock::now();
    const long long attitudeNs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(attitudeEnd - attitudeBegin).count();

    std::cout << "attitude_average_ns " << attitudeNs << "\n";
    std::cout << "total_average_ns " << (starIdNs + attitudeNs) << "\n";

    if (!attitude.IsKnown()) {
        std::cout << "attitude_known 0\n";
        return 0;
    }

    EulerAngles spherical = attitude.ToSpherical();
    std::cout << "attitude_known 1\n";
    std::cout << "attitude_ra " << RadToDeg(spherical.ra) << "\n";
    std::cout << "attitude_de " << RadToDeg(spherical.de) << "\n";
    std::cout << "attitude_roll " << RadToDeg(spherical.roll) << "\n";
    return 0;
}
