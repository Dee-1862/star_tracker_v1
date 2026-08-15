/**
 * Ours: frozen centroids in -> star ID + gated attitude out.
 *
 * Mirrors the interface of bench/adapters/lost_from_centroids so the three
 * solvers can be compared on identical centroid lists with centroiding removed
 * as a variable.
 *
 *   centroid_runner CENTROIDS.tsv WIDTH HEIGHT FOV_DEG [MAX_RESIDUAL_ARCSEC]
 *
 * CENTROIDS.tsv lines: x y [intensity], top-left origin.
 */
#include "star_tracker/attitude_solver.hpp"
#include "star_tracker/centroiding.hpp"
#include "star_tracker/hipparcos_catalog.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

star_tracker::LisGridMatcher matcher;

bool loadCentroids(
    const std::string& path, star_tracker::Centroiding::Result& result) {
    std::ifstream input(path.c_str());
    if (!input) {
        std::cerr << "ERROR: cannot read centroids: " << path << "\n";
        return false;
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream stream(line);
        float x = 0.0F;
        float y = 0.0F;
        long intensity = 1000L;
        if (!(stream >> x >> y)) {
            continue;
        }
        stream >> intensity;
        if (intensity < 0L) {
            intensity = 0L;
        }
        if (result.count >= star_tracker::Centroiding::kMaxCentroids) {
            break;
        }
        result.points[result.count].x = x;
        result.points[result.count].y = y;
        result.points[result.count].intensity =
            static_cast<std::uint32_t>(intensity);
        ++result.count;
    }
    return true;
}

const char* reasonText(star_tracker::AttitudeSolver::Reason reason) {
    switch (reason) {
        case star_tracker::AttitudeSolver::kAccepted:
            return "accepted";
        case star_tracker::AttitudeSolver::kTooFewStars:
            return "too_few_stars";
        case star_tracker::AttitudeSolver::kDegenerateGeometry:
            return "degenerate_geometry";
        case star_tracker::AttitudeSolver::kResidualTooLarge:
            return "residual_too_large";
        default:
            return "unknown";
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: centroid_runner CENTROIDS.tsv WIDTH HEIGHT "
                     "FOV_DEG [MAX_RESIDUAL_ARCSEC=28]\n";
        return 2;
    }

    const std::string centroidsPath = argv[1];
    const int width = std::atoi(argv[2]);
    const int height = std::atoi(argv[3]);
    const double fovDeg = std::atof(argv[4]);
    const float maxResidual = (argc >= 6) ? static_cast<float>(std::atof(argv[5]))
                                          : 28.0F;

    if (width <= 0 || height <= 0 || fovDeg <= 0.0 || maxResidual <= 0.0F) {
        std::cerr << "ERROR: invalid width/height/fov/residual\n";
        return 2;
    }

    star_tracker::Centroiding::Result centroids;
    if (!loadCentroids(centroidsPath, centroids)) {
        return 2;
    }

    const double pi = 3.14159265358979323846;
    const float focal = static_cast<float>(
        (static_cast<double>(width) / 2.0) /
        std::tan((fovDeg * pi / 180.0) / 2.0));

    star_tracker::LisGridMatcher::Config matcherConfig;
    matcherConfig.focal_length_pixels = focal;
    matcherConfig.principal_x = static_cast<float>(width) / 2.0F;
    matcherConfig.principal_y = static_cast<float>(height) / 2.0F;
    matcher = star_tracker::LisGridMatcher(matcherConfig);

    star_tracker::AttitudeSolver::Config solverConfig;
    solverConfig.focal_length_pixels = focal;
    solverConfig.principal_x = static_cast<float>(width) / 2.0F;
    solverConfig.principal_y = static_cast<float>(height) / 2.0F;
    solverConfig.max_residual_arcsec = maxResidual;
    star_tracker::AttitudeSolver solver(solverConfig);

    typedef std::chrono::steady_clock Clock;

    const Clock::time_point indexBegin = Clock::now();
    if (!matcher.buildIndex(
            star_tracker::kHipparcosCatalog,
            star_tracker::kHipparcosCatalogCount)) {
        std::cerr << "ERROR: failed to build the pair index\n";
        return 3;
    }
    const Clock::time_point indexEnd = Clock::now();

    const Clock::time_point matchBegin = Clock::now();
    const star_tracker::LisGridMatcher::MatchResult matches =
        matcher.match(centroids);
    const Clock::time_point matchEnd = Clock::now();

    const Clock::time_point attitudeBegin = Clock::now();
    const star_tracker::AttitudeSolver::Solution solution = solver.solve(
        centroids, matches, star_tracker::kHipparcosCatalog,
        star_tracker::kHipparcosCatalogCount);
    const Clock::time_point attitudeEnd = Clock::now();

    const long long indexNs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            indexEnd - indexBegin).count();
    const long long matchNs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            matchEnd - matchBegin).count();
    const long long attitudeNs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            attitudeEnd - attitudeBegin).count();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "num_centroids " << centroids.count << "\n";
    std::cout << "num_star_ids " << matches.matched_count << "\n";
    std::cout << "num_paired " << solution.matched_count << "\n";
    std::cout << "index_ns " << indexNs << "\n";
    std::cout << "starid_average_ns " << matchNs << "\n";
    std::cout << "attitude_average_ns " << attitudeNs << "\n";
    std::cout << "total_average_ns " << (matchNs + attitudeNs) << "\n";
    std::cout << "clique_size " << matches.clique_size << "\n";
    std::cout << "residual_rms_arcsec " << solution.residual_rms_arcsec << "\n";
    std::cout << "residual_max_arcsec " << solution.residual_max_arcsec << "\n";
    std::cout << "gate_reason " << reasonText(solution.reason) << "\n";

    for (std::size_t index = 0U; index < centroids.count; ++index) {
        if (matches.star_ids[index] == 0U) {
            continue;
        }
        std::cout << "match_" << index << "_id " << matches.star_ids[index]
                  << "\n";
        std::cout << "match_" << index << "_x " << centroids.points[index].x
                  << "\n";
        std::cout << "match_" << index << "_y " << centroids.points[index].y
                  << "\n";
    }

    if (!solution.valid) {
        std::cout << "attitude_known 0\n";
        return 0;
    }

    std::cout << "attitude_known 1\n";
    std::cout << "attitude_ra " << solution.right_ascension_deg << "\n";
    std::cout << "attitude_de " << solution.declination_deg << "\n";
    std::cout << "attitude_roll " << solution.roll_deg << "\n";
    return 0;
}
