/**
 * Ours, sequence mode: many frames through ONE tracker, state preserved.
 *
 * centroid_runner processes a single frame and exits, which cannot exercise
 * anything that depends on history: focal-length tracking, mode transitions, or
 * how refusals cluster in time. This feeds a list of frames through one
 * StarTracker so that behaviour is observable.
 *
 *   sequence_runner LIST.txt WIDTH HEIGHT FOV_DEG [MAX_RESIDUAL] [SEARCH_PCT]
 *
 * LIST.txt holds one centroid TSV path per line, in time order.
 * Emits one CSV row per frame.
 */
#include "star_tracker/centroiding.hpp"
#include "star_tracker/hipparcos_catalog.hpp"
#include "star_tracker/tracker.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

star_tracker::StarTracker tracker;

bool loadCentroids(
    const std::string& path, star_tracker::Centroiding::Result& result) {
    std::ifstream input(path.c_str());
    if (!input) {
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
        if (result.count >= star_tracker::Centroiding::kMaxCentroids) {
            break;
        }
        result.points[result.count].x = x;
        result.points[result.count].y = y;
        result.points[result.count].intensity =
            static_cast<std::uint32_t>((intensity < 0L) ? 0L : intensity);
        ++result.count;
    }
    return true;
}

const char* reasonText(star_tracker::AttitudeSolver::Reason reason) {
    switch (reason) {
        case star_tracker::AttitudeSolver::kAccepted: return "accepted";
        case star_tracker::AttitudeSolver::kTooFewStars: return "too_few_stars";
        case star_tracker::AttitudeSolver::kDegenerateGeometry:
            return "degenerate_geometry";
        case star_tracker::AttitudeSolver::kResidualTooLarge:
            return "residual_too_large";
        default: return "unknown";
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5 || argc > 7) {
        std::cerr << "Usage: sequence_runner LIST.txt WIDTH HEIGHT FOV_DEG "
                     "[MAX_RESIDUAL_ARCSEC=30] [SEARCH_PCT=6]\n";
        return 2;
    }

    const std::string listPath = argv[1];
    const int width = std::atoi(argv[2]);
    const int height = std::atoi(argv[3]);
    const double fovDeg = std::atof(argv[4]);
    const float maxResidual =
        (argc >= 6) ? static_cast<float>(std::atof(argv[5])) : 30.0F;
    const float searchPct =
        (argc >= 7) ? static_cast<float>(std::atof(argv[6])) : 6.0F;

    if (width <= 0 || height <= 0 || fovDeg <= 0.0) {
        std::cerr << "ERROR: invalid width/height/fov\n";
        return 2;
    }

    const double pi = 3.14159265358979323846;
    const float focal = static_cast<float>(
        (static_cast<double>(width) / 2.0) /
        std::tan((fovDeg * pi / 180.0) / 2.0));

    star_tracker::StarTracker::Config config;
    config.nominal_focal_pixels = focal;
    config.principal_x = static_cast<float>(width) / 2.0F;
    config.principal_y = static_cast<float>(height) / 2.0F;
    config.max_residual_arcsec = maxResidual;
    config.search_percent = searchPct;
    tracker = star_tracker::StarTracker(config);

    if (!tracker.buildIndex(star_tracker::kHipparcosCatalog,
                            star_tracker::kHipparcosCatalogCount)) {
        std::cerr << "ERROR: failed to build the pair index\n";
        return 3;
    }

    std::ifstream list(listPath.c_str());
    if (!list) {
        std::cerr << "ERROR: cannot read list: " << listPath << "\n";
        return 2;
    }

    std::cout << "frame,solved,mode,reason,clique,ra_deg,de_deg,roll_deg,"
                 "residual_arcsec,focal_used,focal_measured,focal_tracked,"
                 "focal_innovation,focal_updated,focal_jump_rejected,trials\n";
    std::cout << std::fixed << std::setprecision(6);

    std::string path;
    std::size_t frame = 0U;
    while (std::getline(list, path)) {
        if (path.empty()) {
            continue;
        }
        star_tracker::Centroiding::Result centroids;
        if (!loadCentroids(path, centroids)) {
            std::cerr << "WARN: skipping unreadable frame: " << path << "\n";
            continue;
        }
        const star_tracker::StarTracker::Solution s = tracker.process(centroids);
        std::cout << frame << ','
                  << (s.valid ? 1 : 0) << ','
                  << ((s.mode == star_tracker::StarTracker::kTracking)
                          ? "tracking" : "acquisition") << ','
                  << reasonText(s.reason) << ','
                  << s.clique_size << ','
                  << s.right_ascension_deg << ','
                  << s.declination_deg << ','
                  << s.roll_deg << ','
                  << s.residual_rms_arcsec << ','
                  << s.focal_used_pixels << ','
                  << s.focal_measured_pixels << ','
                  << s.focal_tracked_pixels << ','
                  << s.focal_innovation << ','
                  << (s.focal_updated ? 1 : 0) << ','
                  << (s.focal_jump_rejected ? 1 : 0) << ','
                  << s.trials << '\n';
        ++frame;
    }
    return 0;
}
