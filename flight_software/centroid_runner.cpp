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
    if (argc < 5 || argc > 7) {
        std::cerr << "Usage: centroid_runner CENTROIDS.tsv WIDTH HEIGHT "
                     "FOV_DEG [MAX_RESIDUAL_ARCSEC=30] [SEARCH_PCT=0]\n"
                     "\n"
                     "SEARCH_PCT enables focal-length acquisition: the solver\n"
                     "sweeps focal length over +/- SEARCH_PCT and accepts the\n"
                     "first value that yields a gated solve, then refines it\n"
                     "from the matched pairs. 0 disables the search.\n";
        return 2;
    }

    const std::string centroidsPath = argv[1];
    const int width = std::atoi(argv[2]);
    const int height = std::atoi(argv[3]);
    const double fovDeg = std::atof(argv[4]);
    const float maxResidual = (argc >= 6) ? static_cast<float>(std::atof(argv[5]))
                                          : 30.0F;
    const double searchPct = (argc >= 7) ? std::atof(argv[6]) : 0.0;

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

    // Focal-length acquisition. A wrong focal length yields no clique at all
    // rather than a wrong answer -- measured: 0 correct, 0 wrong, all refused
    // at +/-2% and +/-5%. That unambiguous failure signal is what makes a blind
    // sweep safe here and unsafe for a solver that false-solves.
    const Clock::time_point matchBegin = Clock::now();
    star_tracker::LisGridMatcher::MatchResult matches = matcher.match(centroids);
    star_tracker::AttitudeSolver::Solution solution = solver.solve(
        centroids, matches, star_tracker::kHipparcosCatalog,
        star_tracker::kHipparcosCatalogCount);

    std::size_t trials = 1U;
    float lockedFocal = focal;
    if (!solution.valid && (searchPct > 0.0)) {
        const std::size_t kSteps = 25U;
        const double step_fraction = (searchPct / 100.0) / 12.0;
        for (std::size_t step = 0U; (step < kSteps) && !solution.valid; ++step) {
            // Walk OUTWARD from nominal: 0, +1, -1, +2, -2, ... Thermal drift
            // is usually small, so the common case is found in a few trials
            // instead of at the far end of a linear scan. This also limits the
            // number of hypotheses actually tested, which bounds the
            // multiple-comparisons exposure of the search.
            const long index = static_cast<long>(step + 1U) / 2L;
            const double signed_index =
                ((step % 2U) == 1U) ? static_cast<double>(index)
                                    : -static_cast<double>(index);
            const float trialFocal = static_cast<float>(
                static_cast<double>(focal) * (1.0 + (signed_index * step_fraction)));

            // The catalogue pair index is focal-independent, so it is built
            // once outside this loop and only the bearing scale changes here.
            matcher.setFocalLength(trialFocal);

            star_tracker::AttitudeSolver::Config trialSolver = solverConfig;
            trialSolver.focal_length_pixels = trialFocal;
            star_tracker::AttitudeSolver trial(trialSolver);

            const star_tracker::LisGridMatcher::MatchResult trialMatches =
                matcher.match(centroids);
            const star_tracker::AttitudeSolver::Solution trialSolution =
                trial.solve(centroids, trialMatches,
                            star_tracker::kHipparcosCatalog,
                            star_tracker::kHipparcosCatalogCount);
            ++trials;
            if (trialSolution.valid) {
                matches = trialMatches;
                solution = trialSolution;
                lockedFocal = trialFocal;
            }
        }
    }

    // Turn the coarse lock into a precise value using every matched pair.
    float refinedFocal = lockedFocal;
    if (solution.valid) {
        refinedFocal = star_tracker::AttitudeSolver::refineFocalLength(
            centroids, matches, star_tracker::kHipparcosCatalog,
            star_tracker::kHipparcosCatalogCount,
            solverConfig.principal_x, solverConfig.principal_y, lockedFocal);

        // Re-solve at the refined focal length: the attitude and the residual
        // should both improve, and the residual is what we report as integrity.
        star_tracker::AttitudeSolver::Config finalConfig = solverConfig;
        finalConfig.focal_length_pixels = refinedFocal;
        star_tracker::AttitudeSolver finalSolver(finalConfig);
        const star_tracker::AttitudeSolver::Solution refinedSolution =
            finalSolver.solve(centroids, matches,
                              star_tracker::kHipparcosCatalog,
                              star_tracker::kHipparcosCatalogCount);
        if (refinedSolution.valid) {
            solution = refinedSolution;
        }
    }
    const Clock::time_point matchEnd = Clock::now();
    const Clock::time_point attitudeBegin = matchEnd;
    const Clock::time_point attitudeEnd = matchEnd;

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
    // index_ns is one-time startup work (catalogue pair index), NOT per-frame.
    // solve_ns is the per-frame cost that a frame-rate budget must cover:
    // shortlist, consistency graph, clique search, attitude, integrity gate,
    // and any focal search or refinement.
    std::cout << "index_ns " << indexNs << "\n";
    std::cout << "solve_ns " << matchNs << "\n";
    std::cout << "starid_average_ns " << matchNs << "\n";
    std::cout << "attitude_average_ns " << attitudeNs << "\n";
    std::cout << "total_average_ns " << (matchNs + attitudeNs) << "\n";
    std::cout << "clique_size " << matches.clique_size << "\n";
    std::cout << "node_count " << matches.node_count << "\n";
    std::cout << "expansions " << matches.expansions << "\n";
    std::cout << "focal_nominal_px " << focal << "\n";
    std::cout << "focal_locked_px " << lockedFocal << "\n";
    std::cout << "focal_refined_px " << refinedFocal << "\n";
    std::cout << "focal_trials " << trials << "\n";
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
