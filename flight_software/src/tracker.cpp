#include "star_tracker/tracker.hpp"

namespace star_tracker {

StarTracker::StarTracker()
    : config_(),
      matcher_(),
      focal_(),
      catalog_(0),
      catalog_count_(0U),
      mode_(kAcquisition),
      consecutive_failures_(0U),
      index_ready_(false) {
    focal_.reset(config_.nominal_focal_pixels);
}

StarTracker::StarTracker(const Config& config)
    : config_(config),
      matcher_(),
      focal_(),
      catalog_(0),
      catalog_count_(0U),
      mode_(kAcquisition),
      consecutive_failures_(0U),
      index_ready_(false) {
    focal_.reset(config_.nominal_focal_pixels);
}

bool StarTracker::buildIndex(
    const LisGridMatcher::Catalog& catalog, std::size_t catalog_count) {
    LisGridMatcher::Config matcher_config;
    matcher_config.focal_length_pixels = config_.nominal_focal_pixels;
    matcher_config.principal_x = config_.principal_x;
    matcher_config.principal_y = config_.principal_y;
    matcher_config.minimum_clique_size = config_.minimum_clique_size;
    matcher_ = LisGridMatcher(matcher_config);

    catalog_ = &catalog;
    catalog_count_ = catalog_count;
    index_ready_ = matcher_.buildIndex(catalog, catalog_count);
    reset();
    return index_ready_;
}

void StarTracker::reset() {
    focal_.reset(config_.nominal_focal_pixels);
    mode_ = kAcquisition;
    consecutive_failures_ = 0U;
}

StarTracker::Mode StarTracker::mode() const { return mode_; }

float StarTracker::focalLengthPixels() const {
    return focal_.focalLengthPixels();
}

std::size_t StarTracker::consecutiveFailures() const {
    return consecutive_failures_;
}

bool StarTracker::attempt(
    const Centroiding::Result& centroids,
    float focal_pixels,
    LisGridMatcher::MatchResult& matches,
    AttitudeSolver::Solution& solution) {
    // The pair index is pure catalogue geometry, so only the bearing scale
    // changes between trials. Rebuilding it here would dominate the sweep.
    matcher_.setFocalLength(focal_pixels);

    AttitudeSolver::Config solver_config;
    solver_config.focal_length_pixels = focal_pixels;
    solver_config.principal_x = config_.principal_x;
    solver_config.principal_y = config_.principal_y;
    solver_config.max_residual_arcsec = config_.max_residual_arcsec;
    AttitudeSolver solver(solver_config);

    matches = matcher_.match(centroids);
    solution = solver.solve(centroids, matches, *catalog_, catalog_count_);
    return solution.valid;
}

StarTracker::Solution StarTracker::process(
    const Centroiding::Result& centroids) {
    Solution result;
    result.mode = mode_;
    if (!index_ready_ || (catalog_ == 0)) {
        result.reason = AttitudeSolver::kTooFewStars;
        result.focal_tracked_pixels = focal_.focalLengthPixels();
        result.focal_used_pixels = result.focal_tracked_pixels;
        return result;
    }

    const float tracked = focal_.focalLengthPixels();
    LisGridMatcher::MatchResult matches;
    AttitudeSolver::Solution solution;
    float used = tracked;

    // Always try the tracked focal length first. In tracking mode that is the
    // whole cost; in acquisition it is the cheapest hypothesis and often the
    // right one.
    result.trials = 1U;
    bool solved = attempt(centroids, tracked, matches, solution);

    if (!solved && (mode_ == kAcquisition) && (config_.search_percent > 0.0F)) {
        const std::size_t steps = config_.search_steps;
        const float span = config_.search_percent / 100.0F;
        const float half = static_cast<float>(steps / 2U);
        const float step = (half > 0.0F) ? (span / half) : span;

        for (std::size_t index = 0U; (index < steps) && !solved; ++index) {
            // Walk outward from the current estimate: 0, +1, -1, +2, -2 ...
            // Drift is usually small, so the common case resolves in a few
            // trials rather than at the far end of a linear scan.
            const long magnitude = static_cast<long>(index + 1U) / 2L;
            const float direction = ((index % 2U) == 1U) ? 1.0F : -1.0F;
            const float trial =
                tracked * (1.0F + (direction * static_cast<float>(magnitude) * step));
            if (trial <= 0.0F) {
                continue;
            }
            ++result.trials;
            if (attempt(centroids, trial, matches, solution)) {
                solved = true;
                used = trial;
            }
        }
    }

    result.focal_used_pixels = used;
    result.residual_rms_arcsec = solution.residual_rms_arcsec;
    result.clique_size = matches.clique_size;
    result.matched_count = solution.matched_count;
    result.reason = solution.reason;

    if (!solved) {
        ++consecutive_failures_;
        // A single refusal is normal. A run of them means the calibration or
        // the scene assumption is wrong, so drop back to searching rather than
        // repeating a hypothesis that has stopped working.
        if ((mode_ == kTracking) &&
            (consecutive_failures_ >= config_.failures_before_reacquire)) {
            mode_ = kAcquisition;
        }
        result.focal_tracked_pixels = focal_.focalLengthPixels();
        result.mode = mode_;
        return result;
    }

    // Solved. Recover focal length from the matched pairs and feed the filter.
    const float measured = AttitudeSolver::refineFocalLength(
        centroids, matches, *catalog_, catalog_count_,
        config_.principal_x, config_.principal_y, used);
    result.focal_measured_pixels = measured;

    const FocalTracker::Update update = focal_.update(measured);
    result.focal_tracked_pixels = update.tracked_focal_pixels;
    result.focal_innovation = update.innovation_fraction;
    result.focal_updated = update.accepted;
    result.focal_jump_rejected =
        (update.reason == FocalTracker::kImplausibleJump);

    result.valid = true;
    result.quaternion = solution.quaternion;
    result.right_ascension_deg = solution.right_ascension_deg;
    result.declination_deg = solution.declination_deg;
    result.roll_deg = solution.roll_deg;

    consecutive_failures_ = 0U;
    mode_ = kTracking;
    result.mode = mode_;
    return result;
}

}  // namespace star_tracker
