#ifndef STAR_TRACKER_TRACKER_HPP
#define STAR_TRACKER_TRACKER_HPP

#include "star_tracker/attitude_solver.hpp"
#include "star_tracker/centroiding.hpp"
#include "star_tracker/focal_tracker.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <cstddef>
#include <cstdint>

namespace star_tracker {

/// Top-level pipeline: centroids in, gated attitude out, calibration maintained
/// across frames.
///
/// This is the object that was missing. Matching, attitude and the integrity
/// gate are all single-frame operations, so each could be tested in isolation,
/// but focal-length tracking is inherently multi-frame: it needs somewhere to
/// keep state between images. Without an owner for that state the tracker had
/// to re-derive focal length from scratch every frame, which throws away the
/// noise averaging and makes a step change indistinguishable from measurement
/// scatter.
///
/// Owning the state here also makes the two operating modes real rather than
/// documentary:
///
///   ACQUISITION  focal length unknown or lost. Sweeps focal length, budgeted
///                but not hard real-time. Entered at power-on and after
///                repeated failures.
///   TRACKING     focal length known. One matcher call per frame, fixed cost,
///                suitable for a hard deadline. Each accepted frame refines the
///                focal estimate, so thermal drift is followed rather than
///                fought.
///
/// Centroiding is deliberately NOT owned here: it has no cross-frame state, and
/// keeping it outside lets the same object serve both live images and the
/// frozen-centroid benchmark harness.
class StarTracker {
public:
    enum Mode { kAcquisition = 0, kTracking = 1 };

    struct Config {
        float nominal_focal_pixels;
        float principal_x;
        float principal_y;

        /// Half-width of the acquisition sweep, percent. Covers the plausible
        /// thermal range of the lens.
        float search_percent;
        std::size_t search_steps;

        /// Consecutive refusals in tracking mode before falling back to
        /// acquisition. One bad frame is normal; a run of them means the
        /// calibration or the scene assumption is wrong.
        std::size_t failures_before_reacquire;

        float max_residual_arcsec;
        std::size_t minimum_clique_size;

        Config()
            : nominal_focal_pixels(2904.0F),
              principal_x(512.0F),
              principal_y(512.0F),
              search_percent(6.0F),
              search_steps(25U),
              failures_before_reacquire(5U),
              max_residual_arcsec(30.0F),
              minimum_clique_size(5U) {}
    };

    struct Solution {
        bool valid;
        std::array<float, 4U> quaternion;
        float right_ascension_deg;
        float declination_deg;
        float roll_deg;
        float residual_rms_arcsec;
        std::size_t clique_size;
        std::size_t matched_count;
        AttitudeSolver::Reason reason;

        /// Calibration state, reported every frame as health telemetry.
        Mode mode;
        float focal_used_pixels;
        float focal_measured_pixels;
        float focal_tracked_pixels;
        /// Signed (measured - tracked) / tracked for this frame. Its running
        /// trend is the lens thermal drift rate; a spike is a fault signature.
        float focal_innovation;
        bool focal_updated;
        bool focal_jump_rejected;
        std::size_t trials;

        Solution()
            : valid(false),
              quaternion(),
              right_ascension_deg(0.0F),
              declination_deg(0.0F),
              roll_deg(0.0F),
              residual_rms_arcsec(-1.0F),
              clique_size(0U),
              matched_count(0U),
              reason(AttitudeSolver::kTooFewStars),
              mode(kAcquisition),
              focal_used_pixels(0.0F),
              focal_measured_pixels(0.0F),
              focal_tracked_pixels(0.0F),
              focal_innovation(0.0F),
              focal_updated(false),
              focal_jump_rejected(false),
              trials(0U) {
            quaternion[3U] = 1.0F;
        }
    };

    StarTracker();
    explicit StarTracker(const Config& config);

    /// Build the catalogue pair index. One-time; the index is independent of
    /// focal length, so acquisition never rebuilds it.
    bool buildIndex(
        const LisGridMatcher::Catalog& catalog, std::size_t catalog_count);

    /// Return to power-on state: focal length back to nominal, mode back to
    /// acquisition, tracking history discarded.
    void reset();

    Solution process(const Centroiding::Result& centroids);

    Mode mode() const;
    float focalLengthPixels() const;
    std::size_t consecutiveFailures() const;

private:
    /// One matcher plus solver call at a given focal length.
    bool attempt(
        const Centroiding::Result& centroids,
        float focal_pixels,
        LisGridMatcher::MatchResult& matches,
        AttitudeSolver::Solution& solution);

    Config config_;
    LisGridMatcher matcher_;
    FocalTracker focal_;
    const LisGridMatcher::Catalog* catalog_;
    std::size_t catalog_count_;
    Mode mode_;
    std::size_t consecutive_failures_;
    bool index_ready_;
};

}  // namespace star_tracker

#endif
