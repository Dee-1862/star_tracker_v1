#ifndef STAR_TRACKER_FOCAL_TRACKER_HPP
#define STAR_TRACKER_FOCAL_TRACKER_HPP

#include <cstddef>

namespace star_tracker {

/// Tracks focal length across frames, and rejects implausible jumps.
///
/// AttitudeSolver::refineFocalLength produces an independent estimate from each
/// solved frame. On its own that is a measurement, not a state: it carries no
/// history, so noise appears at full amplitude and a single bad frame is
/// indistinguishable from a real change.
///
/// A lens drifts with temperature slowly and smoothly. So:
///   - low-pass filtering the per-frame estimates turns them into a tracked
///     value whose noise falls as more frames arrive, and
///   - a *step* in the estimate is diagnostic, not physical. Thermal drift
///     cannot move focal length several percent between consecutive frames, so
///     an innovation that large means a bad lock, a corrupted star set, or a
///     genuine anomaly -- and should be refused rather than absorbed.
///
/// The second property is why this is an integrity monitor and not merely a
/// smoother: the same residual that makes the filter work also detects faults,
/// at no extra cost.
class FocalTracker {
public:
    struct Config {
        /// Weight given to each new measurement. 0.02 reaches ~63% of a step
        /// in 50 frames -- far slower than thermal drift, so it averages noise
        /// hard without lagging the physics.
        float smoothing;

        /// Reject an update whose disagreement with the tracked value exceeds
        /// this fraction. Thermal drift between adjacent frames is orders of
        /// magnitude smaller; anything near 1% is a fault signature.
        float maximum_innovation_fraction;

        /// Frames accepted before the innovation test is armed. Until then the
        /// tracked value is still converging and would reject valid data.
        std::size_t warmup_frames;

        Config()
            : smoothing(0.02F),
              maximum_innovation_fraction(0.01F),
              warmup_frames(8U) {}
    };

    enum Reason {
        kAccepted = 0,
        kWarmup = 1,
        kImplausibleJump = 2,
        kNotInitialised = 3
    };

    struct Update {
        bool accepted;
        float tracked_focal_pixels;
        /// Signed (measured - tracked) / tracked. Health telemetry: its
        /// running trend is the lens thermal drift rate.
        float innovation_fraction;
        Reason reason;

        Update()
            : accepted(false),
              tracked_focal_pixels(0.0F),
              innovation_fraction(0.0F),
              reason(kNotInitialised) {}
    };

    FocalTracker();
    explicit FocalTracker(const Config& config);

    void reset(float initial_focal_pixels);
    Update update(float measured_focal_pixels);

    float focalLengthPixels() const;
    bool converged() const;
    std::size_t acceptedCount() const;
    std::size_t rejectedCount() const;

private:
    Config config_;
    float tracked_;
    bool initialised_;
    std::size_t accepted_;
    std::size_t rejected_;
};

}  // namespace star_tracker

#endif
