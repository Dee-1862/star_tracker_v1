#include "star_tracker/focal_tracker.hpp"

#include <cmath>

namespace star_tracker {

FocalTracker::FocalTracker()
    : config_(), tracked_(0.0F), initialised_(false), accepted_(0U),
      rejected_(0U) {}

FocalTracker::FocalTracker(const Config& config)
    : config_(config), tracked_(0.0F), initialised_(false), accepted_(0U),
      rejected_(0U) {}

void FocalTracker::reset(float initial_focal_pixels) {
    tracked_ = initial_focal_pixels;
    initialised_ = (initial_focal_pixels > 0.0F);
    accepted_ = 0U;
    rejected_ = 0U;
}

FocalTracker::Update FocalTracker::update(float measured_focal_pixels) {
    Update result;

    if (!initialised_ || (measured_focal_pixels <= 0.0F)) {
        result.reason = kNotInitialised;
        result.tracked_focal_pixels = tracked_;
        return result;
    }

    const float innovation =
        (measured_focal_pixels - tracked_) / tracked_;
    result.innovation_fraction = innovation;

    const bool warming = accepted_ < config_.warmup_frames;
    if (!warming &&
        (std::fabs(innovation) > config_.maximum_innovation_fraction)) {
        // Thermal drift cannot step this far between frames. Refuse the
        // update and leave the tracked value untouched -- absorbing it would
        // corrupt the very number the projection depends on.
        ++rejected_;
        result.accepted = false;
        result.reason = kImplausibleJump;
        result.tracked_focal_pixels = tracked_;
        return result;
    }

    tracked_ += config_.smoothing * (measured_focal_pixels - tracked_);
    ++accepted_;
    result.accepted = true;
    result.reason = warming ? kWarmup : kAccepted;
    result.tracked_focal_pixels = tracked_;
    return result;
}

float FocalTracker::focalLengthPixels() const { return tracked_; }

bool FocalTracker::converged() const {
    return initialised_ && (accepted_ >= config_.warmup_frames);
}

std::size_t FocalTracker::acceptedCount() const { return accepted_; }

std::size_t FocalTracker::rejectedCount() const { return rejected_; }

}  // namespace star_tracker
