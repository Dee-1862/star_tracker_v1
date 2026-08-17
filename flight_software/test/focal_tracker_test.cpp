#include "star_tracker/focal_tracker.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

int main() {
    const float truth = 2903.7F;

    // --- 1. Noise is averaged down, not tracked ----------------------------
    {
        star_tracker::FocalTracker tracker;
        tracker.reset(truth);
        // Alternating +/-0.3% noise about the truth: a filter that works
        // should sit far closer to truth than any individual sample.
        for (std::size_t frame = 0U; frame < 200U; ++frame) {
            const float sign = ((frame % 2U) == 0U) ? 1.0F : -1.0F;
            tracker.update(truth * (1.0F + (sign * 0.003F)));
        }
        const float error = std::fabs(tracker.focalLengthPixels() - truth) / truth;
        assert(error < 0.0005F);
        assert(tracker.converged());
        assert(tracker.rejectedCount() == 0U);
    }

    // --- 2. Slow thermal drift is followed ---------------------------------
    {
        star_tracker::FocalTracker tracker;
        tracker.reset(truth);
        // 0.5% drift applied gradually, as a warming lens actually behaves.
        const float target = truth * 1.005F;
        for (std::size_t frame = 0U; frame < 600U; ++frame) {
            const float fraction = static_cast<float>(frame) / 600.0F;
            tracker.update(truth + ((target - truth) * fraction));
        }
        const float error = std::fabs(tracker.focalLengthPixels() - target) / target;
        assert(error < 0.002F);
        assert(tracker.rejectedCount() == 0U);
    }

    // --- 3. A step is a fault, and must be refused, not absorbed -----------
    // This is the property that makes the filter an integrity monitor: a bad
    // lock or corrupted star set shows up as a discontinuity, and absorbing it
    // would corrupt the projection every later frame depends on.
    {
        star_tracker::FocalTracker tracker;
        tracker.reset(truth);
        for (std::size_t frame = 0U; frame < 40U; ++frame) {
            tracker.update(truth);
        }
        const float before = tracker.focalLengthPixels();

        const star_tracker::FocalTracker::Update bad =
            tracker.update(truth * 1.08F);
        assert(!bad.accepted);
        assert(bad.reason == star_tracker::FocalTracker::kImplausibleJump);
        assert(bad.innovation_fraction > 0.07F);
        assert(std::fabs(tracker.focalLengthPixels() - before) < 1e-3F);
        assert(tracker.rejectedCount() == 1U);

        // Good data afterwards is still accepted: one bad frame must not
        // wedge the filter.
        const star_tracker::FocalTracker::Update good = tracker.update(truth);
        assert(good.accepted);
    }

    // --- 4. Warmup accepts what a converged filter would reject ------------
    {
        star_tracker::FocalTracker tracker;
        tracker.reset(truth);
        const star_tracker::FocalTracker::Update first =
            tracker.update(truth * 1.05F);
        assert(first.accepted);
        assert(first.reason == star_tracker::FocalTracker::kWarmup);
        assert(!tracker.converged());
    }

    // --- 5. Uninitialised and invalid input are refused, not guessed -------
    {
        star_tracker::FocalTracker tracker;
        const star_tracker::FocalTracker::Update none = tracker.update(truth);
        assert(!none.accepted);
        assert(none.reason == star_tracker::FocalTracker::kNotInitialised);

        tracker.reset(truth);
        const star_tracker::FocalTracker::Update negative = tracker.update(-1.0F);
        assert(!negative.accepted);
    }

    return 0;
}
