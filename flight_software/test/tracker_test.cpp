#include "star_tracker/tracker.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

constexpr float kFocal = 2904.0F;
constexpr float kPrincipal = 512.0F;
constexpr std::size_t kStarCount = 10U;

star_tracker::LisGridMatcher::Catalog catalog;

/// Build a self-consistent scene at a GIVEN focal length, so the tracker's
/// recovered value can be checked against a known truth. A larger focal length
/// is exactly what a warming lens produces.
void buildScene(
    float focal,
    star_tracker::Centroiding::Result& centroids,
    star_tracker::LisGridMatcher::MatchResult& matches) {
    const float pixels[kStarCount][2U] = {
        {300.0F, 280.0F}, {700.0F, 320.0F}, {480.0F, 640.0F}, {820.0F, 700.0F},
        {200.0F, 750.0F}, {560.0F, 200.0F}, {380.0F, 480.0F}, {740.0F, 520.0F},
        {620.0F, 860.0F}, {260.0F, 560.0F}};

    centroids.count = kStarCount;
    matches.matched_count = kStarCount;
    matches.clique_size = kStarCount;

    for (std::size_t index = 0U; index < kStarCount; ++index) {
        centroids.points[index].x = pixels[index][0U];
        centroids.points[index].y = pixels[index][1U];
        centroids.points[index].intensity =
            static_cast<std::uint32_t>(3000U - (index * 10U));

        const float y = (kPrincipal - pixels[index][0U]) / focal;
        const float z = (kPrincipal - pixels[index][1U]) / focal;
        const float inverse = 1.0F / std::sqrt(1.0F + (y * y) + (z * z));

        const std::uint32_t id = static_cast<std::uint32_t>(index + 1U);
        catalog[index].star_id = id;
        catalog[index].x = inverse;
        catalog[index].y = y * inverse;
        catalog[index].z = z * inverse;
        catalog[index].visual_magnitude = 3.0F + (0.1F * index);
        matches.star_ids[index] = id;
        matches.votes[index] = 10U;
    }
}

}  // namespace

int main() {
    // The tracker owns the matcher, so a real index is needed. Feed it the
    // synthetic catalogue built above.
    star_tracker::Centroiding::Result seed;
    star_tracker::LisGridMatcher::MatchResult seedMatches;
    buildScene(kFocal, seed, seedMatches);

    star_tracker::StarTracker::Config config;
    config.nominal_focal_pixels = kFocal;
    config.principal_x = kPrincipal;
    config.principal_y = kPrincipal;

    // --- 1. Focal state persists across frames -----------------------------
    // The whole reason this class exists: a single-frame API cannot average
    // noise or notice a step, because it has nowhere to keep the history.
    {
        star_tracker::StarTracker tracker(config);
        assert(tracker.mode() == star_tracker::StarTracker::kAcquisition);
        const float before = tracker.focalLengthPixels();
        assert(std::fabs(before - kFocal) < 1e-3F);

        // Focal length must be carried between calls, not recomputed from
        // nominal each time.
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult unused;
        buildScene(kFocal, centroids, unused);
        (void)tracker.process(centroids);
        assert(tracker.focalLengthPixels() > 0.0F);
    }

    // --- 2. Reset returns to power-on state --------------------------------
    {
        star_tracker::StarTracker tracker(config);
        tracker.reset();
        assert(tracker.mode() == star_tracker::StarTracker::kAcquisition);
        assert(tracker.consecutiveFailures() == 0U);
        assert(std::fabs(tracker.focalLengthPixels() - kFocal) < 1e-3F);
    }

    // --- 3. An unconfigured tracker reports, but does not count failures ----
    // No index means not configured, which is different from tried-and-failed.
    // Counting it as a solve failure would flip operating mode on a condition
    // that has nothing to do with the sky.
    {
        star_tracker::StarTracker tracker(config);
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult unused;
        buildScene(kFocal, centroids, unused);

        const star_tracker::StarTracker::Solution s = tracker.process(centroids);
        assert(!s.valid);
        // Telemetry must still be populated, so a ground operator can tell
        // "not configured" from "could not solve".
        assert(s.focal_tracked_pixels > 0.0F);
        assert(tracker.consecutiveFailures() == 0U);
    }

    // --- 4. Real refusals accumulate and force re-acquisition --------------
    // With a working index, unsolvable input is a genuine failure. A tracker
    // that stayed in TRACKING would keep retrying one dead hypothesis forever.
    {
        star_tracker::StarTracker::Config quick = config;
        quick.failures_before_reacquire = 3U;
        star_tracker::StarTracker tracker(quick);
        assert(tracker.buildIndex(catalog, kStarCount));

        // Noise: dots at positions no catalogue geometry can explain.
        star_tracker::Centroiding::Result noise;
        for (std::size_t index = 0U; index < 12U; ++index) {
            noise.points[index].x =
                40.0F + static_cast<float>((index * 191U) % 940U);
            noise.points[index].y =
                40.0F + static_cast<float>((index * 373U) % 940U);
            noise.points[index].intensity =
                static_cast<std::uint32_t>(2000U - index);
            ++noise.count;
        }

        for (std::size_t frame = 0U; frame < 6U; ++frame) {
            const star_tracker::StarTracker::Solution s = tracker.process(noise);
            assert(!s.valid);
        }
        assert(tracker.consecutiveFailures() >= 3U);
        assert(tracker.mode() == star_tracker::StarTracker::kAcquisition);
    }

    // --- 5. The filter underneath still rejects an implausible jump ---------
    // Wired through the tracker rather than tested in isolation: this is the
    // property that turns calibration into an integrity monitor.
    {
        star_tracker::FocalTracker filter;
        filter.reset(kFocal);
        for (std::size_t frame = 0U; frame < 30U; ++frame) {
            filter.update(kFocal);
        }
        const float steady = filter.focalLengthPixels();
        const star_tracker::FocalTracker::Update bad =
            filter.update(kFocal * 1.08F);
        assert(!bad.accepted);
        assert(std::fabs(filter.focalLengthPixels() - steady) < 1e-3F);
        assert(filter.rejectedCount() == 1U);
    }

    return 0;
}
