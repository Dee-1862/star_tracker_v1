#include "star_tracker/defect_map.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

/// A boresight that sweeps across the sky, so successive frames look at
/// different stars. Motion is what distinguishes a defect from a star.
void boresightAt(float degrees, float out[3U]) {
    const float radians = degrees * 0.01745329252F;
    out[0U] = std::cos(radians);
    out[1U] = std::sin(radians);
    out[2U] = 0.0F;
}

/// One frame: a stuck pixel always at the same place, plus real stars that
/// move, of which the named count joined the clique.
void makeFrame(
    float defect_x, float defect_y,
    float star_x, float star_y,
    bool star_matched,
    star_tracker::Centroiding::Result& centroids,
    star_tracker::LisGridMatcher::MatchResult& matches) {
    centroids = star_tracker::Centroiding::Result();
    matches = star_tracker::LisGridMatcher::MatchResult();

    centroids.points[0U].x = defect_x;
    centroids.points[0U].y = defect_y;
    centroids.points[0U].intensity = 5000U;
    matches.star_ids[0U] = 0U;  // never joins: it is not a star

    centroids.points[1U].x = star_x;
    centroids.points[1U].y = star_y;
    centroids.points[1U].intensity = 3000U;
    matches.star_ids[1U] = star_matched ? 4242U : 0U;

    centroids.count = 2U;
    matches.clique_size = star_matched ? 8U : 0U;
}

}  // namespace

int main() {
    const float kDefectX = 640.0F;
    const float kDefectY = 200.0F;

    // --- 1. A stuck pixel is learned; a moving star is not ------------------
    {
        star_tracker::DefectMap map;
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        float boresight[3U];

        for (std::size_t frame = 0U; frame < 40U; ++frame) {
            boresightAt(static_cast<float>(frame) * 2.0F, boresight);
            // The star lands somewhere different every frame, as the sky moves.
            const float sx = 100.0F + (static_cast<float>(frame) * 17.0F);
            makeFrame(kDefectX, kDefectY, sx, 400.0F, true, centroids, matches);
            map.observe(centroids, matches, boresight, true);
        }

        assert(map.isDefect(kDefectX, kDefectY));
        assert(map.confirmedCount() == 1U);
        // Stars that moved must not be condemned, even though each was seen at
        // its own location only once.
        assert(!map.isDefect(100.0F, 400.0F));
        assert(!map.isDefect(300.0F, 400.0F));

        // Filtering removes the defect and keeps the star.
        makeFrame(kDefectX, kDefectY, 500.0F, 400.0F, true, centroids, matches);
        const star_tracker::Centroiding::Result kept = map.filter(centroids);
        assert(kept.count == 1U);
        assert(std::fabs(kept.points[0U].x - 500.0F) < 1e-3F);
    }

    // --- 2. A STATIONARY spacecraft must not condemn a real star -----------
    // The critical false-positive case. At a fixed attitude an unindexed real
    // star sits at a fixed pixel and looks exactly like a defect. Scoring it
    // would blind the tracker to a genuine star for the rest of the mission.
    {
        star_tracker::DefectMap map;
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        float boresight[3U];
        boresightAt(0.0F, boresight);  // never moves

        for (std::size_t frame = 0U; frame < 200U; ++frame) {
            makeFrame(kDefectX, kDefectY, 500.0F, 400.0F, true,
                      centroids, matches);
            map.observe(centroids, matches, boresight, true);
        }
        assert(!map.isDefect(kDefectX, kDefectY));
    }

    // --- 3. Unsolved frames carry no evidence ------------------------------
    {
        star_tracker::DefectMap map;
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        float boresight[3U];

        for (std::size_t frame = 0U; frame < 100U; ++frame) {
            boresightAt(static_cast<float>(frame) * 3.0F, boresight);
            makeFrame(kDefectX, kDefectY, 500.0F, 400.0F, false,
                      centroids, matches);
            map.observe(centroids, matches, boresight, false);
        }
        assert(!map.isDefect(kDefectX, kDefectY));
        assert(map.framesObserved() == 100U);
    }

    // --- 4. A healed pixel is forgiven, in bounded time --------------------
    {
        star_tracker::DefectMap map;
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        float boresight[3U];

        for (std::size_t frame = 0U; frame < 40U; ++frame) {
            boresightAt(static_cast<float>(frame) * 2.0F, boresight);
            makeFrame(kDefectX, kDefectY, 100.0F + (frame * 13.0F), 400.0F,
                      true, centroids, matches);
            map.observe(centroids, matches, boresight, true);
        }
        assert(map.isDefect(kDefectX, kDefectY));

        // The pixel recovers and now behaves like a star, joining every clique.
        for (std::size_t frame = 0U; frame < 40U; ++frame) {
            boresightAt(static_cast<float>(frame) * 2.0F, boresight);
            centroids = star_tracker::Centroiding::Result();
            matches = star_tracker::LisGridMatcher::MatchResult();
            centroids.points[0U].x = kDefectX;
            centroids.points[0U].y = kDefectY;
            centroids.points[0U].intensity = 3000U;
            matches.star_ids[0U] = 777U;
            matches.clique_size = 8U;
            centroids.count = 1U;
            map.observe(centroids, matches, boresight, true);
        }
        assert(!map.isDefect(kDefectX, kDefectY));
    }

    // --- 5. Bounded memory: a full table degrades, it does not corrupt ------
    {
        star_tracker::DefectMap map;
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        float boresight[3U];

        for (std::size_t frame = 0U; frame < 600U; ++frame) {
            boresightAt(static_cast<float>(frame) * 1.7F, boresight);
            makeFrame(static_cast<float>((frame * 37U) % 1000U) + 5.0F,
                      static_cast<float>((frame * 53U) % 1000U) + 5.0F,
                      500.0F, 400.0F, true, centroids, matches);
            map.observe(centroids, matches, boresight, true);
        }
        assert(map.trackedCount() <= star_tracker::DefectMap::kCapacity);
    }

    return 0;
}
