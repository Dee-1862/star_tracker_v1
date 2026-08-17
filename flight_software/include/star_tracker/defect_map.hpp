#ifndef STAR_TRACKER_DEFECT_MAP_HPP
#define STAR_TRACKER_DEFECT_MAP_HPP

#include "star_tracker/centroiding.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace star_tracker {

/// Learns which sensor locations are damaged, from data the pipeline already
/// discards.
///
/// Radiation gradually creates stuck-bright pixels. Isolated ones are rejected
/// by the centroider's minimum blob size, and streaks by its elongation limit,
/// but a compact hot *cluster* of a few pixels is indistinguishable from a star
/// in any single frame. No per-frame filter can catch it; the information is
/// simply not there.
///
/// Across frames it is obvious: **the sky moves and a defect does not.** The
/// clique already says which dots are real, so a dot that repeatedly fails to
/// join at the same sensor location, while the spacecraft points elsewhere, is
/// damage rather than an unlucky star.
///
/// Two design constraints worth stating:
///
///   - A full 1024x1024 score array would be 1 MiB, most of the RAM budget, to
///     track at most 30 centroids per frame. This uses a small open-addressed
///     hash table on the quantised coordinate instead: a few kilobytes.
///   - Scores saturate and decay rather than accumulating without bound. An
///     unbounded counter means a pixel bad for 10,000 frames takes 10,000
///     frames to forgive. Bounded scores give a known response time in both
///     directions.
///
/// The critical correctness condition is the attitude gate in observe(): a real
/// star at a fixed attitude also sits at a fixed pixel, so accumulating while
/// the spacecraft is inertially still would condemn perfectly good stars. Score
/// is only added once the boresight has actually moved.
class DefectMap {
public:
    static constexpr std::size_t kCapacity = 256U;
    static constexpr std::uint16_t kEmpty = 0xFFFFU;

    struct Config {
        /// Bucket size in pixels. Centroid position jitters slightly between
        /// frames, so a defect must land in the same bucket despite noise.
        std::uint16_t quantisation_px;
        /// Score at which a location is treated as damaged.
        std::int16_t confirm_threshold;
        std::int16_t score_max;
        std::int16_t score_min;
        /// Penalty subtracted when a dot at this location DOES join a clique.
        /// Larger than the increment, so one confirmed sighting outweighs
        /// several misses and a healed pixel clears quickly.
        std::int16_t clear_bonus;
        /// Frames between decay passes, which pull every score one step toward
        /// zero. Bounds how long a stale verdict can persist.
        std::size_t decay_interval_frames;
        /// Minimum boresight motion, degrees, before misses are scored. Below
        /// this the sky has not moved relative to the sensor and a real star is
        /// indistinguishable from a defect.
        float minimum_motion_deg;

        Config()
            : quantisation_px(4U),
              confirm_threshold(12),
              score_max(48),
              score_min(-8),
              clear_bonus(3),
              decay_interval_frames(512U),
              minimum_motion_deg(0.5F) {}
    };

    struct Entry {
        std::uint16_t key_x;
        std::uint16_t key_y;
        std::int16_t score;
    };

    DefectMap();
    explicit DefectMap(const Config& config);

    void reset();

    /// Fold one solved frame into the map.
    ///
    /// @param centroids  every dot found this frame
    /// @param matches    which of them joined the clique
    /// @param boresight  unit boresight in the inertial frame, used to tell
    ///                   whether the sky has moved since the last update
    /// @param solved     only solved frames are informative: without a clique
    ///                   there is no evidence about which dots were real
    void observe(
        const Centroiding::Result& centroids,
        const LisGridMatcher::MatchResult& matches,
        const float boresight[3U],
        bool solved);

    /// True once a location has been seen to fail repeatedly.
    bool isDefect(float x, float y) const;

    /// Drop centroids at known-bad locations.
    Centroiding::Result filter(const Centroiding::Result& centroids) const;

    std::size_t trackedCount() const;
    std::size_t confirmedCount() const;
    std::size_t framesObserved() const;

private:
    std::size_t slotFor(std::uint16_t key_x, std::uint16_t key_y) const;
    std::size_t findOrInsert(std::uint16_t key_x, std::uint16_t key_y);
    void adjust(float x, float y, std::int16_t delta);
    void decay();

    Config config_;
    std::array<Entry, kCapacity> entries_;
    std::size_t tracked_;
    std::size_t frames_;
    float last_boresight_[3U];
    bool has_last_;
};

}  // namespace star_tracker

#endif
