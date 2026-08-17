#include "star_tracker/defect_map.hpp"

#include <cmath>

namespace star_tracker {

DefectMap::DefectMap() : config_(), entries_(), tracked_(0U), frames_(0U),
                         last_boresight_(), has_last_(false) {
    reset();
}

DefectMap::DefectMap(const Config& config)
    : config_(config), entries_(), tracked_(0U), frames_(0U),
      last_boresight_(), has_last_(false) {
    reset();
}

void DefectMap::reset() {
    for (std::size_t index = 0U; index < kCapacity; ++index) {
        entries_[index].key_x = kEmpty;
        entries_[index].key_y = kEmpty;
        entries_[index].score = 0;
    }
    tracked_ = 0U;
    frames_ = 0U;
    has_last_ = false;
    last_boresight_[0U] = 0.0F;
    last_boresight_[1U] = 0.0F;
    last_boresight_[2U] = 0.0F;
}

std::size_t DefectMap::slotFor(
    std::uint16_t key_x, std::uint16_t key_y) const {
    // Knuth multiplicative mixing, truncated to the table size. The table is a
    // power of two, so the mask is exact.
    const std::uint32_t mixed =
        ((static_cast<std::uint32_t>(key_x) * 2654435761UL) ^
         (static_cast<std::uint32_t>(key_y) * 2246822519UL));
    return static_cast<std::size_t>(mixed) & (kCapacity - 1U);
}

std::size_t DefectMap::findOrInsert(
    std::uint16_t key_x, std::uint16_t key_y) {
    std::size_t slot = slotFor(key_x, key_y);
    for (std::size_t probe = 0U; probe < kCapacity; ++probe) {
        Entry& entry = entries_[slot];
        if ((entry.key_x == key_x) && (entry.key_y == key_y)) {
            return slot;
        }
        if (entry.key_x == kEmpty) {
            entry.key_x = key_x;
            entry.key_y = key_y;
            entry.score = 0;
            ++tracked_;
            return slot;
        }
        slot = (slot + 1U) & (kCapacity - 1U);
    }
    return kCapacity;  // full: silently ignore rather than evict something real
}

void DefectMap::adjust(float x, float y, std::int16_t delta) {
    if ((x < 0.0F) || (y < 0.0F)) {
        return;
    }
    const std::uint16_t key_x = static_cast<std::uint16_t>(
        static_cast<std::uint32_t>(x) / config_.quantisation_px);
    const std::uint16_t key_y = static_cast<std::uint16_t>(
        static_cast<std::uint32_t>(y) / config_.quantisation_px);
    if ((key_x == kEmpty) || (key_y == kEmpty)) {
        return;
    }

    // Only spend a table slot on something that has failed at least once.
    // Clearing an untracked location is a no-op, which keeps the table
    // populated by suspects rather than by every star ever seen.
    if (delta < 0) {
        const std::size_t slot = slotFor(key_x, key_y);
        std::size_t probe_slot = slot;
        for (std::size_t probe = 0U; probe < kCapacity; ++probe) {
            Entry& entry = entries_[probe_slot];
            if (entry.key_x == kEmpty) {
                return;
            }
            if ((entry.key_x == key_x) && (entry.key_y == key_y)) {
                entry.score = static_cast<std::int16_t>(entry.score + delta);
                if (entry.score < config_.score_min) {
                    entry.score = config_.score_min;
                }
                return;
            }
            probe_slot = (probe_slot + 1U) & (kCapacity - 1U);
        }
        return;
    }

    const std::size_t slot = findOrInsert(key_x, key_y);
    if (slot >= kCapacity) {
        return;
    }
    Entry& entry = entries_[slot];
    entry.score = static_cast<std::int16_t>(entry.score + delta);
    if (entry.score > config_.score_max) {
        entry.score = config_.score_max;
    }
}

void DefectMap::decay() {
    for (std::size_t index = 0U; index < kCapacity; ++index) {
        Entry& entry = entries_[index];
        if (entry.key_x == kEmpty) {
            continue;
        }
        if (entry.score > 0) {
            --entry.score;
        } else if (entry.score < 0) {
            ++entry.score;
        }
    }
}

void DefectMap::observe(
    const Centroiding::Result& centroids,
    const LisGridMatcher::MatchResult& matches,
    const float boresight[3U],
    bool solved) {
    ++frames_;
    if ((frames_ % config_.decay_interval_frames) == 0U) {
        decay();
    }
    if (!solved) {
        // No clique means no evidence about which dots were real. Scoring here
        // would punish every star in a frame that failed for unrelated reasons.
        return;
    }

    // A dot that joined the clique is confirmed real: clear any suspicion,
    // regardless of how far the spacecraft has moved.
    const std::size_t count =
        (centroids.count < Centroiding::kMaxCentroids)
            ? centroids.count : Centroiding::kMaxCentroids;
    for (std::size_t index = 0U; index < count; ++index) {
        if (matches.star_ids[index] != 0U) {
            adjust(centroids.points[index].x, centroids.points[index].y,
                   static_cast<std::int16_t>(-config_.clear_bonus));
        }
    }

    // Misses are only informative once the sky has moved. At a fixed attitude
    // a real but unindexed star sits at a fixed pixel and is indistinguishable
    // from a defect, so scoring it would be a false accusation.
    bool moved = !has_last_;
    if (has_last_) {
        const float dot = (boresight[0U] * last_boresight_[0U]) +
                          (boresight[1U] * last_boresight_[1U]) +
                          (boresight[2U] * last_boresight_[2U]);
        const float clamped = (dot > 1.0F) ? 1.0F : ((dot < -1.0F) ? -1.0F : dot);
        const float degrees = std::acos(clamped) * 57.2957795130823F;
        moved = degrees >= config_.minimum_motion_deg;
    }

    if (moved) {
        for (std::size_t index = 0U; index < count; ++index) {
            if (matches.star_ids[index] == 0U) {
                adjust(centroids.points[index].x, centroids.points[index].y, 1);
            }
        }
        last_boresight_[0U] = boresight[0U];
        last_boresight_[1U] = boresight[1U];
        last_boresight_[2U] = boresight[2U];
        has_last_ = true;
    }
}

bool DefectMap::isDefect(float x, float y) const {
    if ((x < 0.0F) || (y < 0.0F)) {
        return false;
    }
    const std::uint16_t key_x = static_cast<std::uint16_t>(
        static_cast<std::uint32_t>(x) / config_.quantisation_px);
    const std::uint16_t key_y = static_cast<std::uint16_t>(
        static_cast<std::uint32_t>(y) / config_.quantisation_px);

    std::size_t slot = slotFor(key_x, key_y);
    for (std::size_t probe = 0U; probe < kCapacity; ++probe) {
        const Entry& entry = entries_[slot];
        if (entry.key_x == kEmpty) {
            return false;
        }
        if ((entry.key_x == key_x) && (entry.key_y == key_y)) {
            return entry.score >= config_.confirm_threshold;
        }
        slot = (slot + 1U) & (kCapacity - 1U);
    }
    return false;
}

Centroiding::Result DefectMap::filter(
    const Centroiding::Result& centroids) const {
    Centroiding::Result kept;
    const std::size_t count =
        (centroids.count < Centroiding::kMaxCentroids)
            ? centroids.count : Centroiding::kMaxCentroids;
    for (std::size_t index = 0U; index < count; ++index) {
        const Point& point = centroids.points[index];
        if (isDefect(point.x, point.y)) {
            continue;
        }
        kept.points[kept.count] = point;
        ++kept.count;
    }
    return kept;
}

std::size_t DefectMap::trackedCount() const { return tracked_; }

std::size_t DefectMap::confirmedCount() const {
    std::size_t total = 0U;
    for (std::size_t index = 0U; index < kCapacity; ++index) {
        if ((entries_[index].key_x != kEmpty) &&
            (entries_[index].score >= config_.confirm_threshold)) {
            ++total;
        }
    }
    return total;
}

std::size_t DefectMap::framesObserved() const { return frames_; }

}  // namespace star_tracker

static_assert(
    sizeof(star_tracker::DefectMap) < 4096U,
    "DefectMap must stay small; a full-frame score array would be 1 MiB");
