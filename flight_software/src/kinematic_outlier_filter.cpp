#include "star_tracker/kinematic_outlier_filter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace star_tracker {

KinematicOutlierFilter::KinematicOutlierFilter()
    : config_(), history_(), next_history_index_(0U), history_count_(0U) {}

KinematicOutlierFilter::KinematicOutlierFilter(const Config& config)
    : config_(config), history_(), next_history_index_(0U), history_count_(0U) {}

void KinematicOutlierFilter::reset() {
    next_history_index_ = 0U;
    history_count_ = 0U;
    for (std::size_t index = 0U; index < kHistoryFrames; ++index) {
        history_[index].count = 0U;
    }
}

void KinematicOutlierFilter::matchFrames(
    const std::array<Point, Centroiding::kMaxCentroids>& current,
    std::size_t current_count,
    const std::array<Point, Centroiding::kMaxCentroids>& previous,
    std::size_t previous_count,
    float maximum_distance,
    MatchArray& matches) {
    std::array<bool, Centroiding::kMaxCentroids> previous_used = {};
    for (std::size_t index = 0U; index < Centroiding::kMaxCentroids; ++index) {
        matches[index] = -1;
    }

    const float maximum_distance_squared =
        maximum_distance * maximum_distance;
    const std::size_t maximum_pairs = std::min(current_count, previous_count);
    for (std::size_t pair = 0U; pair < maximum_pairs; ++pair) {
        float best_distance_squared = maximum_distance_squared;
        std::int16_t best_current = -1;
        std::int16_t best_previous = -1;

        for (std::size_t current_index = 0U;
             current_index < current_count; ++current_index) {
            if (matches[current_index] >= 0) {
                continue;
            }
            for (std::size_t previous_index = 0U;
                 previous_index < previous_count; ++previous_index) {
                if (previous_used[previous_index]) {
                    continue;
                }

                const float delta_x =
                    current[current_index].x - previous[previous_index].x;
                const float delta_y =
                    current[current_index].y - previous[previous_index].y;
                const float distance_squared =
                    delta_x * delta_x + delta_y * delta_y;
                if (distance_squared <= best_distance_squared) {
                    best_distance_squared = distance_squared;
                    best_current = static_cast<std::int16_t>(current_index);
                    best_previous = static_cast<std::int16_t>(previous_index);
                }
            }
        }

        if (best_current < 0 || best_previous < 0) {
            break;
        }
        matches[static_cast<std::size_t>(best_current)] = best_previous;
        previous_used[static_cast<std::size_t>(best_previous)] = true;
    }
}

float KinematicOutlierFilter::median(
    const std::array<float, Centroiding::kMaxCentroids>& values,
    std::size_t count) {
    std::array<float, Centroiding::kMaxCentroids> sorted = values;
    for (std::size_t index = 1U; index < count; ++index) {
        const float value = sorted[index];
        std::size_t insertion = index;
        while (insertion > 0U && sorted[insertion - 1U] > value) {
            sorted[insertion] = sorted[insertion - 1U];
            --insertion;
        }
        sorted[insertion] = value;
    }

    if ((count & 1U) != 0U) {
        return sorted[count / 2U];
    }
    return 0.5F * (sorted[count / 2U - 1U] + sorted[count / 2U]);
}

void KinematicOutlierFilter::pushFrame(
    const Centroiding::Result& centroids) {
    Frame& destination = history_[next_history_index_];
    destination.count =
        std::min(centroids.count, Centroiding::kMaxCentroids);
    for (std::size_t index = 0U; index < destination.count; ++index) {
        destination.points[index] = centroids.points[index];
    }

    next_history_index_ = (next_history_index_ + 1U) % kHistoryFrames;
    history_count_ = std::min(history_count_ + 1U, kHistoryFrames);
}

Centroiding::Result KinematicOutlierFilter::process(
    const Centroiding::Result& centroids) {
    const std::size_t current_count =
        std::min(centroids.count, Centroiding::kMaxCentroids);
    std::array<bool, Centroiding::kMaxCentroids> rejected = {};

    if (history_count_ >= 2U && current_count > 0U) {
        const std::size_t previous_index =
            (next_history_index_ + kHistoryFrames - 1U) % kHistoryFrames;
        const std::size_t older_index =
            (next_history_index_ + kHistoryFrames - 2U) % kHistoryFrames;
        const Frame& previous = history_[previous_index];
        const Frame& older = history_[older_index];

        MatchArray current_to_previous;
        MatchArray previous_to_older;
        matchFrames(
            centroids.points,
            current_count,
            previous.points,
            previous.count,
            config_.maximum_match_distance,
            current_to_previous);
        matchFrames(
            previous.points,
            previous.count,
            older.points,
            older.count,
            config_.maximum_match_distance,
            previous_to_older);

        std::array<Motion, Centroiding::kMaxCentroids> motions = {};
        std::array<float, Centroiding::kMaxCentroids> velocity_x = {};
        std::array<float, Centroiding::kMaxCentroids> velocity_y = {};
        std::size_t motion_count = 0U;

        for (std::size_t current_index = 0U;
             current_index < current_count; ++current_index) {
            const std::int16_t previous_match =
                current_to_previous[current_index];
            if (previous_match < 0) {
                continue;
            }
            const std::size_t matched_previous =
                static_cast<std::size_t>(previous_match);
            const std::int16_t older_match =
                previous_to_older[matched_previous];
            if (older_match < 0) {
                continue;
            }
            const std::size_t matched_older =
                static_cast<std::size_t>(older_match);

            const Point& current_point = centroids.points[current_index];
            const Point& previous_point = previous.points[matched_previous];
            const Point& older_point = older.points[matched_older];
            const float current_velocity_x =
                current_point.x - previous_point.x;
            const float current_velocity_y =
                current_point.y - previous_point.y;
            const float previous_velocity_x =
                previous_point.x - older_point.x;
            const float previous_velocity_y =
                previous_point.y - older_point.y;
            const float acceleration_x =
                current_velocity_x - previous_velocity_x;
            const float acceleration_y =
                current_velocity_y - previous_velocity_y;

            Motion& motion = motions[motion_count];
            motion.velocity_x = 0.5F * (current_point.x - older_point.x);
            motion.velocity_y = 0.5F * (current_point.y - older_point.y);
            motion.acceleration = std::sqrt(
                acceleration_x * acceleration_x +
                acceleration_y * acceleration_y);
            motion.current_index = current_index;
            velocity_x[motion_count] = motion.velocity_x;
            velocity_y[motion_count] = motion.velocity_y;
            ++motion_count;
        }

        if (motion_count > config_.minimum_consensus_stars) {
            const float consensus_x = median(velocity_x, motion_count);
            const float consensus_y = median(velocity_y, motion_count);
            std::array<float, Centroiding::kMaxCentroids> residuals = {};

            for (std::size_t index = 0U; index < motion_count; ++index) {
                const float delta_x = motions[index].velocity_x - consensus_x;
                const float delta_y = motions[index].velocity_y - consensus_y;
                motions[index].residual =
                    std::sqrt(delta_x * delta_x + delta_y * delta_y);
                residuals[index] = motions[index].residual;
            }

            const float typical_residual = median(residuals, motion_count);
            const float rejection_threshold = std::max(
                config_.minimum_outlier_residual,
                config_.residual_scale_multiplier * typical_residual);
            std::size_t consensus_count = 0U;
            for (std::size_t index = 0U; index < motion_count; ++index) {
                if (motions[index].residual <= rejection_threshold) {
                    ++consensus_count;
                }
            }

            if (consensus_count >= config_.minimum_consensus_stars) {
                for (std::size_t index = 0U; index < motion_count; ++index) {
                    if (motions[index].residual > rejection_threshold &&
                        motions[index].acceleration <=
                            config_.maximum_linear_acceleration) {
                        rejected[motions[index].current_index] = true;
                    }
                }
            }
        }
    }

    Centroiding::Result filtered;
    for (std::size_t index = 0U; index < current_count; ++index) {
        if (!rejected[index]) {
            filtered.points[filtered.count] = centroids.points[index];
            ++filtered.count;
        }
    }

    // Keep raw detections privately so a moving object can be tracked for
    // three frames even after it is removed from the downstream result.
    pushFrame(centroids);
    return filtered;
}

}  // namespace star_tracker

static_assert(
    sizeof(star_tracker::KinematicOutlierFilter) < 16U * 1024U,
    "Kinematic filter working memory unexpectedly exceeds 16 KB");
