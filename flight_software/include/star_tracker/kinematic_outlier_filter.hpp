#ifndef STAR_TRACKER_KINEMATIC_OUTLIER_FILTER_HPP
#define STAR_TRACKER_KINEMATIC_OUTLIER_FILTER_HPP

#include "star_tracker/centroiding.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace star_tracker {

class KinematicOutlierFilter {
public:
    static constexpr std::size_t kHistoryFrames = 3U;

    struct Config {
        float maximum_match_distance;
        float minimum_outlier_residual;
        float residual_scale_multiplier;
        float maximum_linear_acceleration;
        std::size_t minimum_consensus_stars;

        Config()
            : maximum_match_distance(20.0F),
              minimum_outlier_residual(2.0F),
              residual_scale_multiplier(4.0F),
              maximum_linear_acceleration(2.0F),
              minimum_consensus_stars(10U) {}
    };

    KinematicOutlierFilter();
    explicit KinematicOutlierFilter(const Config& config);

    Centroiding::Result process(const Centroiding::Result& centroids);
    void reset();

private:
    typedef std::array<std::int16_t, Centroiding::kMaxCentroids> MatchArray;

    struct Frame {
        std::array<Point, Centroiding::kMaxCentroids> points;
        std::size_t count;

        Frame() : points(), count(0U) {}
    };

    struct Motion {
        float velocity_x;
        float velocity_y;
        float acceleration;
        float residual;
        std::size_t current_index;
    };

    static void matchFrames(
        const std::array<Point, Centroiding::kMaxCentroids>& current,
        std::size_t current_count,
        const std::array<Point, Centroiding::kMaxCentroids>& previous,
        std::size_t previous_count,
        float maximum_distance,
        MatchArray& matches);
    static float median(
        const std::array<float, Centroiding::kMaxCentroids>& values,
        std::size_t count);
    void pushFrame(const Centroiding::Result& centroids);

    Config config_;
    std::array<Frame, kHistoryFrames> history_;
    std::size_t next_history_index_;
    std::size_t history_count_;
};

}  // namespace star_tracker

#endif
