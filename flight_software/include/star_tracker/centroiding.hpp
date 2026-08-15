#ifndef STAR_TRACKER_CENTROIDING_HPP
#define STAR_TRACKER_CENTROIDING_HPP

#include <array>
#include <cstddef>
#include <cstdint>

namespace star_tracker {

struct Point {
    float x;
    float y;
    std::uint32_t intensity;
};

class Centroiding {
public:
    static constexpr std::size_t kImageWidth = 1024U;
    static constexpr std::size_t kImageHeight = 1024U;
    static constexpr std::size_t kPixelCount = kImageWidth * kImageHeight;
    static constexpr std::size_t kMaxCentroids = 50U;

    typedef std::array<std::uint8_t, kPixelCount> Image;

    struct Result {
        std::array<Point, kMaxCentroids> points;
        std::size_t count;

        Result() : points(), count(0U) {}
    };

    struct Config {
        std::uint8_t brightness_threshold;
        std::uint16_t minimum_cluster_pixels;
        std::uint16_t maximum_cluster_pixels;
        float maximum_elongation;

        Config()
            : brightness_threshold(200U),
              minimum_cluster_pixels(3U),
              maximum_cluster_pixels(100U),
              maximum_elongation(4.0F) {}
    };

    Centroiding();
    explicit Centroiding(const Config& config);

    Result process(const Image& image);

private:
    static constexpr std::size_t kMaxRunsPerRow = (kImageWidth + 1U) / 2U;
    static constexpr std::size_t kMaxComponents = kMaxRunsPerRow * 2U;
    static constexpr std::uint16_t kInvalidLabel = 0xFFFFU;

    struct Run {
        std::uint16_t start_x;
        std::uint16_t end_x;
        std::uint16_t label;
    };

    struct Component {
        std::uint16_t parent;
        std::uint16_t last_row;
        std::uint16_t min_x;
        std::uint16_t max_x;
        std::uint16_t min_y;
        std::uint16_t max_y;
        std::uint32_t pixel_count;
        std::uint32_t intensity_sum;
        std::uint64_t weighted_x;
        std::uint64_t weighted_y;
        std::uint64_t weighted_x2;
        std::uint64_t weighted_y2;
        std::uint64_t weighted_xy;
    };

    struct RunStatistics {
        std::uint16_t min_x;
        std::uint16_t max_x;
        std::uint16_t y;
        std::uint32_t pixel_count;
        std::uint32_t intensity_sum;
        std::uint64_t weighted_x;
        std::uint64_t weighted_y;
        std::uint64_t weighted_x2;
        std::uint64_t weighted_y2;
        std::uint64_t weighted_xy;
    };

    std::uint16_t findRoot(std::uint16_t label);
    std::uint16_t mergeComponents(std::uint16_t first, std::uint16_t second);
    void addRun(std::uint16_t label, const RunStatistics& statistics);
    void emitComponent(const Component& component, Result& result) const;
    static void addStrongestPoint(const Point& point, Result& result);

    Config config_;
    std::array<Run, kMaxRunsPerRow> previous_runs_;
    std::array<Run, kMaxRunsPerRow> current_runs_;
    std::array<Component, kMaxComponents> components_;
    std::array<Component, kMaxRunsPerRow> next_components_;
    std::array<std::uint16_t, kMaxComponents> compact_labels_;
};

}  // namespace star_tracker

#endif
