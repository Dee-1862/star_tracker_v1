#include "star_tracker/centroiding.hpp"

#include <algorithm>
#include <cmath>

namespace star_tracker {

Centroiding::Centroiding() : config_(), previous_runs_(), current_runs_(),
                             components_(), next_components_(), compact_labels_() {}

Centroiding::Centroiding(const Config& config)
    : config_(config), previous_runs_(), current_runs_(), components_(),
      next_components_(), compact_labels_() {}

std::uint16_t Centroiding::findRoot(std::uint16_t label) {
    std::uint16_t root = label;
    while (components_[root].parent != root) {
        root = components_[root].parent;
    }

    while (components_[label].parent != label) {
        const std::uint16_t parent = components_[label].parent;
        components_[label].parent = root;
        label = parent;
    }
    return root;
}

std::uint16_t Centroiding::mergeComponents(
    std::uint16_t first, std::uint16_t second) {
    std::uint16_t first_root = findRoot(first);
    std::uint16_t second_root = findRoot(second);
    if (first_root == second_root) {
        return first_root;
    }

    if (components_[second_root].pixel_count > components_[first_root].pixel_count) {
        std::swap(first_root, second_root);
    }

    Component& destination = components_[first_root];
    const Component& source = components_[second_root];
    components_[second_root].parent = first_root;

    destination.last_row = std::max(destination.last_row, source.last_row);
    destination.min_x = std::min(destination.min_x, source.min_x);
    destination.max_x = std::max(destination.max_x, source.max_x);
    destination.min_y = std::min(destination.min_y, source.min_y);
    destination.max_y = std::max(destination.max_y, source.max_y);
    destination.pixel_count += source.pixel_count;
    destination.intensity_sum += source.intensity_sum;
    destination.weighted_x += source.weighted_x;
    destination.weighted_y += source.weighted_y;
    destination.weighted_x2 += source.weighted_x2;
    destination.weighted_y2 += source.weighted_y2;
    destination.weighted_xy += source.weighted_xy;
    return first_root;
}

void Centroiding::addRun(
    std::uint16_t label, const RunStatistics& statistics) {
    Component& component = components_[findRoot(label)];
    component.last_row = statistics.y;
    component.min_x = std::min(component.min_x, statistics.min_x);
    component.max_x = std::max(component.max_x, statistics.max_x);
    component.min_y = std::min(component.min_y, statistics.y);
    component.max_y = std::max(component.max_y, statistics.y);
    component.pixel_count += statistics.pixel_count;
    component.intensity_sum += statistics.intensity_sum;
    component.weighted_x += statistics.weighted_x;
    component.weighted_y += statistics.weighted_y;
    component.weighted_x2 += statistics.weighted_x2;
    component.weighted_y2 += statistics.weighted_y2;
    component.weighted_xy += statistics.weighted_xy;
}

void Centroiding::addStrongestPoint(const Point& point, Result& result) {
    if (result.count < kMaxCentroids) {
        result.points[result.count] = point;
        ++result.count;
        return;
    }

    std::size_t weakest_index = 0U;
    for (std::size_t index = 1U; index < kMaxCentroids; ++index) {
        if (result.points[index].intensity <
            result.points[weakest_index].intensity) {
            weakest_index = index;
        }
    }
    if (point.intensity > result.points[weakest_index].intensity) {
        result.points[weakest_index] = point;
    }
}

void Centroiding::emitComponent(
    const Component& component, Result& result) const {
    if (component.pixel_count < config_.minimum_cluster_pixels ||
        component.pixel_count > config_.maximum_cluster_pixels ||
        component.intensity_sum == 0U) {
        return;
    }

    const double inverse_intensity =
        1.0 / static_cast<double>(component.intensity_sum);
    const double center_x =
        static_cast<double>(component.weighted_x) * inverse_intensity;
    const double center_y =
        static_cast<double>(component.weighted_y) * inverse_intensity;
    const double variance_x =
        std::max(0.0, static_cast<double>(component.weighted_x2) *
                          inverse_intensity - center_x * center_x);
    const double variance_y =
        std::max(0.0, static_cast<double>(component.weighted_y2) *
                          inverse_intensity - center_y * center_y);
    const double covariance =
        static_cast<double>(component.weighted_xy) * inverse_intensity -
        center_x * center_y;

    const double trace = variance_x + variance_y;
    const double discriminant = std::sqrt(
        std::max(0.0, (variance_x - variance_y) * (variance_x - variance_y) +
                          4.0 * covariance * covariance));
    const double major_variance = 0.5 * (trace + discriminant);
    const double minor_variance = 0.5 * (trace - discriminant);
    if (minor_variance <= 1.0e-9) {
        return;
    }

    const double elongation = std::sqrt(major_variance / minor_variance);
    if (elongation > static_cast<double>(config_.maximum_elongation)) {
        return;
    }

    const Point point = {
        static_cast<float>(center_x),
        static_cast<float>(center_y),
        component.intensity_sum,
    };
    addStrongestPoint(point, result);
}

Centroiding::Result Centroiding::process(const Image& image) {
    Result result;
    std::size_t previous_run_count = 0U;
    std::size_t active_component_count = 0U;

    for (std::size_t y = 0U; y < kImageHeight; ++y) {
        std::size_t current_run_count = 0U;
        std::size_t component_count = active_component_count;
        std::size_t previous_search_start = 0U;
        std::size_t x = 0U;

        while (x < kImageWidth) {
            const std::size_t pixel_index = y * kImageWidth + x;
            if (image[pixel_index] <= config_.brightness_threshold) {
                ++x;
                continue;
            }

            RunStatistics statistics = {};
            statistics.min_x = static_cast<std::uint16_t>(x);
            statistics.y = static_cast<std::uint16_t>(y);

            while (x < kImageWidth &&
                   image[y * kImageWidth + x] > config_.brightness_threshold) {
                const std::uint64_t intensity = image[y * kImageWidth + x];
                statistics.max_x = static_cast<std::uint16_t>(x);
                ++statistics.pixel_count;
                statistics.intensity_sum += static_cast<std::uint32_t>(intensity);
                statistics.weighted_x += intensity * x;
                statistics.weighted_y += intensity * y;
                statistics.weighted_x2 += intensity * x * x;
                statistics.weighted_y2 += intensity * y * y;
                statistics.weighted_xy += intensity * x * y;
                ++x;
            }

            Run& run = current_runs_[current_run_count];
            run.start_x = statistics.min_x;
            run.end_x = statistics.max_x;
            run.label = kInvalidLabel;

            while (previous_search_start < previous_run_count &&
                   static_cast<std::size_t>(
                       previous_runs_[previous_search_start].end_x) +
                           1U <
                       run.start_x) {
                ++previous_search_start;
            }

            for (std::size_t previous_index = previous_search_start;
                 previous_index < previous_run_count; ++previous_index) {
                const Run& previous = previous_runs_[previous_index];
                if (static_cast<std::size_t>(run.end_x) + 1U <
                    previous.start_x) {
                    break;
                }

                if (run.label == kInvalidLabel) {
                    run.label = findRoot(previous.label);
                } else {
                    run.label = mergeComponents(run.label, previous.label);
                }
            }

            if (run.label == kInvalidLabel) {
                run.label = static_cast<std::uint16_t>(component_count);
                Component& component = components_[component_count];
                component = {};
                component.parent = run.label;
                component.last_row = static_cast<std::uint16_t>(y);
                component.min_x = statistics.min_x;
                component.max_x = statistics.max_x;
                component.min_y = statistics.y;
                component.max_y = statistics.y;
                ++component_count;
            }

            addRun(run.label, statistics);
            ++current_run_count;
        }

        for (std::size_t index = 0U; index < component_count; ++index) {
            if (components_[index].parent == index &&
                components_[index].last_row != y) {
                emitComponent(components_[index], result);
            }
            compact_labels_[index] = kInvalidLabel;
        }

        std::size_t next_component_count = 0U;
        for (std::size_t index = 0U; index < current_run_count; ++index) {
            const std::uint16_t root = findRoot(current_runs_[index].label);
            if (compact_labels_[root] == kInvalidLabel) {
                compact_labels_[root] =
                    static_cast<std::uint16_t>(next_component_count);
                next_components_[next_component_count] = components_[root];
                next_components_[next_component_count].parent =
                    static_cast<std::uint16_t>(next_component_count);
                ++next_component_count;
            }
            current_runs_[index].label = compact_labels_[root];
        }

        for (std::size_t index = 0U; index < next_component_count; ++index) {
            components_[index] = next_components_[index];
        }
        for (std::size_t index = 0U; index < current_run_count; ++index) {
            previous_runs_[index] = current_runs_[index];
        }
        active_component_count = next_component_count;
        previous_run_count = current_run_count;
    }

    for (std::size_t index = 0U; index < active_component_count; ++index) {
        emitComponent(components_[index], result);
    }
    return result;
}

}  // namespace star_tracker

static_assert(
    sizeof(star_tracker::Centroiding) < 500U * 1024U,
    "Centroiding working memory exceeds the 500 KB flight-software budget");
