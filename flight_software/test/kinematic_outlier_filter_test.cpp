#include "star_tracker/kinematic_outlier_filter.hpp"

#include <cassert>
#include <cstddef>

namespace {

star_tracker::Centroiding::Result makeFrame(std::size_t frame_index) {
    star_tracker::Centroiding::Result frame;

    for (std::size_t index = 0U; index < 10U; ++index) {
        star_tracker::Point& point = frame.points[frame.count];
        point.x = 50.0F + static_cast<float>(index * 60U + frame_index);
        point.y =
            100.0F + static_cast<float>((index % 2U) * 100U) +
            0.5F * static_cast<float>(frame_index);
        point.intensity = 2000U + static_cast<unsigned int>(index);
        ++frame.count;
    }

    star_tracker::Point& moving_object = frame.points[frame.count];
    moving_object.x = 900.0F;
    moving_object.y = 100.0F + 5.0F * static_cast<float>(frame_index);
    moving_object.intensity = 2500U;
    ++frame.count;
    return frame;
}

}  // namespace

int main() {
    star_tracker::KinematicOutlierFilter filter;

    const star_tracker::Centroiding::Result first = filter.process(makeFrame(0U));
    const star_tracker::Centroiding::Result second = filter.process(makeFrame(1U));
    const star_tracker::Centroiding::Result third = filter.process(makeFrame(2U));

    assert(first.count == 11U);
    assert(second.count == 11U);
    assert(third.count == 10U);
    for (std::size_t index = 0U; index < third.count; ++index) {
        assert(third.points[index].x < 900.0F);
    }

    filter.reset();
    const star_tracker::Centroiding::Result after_reset =
        filter.process(makeFrame(3U));
    assert(after_reset.count == 11U);
    return 0;
}
