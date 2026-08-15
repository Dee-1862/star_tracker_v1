#include "star_tracker/centroiding.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace {

star_tracker::Centroiding::Image image;

void setPixel(std::size_t x, std::size_t y, std::uint8_t intensity) {
    image[y * star_tracker::Centroiding::kImageWidth + x] = intensity;
}

}  // namespace

int main() {
    image.fill(0U);

    const std::uint8_t star[3][3] = {
        {201U, 205U, 201U},
        {205U, 255U, 220U},
        {201U, 215U, 205U},
    };
    std::uint32_t expected_intensity = 0U;
    double expected_x_sum = 0.0;
    double expected_y_sum = 0.0;
    for (std::size_t row = 0U; row < 3U; ++row) {
        for (std::size_t column = 0U; column < 3U; ++column) {
            const std::size_t x = 100U + column;
            const std::size_t y = 200U + row;
            setPixel(x, y, star[row][column]);
            expected_intensity += star[row][column];
            expected_x_sum += static_cast<double>(x * star[row][column]);
            expected_y_sum += static_cast<double>(y * star[row][column]);
        }
    }

    // Too small: isolated hot pixel.
    setPixel(500U, 500U, 255U);

    // Highly elongated: cosmic-ray-like streak.
    for (std::size_t x = 700U; x < 710U; ++x) {
        setPixel(x, 400U, 255U);
    }

    // Too large: 11 x 11 saturated cluster.
    for (std::size_t y = 800U; y < 811U; ++y) {
        for (std::size_t x = 800U; x < 811U; ++x) {
            setPixel(x, y, 255U);
        }
    }

    star_tracker::Centroiding centroiding;
    const star_tracker::Centroiding::Result result = centroiding.process(image);

    assert(result.count == 1U);
    assert(result.points[0].intensity == expected_intensity);

    const double expected_x = expected_x_sum / expected_intensity;
    const double expected_y = expected_y_sum / expected_intensity;
    assert(std::fabs(result.points[0].x - expected_x) < 1.0e-4);
    assert(std::fabs(result.points[0].y - expected_y) < 1.0e-4);
    return 0;
}
