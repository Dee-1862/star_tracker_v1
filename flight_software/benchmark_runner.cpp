#include "star_tracker/centroiding.hpp"
#include "star_tracker/hipparcos_catalog.hpp"
#include "star_tracker/kinematic_outlier_filter.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

namespace {

star_tracker::Centroiding::Image image;
star_tracker::LisGridMatcher matcher;

std::uint8_t parseThreshold(const char* text) {
    const long value = std::strtol(text, 0, 10);
    if (value < 0L || value > 255L) {
        std::cerr << "Threshold must be in [0, 255]\n";
        std::exit(2);
    }
    return static_cast<std::uint8_t>(value);
}

template <typename Clock>
double elapsedMicroseconds(
    const typename Clock::time_point& begin,
    const typename Clock::time_point& end) {
    return std::chrono::duration<double, std::micro>(end - begin).count();
}

std::uint64_t peakResidentBytes() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS counters = {};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(
            GetCurrentProcess(), &counters, sizeof(counters)) == 0) {
        return 0U;
    }
    return static_cast<std::uint64_t>(counters.PeakWorkingSetSize);
#else
    struct rusage usage = {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0U;
    }
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024U;
#endif
#endif
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: star_tracker_benchmark_runner IMAGE.raw THRESHOLD\n";
        return 2;
    }

    std::ifstream input(argv[1], std::ios::binary);
    if (!input.read(
            reinterpret_cast<char*>(image.data()),
            static_cast<std::streamsize>(image.size()))) {
        std::cerr << "Input must contain exactly one 1024x1024 uint8 image\n";
        return 2;
    }
    if (input.peek() != std::ifstream::traits_type::eof()) {
        std::cerr << "Input contains more than 1024x1024 bytes\n";
        return 2;
    }

    star_tracker::Centroiding::Config centroid_config;
    centroid_config.brightness_threshold = parseThreshold(argv[2]);
    star_tracker::Centroiding centroiding(centroid_config);
    star_tracker::KinematicOutlierFilter outlier_filter;

    typedef std::chrono::steady_clock Clock;
    const Clock::time_point index_begin = Clock::now();
    if (!matcher.buildIndex(
            star_tracker::kHipparcosCatalog,
            star_tracker::kHipparcosCatalogCount)) {
        std::cerr << "Failed to build the Hipparcos pair index\n";
        return 3;
    }
    const Clock::time_point index_end = Clock::now();

    const Clock::time_point centroid_begin = Clock::now();
    const star_tracker::Centroiding::Result centroids =
        centroiding.process(image);
    const Clock::time_point centroid_end = Clock::now();

    const Clock::time_point filter_begin = Clock::now();
    const star_tracker::Centroiding::Result cleaned =
        outlier_filter.process(centroids);
    const Clock::time_point filter_end = Clock::now();

    const Clock::time_point match_begin = Clock::now();
    const star_tracker::LisGridMatcher::MatchResult matches =
        matcher.match(cleaned);
    const Clock::time_point match_end = Clock::now();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "SUMMARY,"
              << centroids.count << ','
              << cleaned.count << ','
              << matches.matched_count << ','
              << matcher.indexedStarCount() << ','
              << matcher.indexedPairCount() << ','
              << elapsedMicroseconds<Clock>(index_begin, index_end) << ','
              << elapsedMicroseconds<Clock>(centroid_begin, centroid_end) << ','
              << elapsedMicroseconds<Clock>(filter_begin, filter_end) << ','
              << elapsedMicroseconds<Clock>(match_begin, match_end) << ','
              << sizeof(star_tracker::Centroiding) << ','
              << sizeof(star_tracker::KinematicOutlierFilter) << ','
              << sizeof(star_tracker::LisGridMatcher) << ','
              << sizeof(star_tracker::LisGridMatcher::Catalog) << ','
              << peakResidentBytes() << '\n';

    for (std::size_t index = 0U; index < cleaned.count; ++index) {
        const star_tracker::Point& point = cleaned.points[index];
        std::cout << "POINT,"
                  << index << ','
                  << point.x << ','
                  << point.y << ','
                  << point.intensity << ','
                  << matches.star_ids[index] << ','
                  << static_cast<unsigned int>(matches.votes[index]) << '\n';
    }
    return 0;
}
