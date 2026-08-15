#include "star_tracker/hipparcos_catalog.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

star_tracker::LisGridMatcher::Catalog catalog;

/// Build a catalogue entry that would be observed at the given pixel. Must use
/// the same pixel-to-bearing convention as the matcher (+x boresight,
/// +y image-left, +z image-up), or the synthetic catalogue is a mirror image of
/// the observations and the chirality test correctly rejects all of it.
star_tracker::LisGridMatcher::CatalogStar makeCatalogStar(
    std::uint32_t id, float pixel_x, float pixel_y, float magnitude) {
    const float focal_length = 2904.0F;
    const float y = (512.0F - pixel_x) / focal_length;
    const float z = (512.0F - pixel_y) / focal_length;
    const float inverse_norm = 1.0F / std::sqrt(1.0F + (y * y) + (z * z));
    return {
        id,
        inverse_norm,
        y * inverse_norm,
        z * inverse_norm,
        magnitude,
    };
}

}  // namespace

int main() {
    const std::array<float, 12U> pixel_x = {
        140.0F, 230.0F, 355.0F, 480.0F, 610.0F, 770.0F,
        880.0F, 180.0F, 430.0F, 690.0F, 820.0F, 300.0F,
    };
    const std::array<float, 12U> pixel_y = {
        170.0F, 330.0F, 140.0F, 500.0F, 250.0F, 420.0F,
        730.0F, 800.0F, 690.0F, 110.0F, 560.0F, 600.0F,
    };

    star_tracker::Centroiding::Result observed;
    for (std::size_t index = 0U; index < pixel_x.size(); ++index) {
        catalog[index] = makeCatalogStar(
            static_cast<std::uint32_t>(1000U + index),
            pixel_x[index],
            pixel_y[index],
            1.0F + 0.1F * static_cast<float>(index));
        observed.points[index] = {
            pixel_x[index],
            pixel_y[index],
            static_cast<std::uint32_t>(3000U - index * 10U),
        };
        ++observed.count;
    }

    star_tracker::LisGridMatcher matcher;
    assert(matcher.buildIndex(catalog, pixel_x.size()));
    assert(matcher.indexedStarCount() == pixel_x.size());
    assert(matcher.indexedPairCount() > 0U);

    const star_tracker::LisGridMatcher::MatchResult result =
        matcher.match(observed);
    assert(result.matched_count == pixel_x.size());
    for (std::size_t index = 0U; index < pixel_x.size(); ++index) {
        assert(result.star_ids[index] == 1000U + index);
        assert(result.votes[index] >= 3U);
    }

    assert(matcher.buildIndex(
        star_tracker::kHipparcosCatalog,
        star_tracker::kHipparcosCatalogCount));
    assert(matcher.indexedStarCount() ==
           star_tracker::LisGridMatcher::kMaxIndexedStars);
    assert(matcher.indexedPairCount() > 10000U);
    assert(matcher.indexedPairCount() <
           star_tracker::LisGridMatcher::kMaxPairEntries);
    return 0;
}
