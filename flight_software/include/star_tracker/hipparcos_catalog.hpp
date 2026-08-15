#ifndef STAR_TRACKER_HIPPARCOS_CATALOG_HPP
#define STAR_TRACKER_HIPPARCOS_CATALOG_HPP

#include "star_tracker/lis_grid_matcher.hpp"

#include <cstddef>

namespace star_tracker {

extern const LisGridMatcher::Catalog kHipparcosCatalog;
static constexpr std::size_t kHipparcosCatalogCount = 5041U;

}  // namespace star_tracker

#endif
