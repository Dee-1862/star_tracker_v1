#include "star_tracker/lis_grid_matcher.hpp"

#include <algorithm>
#include <cmath>

namespace {

const float kPi = 3.14159265358979323846F;

float clampUnit(float value) {
    return std::max(-1.0F, std::min(1.0F, value));
}

}  // namespace

namespace star_tracker {

LisGridMatcher::LisGridMatcher()
    : config_(),
      catalog_(0),
      catalog_count_(0U),
      indexed_star_count_(0U),
      pair_entry_count_(0U),
      maximum_indexed_angle_radians_(0.0F),
      angular_tolerance_radians_(0.0F),
      index_ready_(false),
      indexed_catalog_positions_(),
      bin_heads_(),
      pair_entries_(),
      votes_() {}

LisGridMatcher::LisGridMatcher(const Config& config)
    : config_(config),
      catalog_(0),
      catalog_count_(0U),
      indexed_star_count_(0U),
      pair_entry_count_(0U),
      maximum_indexed_angle_radians_(0.0F),
      angular_tolerance_radians_(0.0F),
      index_ready_(false),
      indexed_catalog_positions_(),
      bin_heads_(),
      pair_entries_(),
      votes_() {}

float LisGridMatcher::catalogPairAngle(
    const CatalogStar& first, const CatalogStar& second) {
    return std::acos(clampUnit(catalogPairCosine(first, second)));
}

float LisGridMatcher::catalogPairCosine(
    const CatalogStar& first, const CatalogStar& second) {
    return (first.x * second.x) + (first.y * second.y) + (first.z * second.z);
}

std::size_t LisGridMatcher::angleToBin(float angle_radians) const {
    if (angle_radians <= 0.0F) {
        return 0U;
    }
    const float scaled =
        angle_radians * static_cast<float>(kAngleBinCount) /
        maximum_indexed_angle_radians_;
    return std::min(
        static_cast<std::size_t>(scaled), kAngleBinCount - 1U);
}

void LisGridMatcher::selectBrightestStars() {
    // Keeping the globally brightest N is the wrong rule: it concentrates the
    // index in already-crowded regions and leaves sparse fields with almost no
    // indexed stars, so the correct identification is simply absent there.
    // Instead thin greedily by angular separation, brightest first, which
    // gives roughly uniform sky density -- the same approach tetra3 uses to
    // hit a target star count per field of view.
    indexed_star_count_ = 0U;

    static std::array<std::uint16_t, kMaxCatalogStars> order;
    std::size_t candidate_count = 0U;
    for (std::size_t catalog_index = 0U;
         catalog_index < catalog_count_; ++catalog_index) {
        const CatalogStar& star = (*catalog_)[catalog_index];
        const float norm_squared =
            (star.x * star.x) + (star.y * star.y) + (star.z * star.z);
        if (star.star_id == 0U || !std::isfinite(star.visual_magnitude) ||
            norm_squared < 0.99F || norm_squared > 1.01F) {
            continue;
        }
        order[candidate_count] = static_cast<std::uint16_t>(catalog_index);
        ++candidate_count;
    }

    const Catalog& catalog = *catalog_;
    std::sort(
        order.begin(), order.begin() + static_cast<std::ptrdiff_t>(candidate_count),
        [&catalog](std::uint16_t left, std::uint16_t right) {
            return catalog[left].visual_magnitude <
                   catalog[right].visual_magnitude;
        });

    // Radius at which kMaxIndexedStars discs tile the sphere: 2*pi*(1-cos r)
    // = 4*pi/N. Backed off slightly so we reach the budget rather than
    // stopping short of it.
    const float target =
        1.0F - (2.0F / static_cast<float>(kMaxIndexedStars));
    const float minimum_cosine = std::cos(0.75F * std::acos(clampUnit(target)));

    for (std::size_t position = 0U;
         (position < candidate_count) &&
         (indexed_star_count_ < kMaxIndexedStars);
         ++position) {
        const CatalogStar& star = catalog[order[position]];
        bool crowded = false;
        for (std::size_t kept = 0U; kept < indexed_star_count_; ++kept) {
            const CatalogStar& other =
                catalog[indexed_catalog_positions_[kept]];
            const float cosine = (star.x * other.x) + (star.y * other.y) +
                                 (star.z * other.z);
            if (cosine > minimum_cosine) {
                crowded = true;
                break;
            }
        }
        if (crowded) {
            continue;
        }
        indexed_catalog_positions_[indexed_star_count_] = order[position];
        ++indexed_star_count_;
    }
}

bool LisGridMatcher::buildIndex(
    const Catalog& catalog, std::size_t catalog_count) {
    index_ready_ = false;
    catalog_ = &catalog;
    catalog_count_ = std::min(catalog_count, kMaxCatalogStars);
    pair_entry_count_ = 0U;
    maximum_indexed_angle_radians_ =
        config_.maximum_indexed_angle_degrees * kPi / 180.0F;
    angular_tolerance_radians_ =
        config_.angular_tolerance_degrees * kPi / 180.0F;

    if (catalog_count_ == 0U || config_.focal_length_pixels <= 0.0F ||
        maximum_indexed_angle_radians_ <= 0.0F ||
        maximum_indexed_angle_radians_ >= kPi ||
        angular_tolerance_radians_ <= 0.0F) {
        return false;
    }

    bin_heads_.fill(kInvalidEntry);
    selectBrightestStars();
    if (indexed_star_count_ < 3U) {
        return false;
    }

    for (std::size_t first = 0U; first < indexed_star_count_; ++first) {
        const CatalogStar& first_star =
            catalog[indexed_catalog_positions_[first]];
        for (std::size_t second = first + 1U;
             second < indexed_star_count_; ++second) {
            const CatalogStar& second_star =
                catalog[indexed_catalog_positions_[second]];
            const float angle = catalogPairAngle(first_star, second_star);
            if (angle > maximum_indexed_angle_radians_) {
                continue;
            }
            if (pair_entry_count_ >= kMaxPairEntries) {
                pair_entry_count_ = 0U;
                return false;
            }

            const std::size_t bin = angleToBin(angle);
            PairEntry& entry = pair_entries_[pair_entry_count_];
            entry.first_star = static_cast<std::uint16_t>(first);
            entry.second_star = static_cast<std::uint16_t>(second);
            entry.next = bin_heads_[bin];
            bin_heads_[bin] =
                static_cast<std::uint16_t>(pair_entry_count_);
            ++pair_entry_count_;
        }
    }

    index_ready_ = pair_entry_count_ > 0U;
    return index_ready_;
}

void LisGridMatcher::selectObservedStars(
    const Centroiding::Result& centroids,
    std::array<std::uint16_t, Centroiding::kMaxCentroids>& selected,
    std::size_t& selected_count) const {
    selected_count =
        std::min(centroids.count, Centroiding::kMaxCentroids);
    for (std::size_t index = 0U; index < selected_count; ++index) {
        selected[index] = static_cast<std::uint16_t>(index);
    }

    for (std::size_t index = 1U; index < selected_count; ++index) {
        const std::uint16_t value = selected[index];
        std::size_t insertion = index;
        while (insertion > 0U &&
               centroids.points[selected[insertion - 1U]].intensity <
                   centroids.points[value].intensity) {
            selected[insertion] = selected[insertion - 1U];
            --insertion;
        }
        selected[insertion] = value;
    }

    selected_count = std::min(
        selected_count,
        std::min(config_.maximum_observed_stars,
                 Centroiding::kMaxCentroids));
}

void LisGridMatcher::addVote(
    std::size_t selection_index, std::size_t indexed_star) {
    if (selection_index >= kMaxCliqueObserved) {
        return;
    }
    std::uint8_t& vote = votes_[selection_index][indexed_star];
    if (vote < 0xFFU) {
        ++vote;
    }
}

void LisGridMatcher::bearingFor(
    float pixel_x, float pixel_y, Bearing& bearing) const {
    // LOST/tetra3/AttitudeSolver convention: +x boresight, +y image-left,
    // +z image-up. Right-handed, and consistent with the catalogue's ICRS
    // handedness -- which the previous (right, up, boresight) ordering was
    // NOT. Pairwise dot products are identical under either, so this changes
    // no existing angle, but it is what makes the chirality test below
    // meaningful rather than systematically inverted.
    const float y = (config_.principal_x - pixel_x) / config_.focal_length_pixels;
    const float z = (config_.principal_y - pixel_y) / config_.focal_length_pixels;
    const float inverse_norm = 1.0F / std::sqrt(1.0F + (y * y) + (z * z));
    bearing[0U] = inverse_norm;
    bearing[1U] = y * inverse_norm;
    bearing[2U] = z * inverse_norm;
}

bool LisGridMatcher::adjacent(std::size_t first, std::size_t second) const {
    return (adjacency_[first][second / 32U] &
            (1UL << (second % 32U))) != 0UL;
}

std::size_t LisGridMatcher::buildNodes(
    const std::array<std::uint16_t, Centroiding::kMaxCentroids>& selected,
    std::size_t selected_count) {
    std::size_t node_count = 0U;
    const std::size_t observed_limit =
        (selected_count < kMaxCliqueObserved) ? selected_count
                                              : kMaxCliqueObserved;

    for (std::size_t selection = 0U; selection < observed_limit; ++selection) {
        const std::size_t centroid_index = selected[selection];

        // Keep the best-voted candidates for this centroid. Voting is only a
        // shortlist now; correctness comes from the clique, not the vote count.
        std::array<std::uint16_t, kCandidatesPerCentroid> best_stars = {};
        std::array<std::uint8_t, kCandidatesPerCentroid> best_votes = {};
        std::size_t kept = 0U;

        for (std::size_t star = 0U; star < indexed_star_count_; ++star) {
            const std::uint8_t vote = votes_[selection][star];
            if (vote < config_.minimum_votes) {
                continue;
            }
            if (kept < kCandidatesPerCentroid) {
                best_stars[kept] = static_cast<std::uint16_t>(star);
                best_votes[kept] = vote;
                ++kept;
            } else {
                std::size_t weakest = 0U;
                for (std::size_t slot = 1U; slot < kept; ++slot) {
                    if (best_votes[slot] < best_votes[weakest]) {
                        weakest = slot;
                    }
                }
                if (vote > best_votes[weakest]) {
                    best_stars[weakest] = static_cast<std::uint16_t>(star);
                    best_votes[weakest] = vote;
                }
            }
        }

        for (std::size_t slot = 0U; slot < kept; ++slot) {
            if (node_count >= kMaxNodes) {
                break;
            }
            nodes_[node_count].centroid =
                static_cast<std::uint16_t>(centroid_index);
            nodes_[node_count].indexed_star = best_stars[slot];
            ++node_count;
        }
    }
    return node_count;
}

void LisGridMatcher::buildAdjacency(
    std::size_t node_count,
    const std::array<Bearing, Centroiding::kMaxCentroids>& bearings) {
    for (std::size_t node = 0U; node < node_count; ++node) {
        adjacency_[node].fill(0UL);
    }

    for (std::size_t first = 0U; first < node_count; ++first) {
        for (std::size_t second = first + 1U; second < node_count; ++second) {
            // A centroid cannot be two stars, and a star cannot be two
            // centroids: these are the one-to-one constraints the previous
            // greedy assignment enforced only after the fact.
            if (nodes_[first].centroid == nodes_[second].centroid) {
                continue;
            }
            if (nodes_[first].indexed_star == nodes_[second].indexed_star) {
                continue;
            }

            const Bearing& a = bearings[nodes_[first].centroid];
            const Bearing& b = bearings[nodes_[second].centroid];
            const float observed_cosine = clampUnit(
                (a[0U] * b[0U]) + (a[1U] * b[1U]) + (a[2U] * b[2U]));

            const CatalogStar& s =
                (*catalog_)[indexed_catalog_positions_[nodes_[first].indexed_star]];
            const CatalogStar& t =
                (*catalog_)[indexed_catalog_positions_[nodes_[second].indexed_star]];
            const float expected_cosine = clampUnit(catalogPairCosine(s, t));

            // Angular error without any acos. Since d(cos w)/dw = -sin w,
            // a difference in cosine maps to a difference in angle by
            // dividing by sin. The residuals of interest are ~10 arcsec
            // (5e-5 rad), where this first-order relation is exact to far
            // better than the measurement itself.
            const float expected_sine = std::sqrt(std::max(
                1.0e-12F, 1.0F - (expected_cosine * expected_cosine)));
            const float error =
                std::fabs(expected_cosine - observed_cosine) / expected_sine;
            if (error > angular_tolerance_radians_) {
                continue;
            }

            const float error_arcsec = error * 206264.806247F;
            edge_error_[first][second] = error_arcsec;
            edge_error_[second][first] = error_arcsec;
            adjacency_[first][second / 32U] |= (1UL << (second % 32U));
            adjacency_[second][first / 32U] |= (1UL << (first % 32U));
        }
    }
}

namespace {

float tripleProduct(
    const std::array<float, 3U>& a,
    const std::array<float, 3U>& b,
    const std::array<float, 3U>& c) {
    const float cross_x = (a[1U] * b[2U]) - (a[2U] * b[1U]);
    const float cross_y = (a[2U] * b[0U]) - (a[0U] * b[2U]);
    const float cross_z = (a[0U] * b[1U]) - (a[1U] * b[0U]);
    return (cross_x * c[0U]) + (cross_y * c[1U]) + (cross_z * c[2U]);
}

}  // namespace

bool LisGridMatcher::chiralityAgrees(
    const std::array<Bearing, Centroiding::kMaxCentroids>& bearings,
    std::size_t first,
    std::size_t second,
    std::size_t candidate) const {
    // Degenerate triples give a triple product near zero, where the sign is
    // noise rather than handedness. Accept those instead of rejecting a valid
    // correspondence on a coin flip.
    const float kEpsilon = 1e-4F;

    const float observed = tripleProduct(
        bearings[nodes_[first].centroid],
        bearings[nodes_[second].centroid],
        bearings[nodes_[candidate].centroid]);

    const CatalogStar& ca =
        (*catalog_)[indexed_catalog_positions_[nodes_[first].indexed_star]];
    const CatalogStar& cb =
        (*catalog_)[indexed_catalog_positions_[nodes_[second].indexed_star]];
    const CatalogStar& cc =
        (*catalog_)[indexed_catalog_positions_[nodes_[candidate].indexed_star]];
    const std::array<float, 3U> va = {ca.x, ca.y, ca.z};
    const std::array<float, 3U> vb = {cb.x, cb.y, cb.z};
    const std::array<float, 3U> vc = {cc.x, cc.y, cc.z};
    const float expected = tripleProduct(va, vb, vc);

    if ((std::fabs(observed) <= kEpsilon) ||
        (std::fabs(expected) <= kEpsilon)) {
        return true;
    }
    return (observed > 0.0F) == (expected > 0.0F);
}

std::size_t LisGridMatcher::findBestClique(
    std::size_t node_count,
    const std::array<Bearing, Centroiding::kMaxCentroids>& bearings,
    std::array<std::uint16_t, kMaxNodes>& best,
    float& best_rms_arcsec,
    std::size_t& expansions) const {
    std::size_t best_size = 0U;
    best_rms_arcsec = -1.0F;
    expansions = 0U;
    std::array<std::uint16_t, kMaxNodes> current = {};

    for (std::size_t seed = 0U; seed < node_count; ++seed) {
        std::size_t size = 0U;
        current[size] = static_cast<std::uint16_t>(seed);
        ++size;
        double sum_squared = 0.0;
        std::size_t pair_count = 0U;

        bool grew = true;
        while (grew && (size < kMaxNodes)) {
            grew = false;

            // Grow by *lowest added error* rather than first-fit. First-fit
            // let a spurious node join early and drag the clique away from the
            // correct correspondence set; the correct star is almost always
            // the closest-agreeing one.
            std::size_t chosen = node_count;
            double chosen_cost = 0.0;
            double chosen_sum_squared = 0.0;

            for (std::size_t candidate = 0U; candidate < node_count;
                 ++candidate) {
                bool compatible = true;
                for (std::size_t member = 0U; member < size; ++member) {
                    ++expansions;
                    if ((candidate == current[member]) ||
                        !adjacent(candidate, current[member])) {
                        compatible = false;
                        break;
                    }
                }
                if (!compatible) {
                    continue;
                }

                // Chirality: a mirrored field preserves every pairwise
                // distance, so the edges above cannot see it. The sign of the
                // scalar triple product can.
                if ((size >= 2U) &&
                    !chiralityAgrees(bearings, current[0U], current[1U],
                                     candidate)) {
                    continue;
                }

                double added = 0.0;
                for (std::size_t member = 0U; member < size; ++member) {
                    const float error = edge_error_[candidate][current[member]];
                    added += static_cast<double>(error) *
                             static_cast<double>(error);
                }
                const double cost = added / static_cast<double>(size);

                if ((chosen >= node_count) || (cost < chosen_cost)) {
                    chosen = candidate;
                    chosen_cost = cost;
                    chosen_sum_squared = added;
                }
            }

            if (chosen < node_count) {
                current[size] = static_cast<std::uint16_t>(chosen);
                sum_squared += chosen_sum_squared;
                pair_count += size;
                ++size;
                grew = true;
            }
        }

        if ((size < config_.minimum_clique_size) || (pair_count == 0U)) {
            continue;
        }

        const float rms = static_cast<float>(
            std::sqrt(sum_squared / static_cast<double>(pair_count)));

        // Verification: a clique whose pairwise angles disagree by far more
        // than centroid noise is a coincidence, however large it is.
        if (rms > config_.maximum_clique_rms_arcsec) {
            continue;
        }

        // Rank hypotheses: more corroboration first, then tighter agreement.
        if ((size > best_size) ||
            ((size == best_size) && (rms < best_rms_arcsec))) {
            best_size = size;
            best_rms_arcsec = rms;
            for (std::size_t index = 0U; index < size; ++index) {
                best[index] = current[index];
            }
        }
    }

    return best_size;
}

LisGridMatcher::MatchResult LisGridMatcher::match(
    const Centroiding::Result& centroids) {
    MatchResult result;
    if (!index_ready_ || centroids.count < 3U) {
        return result;
    }

    for (std::size_t row = 0U; row < kMaxCliqueObserved; ++row) {
        votes_[row].fill(0U);
    }

    std::array<std::uint16_t, Centroiding::kMaxCentroids> selected = {};
    std::size_t selected_count = 0U;
    selectObservedStars(centroids, selected, selected_count);

    std::array<Bearing, Centroiding::kMaxCentroids> rays = {};
    for (std::size_t selection = 0U; selection < selected_count; ++selection) {
        const std::size_t centroid_index = selected[selection];
        const Point& point = centroids.points[centroid_index];
        bearingFor(point.x, point.y, rays[centroid_index]);
    }

    for (std::size_t first_selection = 0U;
         first_selection < selected_count; ++first_selection) {
        const std::size_t first_centroid = selected[first_selection];
        for (std::size_t second_selection = first_selection + 1U;
             second_selection < selected_count; ++second_selection) {
            const std::size_t second_centroid = selected[second_selection];
            const float dot =
                rays[first_centroid][0U] * rays[second_centroid][0U] +
                rays[first_centroid][1U] * rays[second_centroid][1U] +
                rays[first_centroid][2U] * rays[second_centroid][2U];
            const float observed_angle = std::acos(clampUnit(dot));
            if (observed_angle > maximum_indexed_angle_radians_) {
                continue;
            }

            const float lower_angle =
                std::max(0.0F, observed_angle - angular_tolerance_radians_);
            const float upper_angle = std::min(
                maximum_indexed_angle_radians_,
                observed_angle + angular_tolerance_radians_);
            const std::size_t lower_bin = angleToBin(lower_angle);
            const std::size_t upper_bin = angleToBin(upper_angle);

            // Compare cosines, not angles. cos is monotonically decreasing on
            // [0, pi], so |theta_cat - theta_obs| <= tol is exactly
            // cos(theta_obs + tol) <= cos(theta_cat) <= cos(theta_obs - tol).
            // Two cosines computed once per observed pair replace one acos per
            // *catalogue entry scanned* -- and tens of entries are scanned per
            // pair. On a part without hardware transcendentals this was the
            // dominant cost of the whole matcher.
            const float cosine_low =
                std::cos(observed_angle + angular_tolerance_radians_);
            const float cosine_high =
                std::cos(observed_angle - angular_tolerance_radians_);

            for (std::size_t bin = lower_bin; bin <= upper_bin; ++bin) {
                std::uint16_t entry_index = bin_heads_[bin];
                while (entry_index != kInvalidEntry) {
                    const PairEntry& entry = pair_entries_[entry_index];
                    const CatalogStar& first_star =
                        (*catalog_)[indexed_catalog_positions_[entry.first_star]];
                    const CatalogStar& second_star =
                        (*catalog_)[indexed_catalog_positions_[entry.second_star]];
                    const float catalog_cosine =
                        catalogPairCosine(first_star, second_star);
                    if ((catalog_cosine >= cosine_low) &&
                        (catalog_cosine <= cosine_high)) {
                        addVote(first_selection, entry.first_star);
                        addVote(first_selection, entry.second_star);
                        addVote(second_selection, entry.first_star);
                        addVote(second_selection, entry.second_star);
                    }
                    entry_index = entry.next;
                }
            }
        }
    }

    // Votes only shortlist. The accepted set is the largest subset of
    // candidate correspondences in which *every* pair agrees -- so a star can
    // no longer be accepted merely because it scored well in isolation, which
    // is what allowed mutually contradictory assignments to survive before.
    const std::size_t node_count = buildNodes(selected, selected_count);
    result.node_count = node_count;
    if (node_count == 0U) {
        return result;
    }

    buildAdjacency(node_count, rays);

    std::array<std::uint16_t, kMaxNodes> clique = {};
    float clique_rms = -1.0F;
    std::size_t expansions = 0U;
    const std::size_t clique_size =
        findBestClique(node_count, rays, clique, clique_rms, expansions);
    result.node_count = node_count;
    result.expansions = expansions;

    result.clique_size = clique_size;
    if (clique_size < config_.minimum_clique_size) {
        return result;  // refuse: no sufficiently corroborated consistent set
    }

    for (std::size_t member = 0U; member < clique_size; ++member) {
        const Node& node = nodes_[clique[member]];
        const std::size_t centroid_index = node.centroid;
        result.star_ids[centroid_index] =
            (*catalog_)[indexed_catalog_positions_[node.indexed_star]].star_id;
        // Report clique size per member: every one of them is corroborated by
        // every other, which is a far stronger statement than a vote tally.
        result.votes[centroid_index] = static_cast<std::uint8_t>(
            (clique_size > 255U) ? 255U : clique_size);
        ++result.matched_count;
    }

    return result;
}

void LisGridMatcher::setFocalLength(float focal_length_pixels) {
    if (focal_length_pixels > 0.0F) {
        config_.focal_length_pixels = focal_length_pixels;
    }
}

std::size_t LisGridMatcher::indexedStarCount() const {
    return indexed_star_count_;
}

std::size_t LisGridMatcher::indexedPairCount() const {
    return pair_entry_count_;
}

}  // namespace star_tracker

static_assert(
    sizeof(star_tracker::LisGridMatcher::Catalog) +
            sizeof(star_tracker::LisGridMatcher) <
        500U * 1024U,
    "LIS matcher and catalog exceed the 500 KB RAM budget");
