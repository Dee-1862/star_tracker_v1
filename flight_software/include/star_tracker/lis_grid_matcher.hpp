#ifndef STAR_TRACKER_LIS_GRID_MATCHER_HPP
#define STAR_TRACKER_LIS_GRID_MATCHER_HPP

#include "star_tracker/centroiding.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace star_tracker {

class LisGridMatcher {
public:
    static constexpr std::size_t kMaxCatalogStars = 5041U;
    /// A 20 deg field holds ~34 stars from a magnitude-6 catalogue, but only
    /// ~5 of those fall in the 768 brightest of the whole sky -- so the correct
    /// identification was usually absent from the index entirely. 2048 raises
    /// in-field coverage to ~12, enough to support a corroborated clique.
    static constexpr std::size_t kMaxIndexedStars = 2048U;
    static constexpr std::size_t kMaxPairEntries = 32768U;
    static constexpr std::size_t kAngleBinCount = 4096U;

    /// Consistency-graph sizing. Pair voting only shortlists candidates; the
    /// accepted set is the largest mutually consistent clique among them.
    static constexpr std::size_t kMaxCliqueObserved = 20U;
    static constexpr std::size_t kCandidatesPerCentroid = 4U;
    static constexpr std::size_t kMaxNodes =
        kMaxCliqueObserved * kCandidatesPerCentroid;
    static constexpr std::size_t kAdjacencyWords = (kMaxNodes + 31U) / 32U;

    struct CatalogStar {
        std::uint32_t star_id;
        float x;
        float y;
        float z;
        float visual_magnitude;
    };

    typedef std::array<CatalogStar, kMaxCatalogStars> Catalog;

    struct Config {
        float focal_length_pixels;
        float principal_x;
        float principal_y;
        float maximum_indexed_angle_degrees;
        float angular_tolerance_degrees;
        std::size_t maximum_observed_stars;
        std::uint8_t minimum_votes;

        /// Smallest mutually consistent set accepted. Four is the geometric
        /// minimum; five demands one independent corroborating star, which is
        /// what separated true from false solves in calibration.
        std::size_t minimum_clique_size;

        /// Verification threshold. A correct clique's pairwise angles agree to
        /// roughly the centroid noise (~10 arcsec at f=2904); a spurious one
        /// is spread across the whole edge tolerance. Rejecting on this is
        /// what stops a large-but-wrong clique from beating a smaller correct
        /// one.
        float maximum_clique_rms_arcsec;

        Config()
            : focal_length_pixels(2904.0F),
              principal_x(512.0F),
              principal_y(512.0F),
              // Narrowed from 30 deg: pair count grows as N^2*(1-cos(theta)),
              // so indexing 2048 stars only fits the 32768-entry budget with a
              // tighter angle. A 20 deg field still yields ample sub-10 deg
              // pairs.
              maximum_indexed_angle_degrees(10.0F),
              // Was 0.02 deg (72 arcsec) -- roughly ten times the centroid
              // noise, loose enough to admit spurious edges. 0.010 deg is
              // ~3 sigma on a 10 arcsec pair-angle error.
              angular_tolerance_degrees(0.010F),
              maximum_observed_stars(20U),
              minimum_votes(3U),
              minimum_clique_size(5U),
              maximum_clique_rms_arcsec(20.0F) {}
    };

    struct MatchResult {
        std::array<std::uint32_t, Centroiding::kMaxCentroids> star_ids;
        std::array<std::uint8_t, Centroiding::kMaxCentroids> votes;
        std::size_t matched_count;
        /// Size of the accepted clique. Doubles as the confidence measure:
        /// every member agrees pairwise with every other, so a large clique
        /// arising by chance is vanishingly unlikely.
        std::size_t clique_size;

        MatchResult()
            : star_ids(), votes(), matched_count(0U), clique_size(0U) {}
    };

    LisGridMatcher();
    explicit LisGridMatcher(const Config& config);

    bool buildIndex(const Catalog& catalog, std::size_t catalog_count);
    MatchResult match(const Centroiding::Result& centroids);
    std::size_t indexedStarCount() const;
    std::size_t indexedPairCount() const;

private:
    static constexpr std::uint16_t kInvalidEntry = 0xFFFFU;

    struct PairEntry {
        std::uint16_t first_star;
        std::uint16_t second_star;
        std::uint16_t next;
    };

    /// Rows are indexed by *selection* order, not centroid index: only the
    /// brightest kMaxCliqueObserved centroids ever vote, so sizing the table
    /// to all 50 slots wasted 60 KB of a 500 KB budget.
    typedef std::array<
        std::array<std::uint8_t, kMaxIndexedStars>,
        kMaxCliqueObserved>
        VoteTable;

    /// One candidate correspondence: this observed centroid might be this
    /// indexed catalogue star. Nodes of the consistency graph.
    struct Node {
        std::uint16_t centroid;
        std::uint16_t indexed_star;
    };

    typedef std::array<float, 3U> Bearing;

    std::size_t angleToBin(float angle_radians) const;
    static float catalogPairAngle(
        const CatalogStar& first, const CatalogStar& second);
    void selectBrightestStars();
    void selectObservedStars(
        const Centroiding::Result& centroids,
        std::array<std::uint16_t, Centroiding::kMaxCentroids>& selected,
        std::size_t& selected_count) const;
    void addVote(std::size_t centroid_index, std::size_t indexed_star);

    void bearingFor(float pixel_x, float pixel_y, Bearing& bearing) const;
    /// Shortlist the best-voted catalogue stars per centroid into graph nodes.
    std::size_t buildNodes(
        const std::array<std::uint16_t, Centroiding::kMaxCentroids>& selected,
        std::size_t selected_count);
    /// Edge iff the observed and catalogue separations agree within tolerance.
    void buildAdjacency(
        std::size_t node_count,
        const std::array<Bearing, Centroiding::kMaxCentroids>& bearings);
    /// Bounded greedy maximal-clique search, seeded from every node. Grows by
    /// lowest added angular error rather than first-fit, enforces chirality
    /// (which pairwise distances alone cannot see), and verifies each finished
    /// clique against maximum_clique_rms_arcsec. Reports the accepted clique's
    /// RMS so callers can rank hypotheses.
    std::size_t findBestClique(
        std::size_t node_count,
        const std::array<Bearing, Centroiding::kMaxCentroids>& bearings,
        std::array<std::uint16_t, kMaxNodes>& best,
        float& best_rms_arcsec) const;
    bool adjacent(std::size_t first, std::size_t second) const;
    bool chiralityAgrees(
        const std::array<Bearing, Centroiding::kMaxCentroids>& bearings,
        std::size_t first,
        std::size_t second,
        std::size_t candidate) const;

    Config config_;
    const Catalog* catalog_;
    std::size_t catalog_count_;
    std::size_t indexed_star_count_;
    std::size_t pair_entry_count_;
    float maximum_indexed_angle_radians_;
    float angular_tolerance_radians_;
    bool index_ready_;
    std::array<std::uint16_t, kMaxIndexedStars> indexed_catalog_positions_;
    std::array<std::uint16_t, kAngleBinCount> bin_heads_;
    std::array<PairEntry, kMaxPairEntries> pair_entries_;
    VoteTable votes_;
    std::array<Node, kMaxNodes> nodes_;
    std::array<std::array<std::uint32_t, kAdjacencyWords>, kMaxNodes> adjacency_;
    /// |observed - catalogue| separation per edge, arcseconds. Cached from
    /// adjacency construction so clique growth can rank candidates without
    /// recomputing two arc cosines per comparison.
    std::array<std::array<float, kMaxNodes>, kMaxNodes> edge_error_;
};

}  // namespace star_tracker

#endif
