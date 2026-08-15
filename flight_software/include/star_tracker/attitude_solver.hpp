#ifndef STAR_TRACKER_ATTITUDE_SOLVER_HPP
#define STAR_TRACKER_ATTITUDE_SOLVER_HPP

#include "star_tracker/centroiding.hpp"
#include "star_tracker/lis_grid_matcher.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace star_tracker {

/// Stage 04 + 07 + 06: bearing projection, attitude estimation, integrity gate.
///
/// Deliberately owns the camera model rather than burying it inside the matcher,
/// so focal length is an explicit, testable input. Both reference implementations
/// hide this transform inside star-identification, which is why a focal-length
/// error there degrades into a confident wrong answer instead of a refusal.
///
/// Frame convention matches LOST and tetra3 so outputs are directly comparable:
///   +x boresight, +y image-left, +z image-up. Right-handed.
/// (LisGridMatcher uses a different axis order internally, but consumes only
///  dot products, so it is unaffected by the difference.)
class AttitudeSolver {
public:
    static constexpr std::size_t kMaxStars = Centroiding::kMaxCentroids;

    struct Config {
        float focal_length_pixels;
        float principal_x;
        float principal_y;

        /// Integrity gate: reject when reprojection RMS exceeds this.
        ///
        /// Derived, not tuned. Over 190 correct solves spanning uniform sky and
        /// up to 50% bright-outlier contamination, the reprojection RMS ranged
        /// 3.0 to 14.7 arcsec (median ~4.8). 30 arcsec leaves 2x margin on the
        /// worst observed case.
        ///
        /// Note this is now a *backstop*, not the primary discriminator: the
        /// clique consistency check and chirality test eliminate wrong
        /// correspondences upstream, which is why no false solve survived to
        /// this stage in any sweep. Re-derive per optical design -- the value
        /// scales with centroid noise and inverse focal length.
        float max_residual_arcsec;

        /// Integrity gate: minimum matched stars. Four is the geometric
        /// minimum for a verified solve; more is corroboration.
        std::size_t minimum_matched_stars;

        Config()
            : focal_length_pixels(2904.0F),
              principal_x(512.0F),
              principal_y(512.0F),
              max_residual_arcsec(30.0F),
              minimum_matched_stars(5U) {}
    };

    enum Reason {
        kAccepted = 0,
        kTooFewStars = 1,
        kDegenerateGeometry = 2,
        kResidualTooLarge = 3
    };

    struct Solution {
        bool valid;
        /// Inertial -> body rotation. Storage order (x, y, z, w).
        std::array<float, 4U> quaternion;
        float right_ascension_deg;
        float declination_deg;
        float roll_deg;
        float residual_rms_arcsec;
        float residual_max_arcsec;
        std::size_t matched_count;
        Reason reason;

        Solution()
            : valid(false),
              quaternion(),
              right_ascension_deg(0.0F),
              declination_deg(0.0F),
              roll_deg(0.0F),
              residual_rms_arcsec(-1.0F),
              residual_max_arcsec(-1.0F),
              matched_count(0U),
              reason(kTooFewStars) {
            quaternion[3U] = 1.0F;
        }
    };

    AttitudeSolver();
    explicit AttitudeSolver(const Config& config);

    /// Stage 04, exposed so it can be tested and swapped independently.
    /// Returns a unit bearing in the camera frame for a pixel coordinate.
    void pixelToBearing(float pixel_x, float pixel_y, float bearing[3U]) const;

    Solution solve(
        const Centroiding::Result& centroids,
        const LisGridMatcher::MatchResult& matches,
        const LisGridMatcher::Catalog& catalog,
        std::size_t catalog_count);

private:
    /// QUEST (Shuster & Oh). Chosen over Davenport's q-method because it needs
    /// no eigendecomposition, and over Kabsch because it needs no SVD -- both
    /// of which would pull in a matrix library this codebase forbids.
    bool questSolve(
        const std::array<std::array<float, 3U>, kMaxStars>& body,
        const std::array<std::array<float, 3U>, kMaxStars>& reference,
        std::size_t count,
        std::array<float, 4U>& quaternion) const;

    static void quaternionToMatrix(
        const std::array<float, 4U>& quaternion, float matrix[3U][3U]);

    Config config_;
};

}  // namespace star_tracker

#endif
