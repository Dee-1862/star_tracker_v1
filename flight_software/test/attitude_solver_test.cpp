#include "star_tracker/attitude_solver.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

constexpr float kFocalLength = 2904.0F;
constexpr float kPrincipal = 512.0F;
constexpr std::size_t kStarCount = 8U;

star_tracker::LisGridMatcher::Catalog catalog;

/// Rotate a body-frame bearing into the inertial frame by the transpose of a
/// z-y-x rotation, giving us a catalogue vector consistent with a known
/// attitude. Mirrors AttitudeSolver::quaternionToMatrix's convention.
void bodyToInertial(
    const float body[3U],
    float right_ascension_rad,
    float declination_rad,
    float roll_rad,
    float inertial[3U]) {
    const float ca = std::cos(right_ascension_rad);
    const float sa = std::sin(right_ascension_rad);
    const float cd = std::cos(declination_rad);
    const float sd = std::sin(declination_rad);
    const float cr = std::cos(roll_rad);
    const float sr = std::sin(roll_rad);

    // Body-from-inertial matrix rows: boresight, image-left, image-up.
    const float m[3U][3U] = {
        {cd * ca, cd * sa, sd},
        {(-sr * sd * ca) - (cr * sa), (-sr * sd * sa) + (cr * ca), sr * cd},
        {(-cr * sd * ca) + (sr * sa), (-cr * sd * sa) - (sr * ca), cr * cd}};

    for (std::size_t col = 0U; col < 3U; ++col) {
        inertial[col] = (m[0U][col] * body[0U]) + (m[1U][col] * body[1U]) +
                        (m[2U][col] * body[2U]);
    }
}

float angularSeparationDegrees(
    float ra_a_deg, float dec_a_deg, float ra_b_deg, float dec_b_deg) {
    const float to_rad = 3.14159265358979F / 180.0F;
    const float ra_a = ra_a_deg * to_rad;
    const float dec_a = dec_a_deg * to_rad;
    const float ra_b = ra_b_deg * to_rad;
    const float dec_b = dec_b_deg * to_rad;
    const float cosine = (std::sin(dec_a) * std::sin(dec_b)) +
                         (std::cos(dec_a) * std::cos(dec_b) *
                          std::cos(ra_a - ra_b));
    const float clamped = (cosine > 1.0F) ? 1.0F
                                          : ((cosine < -1.0F) ? -1.0F : cosine);
    return std::acos(clamped) * 180.0F / 3.14159265358979F;
}

/// Build a consistent scene: pixel positions, matching catalogue entries for a
/// known attitude, and the centroid/match structures the solver consumes.
void buildScene(
    float ra_deg,
    float dec_deg,
    float roll_deg,
    star_tracker::Centroiding::Result& centroids,
    star_tracker::LisGridMatcher::MatchResult& matches) {
    const float to_rad = 3.14159265358979F / 180.0F;
    const float pixels[kStarCount][2U] = {
        {300.0F, 280.0F}, {700.0F, 320.0F}, {480.0F, 640.0F}, {820.0F, 700.0F},
        {200.0F, 750.0F}, {560.0F, 200.0F}, {380.0F, 480.0F}, {740.0F, 520.0F}};

    star_tracker::AttitudeSolver projector;
    centroids.count = kStarCount;
    matches.matched_count = kStarCount;

    for (std::size_t index = 0U; index < kStarCount; ++index) {
        centroids.points[index].x = pixels[index][0U];
        centroids.points[index].y = pixels[index][1U];
        centroids.points[index].intensity = 1000U;

        float bearing[3U];
        projector.pixelToBearing(pixels[index][0U], pixels[index][1U], bearing);

        float inertial[3U];
        bodyToInertial(bearing, ra_deg * to_rad, dec_deg * to_rad,
                       roll_deg * to_rad, inertial);

        const std::uint32_t star_id = static_cast<std::uint32_t>(index + 1U);
        catalog[index].star_id = star_id;
        catalog[index].x = inertial[0U];
        catalog[index].y = inertial[1U];
        catalog[index].z = inertial[2U];
        catalog[index].visual_magnitude = 4.0F;

        matches.star_ids[index] = star_id;
        matches.votes[index] = 8U;
    }
}

}  // namespace

int main() {
    const float ra_deg = 137.5F;
    const float dec_deg = -22.25F;
    const float roll_deg = 64.0F;

    star_tracker::AttitudeSolver solver;

    // --- 1. Correct identifications must recover the attitude ---------------
    {
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        buildScene(ra_deg, dec_deg, roll_deg, centroids, matches);

        const star_tracker::AttitudeSolver::Solution solution =
            solver.solve(centroids, matches, catalog, kStarCount);

        assert(solution.valid);
        assert(solution.reason == star_tracker::AttitudeSolver::kAccepted);
        assert(solution.matched_count == kStarCount);

        const float pointing_error = angularSeparationDegrees(
            solution.right_ascension_deg, solution.declination_deg,
            ra_deg, dec_deg);
        assert(pointing_error < 0.01F);

        float roll_error = solution.roll_deg - roll_deg;
        while (roll_error > 180.0F) {
            roll_error -= 360.0F;
        }
        while (roll_error < -180.0F) {
            roll_error += 360.0F;
        }
        assert(std::fabs(roll_error) < 0.05F);

        // A consistent scene must produce a near-zero residual.
        assert(solution.residual_rms_arcsec >= 0.0F);
        assert(solution.residual_rms_arcsec < 1.0F);
    }

    // --- 2. The integrity gate must refuse scrambled identifications --------
    // This is the case both baselines get wrong: the matcher is confident, the
    // geometry is self-consistent for the quad it chose, and the attitude is
    // nonsense. The residual is the only thing that knows.
    {
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        buildScene(ra_deg, dec_deg, roll_deg, centroids, matches);

        // Rotate the assignment: every centroid now points at its neighbour's
        // catalogue star. Distances are unchanged; the mapping is wrong.
        const std::uint32_t first = matches.star_ids[0U];
        for (std::size_t index = 0U; index + 1U < kStarCount; ++index) {
            matches.star_ids[index] = matches.star_ids[index + 1U];
        }
        matches.star_ids[kStarCount - 1U] = first;

        const star_tracker::AttitudeSolver::Solution solution =
            solver.solve(centroids, matches, catalog, kStarCount);

        assert(!solution.valid);
        assert(solution.reason ==
               star_tracker::AttitudeSolver::kResidualTooLarge);
        assert(solution.residual_rms_arcsec > 28.0F);
    }

    // --- 3. The residual must MEASURE error, not manufacture it -------------
    // Regression test for a real defect: the residual was computed as
    // acos(dot(predicted, observed)). Near zero separation acos is
    // catastrophically ill-conditioned -- float rounding alone yields
    // sqrt(2*eps) ~ 65 arcsec of phantom residual, which silently vetoed every
    // correct solve. Tests 1 and 2 could not see it, because their catalogue is
    // built from the same floats as the bearings, so dot came out exactly 1.0.
    //
    // Here each catalogue vector is nudged by a known angle through an
    // independent numeric path. A correct implementation reports roughly that
    // angle; the acos version reported ~65 arcsec regardless.
    {
        const float kPerturbArcsec = 10.0F;
        const float kPerturbRadians = kPerturbArcsec / 206264.806247F;

        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        buildScene(ra_deg, dec_deg, roll_deg, centroids, matches);

        for (std::size_t index = 0U; index < kStarCount; ++index) {
            // Displace along an orthogonal direction by a known small angle.
            const double x = static_cast<double>(catalog[index].x);
            const double y = static_cast<double>(catalog[index].y);
            const double z = static_cast<double>(catalog[index].z);
            double ox = -y;
            double oy = x;
            double oz = 0.0;
            const double norm = std::sqrt((ox * ox) + (oy * oy) + (oz * oz));
            if (norm < 1e-9) {
                continue;
            }
            ox /= norm;
            oy /= norm;
            oz /= norm;
            const double t = static_cast<double>(kPerturbRadians);
            double nx = x + (t * ox);
            double ny = y + (t * oy);
            double nz = z + (t * oz);
            const double inverse =
                1.0 / std::sqrt((nx * nx) + (ny * ny) + (nz * nz));
            catalog[index].x = static_cast<float>(nx * inverse);
            catalog[index].y = static_cast<float>(ny * inverse);
            catalog[index].z = static_cast<float>(nz * inverse);
        }

        star_tracker::AttitudeSolver::Config config;
        config.max_residual_arcsec = 1000.0F;  // open the gate; we want the value
        star_tracker::AttitudeSolver measuring(config);

        const star_tracker::AttitudeSolver::Solution solution =
            measuring.solve(centroids, matches, catalog, kStarCount);

        assert(solution.valid);
        // The fit absorbs part of a coherent perturbation, so demand only that
        // the reported residual is the right order of magnitude -- and above
        // all that it is nowhere near the ~65 arcsec floating-point floor the
        // broken metric produced for ANY input.
        assert(solution.residual_rms_arcsec >= 0.0F);
        assert(solution.residual_rms_arcsec < 3.0F * kPerturbArcsec);
    }

    // --- 4. Focal length must be recoverable from the matched stars ---------
    // Each matched pair is a calibration standard: the catalogue fixes their
    // angular separation, the image measures their pixel separation, and only
    // one focal length reconciles the two. This is what lets the tracker
    // recalibrate itself after thermal drift instead of needing a ground pass.
    {
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        buildScene(ra_deg, dec_deg, roll_deg, centroids, matches);

        // The scene was built at kFocalLength. Start the search badly wrong in
        // both directions and require recovery to better than 0.1%.
        const float truth = kFocalLength;
        const float starts[3] = {truth * 0.90F, truth * 1.08F, truth};
        for (std::size_t trial = 0U; trial < 3U; ++trial) {
            const float recovered = star_tracker::AttitudeSolver::refineFocalLength(
                centroids, matches, catalog, kStarCount,
                kPrincipal, kPrincipal, starts[trial]);
            const float relative = std::fabs(recovered - truth) / truth;
            assert(relative < 0.001F);
        }

        // Too few pairs to constrain the fit: return the input untouched
        // rather than inventing a number.
        star_tracker::LisGridMatcher::MatchResult sparse;
        sparse.star_ids[0U] = matches.star_ids[0U];
        sparse.matched_count = 1U;
        const float unchanged = star_tracker::AttitudeSolver::refineFocalLength(
            centroids, sparse, catalog, kStarCount,
            kPrincipal, kPrincipal, 1234.0F);
        assert(std::fabs(unchanged - 1234.0F) < 1e-3F);
    }

    // --- 5. Too few matched stars must refuse, not guess --------------------
    {
        star_tracker::Centroiding::Result centroids;
        star_tracker::LisGridMatcher::MatchResult matches;
        buildScene(ra_deg, dec_deg, roll_deg, centroids, matches);

        for (std::size_t index = 2U; index < kStarCount; ++index) {
            matches.star_ids[index] = 0U;  // unmatched sentinel
        }

        const star_tracker::AttitudeSolver::Solution solution =
            solver.solve(centroids, matches, catalog, kStarCount);

        assert(!solution.valid);
        assert(solution.reason == star_tracker::AttitudeSolver::kTooFewStars);
        assert(solution.matched_count == 2U);
    }

    return 0;
}
