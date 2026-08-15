#include "star_tracker/attitude_solver.hpp"

#include <cmath>

namespace star_tracker {

namespace {

constexpr float kPi = 3.14159265358979323846F;
constexpr float kRadiansToArcseconds = 206264.806247F;
constexpr float kRadiansToDegrees = 57.2957795130823F;

float dot3(const float a[3U], const float b[3U]) {
    return (a[0U] * b[0U]) + (a[1U] * b[1U]) + (a[2U] * b[2U]);
}

float clampUnit(float value) {
    if (value > 1.0F) {
        return 1.0F;
    }
    if (value < -1.0F) {
        return -1.0F;
    }
    return value;
}

}  // namespace

AttitudeSolver::AttitudeSolver() : config_() {}

AttitudeSolver::AttitudeSolver(const Config& config) : config_(config) {}

void AttitudeSolver::pixelToBearing(
    float pixel_x, float pixel_y, float bearing[3U]) const {
    // LOST/tetra3 convention: +y is image-left, +z is image-up, so both
    // offsets are (principal - pixel), not (pixel - principal).
    const float y = (config_.principal_x - pixel_x) / config_.focal_length_pixels;
    const float z = (config_.principal_y - pixel_y) / config_.focal_length_pixels;
    const float inverse_norm = 1.0F / std::sqrt(1.0F + (y * y) + (z * z));
    bearing[0U] = inverse_norm;
    bearing[1U] = y * inverse_norm;
    bearing[2U] = z * inverse_norm;
}

void AttitudeSolver::quaternionToMatrix(
    const std::array<float, 4U>& q, float m[3U][3U]) {
    const float x = q[0U];
    const float y = q[1U];
    const float z = q[2U];
    const float w = q[3U];

    m[0U][0U] = 1.0F - (2.0F * ((y * y) + (z * z)));
    m[0U][1U] = 2.0F * ((x * y) + (z * w));
    m[0U][2U] = 2.0F * ((x * z) - (y * w));
    m[1U][0U] = 2.0F * ((x * y) - (z * w));
    m[1U][1U] = 1.0F - (2.0F * ((x * x) + (z * z)));
    m[1U][2U] = 2.0F * ((y * z) + (x * w));
    m[2U][0U] = 2.0F * ((x * z) + (y * w));
    m[2U][1U] = 2.0F * ((y * z) - (x * w));
    m[2U][2U] = 1.0F - (2.0F * ((x * x) + (y * y)));
}

bool AttitudeSolver::questSolve(
    const std::array<std::array<float, 3U>, kMaxStars>& body,
    const std::array<std::array<float, 3U>, kMaxStars>& reference,
    std::size_t count,
    std::array<float, 4U>& quaternion) const {
    if (count < 2U) {
        return false;
    }

    // B = sum w_i b_i r_i^T, all weights unity.
    double b_matrix[3U][3U] = {{0.0}};
    for (std::size_t index = 0U; index < count; ++index) {
        for (std::size_t row = 0U; row < 3U; ++row) {
            for (std::size_t col = 0U; col < 3U; ++col) {
                b_matrix[row][col] +=
                    static_cast<double>(body[index][row]) *
                    static_cast<double>(reference[index][col]);
            }
        }
    }

    const double sigma = b_matrix[0U][0U] + b_matrix[1U][1U] + b_matrix[2U][2U];

    double s_matrix[3U][3U];
    for (std::size_t row = 0U; row < 3U; ++row) {
        for (std::size_t col = 0U; col < 3U; ++col) {
            s_matrix[row][col] = b_matrix[row][col] + b_matrix[col][row];
        }
    }

    const double z_vector[3U] = {
        b_matrix[1U][2U] - b_matrix[2U][1U],
        b_matrix[2U][0U] - b_matrix[0U][2U],
        b_matrix[0U][1U] - b_matrix[1U][0U]};

    const double det_s =
        (s_matrix[0U][0U] *
         ((s_matrix[1U][1U] * s_matrix[2U][2U]) - (s_matrix[1U][2U] * s_matrix[2U][1U]))) -
        (s_matrix[0U][1U] *
         ((s_matrix[1U][0U] * s_matrix[2U][2U]) - (s_matrix[1U][2U] * s_matrix[2U][0U]))) +
        (s_matrix[0U][2U] *
         ((s_matrix[1U][0U] * s_matrix[2U][1U]) - (s_matrix[1U][1U] * s_matrix[2U][0U])));

    // kappa = trace(adjugate(S))
    const double kappa =
        ((s_matrix[1U][1U] * s_matrix[2U][2U]) - (s_matrix[1U][2U] * s_matrix[2U][1U])) +
        ((s_matrix[0U][0U] * s_matrix[2U][2U]) - (s_matrix[0U][2U] * s_matrix[2U][0U])) +
        ((s_matrix[0U][0U] * s_matrix[1U][1U]) - (s_matrix[0U][1U] * s_matrix[1U][0U]));

    double s_z[3U];
    for (std::size_t row = 0U; row < 3U; ++row) {
        s_z[row] = (s_matrix[row][0U] * z_vector[0U]) +
                   (s_matrix[row][1U] * z_vector[1U]) +
                   (s_matrix[row][2U] * z_vector[2U]);
    }
    double s2_z[3U];
    for (std::size_t row = 0U; row < 3U; ++row) {
        s2_z[row] = (s_matrix[row][0U] * s_z[0U]) +
                    (s_matrix[row][1U] * s_z[1U]) +
                    (s_matrix[row][2U] * s_z[2U]);
    }

    const double z_dot_z = (z_vector[0U] * z_vector[0U]) +
                           (z_vector[1U] * z_vector[1U]) +
                           (z_vector[2U] * z_vector[2U]);
    const double z_s_z = (z_vector[0U] * s_z[0U]) + (z_vector[1U] * s_z[1U]) +
                         (z_vector[2U] * s_z[2U]);
    const double z_s2_z = (z_vector[0U] * s2_z[0U]) + (z_vector[1U] * s2_z[1U]) +
                          (z_vector[2U] * s2_z[2U]);

    const double a_coefficient = (sigma * sigma) - kappa;
    const double b_coefficient = (sigma * sigma) + z_dot_z;
    const double c_coefficient = det_s + z_s_z;
    const double d_coefficient = z_s2_z;

    // Newton-Raphson on the characteristic polynomial; the initial guess is the
    // sum of weights, which is exact for noiseless, consistent observations.
    double lambda = static_cast<double>(count);
    for (std::size_t iteration = 0U; iteration < 24U; ++iteration) {
        const double lambda2 = lambda * lambda;
        const double f = (lambda2 * lambda2) -
                         ((a_coefficient + b_coefficient) * lambda2) -
                         (c_coefficient * lambda) +
                         ((a_coefficient * b_coefficient) +
                          (c_coefficient * sigma) - d_coefficient);
        const double df = (4.0 * lambda2 * lambda) -
                          (2.0 * (a_coefficient + b_coefficient) * lambda) -
                          c_coefficient;
        if ((df < 1e-12) && (df > -1e-12)) {
            break;
        }
        const double step = f / df;
        lambda -= step;
        if ((step < 1e-12) && (step > -1e-12)) {
            break;
        }
    }

    const double alpha = (lambda * lambda) - (sigma * sigma) + kappa;
    const double beta = lambda - sigma;
    const double gamma = ((lambda + sigma) * alpha) - det_s;

    double x_vector[3U];
    for (std::size_t row = 0U; row < 3U; ++row) {
        x_vector[row] = (alpha * z_vector[row]) + (beta * s_z[row]) + s2_z[row];
    }

    const double norm_squared = (gamma * gamma) +
                                (x_vector[0U] * x_vector[0U]) +
                                (x_vector[1U] * x_vector[1U]) +
                                (x_vector[2U] * x_vector[2U]);
    if (norm_squared < 1e-18) {
        return false;  // degenerate: collinear observations or no information
    }

    const double inverse_norm = 1.0 / std::sqrt(norm_squared);
    quaternion[0U] = static_cast<float>(x_vector[0U] * inverse_norm);
    quaternion[1U] = static_cast<float>(x_vector[1U] * inverse_norm);
    quaternion[2U] = static_cast<float>(x_vector[2U] * inverse_norm);
    quaternion[3U] = static_cast<float>(gamma * inverse_norm);
    return true;
}

AttitudeSolver::Solution AttitudeSolver::solve(
    const Centroiding::Result& centroids,
    const LisGridMatcher::MatchResult& matches,
    const LisGridMatcher::Catalog& catalog,
    std::size_t catalog_count) {
    Solution solution;

    std::array<std::array<float, 3U>, kMaxStars> body;
    std::array<std::array<float, 3U>, kMaxStars> reference;
    std::size_t paired = 0U;

    const std::size_t centroid_count =
        (centroids.count < kMaxStars) ? centroids.count : kMaxStars;

    for (std::size_t index = 0U; index < centroid_count; ++index) {
        const std::uint32_t star_id = matches.star_ids[index];
        if (star_id == 0U) {
            continue;
        }
        // Resolve the catalog entry. The matcher reports catalogue star_id
        // rather than an index, deliberately: indices shift whenever the
        // magnitude cut changes, and a stale index fails silently.
        std::size_t catalog_index = catalog_count;
        for (std::size_t scan = 0U; scan < catalog_count; ++scan) {
            if (catalog[scan].star_id == star_id) {
                catalog_index = scan;
                break;
            }
        }
        if (catalog_index >= catalog_count) {
            continue;
        }

        float bearing[3U];
        pixelToBearing(
            centroids.points[index].x, centroids.points[index].y, bearing);

        body[paired][0U] = bearing[0U];
        body[paired][1U] = bearing[1U];
        body[paired][2U] = bearing[2U];
        reference[paired][0U] = catalog[catalog_index].x;
        reference[paired][1U] = catalog[catalog_index].y;
        reference[paired][2U] = catalog[catalog_index].z;
        ++paired;
        if (paired >= kMaxStars) {
            break;
        }
    }

    solution.matched_count = paired;

    if (paired < config_.minimum_matched_stars) {
        solution.reason = kTooFewStars;
        return solution;
    }

    if (!questSolve(body, reference, paired, solution.quaternion)) {
        solution.reason = kDegenerateGeometry;
        return solution;
    }

    // Integrity statistic: rotate each catalogue vector into the body frame and
    // measure the angle to the observed bearing. A quad matched to the wrong
    // catalogue stars fits its own points and nothing else, so this is the
    // posterior evidence that an analytic mismatch prior never examines.
    float matrix[3U][3U];
    quaternionToMatrix(solution.quaternion, matrix);

    double sum_squared = 0.0;
    float worst = 0.0F;
    for (std::size_t index = 0U; index < paired; ++index) {
        float predicted[3U];
        for (std::size_t row = 0U; row < 3U; ++row) {
            predicted[row] = (matrix[row][0U] * reference[index][0U]) +
                             (matrix[row][1U] * reference[index][1U]) +
                             (matrix[row][2U] * reference[index][2U]);
        }
        // Renormalise: catalogue entries are only unit to float precision, and
        // the rotation adds a little more.
        const float inverse_norm =
            1.0F / std::sqrt((predicted[0U] * predicted[0U]) +
                             (predicted[1U] * predicted[1U]) +
                             (predicted[2U] * predicted[2U]));
        predicted[0U] *= inverse_norm;
        predicted[1U] *= inverse_norm;
        predicted[2U] *= inverse_norm;

        // Measure via chord length, NOT acos(dot). Near zero separation acos is
        // catastrophically ill-conditioned: d(acos)/d(x) diverges as x -> 1, so
        // float rounding alone yields sqrt(2*eps) ~ 65 arcsec of phantom
        // residual -- enough to fail the integrity gate on a perfect solve.
        // The chord 2*sin(theta/2) is well conditioned throughout.
        const float dx = predicted[0U] - body[index][0U];
        const float dy = predicted[1U] - body[index][1U];
        const float dz = predicted[2U] - body[index][2U];
        const float chord = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const float angle =
            2.0F * std::asin(clampUnit(0.5F * chord)) * kRadiansToArcseconds;
        sum_squared += static_cast<double>(angle) * static_cast<double>(angle);
        if (angle > worst) {
            worst = angle;
        }
    }

    solution.residual_rms_arcsec = static_cast<float>(
        std::sqrt(sum_squared / static_cast<double>(paired)));
    solution.residual_max_arcsec = worst;

    if (solution.residual_rms_arcsec > config_.max_residual_arcsec) {
        solution.reason = kResidualTooLarge;
        return solution;  // refuse rather than emit a confident wrong attitude
    }

    // Boresight is the body x-axis expressed in the inertial frame, i.e. the
    // first row of the transpose: column 0 of the body-from-inertial matrix.
    const float boresight_x = matrix[0U][0U];
    const float boresight_y = matrix[0U][1U];
    const float boresight_z = matrix[0U][2U];

    float right_ascension = std::atan2(boresight_y, boresight_x) * kRadiansToDegrees;
    if (right_ascension < 0.0F) {
        right_ascension += 360.0F;
    }
    solution.right_ascension_deg = right_ascension;
    solution.declination_deg =
        std::asin(clampUnit(boresight_z)) * kRadiansToDegrees;

    float roll = std::atan2(matrix[1U][2U], matrix[2U][2U]) * kRadiansToDegrees;
    if (roll < 0.0F) {
        roll += 360.0F;
    }
    solution.roll_deg = roll;

    solution.valid = true;
    solution.reason = kAccepted;
    return solution;
}

}  // namespace star_tracker

static_assert(
    sizeof(star_tracker::AttitudeSolver) < 1024U,
    "AttitudeSolver state must stay small; working arrays are stack-local");
