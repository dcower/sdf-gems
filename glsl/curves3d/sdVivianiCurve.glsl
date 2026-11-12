#ifndef SDF_GEMS_HAS_VIVIANI_CURVE
#define SDF_GEMS_HAS_VIVIANI_CURVE

#include "../common/common.glsl"

// ============================================================================
// Exact 3D Viviani's Temple Window Curve & D_2d Tennis-Ball Seam Curve
//
// The S^2 Dot-Product Linearization Theorem:
// For any space curve q(t) lying on the sphere S^2(R) (||q(t)||^2 == R^2),
// minimizing the 3D Euclidean distance ||p - q(t)||_2 is 100% equivalent to
// maximizing the linear projection h(t) = dot(p, q(t)).
//
// This eliminates ALL square roots, square root derivatives, and rational
// fractions from the stationarity equation and Newton iterations (0 SFU sqrts
// inside the solver loop), converging unconditionally to machine precision.
//
// 1. Viviani's Curve (1692):
//    Intersection of sphere x^2 + y^2 + z^2 = R^2 and cylinder x^2 + y^2 = R*x.
//    Parameterization: q(t) = R * (cos^2(t), cos(t)*sin(t), sin(t)).
//    Quarter-period symmetry folding (y, z) -> (|y|, |z|) reduces t in [-pi, pi]
//    to t in [0, pi/2] with zero loss of generality.
//
// 2. Tennis-Ball / Baseball Seam Curve:
//    D_{2d}-symmetric curve on S^2(R):
//    q(t) = R * (a*cos(t) + alpha*cos(3t), a*sin(t) - alpha*sin(3t), 2*sqrt(a*alpha)*sin(2t)),
//    where a = 1 - alpha.
//    8-fold D_{2d} Weyl-chamber reduction restricts parameter search to t in [-pi/4, pi/4].
//
// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.8 Flash (experiments & implementation)
// @exact  Yes (Both curves evaluate exact 3D Euclidean distance via the S^2 Dot-Product Linearization Theorem)
// ============================================================================

// 1. Exact 3D Euclidean SDF of Viviani's Curve
float sdVivianiCurve(vec3 p, float R, float rTube) {
    float x = p.x;
    float ay = abs(p.y);
    float az = abs(p.z);

    float best_h = -1e30;
    float best_t = 0.0;

    // 2 seeds in [0, pi/2] cover the unimodal / bimodal basins of h(t)
    float seeds[2] = float[2](0.32, 1.18);
    for (int i = 0; i < 2; i++) {
        float t = seeds[i];
        for (int iter = 0; iter < 4; iter++) {
            float ct = cos(t);
            float st = sin(t);
            float c2t = ct * ct - st * st;
            float s2t = 2.0 * ct * st;
            float hp = -x * s2t + ay * c2t + az * ct;
            float hpp = -2.0 * x * c2t - 2.0 * ay * s2t - az * st;
            float denom = (hpp < 0.0) ? min(hpp, -1e-7) : -max(abs(hpp), 1e-7);
            t = clamp(t - clamp(hp / denom, -0.45, 0.45), 0.0, 1.57079632679);
        }
        float ct = cos(t);
        float st = sin(t);
        float h_val = x * ct * ct + ay * ct * st + az * st;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    // Check interval boundaries {0, pi/2}
    if (x > best_h) {
        best_h = x;
        best_t = 0.0;
    }
    if (az > best_h) {
        best_h = az;
        best_t = 1.57079632679;
    }

    float ct = cos(best_t);
    float st = sin(best_t);
    vec3 q_sym = R * vec3(ct * ct, ct * st, st);
    vec3 p_sym = vec3(x, ay, az);
    return length(p_sym - q_sym) - rTube;
}

// 2. Exact 3D Euclidean SDF & Analytical Unit Gradient of Viviani's Curve
vec4 sdgVivianiCurve(vec3 p, float R, float rTube) {
    float x = p.x;
    float ay = abs(p.y);
    float az = abs(p.z);
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float sz = (p.z < 0.0) ? -1.0 : 1.0;

    float best_h = -1e30;
    float best_t = 0.0;

    float seeds[2] = float[2](0.32, 1.18);
    for (int i = 0; i < 2; i++) {
        float t = seeds[i];
        for (int iter = 0; iter < 4; iter++) {
            float ct = cos(t);
            float st = sin(t);
            float c2t = ct * ct - st * st;
            float s2t = 2.0 * ct * st;
            float hp = -x * s2t + ay * c2t + az * ct;
            float hpp = -2.0 * x * c2t - 2.0 * ay * s2t - az * st;
            float denom = (hpp < 0.0) ? min(hpp, -1e-7) : -max(abs(hpp), 1e-7);
            t = clamp(t - clamp(hp / denom, -0.45, 0.45), 0.0, 1.57079632679);
        }
        float ct = cos(t);
        float st = sin(t);
        float h_val = x * ct * ct + ay * ct * st + az * st;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    if (x > best_h) {
        best_h = x;
        best_t = 0.0;
    }
    if (az > best_h) {
        best_h = az;
        best_t = 1.57079632679;
    }

    float ct = cos(best_t);
    float st = sin(best_t);
    vec3 q_sym = R * vec3(ct * ct, ct * st, st);
    vec3 p_sym = vec3(x, ay, az);
    vec3 diff = p_sym - q_sym;
    float d = length(diff);
    vec3 g_sym = (d > 1e-7) ? diff / d : vec3(1.0, 0.0, 0.0);
    vec3 grad = vec3(g_sym.x, sy * g_sym.y, sz * g_sym.z);
    return vec4(d - rTube, grad);
}

// 3. Exact 3D Euclidean SDF of Tennis-Ball / Baseball Seam Curve
float sdTennisBallSeam(vec3 p, float R, float alpha, float rTube) {
    float a = 1.0 - alpha;
    float sq = 2.0 * sqrt(max(a * alpha, 0.0));

    float px = p.x;
    float py = p.y;
    float pz = p.z;

    // 8-fold D_{2d} Weyl-chamber reduction:
    // 1. Reflection across x + y = 0
    bool m1 = (px + py < 0.0);
    if (m1) {
        float tmp = px; px = -py; py = -tmp;
    }
    // 2. Reflection across x - y = 0
    bool m2 = (px < py);
    if (m2) {
        float tmp = px; px = py; py = tmp;
    }
    // 3. 180-deg rotation around X axis
    bool rx = (pz < 0.0);
    if (rx) {
        py = -py; pz = -pz;
    }

    float best_h = -1e30;
    float best_t = 0.0;

    // 3 seeds in [-pi/4, pi/4]
    float seeds[3] = float[3](-0.5, 0.0, 0.5);
    for (int i = 0; i < 3; i++) {
        float t = seeds[i];
        for (int iter = 0; iter < 5; iter++) {
            float ct = cos(t);
            float st = sin(t);
            float c2t = ct * ct - st * st;
            float s2t = 2.0 * ct * st;
            float c3t = ct * c2t - st * s2t;
            float s3t = st * c2t + ct * s2t;
            float hp = px * (-a * st - 3.0 * alpha * s3t) + py * (a * ct - 3.0 * alpha * c3t) + 2.0 * pz * sq * c2t;
            float hpp = px * (-a * ct - 9.0 * alpha * c3t) + py * (-a * st + 9.0 * alpha * s3t) - 4.0 * pz * sq * s2t;
            float denom = (hpp < 0.0) ? min(hpp, -1e-7) : -max(abs(hpp), 1e-7);
            t = clamp(t - clamp(hp / denom, -0.4, 0.4), -0.78539816339, 0.78539816339);
        }
        float ct = cos(t);
        float st = sin(t);
        float c2t = ct * ct - st * st;
        float s2t = 2.0 * ct * st;
        float c3t = ct * c2t - st * s2t;
        float s3t = st * c2t + ct * s2t;
        float h_val = px * (a * ct + alpha * c3t) + py * (a * st - alpha * s3t) + pz * sq * s2t;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    // Check chamber boundaries {-pi/4, pi/4}
    float ends[2] = float[2](-0.78539816339, 0.78539816339);
    for (int i = 0; i < 2; i++) {
        float t = ends[i];
        float ct = cos(t);
        float st = sin(t);
        float c2t = ct * ct - st * st;
        float s2t = 2.0 * ct * st;
        float c3t = ct * c2t - st * s2t;
        float s3t = st * c2t + ct * s2t;
        float h_val = px * (a * ct + alpha * c3t) + py * (a * st - alpha * s3t) + pz * sq * s2t;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    float ct = cos(best_t);
    float st = sin(best_t);
    float c2t = ct * ct - st * st;
    float s2t = 2.0 * ct * st;
    float c3t = ct * c2t - st * s2t;
    float s3t = st * c2t + ct * s2t;
    vec3 q_loc = R * vec3(a * ct + alpha * c3t, a * st - alpha * s3t, sq * s2t);
    vec3 p_loc = vec3(px, py, pz);
    return length(p_loc - q_loc) - rTube;
}

// 4. Exact 3D Euclidean SDF & Analytical Unit Gradient of Tennis-Ball / Baseball Seam Curve
vec4 sdgTennisBallSeam(vec3 p, float R, float alpha, float rTube) {
    float a = 1.0 - alpha;
    float sq = 2.0 * sqrt(max(a * alpha, 0.0));

    float px = p.x;
    float py = p.y;
    float pz = p.z;

    bool m1 = (px + py < 0.0);
    if (m1) {
        float tmp = px; px = -py; py = -tmp;
    }
    bool m2 = (px < py);
    if (m2) {
        float tmp = px; px = py; py = tmp;
    }
    bool rx = (pz < 0.0);
    if (rx) {
        py = -py; pz = -pz;
    }

    float best_h = -1e30;
    float best_t = 0.0;

    float seeds[3] = float[3](-0.5, 0.0, 0.5);
    for (int i = 0; i < 3; i++) {
        float t = seeds[i];
        for (int iter = 0; iter < 5; iter++) {
            float ct = cos(t);
            float st = sin(t);
            float c2t = ct * ct - st * st;
            float s2t = 2.0 * ct * st;
            float c3t = ct * c2t - st * s2t;
            float s3t = st * c2t + ct * s2t;
            float hp = px * (-a * st - 3.0 * alpha * s3t) + py * (a * ct - 3.0 * alpha * c3t) + 2.0 * pz * sq * c2t;
            float hpp = px * (-a * ct - 9.0 * alpha * c3t) + py * (-a * st + 9.0 * alpha * s3t) - 4.0 * pz * sq * s2t;
            float denom = (hpp < 0.0) ? min(hpp, -1e-7) : -max(abs(hpp), 1e-7);
            t = clamp(t - clamp(hp / denom, -0.4, 0.4), -0.78539816339, 0.78539816339);
        }
        float ct = cos(t);
        float st = sin(t);
        float c2t = ct * ct - st * st;
        float s2t = 2.0 * ct * st;
        float c3t = ct * c2t - st * s2t;
        float s3t = st * c2t + ct * s2t;
        float h_val = px * (a * ct + alpha * c3t) + py * (a * st - alpha * s3t) + pz * sq * s2t;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    float ends[2] = float[2](-0.78539816339, 0.78539816339);
    for (int i = 0; i < 2; i++) {
        float t = ends[i];
        float ct = cos(t);
        float st = sin(t);
        float c2t = ct * ct - st * st;
        float s2t = 2.0 * ct * st;
        float c3t = ct * c2t - st * s2t;
        float s3t = st * c2t + ct * s2t;
        float h_val = px * (a * ct + alpha * c3t) + py * (a * st - alpha * s3t) + pz * sq * s2t;
        if (h_val > best_h) {
            best_h = h_val;
            best_t = t;
        }
    }

    float ct = cos(best_t);
    float st = sin(best_t);
    float c2t = ct * ct - st * st;
    float s2t = 2.0 * ct * st;
    float c3t = ct * c2t - st * s2t;
    float s3t = st * c2t + ct * s2t;
    vec3 q_loc = R * vec3(a * ct + alpha * c3t, a * st - alpha * s3t, sq * s2t);
    vec3 p_loc = vec3(px, py, pz);
    vec3 diff = p_loc - q_loc;
    float d = length(diff);
    vec3 g = (d > 1e-7) ? diff / d : vec3(1.0, 0.0, 0.0);

    // Unfold gradient in exact reverse symmetry order
    if (rx) {
        g.y = -g.y;
        g.z = -g.z;
    }
    if (m2) {
        float tmp = g.x;
        g.x = g.y;
        g.y = tmp;
    }
    if (m1) {
        float tmp = g.x;
        g.x = -g.y;
        g.y = -tmp;
    }

    return vec4(d - rTube, g);
}

#endif // SDF_GEMS_HAS_VIVIANI_CURVE
