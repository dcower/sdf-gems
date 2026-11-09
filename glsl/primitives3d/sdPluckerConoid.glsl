// ============================================================================
// sdRuledConoid, sdPluckerConoid — Exact 3D Euclidean Ruled Conoids & Plücker's Conoid
// @author Gemini 4 Argon (research & mathematical design), Gemini 3.8 Flash (experiments & implementation)
// @exact  Yes (vectorized vec4 Envelope-Theorem ruling solver with exact spine distance)
// ============================================================================

#ifndef SDF_GEMS_HAS_SD_PLUCKER_CONOID
#define SDF_GEMS_HAS_SD_PLUCKER_CONOID

#include "../common/common.glsl"

// ============================================================================
// Exact 3D Euclidean SDF of the n-Fold Ruled Conoid
// Parameterization: S(r, theta) = (r*cos(theta), r*sin(theta), h*sin(n*theta))
// for radial extension r in [0, R], height amplitude h, thickness th, and fold count nFolds.
// ============================================================================
float sdRuledConoid(vec3 p, float R, float h, float th, float nFolds) {
    float rxy = length(p.xy);
    float z_sp = clamp(p.z, -h, h);
    if (rxy < 1e-7) {
        return abs(p.z - z_sp) - th;
    }

    vec3 diff_sp = vec3(p.xy, p.z - z_sp);
    float best_d2 = dot(diff_sp, diff_sp);

    float phi = atan(p.y, p.x);
    float period = SDF_TAU / nFolds;
    float half_period = SDF_PI / nFolds;
    float quarter_period = 0.5 * half_period;
    float crest = quarter_period + period * round((phi - quarter_period) / period);
    float trough = -quarter_period + period * round((phi + quarter_period) / period);
    float node = half_period * round(phi / half_period);
    vec4 theta = vec4(phi, crest, trough, node);

    for (int i = 0; i < 6; ++i) {
        vec4 ct = cos(theta);
        vec4 st = sin(theta);
        vec4 ntheta = nFolds * theta;
        vec4 cnt = cos(ntheta);
        vec4 snt = sin(ntheta);
        vec4 w = p.x * ct + p.y * st;
        vec4 q = -p.x * st + p.y * ct;
        vec4 ustar = clamp(w, 0.0, R);
        vec4 dz = p.z - h * snt;
        vec4 g = -ustar * q - (nFolds * h) * dz * cnt;
        vec4 mask = step(vec4(0.0), w) * step(w, vec4(R));
        vec4 gp = ustar * w - mask * (q * q) + (nFolds * nFolds * h) * (h * (cnt * cnt - snt * snt) + p.z * snt);
        vec4 step_newton = clamp(g / max(gp, vec4(1e-4)), -0.35, 0.35);
        vec4 step_grad = sign(g) * 0.15;
        vec4 cond = step(vec4(1e-4), gp);
        theta -= mix(step_grad, step_newton, cond);
    }

    vec4 ct = cos(theta);
    vec4 st = sin(theta);
    vec4 snt = sin(nFolds * theta);
    vec4 w = p.x * ct + p.y * st;
    vec4 ustar = clamp(w, 0.0, R);

    vec3 q0 = vec3(ustar.x * ct.x, ustar.x * st.x, h * snt.x);
    vec3 q1 = vec3(ustar.y * ct.y, ustar.y * st.y, h * snt.y);
    vec3 q2 = vec3(ustar.z * ct.z, ustar.z * st.z, h * snt.z);
    vec3 q3 = vec3(ustar.w * ct.w, ustar.w * st.w, h * snt.w);

    vec3 d0 = p - q0;
    vec3 d1 = p - q1;
    vec3 d2 = p - q2;
    vec3 d3 = p - q3;

    float min_d2 = min(best_d2, min(min(dot(d0, d0), dot(d1, d1)), min(dot(d2, d2), dot(d3, d3))));
    return sqrt(min_d2) - th;
}

// ============================================================================
// Exact 3D Euclidean SDF & Analytical Unit Gradient of the n-Fold Ruled Conoid
// Returns vec4(d, grad.xyz) with ||grad|| == 1.0 identically.
// ============================================================================
vec4 sdgRuledConoid(vec3 p, float R, float h, float th, float nFolds) {
    float rxy = length(p.xy);
    float z_sp = clamp(p.z, -h, h);
    if (rxy < 1e-7) {
        float d_sp = abs(p.z - z_sp);
        vec3 g_sp = (abs(p.z) > h) ? vec3(0.0, 0.0, sign(p.z)) : vec3(1.0, 0.0, 0.0);
        return vec4(d_sp - th, g_sp);
    }

    vec3 diff_sp = vec3(p.xy, p.z - z_sp);
    float best_d2 = dot(diff_sp, diff_sp);
    vec3 best_diff = diff_sp;
    float best_theta = 0.0;
    float best_u = 0.0;
    bool found_ruling = false;

    float phi = atan(p.y, p.x);
    float period = SDF_TAU / nFolds;
    float half_period = SDF_PI / nFolds;
    float quarter_period = 0.5 * half_period;
    float crest = quarter_period + period * round((phi - quarter_period) / period);
    float trough = -quarter_period + period * round((phi + quarter_period) / period);
    float node = half_period * round(phi / half_period);
    vec4 theta = vec4(phi, crest, trough, node);

    for (int i = 0; i < 6; ++i) {
        vec4 ct = cos(theta);
        vec4 st = sin(theta);
        vec4 ntheta = nFolds * theta;
        vec4 cnt = cos(ntheta);
        vec4 snt = sin(ntheta);
        vec4 w = p.x * ct + p.y * st;
        vec4 q = -p.x * st + p.y * ct;
        vec4 ustar = clamp(w, 0.0, R);
        vec4 dz = p.z - h * snt;
        vec4 g = -ustar * q - (nFolds * h) * dz * cnt;
        vec4 mask = step(vec4(0.0), w) * step(w, vec4(R));
        vec4 gp = ustar * w - mask * (q * q) + (nFolds * nFolds * h) * (h * (cnt * cnt - snt * snt) + p.z * snt);
        vec4 step_newton = clamp(g / max(gp, vec4(1e-4)), -0.35, 0.35);
        vec4 step_grad = sign(g) * 0.15;
        vec4 cond = step(vec4(1e-4), gp);
        theta -= mix(step_grad, step_newton, cond);
    }

    vec4 ct = cos(theta);
    vec4 st = sin(theta);
    vec4 snt = sin(nFolds * theta);
    vec4 w = p.x * ct + p.y * st;
    vec4 ustar = clamp(w, 0.0, R);

    // Check seed 0
    {
        vec3 q0 = vec3(ustar.x * ct.x, ustar.x * st.x, h * snt.x);
        vec3 diff0 = p - q0;
        float d2_0 = dot(diff0, diff0);
        if (d2_0 < best_d2) {
            best_d2 = d2_0;
            best_diff = diff0;
            best_theta = theta.x;
            best_u = ustar.x;
            found_ruling = true;
        }
    }
    // Check seed 1
    {
        vec3 q1 = vec3(ustar.y * ct.y, ustar.y * st.y, h * snt.y);
        vec3 diff1 = p - q1;
        float d2_1 = dot(diff1, diff1);
        if (d2_1 < best_d2) {
            best_d2 = d2_1;
            best_diff = diff1;
            best_theta = theta.y;
            best_u = ustar.y;
            found_ruling = true;
        }
    }
    // Check seed 2
    {
        vec3 q2 = vec3(ustar.z * ct.z, ustar.z * st.z, h * snt.z);
        vec3 diff2 = p - q2;
        float d2_2 = dot(diff2, diff2);
        if (d2_2 < best_d2) {
            best_d2 = d2_2;
            best_diff = diff2;
            best_theta = theta.z;
            best_u = ustar.z;
            found_ruling = true;
        }
    }
    // Check seed 3
    {
        vec3 q3 = vec3(ustar.w * ct.w, ustar.w * st.w, h * snt.w);
        vec3 diff3 = p - q3;
        float d2_3 = dot(diff3, diff3);
        if (d2_3 < best_d2) {
            best_d2 = d2_3;
            best_diff = diff3;
            best_theta = theta.w;
            best_u = ustar.w;
            found_ruling = true;
        }
    }

    float d = sqrt(best_d2);
    vec3 grad;
    if (d > 1e-7) {
        grad = best_diff / d;
    } else if (found_ruling) {
        float cth = cos(best_theta);
        float sth = sin(best_theta);
        float cnth = cos(nFolds * best_theta);
        vec3 nrm = vec3((nFolds * h) * sth * cnth, -(nFolds * h) * cth * cnth, best_u);
        float nl = length(nrm);
        grad = (nl > 1e-7) ? (nrm / nl) : vec3(1.0, 0.0, 0.0);
    } else {
        float nl_sp = length(best_diff);
        grad = (nl_sp > 1e-7) ? (best_diff / nl_sp) : vec3(1.0, 0.0, 0.0);
    }

    return vec4(d - th, grad);
}

// ============================================================================
// Convenience Wrappers: Julius Plücker's Canonical 2-Fold Ruled Conoid (1846)
// Cartesian form: z * (x^2 + y^2) - 2 * h * x * y = 0
// ============================================================================
float sdPluckerConoid(vec3 p, float R, float h, float th) {
    return sdRuledConoid(p, R, h, th, 2.0);
}

vec4 sdgPluckerConoid(vec3 p, float R, float h, float th) {
    return sdgRuledConoid(p, R, h, th, 2.0);
}

#endif // SDF_GEMS_HAS_SD_PLUCKER_CONOID
