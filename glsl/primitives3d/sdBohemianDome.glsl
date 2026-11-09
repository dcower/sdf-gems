// ============================================================================
// sdBohemianDome, sdgBohemianDome — Exact 3D Euclidean Bohemian Dome & Vaults
// @author Gemini 4 Argon (research & mathematical design), Gemini 3.8 Flash (experiments & implementation)
// @exact  Yes (vectorized vec4 Circular Cross-Section Envelope-Theorem solver)
// ============================================================================

#ifndef SDF_GEMS_HAS_SD_BOHEMIAN_DOME
#define SDF_GEMS_HAS_SD_BOHEMIAN_DOME

#include "../common/common.glsl"

// ============================================================================
// Exact 3D Euclidean Signed Distance of the Quartic Bohemian Dome
// Surface S(u, v) = (a*cos(u), b*sin(u) + c*cos(v), c*sin(v))
// for semi-axes (a, b), circle radius c, and dilation / shell thickness th.
// ============================================================================
float sdBohemianDome(vec3 p, float a, float b, float c, float th) {
    float x = abs(p.x);
    float y = p.y;
    float z = abs(p.z);

    // 4 vectorized seeds: negative flank, center, positive flank, and analytical elevation guess
    float u_guess = asin(clamp(y / max(b + c, 1e-6), -1.0, 1.0));
    vec4 u = vec4(-1.25, 0.0, 1.25, u_guess);

    for (int i = 0; i < 6; ++i) {
        vec4 cu = cos(u);
        vec4 su = sin(u);
        vec4 dx = vec4(x) - a * cu;
        vec4 dy = vec4(y) - b * su;
        vec4 rho2 = dy * dy + vec4(z * z);
        vec4 rho = sqrt(rho2);
        vec4 rho_safe = max(rho, vec4(1e-7));
        vec4 dr = rho - vec4(c);
        vec4 ratio_dy = dy / rho_safe;
        vec4 z2_over_rho3 = vec4(z * z) / (rho_safe * rho2 + vec4(1e-30));

        vec4 g = (a * dx) * su - (b * dr) * ratio_dy * cu;
        vec4 gp = (a * dx) * cu + (a * a) * (su * su)
                + (b * dr) * ratio_dy * su
                + (b * b) * (cu * cu) * (vec4(1.0) - c * z2_over_rho3);

        vec4 step_newton = clamp(g / max(gp, vec4(1e-4)), -0.4, 0.4);
        vec4 step_grad = sign(g) * 0.15;
        vec4 cond = step(vec4(1e-4), gp);
        vec4 step_u = mix(step_grad, step_newton, cond);
        u = clamp(u - step_u, -0.5 * SDF_PI, 0.5 * SDF_PI);
    }

    // Evaluate distance squared for all 4 candidate roots
    vec4 cu = cos(u);
    vec4 su = sin(u);
    vec4 dx = vec4(x) - a * cu;
    vec4 dy = vec4(y) - b * su;
    vec4 rho = sqrt(dy * dy + vec4(z * z));
    vec4 dr = rho - vec4(c);
    vec4 d2_vec = dx * dx + dr * dr;

    // Minimum over the 4 candidate roots
    float min_d2 = min(min(d2_vec.x, d2_vec.y), min(d2_vec.z, d2_vec.w));

    // Also check boundary endpoints u = -pi/2 and u = +pi/2
    float dy_m = y + b;
    float dr_m = sqrt(dy_m * dy_m + z * z) - c;
    float d2_m = x * x + dr_m * dr_m;

    float dy_p = y - b;
    float dr_p = sqrt(dy_p * dy_p + z * z) - c;
    float d2_p = x * x + dr_p * dr_p;

    min_d2 = min(min_d2, min(d2_m, d2_p));
    return sqrt(max(min_d2, 0.0)) - th;
}

// ============================================================================
// Exact 3D Euclidean SDF & Analytical Unit Gradient of the Quartic Bohemian Dome
// Returns vec4(d - th, grad.xyz) with ||grad|| == 1.0 identically.
// ============================================================================
vec4 sdgBohemianDome(vec3 p, float a, float b, float c, float th) {
    float x = abs(p.x);
    float y = p.y;
    float z = abs(p.z);

    // 4 vectorized seeds: negative flank, center, positive flank, and analytical elevation guess
    float u_guess = asin(clamp(y / max(b + c, 1e-6), -1.0, 1.0));
    vec4 u = vec4(-1.25, 0.0, 1.25, u_guess);

    for (int i = 0; i < 6; ++i) {
        vec4 cu = cos(u);
        vec4 su = sin(u);
        vec4 dx = vec4(x) - a * cu;
        vec4 dy = vec4(y) - b * su;
        vec4 rho2 = dy * dy + vec4(z * z);
        vec4 rho = sqrt(rho2);
        vec4 rho_safe = max(rho, vec4(1e-7));
        vec4 dr = rho - vec4(c);
        vec4 ratio_dy = dy / rho_safe;
        vec4 z2_over_rho3 = vec4(z * z) / (rho_safe * rho2 + vec4(1e-30));

        vec4 g = (a * dx) * su - (b * dr) * ratio_dy * cu;
        vec4 gp = (a * dx) * cu + (a * a) * (su * su)
                + (b * dr) * ratio_dy * su
                + (b * b) * (cu * cu) * (vec4(1.0) - c * z2_over_rho3);

        vec4 step_newton = clamp(g / max(gp, vec4(1e-4)), -0.4, 0.4);
        vec4 step_grad = sign(g) * 0.15;
        vec4 cond = step(vec4(1e-4), gp);
        vec4 step_u = mix(step_grad, step_newton, cond);
        u = clamp(u - step_u, -0.5 * SDF_PI, 0.5 * SDF_PI);
    }

    // Evaluate distance squared for all 4 candidate roots
    vec4 cu = cos(u);
    vec4 su = sin(u);
    vec4 dx = vec4(x) - a * cu;
    vec4 dy = vec4(y) - b * su;
    vec4 rho = sqrt(dy * dy + vec4(z * z));
    vec4 dr = rho - vec4(c);
    vec4 d2_vec = dx * dx + dr * dr;

    // Check boundary endpoints u = -pi/2 and u = +pi/2
    float dy_m = y + b;
    float dr_m = sqrt(dy_m * dy_m + z * z) - c;
    float d2_m = x * x + dr_m * dr_m;

    float dy_p = y - b;
    float dr_p = sqrt(dy_p * dy_p + z * z) - c;
    float d2_p = x * x + dr_p * dr_p;

    float best_d2 = d2_vec.x;
    float best_u = u.x;
    if (d2_vec.y < best_d2) { best_d2 = d2_vec.y; best_u = u.y; }
    if (d2_vec.z < best_d2) { best_d2 = d2_vec.z; best_u = u.z; }
    if (d2_vec.w < best_d2) { best_d2 = d2_vec.w; best_u = u.w; }
    if (d2_m < best_d2) { best_d2 = d2_m; best_u = -0.5 * SDF_PI; }
    if (d2_p < best_d2) { best_d2 = d2_p; best_u = 0.5 * SDF_PI; }

    float opt_cu = cos(best_u);
    float opt_su = sin(best_u);
    float opt_dx = x - a * opt_cu;
    float opt_dy = y - b * opt_su;
    float opt_rho = sqrt(opt_dy * opt_dy + z * z);
    float opt_dr = opt_rho - c;
    float d = sqrt(max(opt_dx * opt_dx + opt_dr * opt_dr, 0.0));

    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sz = (p.z < 0.0) ? -1.0 : 1.0;

    vec3 grad;
    if (d < 1e-7) {
        float inv_rho = (opt_rho > 1e-7) ? (1.0 / opt_rho) : 0.0;
        vec3 n = vec3(0.0, opt_dy * inv_rho, sz * z * inv_rho);
        float len_n = length(n);
        grad = (len_n > 1e-7) ? (n / len_n) : vec3(1.0, 0.0, 0.0);
    } else {
        float inv_d = 1.0 / d;
        float dr_over_rho = (opt_rho > 1e-7) ? (opt_dr / opt_rho) : 0.0;
        grad = vec3(sx * opt_dx * inv_d,
                    opt_dy * dr_over_rho * inv_d,
                    sz * z * dr_over_rho * inv_d);
        float len_g = length(grad);
        grad = (len_g > 1e-7) ? (grad / len_g) : vec3(1.0, 0.0, 0.0);
    }

    return vec4(d - th, grad);
}

#endif // SDF_GEMS_HAS_SD_BOHEMIAN_DOME
