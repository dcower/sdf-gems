/**
 * Exact Catenoid, Ruled Helicoid, Hyperbolic-Secant Scherk Saddle Tower,
 * and Exact Bounded Enneper Minimal Surface
 * @author Gemini 3.1 Pro (experiments & implementation), Gemini 4 Argon (research & mathematical design),
 *         and Gemini 4 Carbon (Scherk asymptotic wing bound, unit-gradient fix & exact Enneper minimal surface, 2026)
 */

#ifndef SDF_GEMS_HAS_MINIMAL_SURFACE
#define SDF_GEMS_HAS_MINIMAL_SURFACE

#include "../common/common.glsl"

// --- 1. Exact Bounded Catenoid ---

float sdCatenoid(vec3 p, float c, float h, float th) {
    float rho = length(p.xy);
    float z = abs(p.z);

    float c_safe = max(c, 1e-6);
    float u = clamp(asinh(max(rho, z) / c_safe), 0.0, h / c_safe);

    for (int i = 0; i < 4; ++i) {
        float sinh_u = sinh(u);
        float cosh_u = cosh(u);
        float g = c * (u - z / c_safe) + (c * cosh_u - rho) * sinh_u;
        float g_prime = c * (1.0 + cosh(2.0 * u)) - rho * cosh_u;
        u = clamp(u - g / g_prime, 0.0, h / c_safe);
    }

    float u0 = 0.0;
    float uh = h / c_safe;

    float dx_root = c * cosh(u) - rho;
    float dz_root = c * u - z;
    float d2_root = dx_root * dx_root + dz_root * dz_root;

    float dx_0 = c * cosh(u0) - rho;
    float dz_0 = c * u0 - z;
    float d2_0 = dx_0 * dx_0 + dz_0 * dz_0;

    float dx_h = c * cosh(uh) - rho;
    float dz_h = c * uh - z;
    float d2_h = dx_h * dx_h + dz_h * dz_h;

    float d2 = min(min(d2_root, d2_0), d2_h);
    return sqrt(d2) - th;
}

vec4 sdgCatenoid(vec3 p, float c, float h, float th) {
    float rho = length(p.xy);
    float z = abs(p.z);

    float c_safe = max(c, 1e-6);
    float u = clamp(asinh(max(rho, z) / c_safe), 0.0, h / c_safe);

    for (int i = 0; i < 4; ++i) {
        float sinh_u = sinh(u);
        float cosh_u = cosh(u);
        float g = c * (u - z / c_safe) + (c * cosh_u - rho) * sinh_u;
        float g_prime = c * (1.0 + cosh(2.0 * u)) - rho * cosh_u;
        u = clamp(u - g / g_prime, 0.0, h / c_safe);
    }

    vec3 u_c = vec3(u, 0.0, h / c_safe);
    vec3 dx_c = c * cosh(u_c) - vec3(rho);
    vec3 dz_c = c * u_c - vec3(z);
    vec3 d2_c = dx_c * dx_c + dz_c * dz_c;

    float d2 = d2_c.x;
    float u_best = u_c.x;
    if (d2_c.y < d2) { d2 = d2_c.y; u_best = u_c.y; }
    if (d2_c.z < d2) { d2 = d2_c.z; u_best = u_c.z; }

    float drho = rho - c * cosh(u_best);
    float dz = z - c * u_best;
    float d = sqrt(d2);

    vec2 n_rz;
    if (d > 1e-7) {
        n_rz = vec2(drho, dz) / d;
    } else {
        n_rz = vec2(1.0, -sinh(u_best)) / cosh(u_best);
    }

    vec3 grad = vec3(
        (rho > 1e-7 ? p.x / rho : 1.0) * n_rz.x,
        (rho > 1e-7 ? p.y / rho : 0.0) * n_rz.x,
        (p.z >= 0.0 ? 1.0 : -1.0) * n_rz.y
    );

    return vec4(d - th, grad);
}

// --- 2. Exact Ruled Helicoid ---

float sdHelicoid(vec3 p, float c, float R, float th) {
    float c_safe = max(c, 1e-6);
    float v0 = p.z / c_safe;
    float v_atan = atan(p.y, p.x);
    float pi = 3.14159265359;

    float v_seed1 = v_atan + pi * round((v0 - v_atan) / pi);
    float v_seed2 = v0;
    float v_seed3 = v_seed1 + pi * (v0 > v_seed1 ? 1.0 : -1.0);

    vec3 v = vec3(v_seed1, v_seed2, v_seed3);

    for (int i = 0; i < 4; ++i) {
        vec3 sin_v = sin(v);
        vec3 cos_v = cos(v);
        vec3 w = p.x * cos_v + p.y * sin_v;
        vec3 q = -p.x * sin_v + p.y * cos_v;
        vec3 u_star = clamp(w, -R, R);

        vec3 g = -u_star * q + c * c * (v - v0);
        vec3 du_dv = vec3(
            abs(w.x) < R ? -q.x : 0.0,
            abs(w.y) < R ? -q.y : 0.0,
            abs(w.z) < R ? -q.z : 0.0
        );

        vec3 g_prime = -du_dv * q + u_star * w + c * c;
        v -= g / g_prime;
        v = clamp(v, v0 - pi, v0 + pi);
    }

    vec3 sin_v = sin(v);
    vec3 cos_v = cos(v);
    vec3 w = p.x * cos_v + p.y * sin_v;
    vec3 u_star = clamp(w, -R, R);

    vec3 dx = vec3(p.x) - u_star * cos_v;
    vec3 dy = vec3(p.y) - u_star * sin_v;
    vec3 dz = c * (vec3(v0) - v);

    vec3 d2 = dx*dx + dy*dy + dz*dz;
    return sqrt(min(min(d2.x, d2.y), d2.z)) - th;
}

vec4 sdgHelicoid(vec3 p, float c, float R, float th) {
    float c_safe = max(c, 1e-6);
    float v0 = p.z / c_safe;
    float v_atan = atan(p.y, p.x);
    float pi = 3.14159265359;

    float v_seed1 = v_atan + pi * round((v0 - v_atan) / pi);
    float v_seed2 = v0;
    float v_seed3 = v_seed1 + pi * (v0 > v_seed1 ? 1.0 : -1.0);

    vec3 v = vec3(v_seed1, v_seed2, v_seed3);

    for (int i = 0; i < 4; ++i) {
        vec3 sin_v = sin(v);
        vec3 cos_v = cos(v);
        vec3 w = p.x * cos_v + p.y * sin_v;
        vec3 q = -p.x * sin_v + p.y * cos_v;
        vec3 u_star = clamp(w, -R, R);

        vec3 g = -u_star * q + c * c * (v - v0);
        vec3 du_dv = vec3(
            abs(w.x) < R ? -q.x : 0.0,
            abs(w.y) < R ? -q.y : 0.0,
            abs(w.z) < R ? -q.z : 0.0
        );

        vec3 g_prime = -du_dv * q + u_star * w + c * c;
        v -= g / g_prime;
        v = clamp(v, v0 - pi, v0 + pi);
    }

    vec3 sin_v = sin(v);
    vec3 cos_v = cos(v);
    vec3 w = p.x * cos_v + p.y * sin_v;
    vec3 u_star = clamp(w, -R, R);

    vec3 dx = vec3(p.x) - u_star * cos_v;
    vec3 dy = vec3(p.y) - u_star * sin_v;
    vec3 dz = c * (vec3(v0) - v);

    vec3 d2_c = dx*dx + dy*dy + dz*dz;

    float d2 = d2_c.x;
    float v_best = v.x;
    float u_best = u_star.x;
    float dx_best = dx.x;
    float dy_best = dy.x;
    float dz_best = dz.x;

    if (d2_c.y < d2) {
        d2 = d2_c.y; v_best = v.y; u_best = u_star.y;
        dx_best = dx.y; dy_best = dy.y; dz_best = dz.y;
    }
    if (d2_c.z < d2) {
        d2 = d2_c.z; v_best = v.z; u_best = u_star.z;
        dx_best = dx.z; dy_best = dy.z; dz_best = dz.z;
    }

    float d = sqrt(d2);
    vec3 grad;
    if (d > 1e-7) {
        grad = vec3(dx_best, dy_best, dz_best) / d;
    } else {
        grad = vec3(c * sin(v_best), -c * cos(v_best), u_best);
        grad = normalize(grad);
    }

    return vec4(d - th, grad);
}

// --- 3. Hyperbolic-Secant Scherk Saddle Tower (with Exact Asymptotic Wing Bound) ---
//
// Surface: sin(z) = sinh(x) * sinh(y), existing strictly inside the hyperbolic cross
//          {|sinh(x) * sinh(y)| <= 1} <=> (1 + e^-|x|)(1 + e^-|y|) >= 2.
//
// 1. Interior / Near-Surface Eikonal Bound:
//    g(x,y,z) = sin(z)*sech(x)*sech(y) - tanh(x)*tanh(y) has ||grad g|| == 1 on g = 0
//    and global maximum ||grad g|| = 4*sqrt(6)/9 via Gudermannian substitution, so
//    d1 = (3*sqrt(6)/8) * g is globally 1-Lipschitz across R^3 without overflow.
// 2. Exterior Asymptotic Wing Bound ((1 + e^-|x|)(1 + e^-|y|) < 2):
//    In 45-degree rotated coordinates (delta, sigma) = ((u - v)/2, (u + v)/2), the
//    hyperbolic cross boundary sinh(u)*sinh(v) = 1 collapses to the smooth graph
//    sigma(delta) = asinh(cosh(delta)) with |sigma'(delta)| < 1 everywhere. Combining
//    |d1| with the exact 2D cross distance yields linear O(||p.xy||) sphere-tracing
//    steps in the far field instead of stalling at +-0.919, with 0% overstep.

vec4 sdgScherkTower(vec3 p) {
    const float kLip = 0.91855865; // Exact 3*sqrt(6)/8 = 1 / (4*sqrt(6)/9)
    const float u0   = 0.88137359; // asinh(1) = ln(1 + sqrt(2))

    vec2 a = abs(p.xy);
    vec2 e = exp(-a);
    vec2 e2 = e * e;
    vec2 inv1pe2 = 1.0 / (1.0 + e2);
    vec2 sech_xy = 2.0 * e * inv1pe2;
    vec2 sgn_xy  = vec2(p.x >= 0.0 ? 1.0 : -1.0, p.y >= 0.0 ? 1.0 : -1.0);
    vec2 tanh_xy = sgn_xy * (1.0 - e2) * inv1pe2;

    float sin_z = sin(p.z);
    float g  = sin_z * sech_xy.x * sech_xy.y - tanh_xy.x * tanh_xy.y;
    float d1 = kLip * g;

    if ((1.0 + e.x) * (1.0 + e.y) < 2.0) {
        float U  = max(a.x, a.y);
        float V  = min(a.x, a.y);
        float dq = 0.5 * (U - V);
        float sq = 0.5 * (U + V);

        float pc = 3.0 * (1.0 - 0.70710678 * (sq - u0));
        float qc = 3.0 * dq;
        float q23 = (qc > 1e-12) ? pow(qc * qc, 1.0 / 3.0) : 0.0;
        float d_seed = sqrt(max(-pc, 0.0) + q23);
        d_seed -= (d_seed * (d_seed * d_seed + pc) - qc) / max(3.0 * d_seed * d_seed + pc, 1e-6);
        float del = mix(min(max(d_seed, 0.0), 0.5 * U), 0.5 * U, smoothstep(2.8, 3.6, U));

        for (int i = 0; i < 3; ++i) {
            float ed  = exp(-del);
            float ed2 = ed * ed;
            float R   = sqrt(1.0 + 6.0 * ed2 + ed2 * ed2);
            float sig = del + log(0.5 * (1.0 + ed2 + R));
            float sp  = (1.0 - ed2) / R;
            float spp = 8.0 * ed2 * (1.0 + ed2) / (R * R * R);
            float psi  = (del - dq) + (sig - sq) * sp;
            float psip = 1.0 + sp * sp + (sig - sq) * spp;
            del = max(0.0, del - psi / max(psip, 1e-6));
        }

        float ed  = exp(-del);
        float ed2 = ed * ed;
        float R   = sqrt(1.0 + 6.0 * ed2 + ed2 * ed2);
        float sig = del + log(0.5 * (1.0 + ed2 + R));

        vec2 dv_arm = vec2(U - (sig + del), V - (sig - del));
        vec2 dv_vtx = vec2(U - u0, V - u0);
        float d2_arm = dot(dv_arm, dv_arm);
        float d2_vtx = dot(dv_vtx, dv_vtx);
        bool use_vtx = (d2_vtx < d2_arm);
        vec2 dv_best = use_vtx ? dv_vtx : dv_arm;
        float d_cross = sqrt(min(d2_arm, d2_vtx));

        if (d_cross > abs(d1)) {
            float s_wing = (d1 >= 0.0) ? 1.0 : -1.0;
            vec2 n_uv;
            if (d_cross > 1e-7) {
                n_uv = dv_best / d_cross;
            } else {
                float sp = use_vtx ? 0.0 : (1.0 - ed2) / R;
                n_uv = vec2(1.0 - sp, 1.0 + sp) * inversesqrt(2.0 * (1.0 + sp * sp));
            }
            vec2 n_ab = (a.x >= a.y) ? n_uv : n_uv.yx;
            return vec4(s_wing * d_cross, s_wing * sgn_xy * n_ab, 0.0);
        }
    }

    vec3 raw_g = vec3(
        -sech_xy.x * (sin_z * tanh_xy.x * sech_xy.y + sech_xy.x * tanh_xy.y),
        -sech_xy.y * (sin_z * sech_xy.x * tanh_xy.y + tanh_xy.x * sech_xy.y),
        cos(p.z) * sech_xy.x * sech_xy.y
    );
    float g2 = dot(raw_g, raw_g);
    vec3 grad = (g2 > 1e-20) ? (raw_g * inversesqrt(g2)) : vec3(0.0, 0.0, 1.0);
    return vec4(d1, grad);
}

float sdScherkTower(vec3 p) {
    return sdgScherkTower(p).x;
}

// --- 4. Exact Bounded Enneper Minimal Surface ---
//
// Alfred Enneper's (1864) classical Weierstrass minimal surface on u^2 + v^2 <= R0^2:
//   X(u, v) = s * (u - u^3/3 + u*v^2,  v - v^3/3 + v*u^2,  u^2 - v^2)
// Exploits order-8 D_2d symmetry, conformal isothermal metric ||X_u||^2 = ||X_v||^2 = (1+r^2)^2,
// X_u . X_v = 0, harmonicity X_vv = -X_uu (yielding a 0-trig, 0-sqrt traceless 2x2 Hessian),
// rational stereographic Gauss map N = (-2u, 2v, 1 - r^2)/(1 + r^2), and 0-trig, 0-sqrt
// Cayley SO(2) rim refinement on u^2 + v^2 = R0^2.

vec2 _sdEnneperSolveSeed(vec3 q, vec2 uv, float R0) {
    float R02 = R0 * R0;
    float r2_0 = dot(uv, uv);
    if (r2_0 > R02) uv *= R0 * inversesqrt(r2_0);

    for (int i = 0; i < 6; ++i) {
        float u = uv.x, v = uv.y;
        float u2 = u * u, v2 = v * v;
        vec3 r = vec3(
            u * (1.0 - (1.0 / 3.0) * u2 + v2) - q.x,
            v * (1.0 - (1.0 / 3.0) * v2 + u2) - q.y,
            (u2 - v2) - q.z
        );

        float g1 = r.x * (1.0 - u2 + v2) + r.y * (2.0 * u * v) + r.z * (2.0 * u);
        float g2 = r.x * (2.0 * u * v) + r.y * (1.0 + u2 - v2) - r.z * (2.0 * v);

        float w  = 1.0 + u2 + v2;
        float w2 = w * w;
        float al = -2.0 * u * r.x + 2.0 * v * r.y + 2.0 * r.z;
        float be =  2.0 * v * r.x + 2.0 * u * r.y;
        float ab2 = al * al + be * be;

        float mu  = (ab2 <= 0.78 * w2 * w2) ? w2 : (0.62 * w2 + 0.62 * ab2 / w2);
        float det = mu * mu - ab2;
        vec2 duv  = vec2((mu - al) * g1 - be * g2, -be * g1 + (mu + al) * g2) / det;

        float lim = 1.0 / max(1.0, max(abs(duv.x), abs(duv.y)) * 1.25);
        uv -= duv * lim;

        float rn2 = dot(uv, uv);
        if (rn2 > R02) uv *= R0 * inversesqrt(rn2);
    }

    if (dot(uv, uv) >= 0.995 * R02) {
        vec2 urv = uv * (R0 * inversesqrt(max(dot(uv, uv), 1e-30)));
        float w2r = (1.0 + R02) * (1.0 + R02);
        for (int j = 0; j < 5; ++j) {
            float u = urv.x, v = urv.y;
            float u2 = u * u, v2 = v * v;
            vec3 r = vec3(
                u * (1.0 - (1.0 / 3.0) * u2 + v2) - q.x,
                v * (1.0 - (1.0 / 3.0) * v2 + u2) - q.y,
                (u2 - v2) - q.z
            );
            float g1 = r.x * (1.0 - u2 + v2) + r.y * (2.0 * u * v) + r.z * (2.0 * u);
            float g2 = r.x * (2.0 * u * v) + r.y * (1.0 + u2 - v2) - r.z * (2.0 * v);
            float al = -2.0 * u * r.x + 2.0 * v * r.y + 2.0 * r.z;
            float be =  2.0 * v * r.x + 2.0 * u * r.y;
            float hp  = -v * g1 + u * g2;
            float hpp = R02 * w2r - (u * g1 + v * g2) + (v2 - u2) * al - 2.0 * u * v * be;
            float st  = -0.5 * clamp(hp / max(hpp, 0.15 * R02 * w2r), -0.65, 0.65);
            float st2 = st * st;
            urv = vec2((1.0 - st2) * u - 2.0 * st * v, 2.0 * st * u + (1.0 - st2) * v) / (1.0 + st2);
        }
        float u_o = uv.x, v_o = uv.y, u2_o = u_o * u_o, v2_o = v_o * v_o;
        vec3 d_old = vec3(u_o * (1.0 - (1.0 / 3.0) * u2_o + v2_o), v_o * (1.0 - (1.0 / 3.0) * v2_o + u2_o), u2_o - v2_o) - q;
        float u_n = urv.x, v_n = urv.y, u2_n = u_n * u_n, v2_n = v_n * v_n;
        vec3 d_new = vec3(u_n * (1.0 - (1.0 / 3.0) * u2_n + v2_n), v_n * (1.0 - (1.0 / 3.0) * v2_n + u2_n), u2_n - v2_n) - q;
        if (dot(d_new, d_new) < dot(d_old, d_old)) uv = urv;
    }
    return uv;
}

vec4 sdgEnneperSurface(vec3 p, float s, float R0, float th) {
    if (s <= 1e-6 || R0 <= 1e-6) {
        float l2 = dot(p, p);
        float l  = sqrt(l2);
        vec3 g0  = (l > 1e-12) ? (p / l) : vec3(0.0, 0.0, 1.0);
        return vec4(l - th, g0);
    }

    vec2 sxy = vec2(p.x >= 0.0 ? 1.0 : -1.0, p.y >= 0.0 ? 1.0 : -1.0);
    vec3 q   = vec3(abs(p.x), abs(p.y), p.z) / s;
    bool sw  = (q.x < q.y);
    if (sw) q = vec3(q.y, q.x, -q.z);

    float qz_pos = max(q.z, 0.0);
    float qz_neg = max(-q.z, 0.0);

    vec2 seed1 = vec2(max(min(q.x, 1.0), max(0.65 * R0 * step(0.0, -q.z), sqrt(qz_pos) + 0.05)), q.y / (1.0 + qz_pos));
    vec2 seed2 = vec2(q.x / (1.0 + R0 * R0 + qz_neg), max(max(min(q.y, 1.0), 0.72 * R0), sqrt(qz_neg) + 0.05));
    float a_est = clamp(0.62 * R0, 0.35, 0.70 * R0);
    vec2 seed3 = (q.z >= 0.0) ? vec2(a_est * 1.08, a_est * 0.92) : vec2(a_est * 0.92, a_est * 1.08);

    vec2 uv_best = vec2(0.0);
    vec3 X_best  = vec3(0.0);
    float d2_best = dot(q, q);

    vec2 cands[3] = vec2[3](
        _sdEnneperSolveSeed(q, seed1, R0),
        _sdEnneperSolveSeed(q, seed2, R0),
        _sdEnneperSolveSeed(q, seed3, R0)
    );

    for (int k = 0; k < 3; ++k) {
        vec2 uv = cands[k];
        float u2 = uv.x * uv.x, v2 = uv.y * uv.y;
        float rawX = uv.x * (1.0 - (1.0 / 3.0) * u2 + v2);
        float rawY = uv.y * (1.0 - (1.0 / 3.0) * v2 + u2);
        if (rawX < 0.0) { uv.x = -uv.x; rawX = -rawX; }
        if (rawY < 0.0) { uv.y = -uv.y; rawY = -rawY; }
        vec3 X = vec3(rawX, rawY, u2 - v2);
        vec3 df = q - X;
        float d2 = dot(df, df);
        if (d2 < d2_best) {
            d2_best = d2;
            uv_best = uv;
            X_best  = X;
        }
    }

    if (R0 > 1.7320508) {
        float r_ov = min(0.95 * R0, sqrt(max(3.0 + abs(q.z), 3.05)));
        vec2 seed4 = (q.z >= 0.0) ? vec2(-r_ov, 0.15) : vec2(0.15, -r_ov);
        vec2 uv = _sdEnneperSolveSeed(q, seed4, R0);
        float u2 = uv.x * uv.x, v2 = uv.y * uv.y;
        float rawX = uv.x * (1.0 - (1.0 / 3.0) * u2 + v2);
        float rawY = uv.y * (1.0 - (1.0 / 3.0) * v2 + u2);
        if (rawX < 0.0) { uv.x = -uv.x; rawX = -rawX; }
        if (rawY < 0.0) { uv.y = -uv.y; rawY = -rawY; }
        vec3 X = vec3(rawX, rawY, u2 - v2);
        vec3 df = q - X;
        float d2 = dot(df, df);
        if (d2 < d2_best) {
            d2_best = d2;
            uv_best = uv;
            X_best  = X;
        }
    }

    float dq = sqrt(d2_best);
    vec3 df_best = q - X_best;
    vec3 g_fold;
    if (dq > 1e-6) {
        g_fold = df_best / dq;
    } else {
        float r2 = dot(uv_best, uv_best);
        vec3 N = vec3(-2.0 * uv_best.x, 2.0 * uv_best.y, 1.0 - r2) / (1.0 + r2);
        g_fold = (dot(df_best, N) >= 0.0) ? N : -N;
    }

    if (sw) g_fold = vec3(g_fold.y, g_fold.x, -g_fold.z);
    g_fold.xy *= sxy;

    return vec4(s * dq - th, g_fold);
}

float sdEnneperSurface(vec3 p, float s, float R0, float th) {
    return sdgEnneperSurface(p, s, R0, th).x;
}

#endif // SDF_GEMS_HAS_MINIMAL_SURFACE
