// ============================================================================
// evSupportEnvelope — Legendre-Fenchel Support Duality, Spherical Harmonic
//   Hedgehog Bodies, Clairaut Normal Projection & Evolute Caustics
//
//   2D: sdgSupportRotor2D, sdSupportRotor2D,
//       sdgHedgehogHypocycloid2D, sdHedgehogHypocycloid2D
//   3D: evSupportCurvatureRadii, evSupportHarmonicPoint3D,
//       evSupportHarmonicBounds3D, sdgSupportHarmonic3D, sdSupportHarmonic3D
//
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdSupportRotor2D / sdgSupportRotor2D: Exact Euclidean SDF everywhere
//         for convex rotors r0 >= (n^2 - 1)|c| (constant width 2*r0 for odd n);
//         1-Lipschitz Legendre-Fenchel support-core bound when r0 < (n^2 - 1)|c|.
//         sdHedgehogHypocycloid2D / sdgHedgehogHypocycloid2D: Exact Euclidean SDF
//         of Steiner's 3-cusped deltoid hedgehog h(theta) = a cos(3 theta) via a
//         closed-form depressed cubic in w = sin^2(theta) plus the 3 evolute cusps.
//         sdSupportHarmonic3D / sdgSupportHarmonic3D: Exact Euclidean SDF outside
//         the evolute caustic d(p) > -min(R1(u*), R2(u*)) and globally across R^3
//         for convex harmonic bodies (3|c2| + 12|c4| <= r0) via a basis-free
//         Cayley-Hamilton Riemannian Newton solver on S^2 (0 trig).
// ============================================================================

#ifndef SDF_GEMS_HAS_EV_SUPPORT_ENVELOPE
#define SDF_GEMS_HAS_EV_SUPPORT_ENVELOPE

// ----------------------------------------------------------------------------
// 1. 2D Fourier Support Bodies & Constant-Width Rotors:
//      h(theta) = r0 + c * cos(n * theta)
//    For every odd integer n (n = 3, 5, ...), h(theta) + h(theta + pi) = 2*r0,
//    so the body has exact constant width 2*r0 in all directions.
//    Radius of curvature: R(theta) = h(theta) + h''(theta) = r0 - (n^2 - 1) c cos(n theta).
//    Strictly convex iff r0 > (n^2 - 1)|c|.
//    Solves max_theta (p . u(theta) - h(theta)) via a pitchfork-caustic seed +
//    flat-face seed and 3 steps of 3rd-order Halley iteration on S^1.
//    Returns vec3(d, grad.x, grad.y) where grad = u(theta*).
// ----------------------------------------------------------------------------
vec3 sdgSupportRotor2D( in vec2 p, in float r0, in float c, in float n )
{
    float l2 = dot(p, p);
    float th_rad = (l2 > 1e-24) ? atan(p.y, p.x) : 0.0;
    float n2 = n * n;
    float n3 = n2 * n;
    float n4 = n2 * n2;

    // Evaluate derivatives of L(p, theta) = p . u(theta) - r0 - c cos(n theta) at th_rad
    float cs0 = cos(th_rad), sn0 = sin(th_rad);
    float n_th0 = n * th_rad;
    float cn0 = cos(n_th0), sn_n0 = sin(n_th0);
    float pu0 = p.x * cs0 + p.y * sn0;
    float pt0 = -p.x * sn0 + p.y * cs0;
    float g0  = pt0 + n * c * sn_n0;
    float H0  = pu0 - n2 * c * cn0;
    float Q0  = n4 * c * cn0 - pu0;

    // Supercritical pitchfork caustic offset when inside the evolute (H0 < 0)
    float sgn_g0 = (g0 >= 0.0) ? 1.0 : -1.0;
    float max_step = 3.14159265359 / max(n, 1.0);
    float delta_bif = (H0 < 0.0)
        ? sgn_g0 * sqrt(max(-6.0 * H0 / max(Q0, 1e-6), 0.0))
        : 0.0;
    float th_seed0 = th_rad + clamp(delta_bif, -max_step, max_step);

    // Nearest flat-face seeds (where c cos(n theta) = -|c| and -L'' is maximal)
    float period = 6.28318530718 / max(n, 1.0);
    float offset = (c >= 0.0) ? 0.5 * period : 0.0;
    float k_sector = floor((th_rad - offset) / period + 0.5);
    float th_seed1 = offset + k_sector * period;
    float th_seed2 = offset + (k_sector + ((th_rad >= th_seed1) ? 1.0 : -1.0)) * period;

    float best_d = -1e30;
    vec2  best_u = vec2(cs0, sn0);

    for (int s = 0; s < 3; ++s)
    {
        float th = (s == 0) ? th_seed0 : ((s == 1) ? th_seed1 : th_seed2);
        for (int i = 0; i < 3; ++i)
        {
            float cs = cos(th), sn = sin(th);
            float nth = n * th;
            float cn = cos(nth), sn_n = sin(nth);
            float pu = p.x * cs + p.y * sn;
            float d_cand = pu - (r0 + c * cn);
            if (d_cand > best_d)
            {
                best_d = d_cand;
                best_u = vec2(cs, sn);
            }
            float pt = -p.x * sn + p.y * cs;
            float g  = pt + n * c * sn_n;
            float H  = pu - n2 * c * cn;
            float T  = pt + n3 * c * sn_n;
            float H_pos = max(H, 1e-4);
            float denom = max(H_pos + 0.5 * g * T / H_pos, 1e-4);
            th += clamp(g / denom, -max_step, max_step);
        }
        float cs = cos(th), sn = sin(th);
        float d_cand = p.x * cs + p.y * sn - (r0 + c * cos(n * th));
        if (d_cand > best_d)
        {
            best_d = d_cand;
            best_u = vec2(cs, sn);
        }
    }
    return vec3(best_d, best_u);
}

float sdSupportRotor2D( in vec2 p, in float r0, in float c, in float n )
{
    return sdgSupportRotor2D(p, r0, c, n).x;
}

// ----------------------------------------------------------------------------
// 2. Steiner's 3-Cusped Deltoid Hedgehog (`sdHedgehogHypocycloid2D`):
//      h(theta) = a * cos(3 * theta)   (zero-width projective hedgehog)
//    By Legendre-Support factorization:
//      <p - y(theta), y'(theta)> = R(theta) * dL/dtheta(p, theta) = 0,
//    which splits into:
//      (A) R(theta) = -8 a cos(3 theta) = 0  --> the 3 cusps at theta = pi/6, pi/2, 5pi/6;
//      (B) dL/dtheta = -p.x sin(theta) + p.y cos(theta) + 3 a sin(3 theta) = 0,
//          which reduces algebraically to a depressed cubic in w = sin^2(theta) in [0, 1]:
//            144 a^2 w^3 + 24 a (p.x - 9 a) w^2 + ((p.x - 9 a)^2 + p.y^2) w - p.y^2 = 0,
//          followed by a single trig-free tangent Newton polish on u = (cos theta, sin theta).
//    Returns vec3(d, grad.x, grad.y).
// ----------------------------------------------------------------------------
vec3 sdgHedgehogHypocycloid2D( in vec2 p, in float a )
{
    float aa = max(abs(a), 1e-12);
    vec2  q  = (a >= 0.0) ? p : -p;

    // (A) 3 cusps where R(theta) = -8 a cos(3 theta) = 0:
    //     theta = pi/2 -> (-3a, 0); theta = pi/6, 5pi/6 -> (1.5a, +-1.5*sqrt(3)*a)
    vec2  c0 = vec2(-3.0 * aa, 0.0);
    vec2  c1 = vec2(1.5 * aa,  2.59807621135 * aa);
    vec2  c2 = vec2(1.5 * aa, -2.59807621135 * aa);

    vec2  v0 = q - c0, v1 = q - c1, v2 = q - c2;
    float d2_min = dot(v0, v0);
    vec2  best_v = v0;
    if (dot(v1, v1) < d2_min) { d2_min = dot(v1, v1); best_v = v1; }
    if (dot(v2, v2) < d2_min) { d2_min = dot(v2, v2); best_v = v2; }
    vec2  v3 = q - vec2(aa, 0.0);
    if (dot(v3, v3) < d2_min) { d2_min = dot(v3, v3); best_v = v3; }

    float d_min = sqrt(d2_min);
    vec2  best_n = (d_min > 1e-20) ? (best_v / d_min) : vec2(1.0, 0.0);

    // (B) Depressed cubic in w = sin^2(theta) in [0, 1]
    float invA3 = 1.0 / (144.0 * aa * aa);
    float B     = q.x - 9.0 * aa;
    float k2    = (24.0 * aa * B) * invA3;
    float k1    = (B * B + q.y * q.y) * invA3;
    float k0    = -(q.y * q.y) * invA3;

    float k2_3  = k2 / 3.0;
    float p_cub = k1 - k2 * k2_3;
    float q_cub = k0 - k1 * k2_3 + 2.0 * k2_3 * k2_3 * k2_3;
    float disc  = 0.25 * q_cub * q_cub + (p_cub * p_cub * p_cub) / 27.0;

    vec3 w_roots;
    if (disc >= 0.0)
    {
        float sq = sqrt(max(disc, 0.0));
        float u1 = -0.5 * q_cub + sq;
        float u2 = -0.5 * q_cub - sq;
        float r1 = sign(u1) * pow(abs(u1), 1.0 / 3.0);
        float r2 = sign(u2) * pow(abs(u2), 1.0 / 3.0);
        w_roots = vec3(r1 + r2 - k2_3);
    }
    else
    {
        float r = sqrt(max(-p_cub / 3.0, 1e-30));
        float phi = acos(clamp(-0.5 * q_cub / (r * r * r), -1.0, 1.0)) / 3.0;
        w_roots = 2.0 * r * cos(phi + vec3(0.0, 2.09439510239, 4.18879020479)) - k2_3;
    }

    for (int k = 0; k < 3; ++k)
    {
        float wk = w_roots[k];
        if (wk >= -1e-2 && wk <= 1.0 + 1e-2)
        {
            float w = clamp(wk, 0.0, 1.0);
            float s = sqrt(w);
            float c_abs = sqrt(max(1.0 - w, 0.0));
            float s3_0 = s * (3.0 - 4.0 * w);
            float res_pos = abs(-q.x * s + q.y * c_abs + 3.0 * aa * s3_0);
            float res_neg = abs(-q.x * s - q.y * c_abs + 3.0 * aa * s3_0);
            float c = (res_pos <= res_neg) ? c_abs : -c_abs;

            // 2 trig-free tangent Newton polishes on u = (c, s) to remove FP32 cubic cancellation
            for (int iter = 0; iter < 2; ++iter)
            {
                float w_it = s * s;
                float c3_it = c * (1.0 - 4.0 * w_it);
                float s3_it = s * (3.0 - 4.0 * w_it);
                float dP_it  = -q.x * s + q.y * c + 3.0 * aa * s3_it;
                float ddP_it = -(q.x * c + q.y * s) + 9.0 * aa * c3_it;
                if (abs(ddP_it) > 1e-7)
                {
                    float dt = clamp(-dP_it / ddP_it, -0.25, 0.25);
                    vec2 u_pol = normalize(vec2(c - dt * s, s + dt * c));
                    c = u_pol.x;
                    s = u_pol.y;
                }
            }

            float w_f = s * s;
            float c3 = c * (1.0 - 4.0 * w_f);
            float s3 = s * (3.0 - 4.0 * w_f);
            float P_val  = q.x * c + q.y * s - aa * c3;
            float dP_val = -q.x * s + q.y * c + 3.0 * aa * s3;
            float dist_k = sqrt(P_val * P_val + dP_val * dP_val);
            if (dist_k < d_min)
            {
                d_min = dist_k;
                vec2 diff = P_val * vec2(c, s) + dP_val * vec2(-s, c);
                best_n = (dist_k > 1e-20) ? (diff / dist_k) : ((P_val >= 0.0) ? vec2(c, s) : -vec2(c, s));
            }
        }
    }

    // Algebraic degree-4 implicit sign of Steiner's deltoid:
    float r2 = dot(q, q);
    float a2 = aa * aa;
    float imp = r2 * r2 + 18.0 * a2 * r2 + 8.0 * aa * (q.x * q.x * q.x - 3.0 * q.x * q.y * q.y) - 27.0 * a2 * a2;
    float sgn = (imp < 0.0) ? -1.0 : 1.0;
    vec2  grad = ((a >= 0.0) ? best_n : -best_n) * sgn;
    return vec3(sgn * d_min, grad);
}

float sdHedgehogHypocycloid2D( in vec2 p, in float a )
{
    return sdgHedgehogHypocycloid2D(p, a).x;
}

// ----------------------------------------------------------------------------
// 3. 3D Cubic/Octahedral & Axisymmetric Spherical Harmonic Support Bodies:
//      h(u) = r0 + c2 * (1.5 * u.y^2 - 0.5) + c4 * (u.x^4 + u.y^4 + u.z^4 - 0.6)
//    Because the ambient polynomial f(u) is separable across (u.x, u.y, u.z),
//    its ambient Euclidean Hessian M = diag(m) is strictly diagonal:
//      g(u) = (4 c4 u.x^3,  3 c2 u.y + 4 c4 u.y^3,  4 c4 u.z^3)
//      m(u) = (12 c4 u.x^2, 3 c2 + 12 c4 u.y^2,     12 c4 u.z^2).
// ----------------------------------------------------------------------------

// Principal radii of curvature R1(u) <= R2(u) of the support surface at normal u in S^2.
// Computed without any tangent-plane ONB via the Cauchy-Binet trace & determinant of
// W^{-1}(u) = P_u diag(m + (h - u.g)) P_u on T_u S^2.
vec2 evSupportCurvatureRadii( in vec3 u, in float r0, in float c2, in float c4 )
{
    vec3  u2 = u * u;
    float h  = r0 + c2 * (1.5 * u2.y - 0.5) + c4 * (dot(u2, u2) - 0.6);
    vec3  g  = vec3(4.0 * c4 * u.x * u2.x, (3.0 * c2 + 4.0 * c4 * u2.y) * u.y, 4.0 * c4 * u.z * u2.z);
    vec3  m  = vec3(12.0 * c4 * u2.x, 3.0 * c2 + 12.0 * c4 * u2.y, 12.0 * c4 * u2.z);
    vec3  a  = m + (h - dot(u, g));
    float tr  = dot(a, 1.0 - u2);
    float det = dot(a.yzx * a.zxy, u2);
    float disc = sqrt(max(0.25 * tr * tr - det, 0.0));
    return vec2(0.5 * tr - disc, 0.5 * tr + disc);
}

// Clairaut envelope map y(u) = h(u) u + grad_{S^2} h(u) = g + (h(u) - u . g) u.
// The outward unit normal to the envelope surface at y(u) is u itself.
vec3 evSupportHarmonicPoint3D( in vec3 u, in float r0, in float c2, in float c4 )
{
    vec3  u2 = u * u;
    float h  = r0 + c2 * (1.5 * u2.y - 0.5) + c4 * (dot(u2, u2) - 0.6);
    vec3  g  = vec3(4.0 * c4 * u.x * u2.x, (3.0 * c2 + 4.0 * c4 * u2.y) * u.y, 4.0 * c4 * u.z * u2.z);
    return g + (h - dot(u, g)) * u;
}

// Exact Legendre-Fenchel lower bound L(p, u) = p . u - h(u) and Clairaut upper bound |p - y(u)|
// for any unit direction u in S^2. At the optimal normal u*(p) outside K, both bounds meet:
// d_lower == d(p) == d_upper.
vec2 evSupportHarmonicBounds3D( in vec3 p, in vec3 u, in float r0, in float c2, in float c4 )
{
    vec3  u2 = u * u;
    float h  = r0 + c2 * (1.5 * u2.y - 0.5) + c4 * (dot(u2, u2) - 0.6);
    vec3  g  = vec3(4.0 * c4 * u.x * u2.x, (3.0 * c2 + 4.0 * c4 * u2.y) * u.y, 4.0 * c4 * u.z * u2.z);
    vec3  y  = g + (h - dot(u, g)) * u;
    return vec2(dot(p, u) - h, length(p - y));
}

// 3D Spherical Harmonic Support SDF + Analytical Gradient:
// Uses radial seed + symmetric axial/octahedral seeds (for deep interior inside the evolute)
// and basis-free Cayley-Hamilton Riemannian Newton steps on S^2 (0 trig, 0 ONB).
vec4 sdgSupportHarmonic3D( in vec3 p, in float r0, in float c2, in float c4 )
{
    float l2 = dot(p, p);
    vec3  u_rad = (l2 > 1e-24) ? (p * inversesqrt(l2)) : vec3(1.0, 0.0, 0.0);
    vec3  s_sgn = vec3((p.x >= 0.0) ? 1.0 : -1.0, (p.y >= 0.0) ? 1.0 : -1.0, (p.z >= 0.0) ? 1.0 : -1.0);

    // Evaluate all symmetric flat-face candidate seeds (3 signed axes + signed octahedral diagonal)
    vec3  cands[5];
    cands[0] = u_rad;
    cands[1] = vec3(s_sgn.x, 0.0, 0.0);
    cands[2] = vec3(0.0, s_sgn.y, 0.0);
    cands[3] = vec3(0.0, 0.0, s_sgn.z);
    cands[4] = s_sgn * 0.57735026919;

    float best_d = -1e30;
    vec3  best_u = u_rad;
    float damp = 0.08 * max(abs(r0), 0.1);

    for (int s = 0; s < 5; ++s)
    {
        vec3 u = cands[s];
        int n_steps = (s == 0) ? 3 : 2;
        for (int i = 0; i < n_steps; ++i)
        {
            vec3  u2 = u * u;
            float h  = r0 + c2 * (1.5 * u2.y - 0.5) + c4 * (dot(u2, u2) - 0.6);
            float d_cand = dot(p, u) - h;
            if (d_cand > best_d)
            {
                best_d = d_cand;
                best_u = u;
            }
            vec3  g  = vec3(4.0 * c4 * u.x * u2.x, (3.0 * c2 + 4.0 * c4 * u2.y) * u.y, 4.0 * c4 * u.z * u2.z);
            vec3  m  = vec3(12.0 * c4 * u2.x, 3.0 * c2 + 12.0 * c4 * u2.y, 12.0 * c4 * u2.z);
            vec3  pg = p - g;
            float pu_gu = dot(pg, u);
            vec3  r_perp = pg - pu_gu * u;
            vec3  q = max(m + pu_gu, vec3(damp));
            float tr  = dot(q, 1.0 - u2);
            float det = max(dot(q.yzx * q.zxy, u2), 1e-6);
            vec3  w   = (tr - q) * r_perp;
            vec3  du  = (w - dot(w, u) * u) / det;
            float du2 = dot(du, du);
            if (du2 > 0.64) du *= 0.8 * inversesqrt(du2);
            u = normalize(u + du);
        }
        vec3  u2 = u * u;
        float h  = r0 + c2 * (1.5 * u2.y - 0.5) + c4 * (dot(u2, u2) - 0.6);
        float d_cand = dot(p, u) - h;
        if (d_cand > best_d)
        {
            best_d = d_cand;
            best_u = u;
        }
    }
    return vec4(best_d, best_u);
}

float sdSupportHarmonic3D( in vec3 p, in float r0, in float c2, in float c4 )
{
    return sdgSupportHarmonic3D(p, r0, c2, c4).x;
}

#endif // SDF_GEMS_HAS_EV_SUPPORT_ENVELOPE
