// ============================================================================
// syHamiltonianFlow — Symplectic Phase-Space Pullbacks, Singular-Value Bounds &
//                     Hamiltonian Whorl SDFs
//   sySingularValue2D, sySingularValueShear2D,
//   opSymplecticHenonMap, opSymplecticHenonInv, opSymplecticHenonPointwise, opSymplecticHenon,
//   opSymplecticChirikovMap, opSymplecticChirikovInv, opSymplecticChirikovPointwise, opSymplecticChirikov,
//   sdgDuffingSeparatrix2D, sdDuffingSeparatrix2D,
//   sdgPhaseSpaceVortex3D, sdPhaseSpaceVortex3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sySingularValue2D / sySingularValueShear2D: Exact spectral norm sigma_max(J) of any 2x2
//           unimodular/symplectic matrix (|det(J)| = 1) and its 3D extrusion diag(J, 1).
//         sdDuffingSeparatrix2D / sdgDuffingSeparatrix2D: Exact 2D Euclidean SDF and unit gradient
//           to the double-well Duffing homoclinic figure-eight separatrix H(q,p) = p^2/2 - q^2/2 + q^4/4 = 0.
//         opSymplecticHenon / opSymplecticChirikov / sdPhaseSpaceVortex3D: Certified conservative
//           Euclidean distance bounds (0% raymarch overstep everywhere; exact at a = 0 / k = 0).
//         opSymplecticHenonPointwise / opSymplecticChirikovPointwise: First-order surface-normalized
//           pullbacks (exact local metric at d = 0, but can overstep far from the surface when curvature != 0).
// ============================================================================

#ifndef SDF_GEMS_HAS_SY_HAMILTONIAN_FLOW
#define SDF_GEMS_HAS_SY_HAMILTONIAN_FLOW

// ----------------------------------------------------------------------------
// 1. Exact Symplectic Singular-Value Theorem (2D & 3D Extruded Phase Space)
//    For any 2x2 matrix J with |det(J)| = 1, its singular values are reciprocals
//    (sigma, 1/sigma) and ||J||_F^2 = sigma^2 + sigma^-2 >= 2. Hence:
//      sigma_max(J) = 0.5 * (sqrt(||J||_F^2 + 2) + sqrt(||J||_F^2 - 2)).
//    Moreover, for 3D orthogonal phase-space extrusion Psi(q, p, z) = (Phi(q, p), z),
//    J_Psi = diag(J_Phi, 1) has singular values {sigma, 1, 1/sigma}, so
//    sigma_max(J_Psi) == sigma_max(J_Phi) >= 1.
// ----------------------------------------------------------------------------

// Exact maximum singular value from squared Frobenius norm ||J||_F^2 of a |det(J)| = 1 matrix.
float sySingularValue2D( in float frob2 )
{
    return 0.5 * (sqrt(max(frob2 + 2.0, 4.0)) + sqrt(max(frob2 - 2.0, 0.0)));
}

// Single-sqrt specialization for symplectic shears where ||J||_F^2 = 2 + tau^2.
float sySingularValueShear2D( in float tau )
{
    return 0.5 * (sqrt(4.0 + tau * tau) + abs(tau));
}

// Cancellation-free exact maximum singular value of a 2x2 matrix J.
// For |det(J)| = 1, the two hypot terms equal sqrt(||J||_F^2 + 2) and sqrt(||J||_F^2 - 2)
// without FP32 subtractive cancellation near J = I (where ||J||_F^2 - 2 = O(eps^2)).
float sySingularValue2D( in mat2 J )
{
    vec2 u = vec2(J[0][0] + J[1][1], J[0][1] - J[1][0]);
    vec2 v = vec2(J[0][0] - J[1][1], J[0][1] + J[1][0]);
    return 0.5 * (length(u) + length(v));
}

// ----------------------------------------------------------------------------
// 2. Reversible Hénon Quadratic Symplectomorphism
//      Phi(q, p) = (p - k - a q^2, q),   det(J_Phi) = -1,   ||J_Phi||_F^2 = 2 + 4 a^2 q^2
//      sigma_max(J_Phi(q, p)) = sqrt(1 + a^2 q^2) + |a q|
// ----------------------------------------------------------------------------

// Forward Hénon area-preserving quadratic map.
vec2 opSymplecticHenonMap( in vec2 qp, in float a, in float k )
{
    return vec2(qp.y - k - a * qp.x * qp.x, qp.x);
}

// Exact degree-2 polynomial inverse of opSymplecticHenonMap.
vec2 opSymplecticHenonInv( in vec2 qp, in float a, in float k )
{
    return vec2(qp.y, qp.x + k + a * qp.y * qp.y);
}

// (a) Pointwise Symplectic Metric Normalization: divides d0(Phi(qp)) by local sigma_max(J(qp)).
// Exact local Lipschitz factor on the zero surface d0 = 0, but NOT globally conservative when a != 0.
float opSymplecticHenonPointwise( in vec2 qp, in float d0, in float a )
{
    float aq = abs(a * qp.x);
    return d0 / (sqrt(1.0 + aq * aq) + aq);
}

// (b) Certified Conservative Symplectic Bound over the self-bounding ball B(qp, |d0|).
// Because sigma_max(J_Phi) >= 1 everywhere, any nearest surface point at distance R < |d0| lies
// strictly inside B(qp, |d0|), where max_{|y - qp| <= |d0|} |y_x| = |q| + |d0|.
// Eliminates 100% of sphere-tracing oversteps while converging to opSymplecticHenonPointwise as d0 -> 0.
float opSymplecticHenon( in vec2 qp, in float d0, in float a )
{
    float q_max = abs(qp.x) + abs(d0);
    float aq = abs(a) * q_max;
    return d0 / (sqrt(1.0 + aq * aq) + aq);
}

// Exact unit surface normal and conservative distance for a Hénon-pulled-back 2D SDF dg0 = (d0, gx, gy).
// Uses directional Jacobian pullback J^T * g0 so ||grad d|| == 1.0 identically on the zero set d0 = 0.
vec3 sdgSymplecticHenonSurface( in vec2 qp, in vec3 dg0, in float a )
{
    // J = [[-2aq, 1], [1, 0]] => J^T * g0 = (-2aq * g0.x + g0.y, g0.x)
    vec2 g_pull = vec2(-2.0 * a * qp.x * dg0.y + dg0.z, dg0.y);
    float l = length(g_pull);
    vec2 n = (l > 1e-20) ? (g_pull / l) : vec2(1.0, 0.0);
    return vec3(opSymplecticHenon(qp, dg0.x, a), n);
}

// ----------------------------------------------------------------------------
// 3. Symplectic Kick-Drift / Chirikov Standard Map
//      p' = p + k sin(w q),   q' = q + tau p'
//      J = [[1 + tau k w cos(w q), tau], [k w cos(w q), 1]],   det(J) == 1
//      ||J||_F^2 = 2 + (tau + c)^2 + (tau c)^2   where c = k w cos(w q).
// ----------------------------------------------------------------------------

// Forward symplectic kick-drift map.
vec2 opSymplecticChirikovMap( in vec2 qp, in float k, in float w, in float tau )
{
    float p1 = qp.y + k * sin(w * qp.x);
    float q1 = qp.x + tau * p1;
    return vec2(q1, p1);
}

// Exact closed-form inverse (drift back, then un-kick).
vec2 opSymplecticChirikovInv( in vec2 qp, in float k, in float w, in float tau )
{
    float q0 = qp.x - tau * qp.y;
    float p0 = qp.y - k * sin(w * q0);
    return vec2(q0, p0);
}

// (a) Pointwise Symplectic Metric Normalization for the Chirikov kick-drift map.
float opSymplecticChirikovPointwise( in vec2 qp, in float d0, in float k, in float w, in float tau )
{
    float c = k * w * cos(w * qp.x);
    float tc = tau + c;
    float p  = tau * c;
    float frob2 = 2.0 + tc * tc + p * p;
    return d0 / sySingularValue2D(frob2);
}

// (b) Certified Conservative Symplectic Bound for the Chirikov kick-drift map over B(qp, |d0|).
// Uses the exact interval maximum of |cos(w y)| on [q - |d0|, q + |d0|]:
//   max_{|y - q| <= r} |cos(w y)| = cos(max(|w q - pi * round(w q / pi)| - |w| r, 0)).
float opSymplecticChirikov( in vec2 qp, in float d0, in float k, in float w, in float tau )
{
    float r = abs(d0);
    float th = w * qp.x;
    float d_pi = abs(th - SDF_PI * round(th / SDF_PI));
    float cos_max = cos(max(d_pi - abs(w) * r, 0.0));
    float c_max = abs(k * w) * cos_max;
    float at = abs(tau);
    float tc = at + c_max;
    float p  = at * c_max;
    float frob2 = 2.0 + tc * tc + p * p;
    return d0 / sySingularValue2D(frob2);
}

// ----------------------------------------------------------------------------
// 4. Exact 2D Duffing Figure-Eight Homoclinic Separatrix SDF
//    Hamiltonian: H(q, p) = 0.5 p^2 - 0.5 q^2 + 0.25 q^4 = 0.
//    Rational Weierstrass parameterization in first quadrant (u = tan(theta/2) in [0, 1]):
//      q(u) = 2 sqrt(2) u / (1 + u^2),   p(u) = 2 sqrt(2) u (1 - u^2) / (1 + u^2)^2.
//    0 trig calls, 1 sqrt total, 2-branch safeguarded Newton on u in [0, 1].
// ----------------------------------------------------------------------------
float sdgDuffingRefineBranch( in vec2 x, in float u0, in float u_lo_in, in float u_hi_in )
{
    const float S2 = 2.8284271247461903;
    float u = u0, u_lo = u_lo_in, u_hi = u_hi_in;
    for (int i = 0; i < 4; i++)
    {
        float u2  = u * u;
        float w   = 1.0 / (1.0 + u2);
        float w2  = w * w;
        float w3  = w2 * w;
        float q   = S2 * u * w;
        float p   = q * (1.0 - u2) * w;
        float dq  = S2 * (1.0 - u2) * w2;
        float dp  = S2 * (1.0 - 6.0 * u2 + u2 * u2) * w3;
        float ddq = 2.0 * S2 * u * (u2 - 3.0) * w3;
        float ddp = -4.0 * S2 * u * (u2 * u2 - 10.0 * u2 + 5.0) * w2 * w2;

        vec2  r  = vec2(q, p) - x;
        float f  = r.x * dq + r.y * dp;
        float df = dq * dq + dp * dp + r.x * ddq + r.y * ddp;

        if (f < 0.0) u_lo = max(u_lo, u);
        else         u_hi = min(u_hi, u);

        float u_new = (df > 1e-12) ? (u - f / df) : (0.5 * (u_lo + u_hi));
        u = (df > 1e-12 && u_new >= u_lo && u_new <= u_hi) ? u_new : (0.5 * (u_lo + u_hi));
    }
    return u;
}

vec3 sdgDuffingSeparatrix2D( in vec2 qp )
{
    vec2  x  = abs(qp);
    const float S2 = 2.8284271247461903; // 2 * sqrt(2)

    // The lobe peak is at u_peak = sqrt(2) - 1 = 0.41421356.
    // On [0, u_peak] and [u_peak, 1], D^2(u) is unimodal; refine the best candidate in each half!
    float best_u1 = 0.0, best_d2_1 = dot(x, x);
    for (int i = 1; i <= 4; i++)
    {
        float u = float(i) * 0.10355339;
        float w = 1.0 / (1.0 + u * u);
        vec2  d = vec2(S2 * u * w, S2 * u * (1.0 - u * u) * w * w) - x;
        float d2 = dot(d, d);
        if (d2 < best_d2_1) { best_d2_1 = d2; best_u1 = u; }
    }
    float u1 = sdgDuffingRefineBranch(x, best_u1, max(best_u1 - 0.10355339, 0.0), min(best_u1 + 0.10355339, 0.41421356));

    float best_u2 = 0.41421356, best_d2_2 = best_d2_1;
    for (int i = 1; i <= 4; i++)
    {
        float u = 0.41421356 + float(i) * 0.14644661;
        float w = 1.0 / (1.0 + u * u);
        vec2  d = vec2(S2 * u * w, S2 * u * (1.0 - u * u) * w * w) - x;
        float d2 = dot(d, d);
        if (d2 < best_d2_2) { best_d2_2 = d2; best_u2 = u; }
    }
    float u2_cand = sdgDuffingRefineBranch(x, best_u2, max(best_u2 - 0.14644661, 0.41421356), min(best_u2 + 0.14644661, 1.0));

    float w1 = 1.0 / (1.0 + u1 * u1);
    vec2  c1 = vec2(S2 * u1 * w1, S2 * u1 * (1.0 - u1 * u1) * w1 * w1);
    float d2 = dot(x - c1, x - c1);
    float u  = u1;
    vec2  c  = c1;

    float w2 = 1.0 / (1.0 + u2_cand * u2_cand);
    vec2  c2 = vec2(S2 * u2_cand * w2, S2 * u2_cand * (1.0 - u2_cand * u2_cand) * w2 * w2);
    float d2_b = dot(x - c2, x - c2);
    if (d2_b < d2) { d2 = d2_b; u = u2_cand; c = c2; }

    // Also check exact endpoints u = 0 -> (0, 0) and u = 1 -> (sqrt(2), 0)
    float d2_0 = dot(x, x);
    if (d2_0 < d2) { d2 = d2_0; c = vec2(0.0); u = 0.0; }
    vec2  ce1  = vec2(1.4142135623730951, 0.0);
    float d2_1 = dot(x - ce1, x - ce1);
    if (d2_1 < d2) { d2 = d2_1; c = ce1; u = 1.0; }

    vec2  diff = x - c;
    float dist = sqrt(d2);
    // Inside the two homoclinic lobes iff H(q, p) = 0.5 p^2 - 0.5 q^2 + 0.25 q^4 < 0
    float H = 0.5 * qp.y * qp.y - 0.5 * qp.x * qp.x * (1.0 - 0.5 * qp.x * qp.x);
    float sgn = (H < 0.0) ? -1.0 : 1.0;

    vec2 n;
    if (dist > 1e-6)
    {
        n = sgn * (diff / dist);
    }
    else
    {
        float u2 = u * u;
        float w  = 1.0 / (1.0 + u2);
        float dq = S2 * (1.0 - u2) * w * w;
        float dp = S2 * (1.0 - 6.0 * u2 + u2 * u2) * w * w * w;
        n = normalize(vec2(-dp, dq));
    }
    n *= vec2((qp.x >= 0.0) ? 1.0 : -1.0, (qp.y >= 0.0) ? 1.0 : -1.0);
    return vec3(sgn * dist, n);
}

float sdDuffingSeparatrix2D( in vec2 qp )
{
    return sdgDuffingSeparatrix2D(qp).x;
}

// ----------------------------------------------------------------------------
// 5. 3D Symplectic Phase-Space Whorl (sdPhaseSpaceVortex3D / sdgPhaseSpaceVortex3D)
//    Extrudes a Hénon-sheared Duffing homoclinic separatrix (offset by r_offset in
//    canonical phase space) in the horizontal (x, z) plane with vertical half-height h.
//    - At a == 0: Exact 3D Euclidean SDF with exact unit analytical gradient.
//    - At a != 0: Certified conservative 3D Euclidean bound (0% overstep everywhere).
// ----------------------------------------------------------------------------
vec4 sdgPhaseSpaceVortex3D( in vec3 p, in float a, in float h, in float r )
{
    vec2 qp = p.xz;
    vec2 w_qp = (abs(a) > 1e-8) ? vec2(qp.x, qp.y - a * qp.x * qp.x) : qp;
    vec3 dg2 = sdgDuffingSeparatrix2D(w_qp);
    float d0_2d = dg2.x - r;

    // Apply certified symplectic shear bound in 2D phase space and pullback normal
    float d2;
    vec2  g2;
    if (abs(a) > 1e-8)
    {
        float q_max = abs(qp.x) + abs(d0_2d);
        float aq = abs(a) * q_max;
        d2 = d0_2d / (sqrt(1.0 + aq * aq) + aq);
        // J = [[1, 0], [-2aq, 1]] => J^T g = (g.x - 2aq g.y, g.y)
        vec2 gp = vec2(dg2.y - 2.0 * a * qp.x * dg2.z, dg2.z);
        float lg = length(gp);
        g2 = (lg > 1e-20) ? (gp / lg) : vec2(1.0, 0.0);
    }
    else
    {
        d2 = d0_2d;
        g2 = dg2.yz;
    }

    // Orthogonal y-extrusion with half-height h
    float dy = abs(p.y) - h;
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;
    vec2  w  = vec2(d2, dy);
    vec2  wp = max(w, 0.0);
    float l  = length(wp);
    float d  = min(max(w.x, w.y), 0.0) + l;

    vec3 n;
    if (w.x > 0.0 && w.y > 0.0 && l > 1e-12)
    {
        vec2 wn = wp / l;
        n = vec3(wn.x * g2.x, wn.y * sy, wn.x * g2.y);
    }
    else if (w.x > w.y)
    {
        n = vec3(g2.x, 0.0, g2.y);
    }
    else
    {
        n = vec3(0.0, sy, 0.0);
    }
    return vec4(d, normalize(n));
}

float sdPhaseSpaceVortex3D( in vec3 p, in float a, in float h, in float r )
{
    return sdgPhaseSpaceVortex3D(p, a, h, r).x;
}

#endif // SDF_GEMS_HAS_SY_HAMILTONIAN_FLOW
