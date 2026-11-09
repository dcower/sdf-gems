#ifndef SDF_GEMS_OP_DIRAC_BELT_GLSL
#define SDF_GEMS_OP_DIRAC_BELT_GLSL

#include "../common/common.glsl"

/**
 * @name opDiracBelt
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.8 Flash (experiments & implementation)
 * @summary Dirac's 720-degree SU(2) spinor belt trick & volume-preserving spherical-shell deformation.
 *
 * In SO(3), a 360-degree (2*pi) rotation is topologically non-trivial (pi_1(SO(3)) = Z_2),
 * whereas a 720-degree (4*pi) rotation lifts to a contractible loop in SU(2) = S^3.
 *
 * This operator implements Dirac's belt trick / Feynman's plate trick in closed algebraic form.
 * For an inner core radius r1 and outer stationary radius r2, space is smoothly deformed by
 * a quaternion spinor homotopy q(r, t). At phase t = 0.5 (360 degrees), the inner core has
 * rotated by 360 degrees while the outer boundary r >= r2 is strictly stationary, producing a
 * twisted intermediate shell. At phase t = 1.0 (720 degrees), the entire deformation smoothly
 * untangles back to the identity map without any seam or discontinuity.
 *
 * The transformation preserves concentric spheres ||Phi_t(p)|| = ||p||, yielding det(J) = 1.000000
 * (strictly incompressible 3D shear) with an exact SVD singular-value Lipschitz bound.
 */

// Rotate vector v by unit quaternion q = vec4(xyz, w)
vec3 diracQuatRotate( in vec3 v, in vec4 q )
{
    vec3 u = q.xyz;
    float w = q.w;
    return v + 2.0 * w * cross(u, v) + 2.0 * cross(u, cross(u, v));
}

// Rotate vector v by conjugate quaternion q* = vec4(-xyz, w) (inverse rotation R^T)
vec3 diracQuatRotateInv( in vec3 v, in vec4 q )
{
    vec3 u = -q.xyz;
    float w = q.w;
    return v + 2.0 * w * cross(u, v) + 2.0 * cross(u, cross(u, v));
}

// 1. Dirac Belt-Trick SU(2) Spinor Coordinate Deformation
//    Maps p -> p' = R(Q(r, phase)) * p
vec3 opDiracBelt( in vec3 p, in float r1, in float r2, in float phase )
{
    if (r1 >= r2) return p;
    float r = length(p);
    if (r >= r2 || r < 1e-12) return p;

    float sin_pt = sin(3.141592653589793 * phase);
    if (abs(sin_pt) < 1e-7) return p;

    float sin_2pt = sin(6.283185307179586 * phase);
    float cos_2pt = cos(6.283185307179586 * phase);

    if (r <= r1) {
        // Rigid core rotation around y-axis by 4 * pi * phase
        vec4 q_core = vec4(0.0, sin_2pt, 0.0, cos_2pt);
        return diracQuatRotate(p, q_core);
    }

    float u = (r - r1) / (r2 - r1);
    float s = u * u * (3.0 - 2.0 * u);
    float beta = 1.5707963267948966 * (1.0 - s);
    float sin_b = sin(beta);
    float cos_b = cos(beta);
    float sin_2b = 2.0 * sin_b * cos_b;
    float sin_pt2 = sin_pt * sin_pt;

    vec4 q = vec4(sin_2b * sin_pt2,
                  sin_b * sin_b * sin_2pt,
                 -0.5 * sin_2b * sin_2pt,
                  1.0 - 2.0 * sin_b * sin_b * sin_pt2);

    return diracQuatRotate(p, q);
}

// 2. Exact Shear-SVD Lipschitz Distance Bound Scaling
//    Scales a raw pulled-back distance dRaw to guarantee a 1-Lipschitz conservative bound
//    (0.0% raymarching overstep across the entire shell and its boundaries).
float opDiracBeltScale( in vec3 p, in float dRaw, in float r1, in float r2, in float phase )
{
    if (r1 >= r2) return dRaw;
    float sin_pt = abs(sin(3.141592653589793 * phase));
    if (sin_pt < 1e-7) return dRaw;

    // Exact maximum shear rate k_max over the shell r in [r1, r2]
    float dr = r2 - r1;
    float disc = r2 * r2 - r1 * r2 + r1 * r1;
    float u_star = clamp((r2 - 2.0 * r1 + sqrt(disc)) / (3.0 * dr), 0.0, 1.0);
    float r_star = r1 + u_star * dr;
    float dbeta_star = 9.42477796076938 * u_star * (1.0 - u_star) / dr;
    float k_max = 4.0 * r_star * dbeta_star * sin_pt;

    // Reciprocal singular value 1 / sigma_max = 0.5 * (sqrt(4 + k^2) - k) (0 divisions)
    float h = sqrt(4.0 + k_max * k_max);
    float inv_sig = 0.5 * (h - k_max);

    float r = length(p);
    if (r >= r2) {
        float db = r - r2;
        return (dRaw <= db) ? dRaw : (db + (dRaw - db) * inv_sig);
    } else if (r <= r1) {
        float db = r1 - r;
        return (dRaw <= db) ? dRaw : (db + (dRaw - db) * inv_sig);
    } else {
        return dRaw * inv_sig;
    }
}

// 3. Exact 1-Jet Chain-Rule Gradient Propagation (sdgDiracBelt)
//    Transforms dgWarped = (d, grad_warped) by the transpose Jacobian J^T:
//      g = J^T * g_warped = R^T * g_warped + p_hat * (p'_r . g_warped)
//    and normalizes to unit gradient length.
vec4 sdgDiracBelt( in vec3 p, in vec4 dgWarped, in float r1, in float r2, in float phase )
{
    if (r1 >= r2) return dgWarped;
    float sin_pt = sin(3.141592653589793 * phase);
    if (abs(sin_pt) < 1e-7) return dgWarped;

    float r = length(p);
    if (r >= r2 || r < 1e-12) return dgWarped;

    float sin_2pt = sin(6.283185307179586 * phase);
    float cos_2pt = cos(6.283185307179586 * phase);

    if (r <= r1) {
        // Rigid core rotation: J = R, so J^T = R^T
        vec4 q_core = vec4(0.0, sin_2pt, 0.0, cos_2pt);
        vec3 g = diracQuatRotateInv(dgWarped.yzw, q_core);
        float invG = inversesqrt(max(dot(g, g), 1e-30));
        return vec4(dgWarped.x * invG, g * invG);
    }

    float dr = r2 - r1;
    float u = (r - r1) / dr;
    float s = u * u * (3.0 - 2.0 * u);
    float beta = 1.5707963267948966 * (1.0 - s);
    float sin_b = sin(beta);
    float cos_b = cos(beta);
    float sin_2b = 2.0 * sin_b * cos_b;
    float cos_2b = cos_b * cos_b - sin_b * sin_b;
    float sin_pt2 = sin_pt * sin_pt;

    vec4 q = vec4(sin_2b * sin_pt2,
                  sin_b * sin_b * sin_2pt,
                 -0.5 * sin_2b * sin_2pt,
                  1.0 - 2.0 * sin_b * sin_b * sin_pt2);

    vec4 dq_db = vec4(2.0 * cos_2b * sin_pt2,
                      sin_2b * sin_2pt,
                     -cos_2b * sin_2pt,
                     -2.0 * sin_2b * sin_pt2);

    float dbeta_dr = -9.42477796076938 * u * (1.0 - u) / dr;
    vec4 dq_dr = dq_db * dbeta_dr;

    // Radial velocity p'_r = d/dr(R(r) * p)
    vec3 cx = cross(q.xyz, p);
    vec3 cdx = cross(dq_dr.xyz, p);
    vec3 p_prime_r = 2.0 * (dq_dr.w * cx + q.w * cdx + cross(dq_dr.xyz, cx) + cross(q.xyz, cdx));

    // J^T * g_warped = R^T * g_warped + p_hat * (p'_r . g_warped)
    vec3 g_rot = diracQuatRotateInv(dgWarped.yzw, q);
    vec3 g = g_rot + (p / r) * dot(p_prime_r, dgWarped.yzw);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dgWarped.x * invG, g * invG);
}

// 4. 3D Demonstration Primitive (sdDiracBeltDemo)
//    Stationary outer anchor ring at r2 = 1.95, connected by 4 rectangular straps to a
//    central hub at r1 = 0.55 twisted by 360 degrees (phase = 0.5), illustrating the
//    SU(2) spinor deformation and untangling.
float sdDiracBeltDemo( in vec3 p )
{
    const float r1 = 0.55;
    const float r2 = 1.95;
    const float phase = 0.5;

    // Warped coordinates
    vec3 q = opDiracBelt(p, r1, r2, phase);

    // Central hub (rounded box inside r <= r1)
    vec3 q_hub = abs(q) - vec3(0.26);
    float d_hub = length(max(q_hub, 0.0)) + min(max(q_hub.x, max(q_hub.y, q_hub.z)), 0.0) - 0.04;

    // Stationary outer anchor ring at r = 1.95
    float r_xz = length(q.xz);
    vec2 q_ring = abs(vec2(r_xz - 1.95, q.y)) - vec2(0.08, 0.05);
    float d_ring = length(max(q_ring, 0.0)) + min(max(q_ring.x, q_ring.y), 0.0) - 0.02;

    // 4 connecting straps in X and Z
    vec3 aq = abs(q);
    vec3 qx = max(vec3(aq.x - 1.95, aq.y - 0.02, aq.z - 0.12), 0.0);
    float dx = length(qx) - 0.01;
    vec3 qz = max(vec3(aq.z - 1.95, aq.y - 0.02, aq.x - 0.12), 0.0);
    float dz = length(qz) - 0.01;
    float d_straps = min(dx, dz);

    float d_raw = min(d_hub, min(d_ring, d_straps));
    return opDiracBeltScale(p, d_raw, r1, r2, phase);
}

#endif // SDF_GEMS_OP_DIRAC_BELT_GLSL
