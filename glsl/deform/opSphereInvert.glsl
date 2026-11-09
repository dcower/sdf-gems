// ============================================================================
// opSphereInvert, sdSphereInvertMetric, sdgSphereInvertFinalize,
// sdInvertPlane, sdInvertSphere, sdDupinCyclide & sdgDupinCyclide
// — Conformal Kelvin-Mobius Spherical Inversion & Dupin Cyclides
// @author Discovered and derived using Claude Opus 5.5; audited and refined
//         (compact-support origin scaling & Householder normal notes) using Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdInvertPlane & sdInvertSphere: Exact Euclidean SDFs (`||grad d|| = 1`);
//         sdSphereInvertMetric & sdDupinCyclide: Globally 1-Lipschitz (`||grad d|| <= 1`,
//         0% overstep) and 1st-Order Eikonal-Exact (`||grad d|| = 1`) on the zero set;
//         sdgSphereInvertFinalize: Exact 0-SFU Householder surface normal (`||g|| = 1`).
// ============================================================================

#ifndef SDF_GEMS_OP_SPHERE_INVERT_GLSL
#define SDF_GEMS_OP_SPHERE_INVERT_GLSL

// 1. Conformal Kelvin-Mobius Spherical Inversion:
//    Maps p in R^3 \ {0} to q = (R^2 / ||p||^2) p in the sphere of radius R.
//    Clamping ||p||^2 >= 1e-16 bounds ||q|| <= 1e8 * R^2 so downstream dot(q, q)
//    evaluations never overflow FP32 (FLT_MAX ~ 3.4e38) at or near the origin.
vec3 opSphereInvert( in vec3 p, in float R )
{
    return (R * R / max(dot(p, p), 1e-16)) * p;
}

// 2. Apollonius-Kelvin Spherical Inversion Metric:
//    Naive Jacobian scaling d_naive = (r^2 / R^2) * dq strictly overestimates Euclidean
//    distance everywhere off the surface and diverges quadratically (~ c_0 r^2 / R^2)
//    as r -> inf, causing 100% raymarching overstep from any camera outside r > 2 R^2 / c_0.
//    By the Kelvin-Mobius chordal identity ||x - p|| = r ||y - q|| / ||y|| and triangle
//    inequality ||y|| <= R^2 / r + |dq|, the exact distance to the inverted Apollonius
//    empty ball Phi(B(q, |dq|)) is achieved radially inward toward the origin:
//      d_Kelvin(dq, r, R) = (r^2 * dq) / (R^2 + r * |dq|).
//    Satisfies ||grad_p d_Kelvin|| <= 1 identically across R^3 (0% overstep),
//    ||grad_p d_Kelvin|| = ||grad_q dq|| = 1 on dq = 0, and d_Kelvin ~ r - R^2 / c_0 as r -> inf.
float sdSphereInvertMetric( in float dq, in float r, in float R )
{
    return (r * r * dq) / max(R * R + r * abs(dq), 1e-30);
}

float sdSphereInvertMetric( in float dq, in vec3 p, in float R )
{
    return sdSphereInvertMetric(dq, length(p), R);
}

// Compact-Support Overload: when the q-space surface S_q is enclosed in ||y|| <= rMaxQ,
// bounding ||y|| <= min(R^2 / r + |dq|, rMaxQ) eliminates the artificial zero at p = 0
// while preserving ||grad_p d|| <= 1 everywhere and exactness on the surface.
float sdSphereInvertMetric( in float dq, in float r, in float R, in float rMaxQ )
{
    return (r > 1e-7) ? (r * dq) / max(min(R * R / r + abs(dq), rMaxQ), 1e-30)
                      : (R * R) / max(rMaxQ, 1e-30);
}

// 3. 0-SFU Householder Conformal Surface Normal Pullback:
//    The Jacobian of q(p) = (R^2 / r^2) p is J(p) = (R^2 / r^2) (I - 2 p_hat p_hat^T).
//    Because H(p_hat) = I - 2 p_hat p_hat^T is an orthogonal Householder reflection
//    (H^T H = I), any unit gradient g_q in q-space pulls back to the exact unit surface
//    normal g_p = g_q - (2 dot(p, g_q) / r^2) p on dq = 0 with 0 sqrt and 0 inversesqrt calls!
//    (Note: off the surface dq != 0, g_p is the conformal Householder pullback of g_q,
//    whereas the spatial gradient of d_Kelvin is s*(1 - w)*p_hat + w*g_p with w = (R^2/(R^2+r*|dq|))^2.)
vec4 sdgSphereInvertFinalize( in vec3 p, in float R, in vec4 dg )
{
    float r2 = dot(p, p);
    float r  = sqrt(r2);
    float d  = sdSphereInvertMetric(dg.x, r, R);
    vec3  gp = (r2 > 1e-16) ? (dg.yzw - (2.0 * dot(p, dg.yzw) / r2) * p) : vec3(-1.0, 0.0, 0.0);
    return vec4(d, gp);
}

vec4 sdgSphereInvertFinalize( in vec3 p, in float R, in float rMaxQ, in vec4 dg )
{
    float r2 = dot(p, p);
    float r  = sqrt(r2);
    float d  = sdSphereInvertMetric(dg.x, r, R, rMaxQ);
    vec3  gp = (r2 > 1e-16) ? (dg.yzw - (2.0 * dot(p, dg.yzw) / r2) * p) : vec3(-1.0, 0.0, 0.0);
    return vec4(d, gp);
}

// 4. Exact Closed-Form Euclidean SDF of an Inverted Plane (n . q - h = 0, ||n|| = 1):
//    Conjugate rationalization of (R^2 - ||2 h p - R^2 n||) / (2 h) cancels 2 h
//    analytically via R^4 - ||2 h p - R^2 n||^2 = 4 h (R^2 (n . p) - h ||p||^2),
//    eliminating 0/0 catastrophic cancellation as curvature h -> 0.
float sdInvertPlane( in vec3 p, in vec3 n, in float h, in float R )
{
    float r2 = dot(p, p);
    float R2 = R * R;
    float w  = length(2.0 * h * p - R2 * n);
    return 2.0 * (R2 * dot(p, n) - h * r2) / max(R2 + w, 1e-30);
}

vec4 sdgInvertPlane( in vec3 p, in vec3 n, in float h, in float R )
{
    float r2 = dot(p, p);
    float R2 = R * R;
    vec3  wv = 2.0 * h * p - R2 * n;
    float w2 = dot(wv, wv);
    float w  = sqrt(w2);
    float d  = 2.0 * (R2 * dot(p, n) - h * r2) / max(R2 + w, 1e-30);
    vec3  g  = (w > 0.0) ? (-wv / w) : n;
    return vec4(d, g);
}

// 5. Exact Closed-Form Euclidean SDF of an Inverted Offset Sphere (||q - c|| - a = 0):
//    With k = ||c||^2 - a^2 (power of the origin with respect to B(c, a)), conjugate
//    rationalization of (||k p - R^2 c|| - R^2 a) / k via the identity
//    ||k p - R^2 c||^2 - R^4 a^2 = k (k ||p||^2 - 2 R^2 (c . p) + R^4) eliminates
//    0/0 cancellation at k -> 0 (where the inverted sphere degenerates into a plane).
float sdInvertSphere( in vec3 p, in vec3 c, in float a, in float R )
{
    float r2 = dot(p, p);
    float R2 = R * R;
    float k  = dot(c, c) - a * a;
    float w  = length(k * p - R2 * c);
    return (k * r2 - 2.0 * R2 * dot(p, c) + R2 * R2) / max(w + R2 * a, 1e-30);
}

vec4 sdgInvertSphere( in vec3 p, in vec3 c, in float a, in float R )
{
    float r2 = dot(p, p);
    float R2 = R * R;
    float k  = dot(c, c) - a * a;
    vec3  wv = k * p - R2 * c;
    float w2 = dot(wv, wv);
    float w  = sqrt(w2);
    float d  = (k * r2 - 2.0 * R2 * dot(p, c) + R2 * R2) / max(w + R2 * a, 1e-30);
    vec3  g  = (w > 0.0) ? (wv / w) : vec3(1.0, 0.0, 0.0);
    return vec4(d, g);
}

// 6. Dupin Cyclide (Ring ra > rb, Horn ra = rb, & Spindle/Apple ra < rb) via
//    Conformal Spherical Inversion of a Translated All-Regime Torus:
float sdDupinCyclide( in vec3 p, in float R, in float ra, in float rb, in vec3 offset )
{
    float r2 = dot(p, p);
    float r  = sqrt(r2);
    vec3  q  = (R * R / max(r2, 1e-16)) * p - offset;
    float x2 = dot(q.xz, q.xz);
    float ay = abs(q.y);
    float b2   = (rb - ra) * (rb + ra);
    float rMaxQ = length(offset) + ra + rb;
    if( b2 > 0.0 )
    {
        float b  = sqrt(b2);
        float dy = ay - b;
        if( dy < 0.0 && b2 * x2 < ra * ra * dy * dy )
        {
            float l2 = x2 + dy * dy;
            return sdSphereInvertMetric(-l2 * inversesqrt(max(l2, 1e-30)), r, R, rMaxQ);
        }
    }
    float invH = inversesqrt(max(x2, 1e-30));
    float h    = (x2 > 0.0) ? (x2 * invH) : 0.0;
    float u    = h - ra, l2 = u * u + q.y * q.y;
    float dq   = l2 * inversesqrt(max(l2, 1e-30)) - rb;
    return sdSphereInvertMetric(dq, r, R, rMaxQ);
}

// 7. Dupin Cyclide Distance + Exact 0-SFU Householder Surface Normal:
vec4 sdgDupinCyclide( in vec3 p, in float R, in float ra, in float rb, in vec3 offset )
{
    float r2    = dot(p, p);
    vec3  q     = (R * R / max(r2, 1e-16)) * p - offset;
    float x2    = dot(q.xz, q.xz);
    float ay    = abs(q.y);
    float sy    = (q.y < 0.0) ? -1.0 : 1.0;
    float b2    = (rb - ra) * (rb + ra);
    float rMaxQ = length(offset) + ra + rb;
    if( b2 > 0.0 )
    {
        float b  = sqrt(b2);
        float dy = ay - b;
        if( dy < 0.0 && b2 * x2 < ra * ra * dy * dy )
        {
            float l2   = x2 + dy * dy;
            float invL = inversesqrt(max(l2, 1e-30));
            vec4  dg0  = vec4(-l2 * invL, (l2 > 0.0) ? (-vec3(q.x, sy * dy, q.z) * invL) : vec3(0.0, sy, 0.0));
            return sdgSphereInvertFinalize(p, R, rMaxQ, dg0);
        }
    }
    float invH = inversesqrt(max(x2, 1e-30));
    float h    = (x2 > 0.0) ? (x2 * invH) : 0.0;
    vec2  er   = (x2 > 0.0) ? (q.xz * invH) : vec2(1.0, 0.0);
    float u    = h - ra, l2 = u * u + q.y * q.y;
    float invL = inversesqrt(max(l2, 1e-30));
    vec4  dg   = vec4(l2 * invL - rb, (l2 > 0.0) ? (vec3(u * er.x, q.y, u * er.y) * invL) : vec3(er.x, 0.0, er.y));
    return sdgSphereInvertFinalize(p, R, rMaxQ, dg);
}

#endif // SDF_GEMS_OP_SPHERE_INVERT_GLSL
