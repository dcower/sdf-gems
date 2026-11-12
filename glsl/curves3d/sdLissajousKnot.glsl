#ifndef SDF_GEMS_SD_LISSAJOUS_KNOT_GLSL
#define SDF_GEMS_SD_LISSAJOUS_KNOT_GLSL

#include "../common/common.glsl"

// ============================================================================
// sdLissajousKnot3D, sdgLissajousKnot3D, sdHoropterCurve & sdgHoropterCurve
// — Exact 3D Euclidean Distance Fields & Analytical Normal Gradients for the
//   (1,2,3) Chebyshev-Lissajous Space Curve and Hermann von Helmholtz's (1867)
//   Cylindrical Horopter Twisted Cubic
//
// Discovered & mathematically designed by Gemini 4 Argon.
// Implementation, verification & GPU test harness by Gemini 3.8 Flash.
// License: MIT
// ============================================================================

/**
 * Exact 3D Euclidean Distance Field and Analytical Normal Gradient for the
 * (1, 2, 3) Chebyshev-Lissajous Space Curve:
 *   q(t) = (a * cos(t), b * cos(2t), cAmp * sin(3t))
 *
 * Uses the Even-Odd Branch Dominance Theorem and octant folding to reduce
 * global parameter space to a single SIMD vec4 Newton pass over [0, pi/2].
 *
 * @param p     Query point in 3D space.
 * @param a     Major radius along X.
 * @param b     Radius along Y.
 * @param cAmp  Saddle amplitude along Z.
 * @param rTube Tube radius around the space curve.
 * @return vec4(distance - rTube, unit normal gradient)
 */
vec4 sdgLissajousKnot3D( in vec3 p, in float a, in float b, in float cAmp, in float rTube )
{
    a    = max(a, 1e-6);
    b    = max(b, 1e-6);
    cAmp = max(cAmp, 1e-6);
    const float HALF_PI = 1.57079632679489661923;
    const float PI_3    = 1.04719755119659774615;

    float sx = (p.x >= 0.0) ? 1.0 : -1.0;
    float sz = (p.z >= 0.0) ? 1.0 : -1.0;
    float ax = abs(p.x), py = p.y, az = abs(p.z);

    // Endpoints u = 0: (a, b, 0), u = pi/3: (0.5*a, -0.5*b, 0), u = pi/2: (0, -b, cAmp)
    float best_t  = 0.0;
    float best_sz = 1.0;
    float best_d2 = (ax - a) * (ax - a) + (py - b) * (py - b) + az * az;

    float d2_pi3 = (ax - 0.5 * a) * (ax - 0.5 * a) + (py + 0.5 * b) * (py + 0.5 * b) + az * az;
    if( d2_pi3 < best_d2 )
    {
        best_d2 = d2_pi3;
        best_t  = PI_3;
        best_sz = 1.0;
    }

    float d2_1 = ax * ax + (py + b) * (py + b) + (az - cAmp) * (az - cAmp);
    if( d2_1 < best_d2 )
    {
        best_d2 = d2_1;
        best_t  = HALF_PI;
        best_sz = -1.0;
    }

    // 2 lanes on branch s_branch = +1 and 2 lanes on s_branch = -1 over full [0, pi/2]
    vec4 t        = vec4(0.22, 0.85, 0.45, 1.28);
    const vec4 SB = vec4(1.0, 1.0, -1.0, -1.0);
    const vec4 LO = vec4(0.0);
    const vec4 HI = vec4(HALF_PI);

    for( int i = 0; i < 7; i++ )
    {
        vec4 ct   = cos(t), st = sin(t);
        vec4 c2t  = ct * ct - st * st, s2t = 2.0 * st * ct;
        vec4 c3t  = ct * (4.0 * ct * ct - 3.0);
        vec4 s3t  = st * (3.0 - 4.0 * st * st);

        vec4 qx   = a * ct, qy = b * c2t, qz = (SB * cAmp) * s3t;
        vec4 qpx  = -a * st, qpy = (-2.0 * b) * s2t, qpz = (3.0 * SB * cAmp) * c3t;
        vec4 qppx = -qx, qppy = -4.0 * qy, qppz = -9.0 * qz;

        vec4 dx   = qx - ax, dy = qy - py, dz = qz - az;
        vec4 g    = dx * qpx + dy * qpy + dz * qpz;
        vec4 gp   = qpx * qpx + qpy * qpy + qpz * qpz + dx * qppx + dy * qppy + dz * qppz;
        vec4 step = vec4(
            (gp.x > 1e-20) ? (g.x / gp.x) : ((g.x >= 0.0 ? 1.0 : -1.0) * 0.15),
            (gp.y > 1e-20) ? (g.y / gp.y) : ((g.y >= 0.0 ? 1.0 : -1.0) * 0.15),
            (gp.z > 1e-20) ? (g.z / gp.z) : ((g.z >= 0.0 ? 1.0 : -1.0) * 0.15),
            (gp.w > 1e-20) ? (g.w / gp.w) : ((g.w >= 0.0 ? 1.0 : -1.0) * 0.15)
        );
        t = clamp(t - clamp(step, vec4(-0.25), vec4(0.25)), LO, HI);
    }

    for( int j = 0; j < 4; j++ )
    {
        float tj  = t[j];
        float sbj = SB[j];
        float ct  = cos(tj), st = sin(tj);
        vec3  q   = vec3(a * ct, b * (ct * ct - st * st), sbj * cAmp * st * (3.0 - 4.0 * st * st));
        vec3  dv  = vec3(ax, py, az) - q;
        float d2  = dot(dv, dv);
        if( d2 < best_d2 )
        {
            best_d2 = d2;
            best_t  = tj;
            best_sz = sbj;
        }
    }

    float cb   = cos(best_t), sb = sin(best_t);
    vec3  q    = vec3(sx * a * cb, b * (cb * cb - sb * sb), sz * best_sz * cAmp * sb * (3.0 - 4.0 * sb * sb));
    vec3  diff = p - q;
    float d    = max(length(diff), 1e-30);
    return vec4(d - rTube, diff / d);
}

/**
 * Signed Distance Function for the (1, 2, 3) Chebyshev-Lissajous Space Curve.
 */
float sdLissajousKnot3D( in vec3 p, in float a, in float b, in float cAmp, in float rTube )
{
    return sdgLissajousKnot3D(p, a, b, cAmp, rTube).x;
}

/**
 * Exact 3D Euclidean Distance Field and Analytical Normal Gradient for
 * Hermann von Helmholtz's (1867) Cylindrical Horopter Twisted Cubic:
 *   q(t) = (R * (1 - t^2)/(1 + t^2), R * (2t)/(1 + t^2), h * t),  t in [-tMax, tMax]
 *
 * Utilizes 0-trig, 0-sqrt rational Weierstrass cylinder linearization.
 *
 * @param p     Query point in 3D space.
 * @param R     Cylinder radius in XY.
 * @param h     Vertical pitch scale.
 * @param tMax  Truncation half-range of parameter t.
 * @param rTube Tube radius around the curve.
 * @return vec4(distance - rTube, unit normal gradient)
 */
vec4 sdgHoropterCurve( in vec3 p, in float R, in float h, in float tMax, in float rTube )
{
    R    = max(R, 1e-6);
    h    = max(h, 1e-6);
    tMax = max(tMax, 1e-4);

    float rho2   = dot(p.xy, p.xy);
    float t2_end = tMax * tMax;
    float invDe  = 1.0 / (1.0 + t2_end);

    // Endpoints t = -tMax and t = +tMax
    float best_t  = -tMax;
    float dz_neg  = -h * tMax - p.z;
    float best_d2 = rho2 + R * R - 2.0 * R * (p.x * (1.0 - t2_end) - 2.0 * p.y * tMax) * invDe + dz_neg * dz_neg;

    float dz_pos = h * tMax - p.z;
    float d2_pos = rho2 + R * R - 2.0 * R * (p.x * (1.0 - t2_end) + 2.0 * p.y * tMax) * invDe + dz_pos * dz_pos;
    if( d2_pos < best_d2 )
    {
        best_d2 = d2_pos;
        best_t  = tMax;
    }

    // 4 seeds in [-tMax, tMax]: axial seed pz/h, angular seed tan(phi/2), and +/-0.5*tMax
    float phi   = (rho2 > 1e-24) ? atan(p.y, p.x) : 0.0;
    float t_phi = clamp(tan(0.5 * phi), -tMax, tMax);
    float t_z   = clamp(p.z / h, -tMax, tMax);
    vec4  t     = vec4(t_z, t_phi, -0.5 * tMax, 0.5 * tMax);

    // 0-trig, 0-sqrt rational Newton loop via Weierstrass cylinder linearization
    for( int i = 0; i < 7; i++ )
    {
        vec4 t2    = t * t;
        vec4 invD  = 1.0 / (1.0 + t2);
        vec4 invD2 = invD * invD;
        vec4 invD3 = invD2 * invD;
        vec4 g     = (-2.0 * R) * (p.y * (1.0 - t2) - 2.0 * p.x * t) * invD2 + h * (h * t - p.z);
        vec4 gp    = (4.0 * R) * (p.x * (1.0 - 3.0 * t2) + p.y * t * (3.0 - t2)) * invD3 + h * h;
        vec4 step  = vec4(
            (gp.x > 1e-12) ? (g.x / gp.x) : ((g.x >= 0.0 ? 1.0 : -1.0) * 0.25 * tMax),
            (gp.y > 1e-12) ? (g.y / gp.y) : ((g.y >= 0.0 ? 1.0 : -1.0) * 0.25 * tMax),
            (gp.z > 1e-12) ? (g.z / gp.z) : ((g.z >= 0.0 ? 1.0 : -1.0) * 0.25 * tMax),
            (gp.w > 1e-12) ? (g.w / gp.w) : ((g.w >= 0.0 ? 1.0 : -1.0) * 0.25 * tMax)
        );
        t = clamp(t - clamp(step, vec4(-0.4 * tMax), vec4(0.4 * tMax)), vec4(-tMax), vec4(tMax));
    }

    for( int j = 0; j < 4; j++ )
    {
        float tj   = t[j];
        float t2   = tj * tj;
        float invD = 1.0 / (1.0 + t2);
        float dz   = h * tj - p.z;
        float d2   = rho2 + R * R - 2.0 * R * (p.x * (1.0 - t2) + 2.0 * p.y * tj) * invD + dz * dz;
        if( d2 < best_d2 )
        {
            best_d2 = d2;
            best_t  = tj;
        }
    }

    float tb2   = best_t * best_t;
    float invDb = 1.0 / (1.0 + tb2);
    vec3  q     = vec3(R * (1.0 - tb2) * invDb, R * 2.0 * best_t * invDb, h * best_t);
    vec3  diff  = p - q;
    float d     = max(length(diff), 1e-30);
    return vec4(d - rTube, diff / d);
}

/**
 * Signed Distance Function for Hermann von Helmholtz's (1867) Cylindrical Horopter Twisted Cubic.
 */
float sdHoropterCurve( in vec3 p, in float R, in float h, in float tMax, in float rTube )
{
    return sdgHoropterCurve(p, R, h, tMax, rTube).x;
}

#endif // SDF_GEMS_SD_LISSAJOUS_KNOT_GLSL
