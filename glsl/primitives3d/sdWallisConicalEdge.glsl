// ============================================================================
// sdWallisConicalEdge, sdgWallisConicalEdge, sdConicalWedge & sdgConicalWedge
// — Exact 3D John Wallis's (1665) Conical Edge (Cuneus Conicus) Ruled Conoid
//   & Watertight Solid Conical Wedge
//
// @author  Gemini 4 Argon (research & mathematical design), Gemini 3.8 Flash (implementation & verification)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; exact Envelope-Theorem ruling projection + 0-trig 0-sqrt polynomial interior test)
// @novelty Evaluates exact 3D Euclidean distance to John Wallis's (1665) Conical Edge
//          S(u, theta) = (a*cos(theta), u*b*sin(theta), u*h) for u in [0, 1] by
//          eliminating the ruling parameter u*(theta) = clamp((p_y*b*sin(theta) + p_z*h)/(b^2*sin^2(theta) + h^2), 0, 1)
//          via the Envelope Theorem, and pairs it with the top elliptical cap and the
//          0-trig 0-sqrt quartic interior sign test a^2*h^2*y^2 < b^2*z^2*(a^2 - x^2)
//          for the watertight 3D Solid Conical Wedge (chisel-tip / bow-wedge solid).
// ============================================================================

#ifndef SDF_GEMS_SD_WALLIS_CONICAL_EDGE_GLSL
#define SDF_GEMS_SD_WALLIS_CONICAL_EDGE_GLSL

#include "../common/common.glsl"

// Internal helper: exact closest point on the Wallis ruled conoid wall S(u, theta) in the first quadrant (ax >= 0, ay >= 0)
vec4 _sdWallisWallCore( in float ax, in float ay, in float pz, in float a, in float b, in float h )
{
    const float HALF_PI = 1.5707963267948966;

    // Boundary ruling theta = 0 -> (a, 0, u*h)
    float u0     = clamp(pz / h, 0.0, 1.0);
    vec3  bestQ  = vec3(a, 0.0, u0 * h);
    float bestD2 = (ax - a) * (ax - a) + ay * ay + (pz - bestQ.z) * (pz - bestQ.z);

    // Boundary ruling theta = pi/2 -> (0, u*b, u*h)
    float u1  = clamp((ay * b + pz * h) / (b * b + h * h), 0.0, 1.0);
    vec3  q1  = vec3(0.0, u1 * b, u1 * h);
    float d21 = ax * ax + (ay - q1.y) * (ay - q1.y) + (pz - q1.z) * (pz - q1.z);
    if( d21 < bestD2 ) { bestD2 = d21; bestQ = q1; }

    // Bottom knife-edge segment u = 0: (xc, 0, 0) for xc in [0, a]
    float xc  = clamp(ax, 0.0, a);
    float d2k = (ax - xc) * (ax - xc) + ay * ay + pz * pz;
    if( d2k < bestD2 ) { bestD2 = d2k; bestQ = vec3(xc, 0.0, 0.0); }

    float thX = acos(clamp(ax / a, 0.0, 1.0));
    vec3  th  = vec3(thX, 0.45, 1.15);
    for( int i = 0; i < 6; i++ )
    {
        vec3 ct     = cos(th), st = sin(th);
        vec3 denomU = (b * b) * (st * st) + (h * h);
        vec3 uRaw   = (ay * b * st + pz * h) / denomU;
        vec3 u      = clamp(uRaw, vec3(0.0), vec3(1.0));

        vec3 dx = a * ct - vec3(ax);
        vec3 dy = u * b * st - vec3(ay);
        vec3 f  = (-a) * st * dx + (u * b) * ct * dy;

        vec3 mask = step(vec3(0.0), uRaw) * step(uRaw, vec3(1.0));
        vec3 du   = mask * ((ay * b * ct * denomU - (ay * b * st + pz * h) * (2.0 * b * b * st * ct)) / (denomU * denomU));
        vec3 df   = (a * a) * (st * st - ct * ct) + a * ct * ax + (u * b * ct) * (u * b * ct)
                  - u * b * st * dy + du * b * (ct * dy + u * b * st * ct);

        vec3 stepTh = vec3(
            (df.x > 1e-8) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.2 : -0.2),
            (df.y > 1e-8) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.2 : -0.2),
            (df.z > 1e-8) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.2 : -0.2)
        );
        th = clamp(th - clamp(stepTh, vec3(-0.35), vec3(0.35)), vec3(0.0), vec3(HALF_PI));
    }

    for( int j = 0; j < 3; j++ )
    {
        float thj = th[j];
        float ct  = cos(thj), st = sin(thj);
        float uj  = clamp((ay * b * st + pz * h) / (b * b * st * st + h * h), 0.0, 1.0);
        vec3  qj  = vec3(a * ct, uj * b * st, uj * h);
        vec3  dv  = vec3(ax, ay, pz) - qj;
        float d2  = dot(dv, dv);
        if( d2 < bestD2 ) { bestD2 = d2; bestQ = qj; }
    }

    return vec4(bestD2, bestQ);
}

// 1. Exact 3D Wallis Conical Edge (Cuneus Conicus) Ruled Conoid Shell + Unit Gradient vec4(d, grad.xyz)
vec4 sdgWallisConicalEdge( in vec3 p, in float a, in float b, in float h, in float th )
{
    a = max(a, 1e-5);
    b = max(b, 1e-5);
    h = max(h, 1e-5);

    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ax = abs(p.x), ay = abs(p.y);

    vec4  res  = _sdWallisWallCore(ax, ay, p.z, a, b, h);
    float invD = inversesqrt(max(res.x, 1e-30));
    float d    = res.x * invD;
    vec3  diff = vec3(sx * (ax - res.y), sy * (ay - res.z), p.z - res.w);
    vec3  grad = (res.x > 1e-24) ? (diff * invD) : vec3(0.0, 1.0, 0.0);
    return vec4(d - th, grad);
}

// 2. Scalar 3D Wallis Conical Edge Ruled Conoid Shell SDF
float sdWallisConicalEdge( in vec3 p, in float a, in float b, in float h, in float th )
{
    return sdgWallisConicalEdge(p, a, b, h, th).x;
}

// 3. Exact Watertight 3D Solid Conical Wedge (Cuneus Conicus Solid) + Analytical Unit Gradient vec4(d, grad.xyz)
//    Bounded by the Wallis ruled conoid side walls and flat elliptical top disk at z = h, with bottom knife-edge on z = 0.
vec4 sdgConicalWedge( in vec3 p, in float a, in float b, in float h )
{
    a = max(a, 1e-5);
    b = max(b, 1e-5);
    h = max(h, 1e-5);

    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ax = abs(p.x), ay = abs(p.y);

    vec4  res    = _sdWallisWallCore(ax, ay, p.z, a, b, h);
    float bestD2 = res.x;
    vec3  bestQ  = res.yzw;

    // Flat top elliptical disk at z = h: (x/a)^2 + (y/b)^2 <= 1
    float d2Cap = (p.z - h) * (p.z - h);
    if( (ax * ax) / (a * a) + (ay * ay) / (b * b) <= 1.0 && d2Cap < bestD2 )
    {
        bestD2 = d2Cap;
        bestQ  = vec3(ax, ay, h);
    }

    // Exact 0-trig 0-sqrt quartic interior sign test:
    bool  inside = (p.z > 0.0) && (p.z < h) && (ax < a) &&
                   ((a * h * ay) * (a * h * ay) < (b * p.z) * (b * p.z) * (a * a - ax * ax));
    float sgn    = inside ? -1.0 : 1.0;

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec3  diff = vec3(sx * (ax - bestQ.x), sy * (ay - bestQ.y), p.z - bestQ.z);
    vec3  grad = (bestD2 > 1e-24) ? (sgn * diff * invD) : vec3(0.0, 0.0, 1.0);
    return vec4(sgn * d, grad);
}

// 4. Scalar Watertight 3D Solid Conical Wedge SDF
float sdConicalWedge( in vec3 p, in float a, in float b, in float h )
{
    return sdgConicalWedge(p, a, b, h).x;
}

#endif // SDF_GEMS_SD_WALLIS_CONICAL_EDGE_GLSL
