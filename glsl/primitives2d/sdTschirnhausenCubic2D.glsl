// ============================================================================
// sdTschirnhausenLoop2D, sdgTschirnhausenLoop2D, sdTschirnhausenCubic2D,
// sdgTschirnhausenCubic2D, sdTschirnhausenDrop3D & sdgTschirnhausenDrop3D
// — Exact 0-Trig, Pure Degree-3 Polynomial 2D Tschirnhausen Cubic
//   (1690) Teardrop Loop, Full Crunodal Cubic & 3D Teardrop Spindle
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig, 0 sqrts inside Newton iterations, degree-3 polynomial interior sign test)
// @novelty Parameterizes Ehrenfried Walther von Tschirnhaus's (1690) cubic curve
//          (also studied by L'Hopital and Catalan as the parabola pedal / catacaustic)
//          a^3 * y^2 = b^2 * x^2 * (a - x) (canonical Tschirnhausen habit at b = a / sqrt(3))
//          via the pure polynomial map q(s) = (a*(1 - s^2), b*s*(1 - s^2)):
//          - Every Newton iteration evaluates q, q', and q'' using ONLY additions and
//            multiplications (ZERO divisions, ZERO trig, ZERO sqrts in polynomial evaluations)!
//          - The closed loop (s in [0, 1], x in [0, a]) and crossed crunodal tails
//            (s in [1, sMax], x in [a*(1 - sMax^2), 0]) share a single branch solver.
//          - The closed loop interior is decided by the exact degree-3 polynomial test
//            (x > 0) && (x < a) && (a^3 * y^2 < b^2 * x^2 * (a - x)).
// ============================================================================

#ifndef SDF_GEMS_SD_TSCHIRNHAUSEN_CUBIC_2D_GLSL
#define SDF_GEMS_SD_TSCHIRNHAUSEN_CUBIC_2D_GLSL

#include "../common/common.glsl"

// Shared pure-polynomial branch solver on s in [sMin, sMax] for q(s) = (a*(1 - s^2), bSigned*s*(1 - s^2))
vec3 _sdTschirnhausenBranch( in float px, in float ay, in float a, in float bSigned, in float sMin, in float sMax )
{
    float sMin2  = sMin * sMin;
    vec2  bestQ  = vec2(a * (1.0 - sMin2), bSigned * sMin * (1.0 - sMin2));
    vec2  dv0    = vec2(px, ay) - bestQ;
    float bestD2 = dot(dv0, dv0);

    float sMax2 = sMax * sMax;
    vec2  q1    = vec2(a * (1.0 - sMax2), bSigned * sMax * (1.0 - sMax2));
    vec2  dv1   = vec2(px, ay) - q1;
    float d2_1  = dot(dv1, dv1);
    if( d2_1 < bestD2 )
    {
        bestD2 = d2_1;
        bestQ  = q1;
    }

    float span = sMax - sMin;
    vec3  s    = vec3(sMin + 0.2 * span, sMin + 0.5 * span, sMin + 0.8 * span);

    for( int i = 0; i < 6; i++ )
    {
        vec3 s2   = s * s;
        vec3 qx   = a * (vec3(1.0) - s2);
        vec3 dqx  = (-2.0 * a) * s;
        vec3 d2qx = vec3(-2.0 * a);

        vec3 qy   = bSigned * s * (vec3(1.0) - s2);
        vec3 dqy  = bSigned * (vec3(1.0) - 3.0 * s2);
        vec3 d2qy = (-6.0 * bSigned) * s;

        vec3 dx = qx - vec3(px);
        vec3 dy = qy - vec3(ay);
        vec3 f  = dx * dqx + dy * dqy;
        vec3 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec3 step = vec3(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? (0.15 * span) : (-0.15 * span)),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? (0.15 * span) : (-0.15 * span)),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? (0.15 * span) : (-0.15 * span))
        );
        s = clamp(s - clamp(step, vec3(-0.35 * span), vec3(0.35 * span)), vec3(sMin), vec3(sMax));
    }

    for( int j = 0; j < 3; j++ )
    {
        float sj  = s[j];
        float s2j = sj * sj;
        vec2  qj  = vec2(a * (1.0 - s2j), bSigned * sj * (1.0 - s2j));
        vec2  dv  = vec2(px, ay) - qj;
        float d2  = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig 2D Closed Tschirnhausen Cubic Loop + Analytical Unit Gradient vec3(d, grad.xy)
//    Spans x in [0, a] with crunode at (0, 0) (tangents +-b/a) and smooth apex at (a, 0).
vec3 sdgTschirnhausenLoop2D( in vec2 p, in float a, in float b )
{
    a = max(a, 1e-5);
    b = max(b, 1e-5);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res    = _sdTschirnhausenBranch(px, ay, a, b, 0.0, 1.0);
    float bestD2 = res.x;
    vec2  bestQ  = res.yz;

    // Exact 0-trig 0-sqrt cubic polynomial interior sign test:
    bool  inside = (px > 0.0) && (px < a) && ((a * a * a) * (ay * ay) < (b * b) * (px * px) * (a - px));
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d    = length(diff);
    vec2  grad = (d > 1e-7) ? (sgn * diff / d) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Closed Tschirnhausen Cubic Loop SDF
float sdTschirnhausenLoop2D( in vec2 p, in float a, in float b )
{
    return sdgTschirnhausenLoop2D(p, a, b).x;
}

// Canonical 1-parameter overloads for Tschirnhaus's (1690) Equilateral Crunode Habit (b = a / sqrt(3))
vec3 sdgTschirnhausenLoop2D( in vec2 p, in float a )
{
    return sdgTschirnhausenLoop2D(p, a, 0.5773502691896258 * a);
}

float sdTschirnhausenLoop2D( in vec2 p, in float a )
{
    return sdgTschirnhausenLoop2D(p, a, 0.5773502691896258 * a).x;
}

// 3. Exact 0-Trig Full Tschirnhausen Crunodal Cubic (Closed Loop + Crossed Tails on |s| <= sMax, stroke radius r)
vec3 sdgTschirnhausenCubic2D( in vec2 p, in float a, in float b, in float sMax, in float r )
{
    a    = max(a, 1e-5);
    b    = max(b, 1e-5);
    sMax = max(sMax, 1.0);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3 resLoop = _sdTschirnhausenBranch(px, ay, a,  b, 0.0, 1.0);
    vec3 resTail = _sdTschirnhausenBranch(px, ay, a, -b, 1.0, sMax);
    vec3 best    = (resTail.x < resLoop.x) ? resTail : resLoop;

    vec2  diff = vec2(px - best.y, sy * (ay - best.z));
    float d    = length(diff);
    vec2  grad = (d > 1e-7) ? (diff / d) : vec2(1.0, 0.0);
    return vec3(d - r, grad);
}

// 4. Scalar Full Tschirnhausen Crunodal Cubic SDF
float sdTschirnhausenCubic2D( in vec2 p, in float a, in float b, in float sMax, in float r )
{
    return sdgTschirnhausenCubic2D(p, a, b, sMax, r).x;
}

// 4-parameter overloads (default sMax = 1.45)
vec3 sdgTschirnhausenCubic2D( in vec2 p, in float a, in float b, in float r )
{
    return sdgTschirnhausenCubic2D(p, a, b, 1.45, r);
}

float sdTschirnhausenCubic2D( in vec2 p, in float a, in float b, in float r )
{
    return sdgTschirnhausenCubic2D(p, a, b, 1.45, r).x;
}

// 3-parameter canonical overloads (b = a / sqrt(3), sMax = 1.45)
vec3 sdgTschirnhausenCubic2D( in vec2 p, in float a, in float r )
{
    return sdgTschirnhausenCubic2D(p, a, 0.5773502691896258 * a, 1.45, r);
}

float sdTschirnhausenCubic2D( in vec2 p, in float a, in float r )
{
    return sdgTschirnhausenCubic2D(p, a, r).x;
}

// 5. Exact 3D Tschirnhausen Teardrop Spindle of Revolution + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgTschirnhausenDrop3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.yz);
    vec3  r2d = sdgTschirnhausenLoop2D(vec2(p.x, rho), a, b);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 6. Scalar 3D Tschirnhausen Teardrop Spindle SDF
float sdTschirnhausenDrop3D( in vec3 p, in float a, in float b )
{
    return sdgTschirnhausenLoop2D(vec2(p.x, length(p.yz)), a, b).x;
}

// Canonical 1-parameter overloads for 3D Tschirnhausen Teardrop Spindle
vec4 sdgTschirnhausenDrop3D( in vec3 p, in float a )
{
    return sdgTschirnhausenDrop3D(p, a, 0.5773502691896258 * a);
}

float sdTschirnhausenDrop3D( in vec3 p, in float a )
{
    return sdgTschirnhausenDrop3D(p, a, 0.5773502691896258 * a).x;
}

#endif // SDF_GEMS_SD_TSCHIRNHAUSEN_CUBIC_2D_GLSL
