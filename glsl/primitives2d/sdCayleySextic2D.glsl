// ============================================================================
// sdCayleyOuter2D, sdgCayleyOuter2D, sdCayleySextic2D, sdgCayleySextic2D,
// sdCayleyShell3D & sdgCayleyShell3D
// — Exact 0-Trig 2D Arthur Cayley's (1859) Sextic (Cardioid Pedal Curve) & 3D Shell
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & 0-trig 0-sqrt polynomial interior test)
// @novelty Eliminates all trigonometric functions from Cayley's (1859) sextic
//          r(theta) = 4*a*cos^3(theta/3), 4*(x^2 + y^2 - a*x)^3 = 27*a^2*(x^2 + y^2)^2,
//          via the rational tangent substitution u = tan(theta/3):
//            q(u) = 4*a / (1 + u^2)^3 * (1 - 3*u^2, u*(3 - u^2)),
//          with exact arc-length derivative ||q'(u)|| = 12*a / (1 + u^2)^2 and the
//          algebraic circle-separator identity r^2 - a*x = 12*a^2 / (1 + u^2)^2, which
//          cleanly partitions the outer kidney loop u in [0, sqrt(3)] from the compact
//          inverted inner teardrop loop v = sqrt(3)/u in [0, 1].
// ============================================================================

#ifndef SDF_GEMS_SD_CAYLEY_SEXTIC_2D_GLSL
#define SDF_GEMS_SD_CAYLEY_SEXTIC_2D_GLSL

#include "../common/common.glsl"

// Internal 0-trig solver for the Outer Loop u in [0, sqrt(3)]
vec3 _sdCayleyOuterCore( in float px, in float ay, in float a )
{
    const float s3 = 1.7320508075688772;
    vec2  bestQ  = vec2(4.0 * a, 0.0);
    float bestD2 = (px - 4.0 * a) * (px - 4.0 * a) + ay * ay;

    float d2Crunode = (px + 0.5 * a) * (px + 0.5 * a) + ay * ay;
    if( d2Crunode < bestD2 )
    {
        bestD2 = d2Crunode;
        bestQ  = vec2(-0.5 * a, 0.0);
    }

    // Exact algebraic seed from r^2 - a*x = 12*a^2 / (1 + u^2)^2
    float r2p  = px * px + ay * ay;
    float wp   = max(r2p - a * px, 0.75 * a * a);
    float uAlg = clamp(sqrt(max(2.0 * s3 * a * inversesqrt(wp) - 1.0, 0.0)), 0.0, s3);

    vec4 u = vec4(uAlg, 0.3, 0.9, 1.5);
    for( int i = 0; i < 6; i++ )
    {
        vec4 u2    = u * u;
        vec4 invW  = 1.0 / (1.0 + u2);
        vec4 invW3 = invW * invW * invW;
        vec4 invW4 = invW3 * invW;
        vec4 invW5 = invW4 * invW;

        vec4 qx   = (4.0 * a) * (1.0 - 3.0 * u2) * invW3;
        vec4 qy   = (4.0 * a) * u * (3.0 - u2) * invW3;
        vec4 dqx  = (48.0 * a) * u * (u2 - 1.0) * invW4;
        vec4 dqy  = (12.0 * a) * (1.0 - 6.0 * u2 + u2 * u2) * invW4;
        vec4 d2qx = (48.0 * a) * (-1.0 + 10.0 * u2 - 5.0 * u2 * u2) * invW5;
        vec4 d2qy = (48.0 * a) * u * (-5.0 + 10.0 * u2 - u2 * u2) * invW5;

        vec4 dx = qx - vec4(px);
        vec4 dy = qy - vec4(ay);
        vec4 f  = dx * dqx + dy * dqy;
        vec4 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec4 step = vec4(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.2 : -0.2),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.2 : -0.2),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.2 : -0.2),
            (df.w > 1e-12) ? (f.w / df.w) : ((f.w >= 0.0) ? 0.2 : -0.2)
        );
        u = clamp(u - clamp(step, vec4(-0.35), vec4(0.35)), vec4(0.0), vec4(s3));
    }

    for( int j = 0; j < 4; j++ )
    {
        float uj    = u[j];
        float u2j   = uj * uj;
        float invW  = 1.0 / (1.0 + u2j);
        float invW3 = invW * invW * invW;
        vec2  qj    = vec2(4.0 * a * (1.0 - 3.0 * u2j) * invW3, 4.0 * a * uj * (3.0 - u2j) * invW3);
        vec2  dv    = vec2(px, ay) - qj;
        float d2    = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// Internal 0-trig solver for the Inner Nodal Teardrop Loop v = sqrt(3)/u in [0, 1]
vec3 _sdCayleyInnerCore( in float px, in float ay, in float a )
{
    const float s3 = 1.7320508075688772;
    vec2  bestQ  = vec2(0.0, 0.0);
    float bestD2 = px * px + ay * ay;

    float d2Crunode = (px + 0.5 * a) * (px + 0.5 * a) + ay * ay;
    if( d2Crunode < bestD2 )
    {
        bestD2 = d2Crunode;
        bestQ  = vec2(-0.5 * a, 0.0);
    }

    vec2 v = vec2(0.35, 0.75);
    for( int i = 0; i < 6; i++ )
    {
        vec2 v2    = v * v;
        vec2 v3    = v2 * v;
        vec2 v4    = v2 * v2;
        vec2 v6    = v4 * v2;
        vec2 invD  = 1.0 / (v2 + 3.0);
        vec2 invD3 = invD * invD * invD;
        vec2 invD4 = invD3 * invD;
        vec2 invD5 = invD4 * invD;

        vec2 qx   = (4.0 * a) * v4 * (v2 - 9.0) * invD3;
        vec2 qy   = (12.0 * s3 * a) * v3 * (1.0 - v2) * invD3;
        vec2 dqx  = (144.0 * a) * v3 * (v2 - 3.0) * invD4;
        vec2 dqy  = (12.0 * s3 * a) * v2 * (9.0 - 18.0 * v2 + v4) * invD4;
        vec2 d2qx = (432.0 * a) * v2 * (-9.0 + 10.0 * v2 - v4) * invD5;
        vec2 d2qy = (24.0 * s3 * a) * v * (27.0 - 135.0 * v2 + 45.0 * v4 - v6) * invD5;

        vec2 dx = qx - vec2(px);
        vec2 dy = qy - vec2(ay);
        vec2 f  = dx * dqx + dy * dqy;
        vec2 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec2 step = vec2(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.15 : -0.15),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.15 : -0.15)
        );
        v = clamp(v - clamp(step, vec2(-0.25), vec2(0.25)), vec2(0.0), vec2(1.0));
    }

    for( int j = 0; j < 2; j++ )
    {
        float vj    = v[j];
        float v2j   = vj * vj;
        float v3j   = v2j * vj;
        float v4j   = v2j * v2j;
        float invD  = 1.0 / (v2j + 3.0);
        float invD3 = invD * invD * invD;
        vec2  qj    = vec2(4.0 * a * v4j * (v2j - 9.0) * invD3, 12.0 * s3 * a * v3j * (1.0 - v2j) * invD3);
        vec2  dv    = vec2(px, ay) - qj;
        float d2    = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig Closed Outer Cayley Sextic Solid + Analytical Unit Gradient vec3(d, grad.xy)
vec3 sdgCayleyOuter2D( in vec2 p, in float a )
{
    a = max(a, 1e-5);
    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res    = _sdCayleyOuterCore(px, ay, a);
    float bestD2 = res.x;
    vec2  bestQ  = res.yz;

    // Exact 0-trig 0-sqrt polynomial interior test via r^2 - a*x = 12*a^2/(1+u^2)^2 circle separator:
    float r2     = px * px + ay * ay;
    float w      = r2 - a * px;
    bool  inside = (4.0 * w * w * w < 27.0 * a * a * r2 * r2) || ((px > -0.5 * a) && (w <= 0.75 * a * a));
    float sgn    = inside ? -1.0 : 1.0;

    vec2 diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d   = length(diff);
    vec2 grad = (d > 1e-12) ? (sgn * (diff / d)) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Closed Outer Cayley Sextic SDF
float sdCayleyOuter2D( in vec2 p, in float a )
{
    return sdgCayleyOuter2D(p, a).x;
}

// 3. Exact 0-Trig Full Self-Intersecting Cayley Sextic Curve (Outer Loop + Inner Teardrop Loop, stroke radius r)
vec3 sdgCayleySextic2D( in vec2 p, in float a, in float r )
{
    a = max(a, 1e-5);
    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3 resOuter = _sdCayleyOuterCore(px, ay, a);
    vec3 resInner = _sdCayleyInnerCore(px, ay, a);
    vec3 best     = (resInner.x < resOuter.x) ? resInner : resOuter;

    vec2 diff = vec2(px - best.y, sy * (ay - best.z));
    float d   = length(diff);
    vec2 grad = (d > 1e-12) ? (diff / d) : vec2(1.0, 0.0);
    return vec3(d - r, grad);
}

// 4. Scalar Full Cayley Sextic Curve SDF (with stroke radius r)
float sdCayleySextic2D( in vec2 p, in float a, in float r )
{
    return sdgCayleySextic2D(p, a, r).x;
}

// 5. Zero-stroke Overload for Full Cayley Sextic Curve SDF
vec3 sdgCayleySextic2D( in vec2 p, in float a )
{
    return sdgCayleySextic2D(p, a, 0.0);
}

float sdCayleySextic2D( in vec2 p, in float a )
{
    return sdCayleySextic2D(p, a, 0.0);
}

// 6. Exact 3D Revolved Cayley Shell + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgCayleyShell3D( in vec3 p, in float a )
{
    float rho = length(p.yz);
    vec3  r2d = sdgCayleyOuter2D(vec2(p.x, rho), a);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 7. Scalar 3D Revolved Cayley Shell SDF
float sdCayleyShell3D( in vec3 p, in float a )
{
    return sdgCayleyOuter2D(vec2(p.x, length(p.yz)), a).x;
}

#endif // SDF_GEMS_SD_CAYLEY_SEXTIC_2D_GLSL
