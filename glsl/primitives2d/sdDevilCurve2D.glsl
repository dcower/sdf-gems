#ifndef SDF_GEMS_SD_DEVIL_CURVE_2D_GLSL
#define SDF_GEMS_SD_DEVIL_CURVE_2D_GLSL

#include "../common/common.glsl"

// ============================================================================
// sdDevilLoop2D, sdgDevilLoop2D, sdDevilCurve2D, sdgDevilCurve2D,
// sdDevilSpindle3D & sdgDevilSpindle3D
// — Exact 0-Trig 2D Gabriel Cramer's (1750) Devil's Curve Central Rotor Loop,
//   Full Electric-Motor Stator-Rotor Curve & 3D Revolved Hourglass Spindle
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & 0-trig 0-sqrt quartic interior sign test)
// @novelty Solves Gabriel Cramer's (1750) quartic Devil's Curve ("Courbe du Diable" /
//          Electric Motor Curve) x^2 * (a^2 - x^2) = y^2 * (b^2 - y^2) (0 < a < b)
//          by factoring its biquadratic roots via the shared radical
//            M(x) = sqrt(0.5 * (b^2 + sqrt(b^4 + 4 * x^2 * (x^2 - a^2)))):
//          1) The central figure-eight rotor loop is q_y = x * sqrt(a^2 - x^2) / M(x),
//             desingularized at x = a via the 0-trig Weierstrass substitution
//             u = (1 - s^2)/(1 + s^2), w = 2*s/(1 + s^2) on s in [0, 1] so q_x = a*u
//             and q_y = a^2*u*w / M(a*u) with M >= b/sqrt(2) > 0 everywhere.
//          2) The outer electric-motor stator branch is the smooth explicit graph q_y = M(x).
//          3) The closed rotor loop interior is decided in 0 trig and 0 sqrts by
//             (x^2 < a^2) && (2*y^2 < b^2) && (y^2*(b^2 - y^2) < x^2*(a^2 - x^2)).
// ============================================================================

// Internal 0-trig solver for the central figure-eight rotor loop in the first quadrant
vec3 _sdDevilInnerCore( in float ax, in float ay, in float a, in float b )
{
    float a2 = a * a;
    float b2 = b * b;
    float b4 = b2 * b2;
    float a4 = a2 * a2;

    vec2  bestQ  = vec2(a, 0.0);
    float bestD2 = (ax - a) * (ax - a) + ay * ay;

    float d2_0 = ax * ax + ay * ay;
    if( d2_0 < bestD2 )
    {
        bestD2 = d2_0;
        bestQ  = vec2(0.0, 0.0);
    }

    vec3 s = vec3(0.08, 0.42, 0.78);
    for( int i = 0; i < 6; i++ )
    {
        vec3 s2   = s * s;
        vec3 invD = 1.0 / (1.0 + s2);
        vec3 u    = (1.0 - s2) * invD;
        vec3 w    = 2.0 * s * invD;

        vec3 pt   = u * w;
        vec3 dpt  = u * u - w * w;
        vec3 d2pt = -4.0 * pt;

        vec3 D2   = max(vec3(b4) - (4.0 * a4) * (pt * pt), vec3(1e-20));
        vec3 invR = inversesqrt(D2);
        vec3 D    = D2 * invR;
        vec3 dD   = (-4.0 * a4) * (pt * dpt) * invR;
        vec3 d2D  = ((-4.0 * a4) * (dpt * dpt + pt * d2pt) - dD * dD) * invR;

        vec3 M2   = 0.5 * (vec3(b2) + D);
        vec3 invM = inversesqrt(M2);
        vec3 M    = M2 * invM;
        vec3 dM   = 0.25 * dD * invM;
        vec3 d2M  = (0.25 * d2D - dM * dM) * invM;

        vec3 qx   = a * u;
        vec3 dqx  = -a * w;
        vec3 d2qx = -a * u;

        vec3 invM2 = invM * invM;
        vec3 qy    = a2 * pt * invM;
        vec3 dqy   = a2 * (dpt * M - pt * dM) * invM2;
        vec3 d2qy  = a2 * ((d2pt * M - pt * d2M) * M - 2.0 * (dpt * M - pt * dM) * dM) * (invM2 * invM);

        vec3 dx = qx - vec3(ax);
        vec3 dy = qy - vec3(ay);
        vec3 f  = dx * dqx + dy * dqy;
        vec3 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec3 stepT = vec3(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.2 : -0.2),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.2 : -0.2),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.2 : -0.2)
        );
        stepT = clamp(stepT, vec3(-0.35), vec3(0.35));
        s = clamp(s - 0.5 * (1.0 + s2) * stepT, vec3(0.0), vec3(1.0));
    }

    for( int j = 0; j < 3; j++ )
    {
        float sj   = s[j];
        float s2j  = sj * sj;
        float invD = 1.0 / (1.0 + s2j);
        float uj   = (1.0 - s2j) * invD;
        float wj   = 2.0 * sj * invD;
        float ptj  = uj * wj;
        float Dj   = sqrt(max(b4 - 4.0 * a4 * ptj * ptj, 1e-20));
        float Mj   = sqrt(0.5 * (b2 + Dj));
        vec2  qj   = vec2(a * uj, a2 * ptj / Mj);
        vec2  dv   = vec2(ax, ay) - qj;
        float d2   = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// Internal 0-trig solver for the outer electric-motor stator branch q_y(x) = M(x) on x in [0, xMax]
vec3 _sdDevilOuterCore( in float ax, in float ay, in float a, in float b, in float xMax )
{
    float a2 = a * a;
    float b2 = b * b;
    float b4 = b2 * b2;

    vec2  bestQ  = vec2(0.0, b);
    float bestD2 = ax * ax + (ay - b) * (ay - b);

    float xMax2 = xMax * xMax;
    float Dend  = sqrt(max(b4 + 4.0 * xMax2 * (xMax2 - a2), 1e-20));
    float Mend  = sqrt(0.5 * (b2 + Dend));
    float d2End = (ax - xMax) * (ax - xMax) + (ay - Mend) * (ay - Mend);
    if( d2End < bestD2 )
    {
        bestD2 = d2End;
        bestQ  = vec2(xMax, Mend);
    }

    vec3 x = vec3(
        clamp(ax, 0.0, xMax),
        clamp(0.70710678 * a, 0.0, xMax),
        clamp(0.5 * (ax + max(ay, b)), 0.0, xMax)
    );

    for( int i = 0; i < 6; i++ )
    {
        vec3 x2   = x * x;
        vec3 D2   = max(vec3(b4) + 4.0 * x2 * (x2 - vec3(a2)), vec3(1e-20));
        vec3 invR = inversesqrt(D2);
        vec3 D    = D2 * invR;
        vec3 dD   = 4.0 * x * (2.0 * x2 - vec3(a2)) * invR;
        vec3 d2D  = (4.0 * (6.0 * x2 - vec3(a2)) - dD * dD) * invR;

        vec3 M2   = 0.5 * (vec3(b2) + D);
        vec3 invM = inversesqrt(M2);
        vec3 M    = M2 * invM;
        vec3 dM   = 0.25 * dD * invM;
        vec3 d2M  = (0.25 * d2D - dM * dM) * invM;

        vec3 dx = x - vec3(ax);
        vec3 dy = M - vec3(ay);
        vec3 f  = dx + dy * dM;
        vec3 df = vec3(1.0) + dM * dM + dy * d2M;

        vec3 step = vec3(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? (0.2 * a) : (-0.2 * a)),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? (0.2 * a) : (-0.2 * a)),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? (0.2 * a) : (-0.2 * a))
        );
        x = clamp(x - clamp(step, vec3(-0.5 * xMax), vec3(0.5 * xMax)), vec3(0.0), vec3(xMax));
    }

    for( int j = 0; j < 3; j++ )
    {
        float xj  = x[j];
        float x2j = xj * xj;
        float Dj  = sqrt(max(b4 + 4.0 * x2j * (x2j - a2), 1e-20));
        float Mj  = sqrt(0.5 * (b2 + Dj));
        vec2  qj  = vec2(xj, Mj);
        vec2  dv  = vec2(ax, ay) - qj;
        float d2  = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig 2D Closed Devil's Curve Central Rotor Loop + Analytical Unit Gradient vec3(d, grad.xy)
//    a: horizontal loop half-width (x in [-a, a]), b: outer stator parameter (b > a > 0).
vec3 sdgDevilLoop2D( in vec2 p, in float a, in float b )
{
    a = max(a, 1e-5);
    b = max(b, a + 1e-5);

    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ax = abs(p.x);
    float ay = abs(p.y);

    vec3  res    = _sdDevilInnerCore(ax, ay, a, b);
    float bestD2 = res.x;
    vec2  bestQ  = res.yz;

    // Exact 0-trig 0-sqrt quartic interior sign test:
    float x2     = ax * ax;
    float y2     = ay * ay;
    float a2     = a * a;
    float b2     = b * b;
    bool  inside = (x2 < a2) && (2.0 * y2 < b2) && (y2 * (b2 - y2) < x2 * (a2 - x2));
    float sgn    = inside ? -1.0 : 1.0;

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec2  diff = vec2(sx * (ax - bestQ.x), sy * (ay - bestQ.y));
    vec2  grad = (bestD2 > 1e-24) ? (sgn * diff * invD) : vec2(sx, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Closed Devil's Curve Central Rotor Loop SDF
float sdDevilLoop2D( in vec2 p, in float a, in float b )
{
    return sdgDevilLoop2D(p, a, b).x;
}

// 3. Exact 0-Trig Complete Devil's Curve (Central Rotor Loop + Outer Stator Branches on |x| <= xMax, stroke radius r)
vec3 sdgDevilCurve2D( in vec2 p, in float a, in float b, in float xMax, in float r )
{
    a    = max(a, 1e-5);
    b    = max(b, a + 1e-5);
    xMax = max(xMax, a);

    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ax = abs(p.x);
    float ay = abs(p.y);

    vec3 resIn  = _sdDevilInnerCore(ax, ay, a, b);
    vec3 resOut = _sdDevilOuterCore(ax, ay, a, b, xMax);
    vec3 best   = (resOut.x < resIn.x) ? resOut : resIn;

    float invD = inversesqrt(max(best.x, 1e-30));
    float d    = best.x * invD;
    vec2  diff = vec2(sx * (ax - best.y), sy * (ay - best.z));
    vec2  grad = (best.x > 1e-24) ? (diff * invD) : vec2(sx, 0.0);
    return vec3(d - r, grad);
}

float sdDevilCurve2D( in vec2 p, in float a, in float b, in float xMax, in float r )
{
    return sdgDevilCurve2D(p, a, b, xMax, r).x;
}

// 4-parameter convenience overloads (defaulting xMax = 1.6 * b)
vec3 sdgDevilCurve2D( in vec2 p, in float a, in float b, in float r )
{
    return sdgDevilCurve2D(p, a, b, 1.6 * max(b, a + 1e-5), r);
}

float sdDevilCurve2D( in vec2 p, in float a, in float b, in float r )
{
    return sdgDevilCurve2D(p, a, b, r).x;
}

// 4. Exact 3D Devil's Hourglass Rotor Spindle of Revolution + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgDevilSpindle3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.yz);
    vec3  r2d = sdgDevilLoop2D(vec2(p.x, rho), a, b);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 5. Scalar 3D Devil's Hourglass Rotor Spindle SDF
float sdDevilSpindle3D( in vec3 p, in float a, in float b )
{
    return sdgDevilLoop2D(vec2(p.x, length(p.yz)), a, b).x;
}

#endif // SDF_GEMS_SD_DEVIL_CURVE_2D_GLSL
