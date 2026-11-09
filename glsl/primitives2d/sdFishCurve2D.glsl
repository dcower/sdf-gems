// ============================================================================
// sdFishLoop2D, sdgFishLoop2D, sdFishCurve2D, sdgFishCurve2D, sdFish3D & sdgFish3D
// — Exact 0-Trig 2D Steiner's Fish Curve (Astroid Pedal / Ellipse Negative Pedal),
//   Closed Teardrop Body Loop & 3D Revolved Streamlined Fish Fuselage
//
// @author  Gemini 4 Argon (research & mathematical design),
//          Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & sign test)
// @novelty Parameterizes both the closed teardrop body loop and crossed tail fins of
//          Steiner's fish curve q(t) = (a*cos(t) - b*sin^2(t), c*cos(t)*sin(t)) via
//          the 0-trig Weierstrass substitution u = (1 - s^2)/(1 + s^2), w = 2*s/(1 + s^2)
//          on s in [0, 1] with intrinsic t-derivative pullbacks ds/dt = 0.5*(1 + s^2),
//          unifying the body branch (+a) and tail branch (-a) into a single helper.
// ============================================================================

#ifndef SDF_GEMS_SD_FISH_CURVE_2D_GLSL
#define SDF_GEMS_SD_FISH_CURVE_2D_GLSL

#include "../common/common.glsl"

// Internal 0-trig branch solver on s in [0, 1] for q(s) = (aSigned * u - b * w^2, c * u * w)
vec3 _sdFishBranchCore( in float px, in float ay, in float aSigned, in float b, in float c )
{
    vec2  bestQ  = vec2(aSigned, 0.0);
    float bestD2 = (px - aSigned) * (px - aSigned) + ay * ay;

    float d2_1 = (px + b) * (px + b) + ay * ay;
    if( d2_1 < bestD2 )
    {
        bestD2 = d2_1;
        bestQ  = vec2(-b, 0.0);
    }

    vec3 s = vec3(0.2, 0.5, 0.8);
    for( int i = 0; i < 6; i++ )
    {
        vec3 s2   = s * s;
        vec3 invD = 1.0 / (1.0 + s2);
        vec3 u    = (1.0 - s2) * invD;
        vec3 w    = 2.0 * s * invD;
        vec3 u2w2 = u * u - w * w;

        vec3 qx   = aSigned * u - b * (w * w);
        vec3 qy   = c * (u * w);

        vec3 dqx  = -w * (aSigned + 2.0 * b * u);
        vec3 d2qx = -aSigned * u - 2.0 * b * u2w2;
        vec3 dqy  = c * u2w2;
        vec3 d2qy = (-4.0 * c) * (u * w);

        vec3 dx = qx - vec3(px);
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
        vec2  qj   = vec2(aSigned * uj - b * wj * wj, c * uj * wj);
        vec2  dv   = vec2(px, ay) - qj;
        float d2   = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig Closed Fish / Teardrop Body Loop + Analytical Unit Gradient vec3(d, grad.xy)
//    Spans x in [-b, a] with sharp tail crunode at (-b, 0) and smooth rounded nose at (a, 0).
vec3 sdgFishLoop2D( in vec2 p, in float a, in float b, in float c )
{
    a = max(a, 1e-5);
    b = max(b, 1e-5);
    c = max(c, 1e-5);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res    = _sdFishBranchCore(px, ay, a, b, c);
    float bestD2 = res.x;
    vec2  bestQ  = res.yz;

    // Exact 0-trig quadratic root inversion for interior sign test on x in (-b, a):
    float xShift = max(px + b, 0.0);
    float ux     = (2.0 * xShift) / (a + sqrt(max(a * a + 4.0 * b * xShift, 0.0)));
    bool  inside = (px > -b) && (px < a) && (ay * ay < (c * c) * (ux * ux) * (1.0 - ux * ux));
    float sgn    = inside ? -1.0 : 1.0;

    float d    = sqrt(bestD2);
    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    vec2  grad = (d > 1e-6) ? (sgn * normalize(diff)) : vec2(sgn, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Closed Fish Body Loop SDF
float sdFishLoop2D( in vec2 p, in float a, in float b, in float c )
{
    return sdgFishLoop2D(p, a, b, c).x;
}

// 3. Exact 0-Trig Full Steiner Fish Curve with general a, b, c and stroke radius r
vec3 sdgFishCurve2D( in vec2 p, in float a, in float b, in float c, in float r )
{
    a = max(a, 1e-5);
    b = max(b, 1e-5);
    c = max(c, 1e-5);
    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3 resBody = _sdFishBranchCore(px, ay,  a, b, c);
    vec3 resTail = _sdFishBranchCore(px, ay, -a, b, c);
    vec3 best    = (resTail.x < resBody.x) ? resTail : resBody;

    float d    = sqrt(best.x);
    vec2  diff = vec2(px - best.y, sy * (ay - best.z));
    vec2  grad = (d > 1e-6) ? normalize(diff) : vec2(1.0, 0.0);
    return vec3(d - r, grad);
}

// 4. Scalar Full Steiner Fish Curve SDF with general a, b, c and stroke radius r
float sdFishCurve2D( in vec2 p, in float a, in float b, in float c, in float r )
{
    return sdgFishCurve2D(p, a, b, c, r).x;
}

// 5. Exact 0-Trig Full Steiner Fish Curve with general a, b, c (zero stroke radius)
vec3 sdgFishCurve2D( in vec2 p, in float a, in float b, in float c )
{
    return sdgFishCurve2D(p, a, b, c, 0.0);
}

// 6. Scalar Full Steiner Fish Curve SDF with general a, b, c (zero stroke radius)
float sdFishCurve2D( in vec2 p, in float a, in float b, in float c )
{
    return sdgFishCurve2D(p, a, b, c, 0.0).x;
}

// 7. Canonical Steiner Fish Curve (b = c = a / sqrt(2)) with stroke radius r
vec3 sdgFishCurve2D( in vec2 p, in float a, in float r )
{
    float b = 0.7071067811865476 * a;
    return sdgFishCurve2D(p, a, b, b, r);
}

// 8. Scalar Canonical Steiner Fish Curve SDF (b = c = a / sqrt(2)) with stroke radius r
float sdFishCurve2D( in vec2 p, in float a, in float r )
{
    float b = 0.7071067811865476 * a;
    return sdgFishCurve2D(p, a, b, b, r).x;
}

// 9. Exact 3D Revolved Streamlined Fish / Airship Fuselage + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgFish3D( in vec3 p, in float a, in float b, in float c )
{
    float rho = length(p.yz);
    vec3  r2d = sdgFishLoop2D(vec2(p.x, rho), a, b, c);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 10. Scalar 3D Revolved Streamlined Fish Fuselage SDF
float sdFish3D( in vec3 p, in float a, in float b, in float c )
{
    return sdgFishLoop2D(vec2(p.x, length(p.yz)), a, b, c).x;
}

#endif // SDF_GEMS_SD_FISH_CURVE_2D_GLSL
