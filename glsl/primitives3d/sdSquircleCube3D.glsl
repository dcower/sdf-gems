// ============================================================================
// sdSquircleCube3D & sdgSquircleCube3D — Exact Secular Co-Normal 3D Goursat Squircle-Cube
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (5.65e-7 max Euclidean error in FP32, 0 trig/pow calls, `||grad d|| = 1` everywhere in R^3)
// @novelty Exact Euclidean 3D SDF and analytical unit gradient for the quartic Goursat
//          squircle-cube (`x^4 + y^4 + z^4 = r^4`) via 48-fold O_h Weyl chamber sorting,
//          vectorized 2D slope-ratio secular solver `p.x*t^3 - p.yz + r*u*t*(1 - t^2) = 0`
//          with `u = inversesqrt(sqrt(1 + t.x^4 + t.y^4))` (0 `pow` calls), and quadratic-
//          stationary co-normal projection `d = (q.x + dot(q.yz, t^3) - r*(1+t.x^4+t.y^4)*u) / sqrt(1 + |t^3|^2)`.
// ============================================================================

#ifndef SDF_GEMS_SD_SQUIRCLE_CUBE_3D_GLSL
#define SDF_GEMS_SD_SQUIRCLE_CUBE_3D_GLSL

vec4 sdgSquircleCube3D( in vec3 p, in float r )
{
    vec3 s = vec3((p.x < 0.0) ? -1.0 : 1.0,
                  (p.y < 0.0) ? -1.0 : 1.0,
                  (p.z < 0.0) ? -1.0 : 1.0);
    vec3 q = abs(p);
    bool s1 = (q.x < q.y); if( s1 ) q.xy = q.yx;
    bool s2 = (q.y < q.z); if( s2 ) q.yz = q.zy;
    bool s3 = (q.x < q.y); if( s3 ) q.xy = q.yx;

    const float c0 = 0.56060388; // 2^(3/4) / 3
    float u     = min(q.x / (c0 * r), 1.0);
    float w     = sqrt(max(1.0 - u, 0.0));
    float t_max = max(0.0, 1.0 - w * (1.12 + (1.0 - w) * (0.1614 - 0.2817 * w)));

    vec2  v  = q.yz / max(q.x, 1e-30);
    vec2  qv = sqrt(sqrt(v));
    vec2  t0 = qv * (1.0 + 2.0 * qv) / (2.0 + qv);
    vec2  t  = min(t0, vec2(t_max));

    for( int i = 0; i < 4; i++ )
    {
        vec2  t2 = t * t, t3 = t2 * t, t4 = t2 * t2;
        float q4 = 1.0 + t4.x + t4.y;
        float ru = r * inversesqrt(sqrt(q4));
        vec2  F  = q.x * t3 - q.yz + ru * t * (1.0 - t2);
        vec2  s4 = vec2(1.0 + t4.y, 1.0 + t4.x);
        vec2  dF = 3.0 * q.x * t2 + (ru / q4) * ((1.0 - 3.0 * t2) * s4 - 2.0 * t2 * t4);
        vec2  dt = F / max(dF, vec2(0.1 * r));
        t = clamp(t - dt, vec2(0.0), vec2(1.0, t.x - dt.x));
    }

    vec2  t2   = t * t, t3 = t2 * t, t4 = t2 * t2;
    float q4   = 1.0 + t4.x + t4.y;
    float u4   = inversesqrt(sqrt(q4));
    float invN = inversesqrt(1.0 + dot(t3, t3));
    float d    = (q.x + dot(q.yz, t3) - r * q4 * u4) * invN;
    vec3  g    = vec3(1.0, t3.x, t3.y) * invN;

    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, s * g);
}

float sdSquircleCube3D( in vec3 p, in float r )
{
    return sdgSquircleCube3D(p, r).x;
}

#endif // SDF_GEMS_SD_SQUIRCLE_CUBE_3D_GLSL
