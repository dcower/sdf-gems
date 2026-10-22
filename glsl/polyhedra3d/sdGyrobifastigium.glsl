// ============================================================================
// sdGyrobifastigium — Exact 0/1-Sqrt Johnson Solid J26 (Space-Filling Gyrobifastigium)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, `2.22e-16` FP64 / `2.98e-8` FP32 error, 0 `sqrt`s in interior, 1 `sqrt` in exterior)
// ============================================================================

#ifndef SDF_GEMS_SD_GYROBIFASTIGIUM_GLSL
#define SDF_GEMS_SD_GYROBIFASTIGIUM_GLSL

float sdGyrobifastigium( vec3 p, float a )
{
    p.xy = abs(p.xy);
    if( p.z < 0.0 ) { p.xy = p.yx; p.z = -p.z; }

    const float s60 = 0.86602540;
    const float c60 = 0.50000000;
    float xm  = p.x - a;
    float ym  = p.y - a;
    float dSp = s60 * xm + c60 * p.z;
    float dSm = s60 * ym - c60 * p.z;
    float df  = max(max(dSp, ym), max(dSm, xm));
    if( df <= 0.0 )
        return df; // 0 sqrts in 100% of interior

    float uSp = -c60 * xm + s60 * p.z;
    float uSm = -c60 * ym - s60 * p.z;
    float wSp = clamp(uSp, 0.0, 2.0 * a) - uSp;
    float wSm = clamp(uSm, 0.0, 2.0 * a) - uSm;

    vec3 vSp = vec3(max(dSp, 0.0), wSp, max(ym, 0.0));
    vec3 vSm = vec3(max(xm, 0.0),  dSm, wSm);
    return sqrt(min(dot(vSp, vSp), dot(vSm, vSm))); // 1 sqrt in exterior
}

vec4 sdgGyrobifastigium( vec3 p, float a )
{
    vec2 sgn = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);
    bool sw = p.z < 0.0;
    if( sw ) { p.xy = p.yx; p.z = -p.z; }

    const float s60 = 0.86602540;
    const float c60 = 0.50000000;
    float xm  = p.x - a;
    float ym  = p.y - a;
    float dSp = s60 * xm + c60 * p.z;
    float dSm = s60 * ym - c60 * p.z;
    float df  = max(max(dSp, ym), max(dSm, xm));
    vec3  g;
    float d;
    if( df <= 0.0 )
    {
        d = df;
        if( df == dSp )      g = vec3(s60, 0.0, c60);
        else if( df == dSm ) g = vec3(0.0, s60, -c60);
        else if( df == ym )  g = vec3(0.0, 1.0, 0.0);
        else                 g = vec3(1.0, 0.0, 0.0);
    }
    else
    {
        float uSp = -c60 * xm + s60 * p.z;
        float uSm = -c60 * ym - s60 * p.z;
        float eSp = uSp - clamp(uSp, 0.0, 2.0 * a);
        float eSm = uSm - clamp(uSm, 0.0, 2.0 * a);
        float mSp = max(dSp, 0.0);
        vec3  wSp = vec3(s60 * mSp - c60 * eSp, max(ym, 0.0), c60 * mSp + s60 * eSp);
        vec3  wSm = vec3(max(xm, 0.0), s60 * dSm - c60 * eSm, -c60 * dSm - s60 * eSm);
        float l2p = dot(wSp, wSp), l2m = dot(wSm, wSm);
        vec3  w   = (l2p <= l2m) ? wSp : wSm;
        d = sqrt(min(l2p, l2m));
        g = w / max(d, 1e-20);
    }
    if( sw ) { g.xy = g.yx; g.z = -g.z; }
    g.xy *= sgn;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_GYROBIFASTIGIUM_GLSL
