// ============================================================================
// sdVesicaSegment & sdgVesicaSegment — Exact 1-to-2 Sqrt 3D Revolved Vesica / Spindle + Gradient
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt`/`rsq` in tip cones & 2 in body, down from 3 `sqrt`s)
// @requires dot2 (from common.glsl)
// ============================================================================

#ifndef SDF_GEMS_SD_VESICA_SEGMENT_GLSL
#define SDF_GEMS_SD_VESICA_SEGMENT_GLSL

// 1. Exact 1-to-2 Sqrt 3D Revolved Vesica / Spindle
float sdVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3  ba = b - a;
    vec3  pc = p - (a + b) * 0.5;
    float l2 = dot(ba, ba);
    float u  = dot(pc, ba);
    float d  = (l2 - 4.0 * w * w) / (8.0 * w);
    float x2 = dot2(pc - ba * (u / l2));
    float z  = d * (2.0 * abs(u) - l2);
    return (z > 0.0 && l2 * l2 * x2 < z * z)
         ? length(pc - ba * (u >= 0.0 ? 0.5 : -0.5))
         : sqrt(dot2(pc) + d * d + 2.0 * d * sqrt(x2)) - (d + w);
}

// 2. Exact 1-to-2 SFU 3D Revolved Vesica / Spindle Distance + Analytical Gradient (vec4(d, grad.xyz))
vec4 sdgVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3  ba = b - a, pc = p - 0.5 * (a + b);
    float l2 = dot(ba, ba), u = dot(pc, ba);
    float d  = (l2 - 4.0 * w * w) / (8.0 * w);
    vec3  up = pc - ba * (u / l2);
    float x2 = dot(up, up), z = d * (2.0 * abs(u) - l2);
    if( z > 0.0 && l2 * l2 * x2 < z * z )
    {
        vec3  vt   = pc - ba * ((u >= 0.0) ? 0.5 : -0.5);
        float r2   = dot(vt, vt), invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR, (r2 > 0.0) ? (vt * invR) : vec3(0.0, 1.0, 0.0));
    }
    float invX = inversesqrt(max(x2, 1e-24));
    float qx   = (x2 > 1e-24) ? (x2 * invX) : 0.0;
    vec3  er   = (x2 > 1e-24) ? (up * invX)
               : normalize((abs(ba.x) < abs(ba.z)) ? vec3(-ba.y, ba.x, 0.0) : vec3(0.0, -ba.z, ba.y));
    vec3  wv   = pc + d * er;
    float w2   = dot(pc, pc) + d * d + 2.0 * d * qx;
    float invW = inversesqrt(max(w2, 1e-30));
    return vec4(w2 * invW - (d + w), wv * invW);
}

#endif // SDF_GEMS_SD_VESICA_SEGMENT_GLSL
