// ============================================================================
// sdgLink & sdgCapsule — 1/2-Sqrt Chain Link & 1-Sqrt Capsule with 3D Gradients
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, fixes 100% 2D planar strip `0.0 / 0.0 = NaN` in IQ's `wc3yD7`)
// ============================================================================

// 1. 1/2-Sqrt Chain Link + Gradient (1 sqrt in slab |p.y| <= le, 0 NaNs on p.x = 0 strip)
vec4 sdgLink( vec3 p, float le, float r1, float r2 )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0, sy = (p.y < 0.0) ? -1.0 : 1.0;
    float x  = abs(p.x), qy = abs(p.y) - le;
    if( qy <= 0.0 ) {
        float u  = x - r1, d0 = length(vec2(u, p.z));
        vec3  g  = (d0 > 0.0) ? (vec3(sx * u, 0.0, p.z) / d0) : vec3(sx, 0.0, 0.0);
        return vec4(d0 - r2, g);
    }
    float w  = length(vec2(x, qy)), u = w - r1, d0 = length(vec2(u, p.z));
    vec2  du = (w > 0.0) ? (vec2(sx * x, sy * qy) / w) : vec2(0.0, sy);
    vec3  g  = (d0 > 0.0) ? (vec3(u * du, p.z) / d0) : vec3(du, 0.0);
    return vec4(d0 - r2, g);
}

// 2. 1-Sqrt Arbitrary-Segment & Vertical Capsule + Gradient
vec4 sdgCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float h  = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    vec3  q  = pa - h * ba;
    float d  = length(q);
    return vec4(d - r, (d > 0.0) ? (q / d) : vec3(0.0, 1.0, 0.0));
}

vec4 sdgVerticalCapsule( vec3 p, float h, float r )
{
    vec3  q = vec3(p.x, p.y - clamp(p.y, 0.0, h), p.z);
    float d = length(q);
    return vec4(d - r, (d > 0.0) ? (q / d) : vec3(0.0, 1.0, 0.0));
}
