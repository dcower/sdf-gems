// ============================================================================
// sdCutDisk2D, sdCutHollowDisk2D & sdgCutDisk2D — Unified 0-Sqrt Cap-Voronoi 2D Cut Disk
// @author Inigo Quilez — upgraded via Unified Cap-Voronoi & Signed-Square Reduction using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1` everywhere; 0 sqrts in chord strip & interior parabola, 1 sqrt in disk & rim)
// ============================================================================

#ifndef SDF_GEMS_SD_CUT_DISK2D_GLSL
#define SDF_GEMS_SD_CUT_DISK2D_GLSL

// 1. 0/1-Sqrt 2D Cut Disk (0 sqrts in exterior chord strip + interior parabola, 1 sqrt in disk & corner)
float sdCutDisk2D( in vec2 p, in float r, in float h )
{
    p.x = abs(p.x);
    float x2 = p.x * p.x, dy = h - p.y;
    if( x2 <= (r + h) * (r - h + 2.0 * min(dy, 0.0)) ) return dy;
    float w2 = r * r - h * h;
    if( h * abs(h) * x2 < w2 * p.y * abs(p.y) ) return sqrt(x2 + p.y * p.y) - r;
    return length(vec2(p.x - sqrt(w2), dy));
}

// 2. 1-Sqrt-Shell 2D Cut Hollow Disk
float sdCutHollowDisk2D( in vec2 p, in float r, in float h, in float t )
{
    float w2 = r * r - h * h;
    p.x = abs(p.x);
    float x2 = p.x * p.x;
    return ((h * abs(h) * x2 < w2 * p.y * abs(p.y))
        ? length(vec2(p.x - sqrt(w2), p.y - h))
        : abs(sqrt(x2 + p.y * p.y) - r)) - t;
}

// 3. Analytical 2D Distance and Unit Gradient vec3(d, grad.xy) (0/1-SFU)
vec3 sdgCutDisk2D( in vec2 p, in float r, in float h )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float x2 = p.x * p.x, dy = h - p.y;
    if( x2 <= (r + h) * (r - h + 2.0 * min(dy, 0.0)) )
        return vec3(dy, 0.0, -1.0);
    float w2 = r * r - h * h;
    if( h * abs(h) * x2 < w2 * p.y * abs(p.y) )
    {
        float r2 = x2 + p.y * p.y, invR = inversesqrt(max(r2, 1e-30));
        return vec3(r2 * invR - r, sx * p.x * invR, p.y * invR);
    }
    vec2  v  = vec2(p.x - sqrt(w2), -dy);
    float v2 = dot(v, v), invV = inversesqrt(max(v2, 1e-30));
    return vec3(v2 * invV, sx * v.x * invV, v.y * invV);
}

#endif // SDF_GEMS_SD_CUT_DISK2D_GLSL
