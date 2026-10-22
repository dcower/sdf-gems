// ============================================================================
// sdSphericalSegmentSym, sdgSphericalSegmentSym, sdSphericalSegment & sdgSphericalSegment
// — Exact 0/1/2-SFU Symmetric & Asymmetric Two-Plane Spherical Segments + Analytical Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, `0` SFUs in exterior + interior caps, `1` SFU in sphere, `2` SFUs in rim)
// ============================================================================

#ifndef SDF_GEMS_SD_SPHERICAL_SEGMENT_SYM_GLSL
#define SDF_GEMS_SD_SPHERICAL_SEGMENT_SYM_GLSL

float sdSphericalSegmentSym( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), y = abs(p.y), dy = y - h;
    if( x2 <= (r - h) * (r + h + 2.0 * min(dy, 0.0)) ) return dy;
    float w2 = r * r - h * h;
    if( w2 * y * y < h * h * x2 ) return sqrt(x2 + y * y) - r;
    float d = x2 - w2;
    return sqrt(d * d / (x2 + w2 + 2.0 * sqrt(w2 * x2)) + dy * dy);
}

vec4 sdgSphericalSegmentSym( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), sy = (p.y < 0.0) ? -1.0 : 1.0, y = abs(p.y), dy = y - h;
    if( x2 <= (r - h) * (r + h + 2.0 * min(dy, 0.0)) )
        return vec4(dy, 0.0, sy, 0.0);
    float w2 = r * r - h * h;
    if( w2 * y * y < h * h * x2 )
    {
        float r2 = x2 + y * y, invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR - r, p * invR);
    }
    float d  = x2 - w2, wx = sqrt(w2 * x2);
    float sc = d / max(x2 + wx, 1e-30);
    vec3  v  = vec3(sc * p.x, sy * dy, sc * p.z);
    float v2 = dot(v, v), invV = inversesqrt(max(v2, 1e-30));
    return vec4(v2 * invV, v * invV);
}

float sdSphericalSegment( vec3 p, float r, float h1, float h2 )
{
    float x2 = dot(p.xz, p.xz), sy = (2.0 * p.y < h1 + h2) ? -1.0 : 1.0;
    float hc = (sy < 0.0) ? -h1 : h2, dy = sy * p.y - hc;
    if( x2 <= (r - hc) * (r + hc + 2.0 * min(dy, 0.0)) ) return dy;
    float w1 = r * r - h1 * h1, w2 = r * r - h2 * h2, ys = p.y * abs(p.y);
    bool  c1 = h1 * abs(h1) * x2 > w1 * ys;
    bool  c2 = !c1 && (w2 * ys > h2 * abs(h2) * x2);
    if( !c1 && !c2 ) return sqrt(x2 + p.y * p.y) - r;
    float w = c2 ? w2 : w1, h = c2 ? h2 : h1, d = x2 - w, pyh = p.y - h;
    return sqrt(d * d / max(x2 + w + 2.0 * sqrt(w * x2), 1e-30) + pyh * pyh);
}

vec4 sdgSphericalSegment( vec3 p, float r, float h1, float h2 )
{
    float x2 = dot(p.xz, p.xz), sy = (2.0 * p.y < h1 + h2) ? -1.0 : 1.0;
    float hc = (sy < 0.0) ? -h1 : h2, dy = sy * p.y - hc;
    if( x2 <= (r - hc) * (r + hc + 2.0 * min(dy, 0.0)) )
        return vec4(dy, 0.0, sy, 0.0);
    float w1 = r * r - h1 * h1, w2 = r * r - h2 * h2, ys = p.y * abs(p.y);
    bool  c1 = h1 * abs(h1) * x2 > w1 * ys;
    bool  c2 = !c1 && (w2 * ys > h2 * abs(h2) * x2);
    if( !c1 && !c2 )
    {
        float r2 = x2 + p.y * p.y, invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR - r, p * invR);
    }
    float w  = c2 ? w2 : w1, h = c2 ? h2 : h1, d = x2 - w, wx = sqrt(w * x2);
    float sc = d / max(x2 + wx, 1e-30);
    vec3  v  = vec3(sc * p.x, p.y - h, sc * p.z);
    float v2 = dot(v, v), invV = inversesqrt(max(v2, 1e-30));
    return vec4(v2 * invV, v * invV);
}

#endif // SDF_GEMS_SD_SPHERICAL_SEGMENT_SYM_GLSL
