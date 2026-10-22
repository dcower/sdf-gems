// ============================================================================
// sdPyramidFrustum & sdgPyramidFrustum — Exact 0/1-Sqrt 3D Square Frustum & Analytical Gradient
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 per-pixel sqrts in interior, caps, & slanted faces)
// ============================================================================

#ifndef SDF_GEMS_SD_PYRAMID_FRUSTUM_GLSL
#define SDF_GEMS_SD_PYRAMID_FRUSTUM_GLSL

float sdPyramidFrustum( vec3 p, float b1, float b2, float h )
{
    p.xz = abs(p.xz);
    p.xz = (p.z > p.x) ? p.zx : p.xz;

    float bm = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q  = p - vec3(bm, 0.0, bm);
    float m2 = h * h + k * k;
    float ds = (h * q.x + k * q.y) * inversesqrt(m2);
    float dc = abs(q.y) - h;
    float u1 = clamp((h * q.y - k * q.x) / m2, -1.0, 1.0);
    if( max(ds, dc) <= 0.0 || (dc > 0.0 && p.x <= ((p.y < 0.0) ? b1 : b2)) ||
        (ds > 0.0 && abs(h * q.y - k * q.x) <= m2 && q.z + k * u1 <= 0.0) )
        return max(ds, dc);

    vec3 d1 = vec3(q.x + k * u1, q.y - h * u1, max(q.z + k * u1, 0.0));
    if( q.z + k * u1 <= 0.0 )
        return length(d1.xy);

    float u2 = clamp((h * q.y - k * (q.x + q.z)) / (m2 + k * k), -1.0, 1.0);
    vec3  d2 = vec3(q.x + k * u2, q.y - h * u2, q.z + k * u2);
    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}

vec4 sdgPyramidFrustum( vec3 p, float b1, float b2, float h )
{
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p.xz = abs(p.xz);
    bool sw = (p.z > p.x);
    if( sw ) p.xz = p.zx;

    float bm   = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q    = p - vec3(bm, 0.0, bm);
    float m2   = h * h + k * k;
    float invM = inversesqrt(m2);
    float ds   = (h * q.x + k * q.y) * invM;
    float dc   = abs(q.y) - h;
    float df   = max(ds, dc);
    float u1   = clamp((h * q.y - k * q.x) / m2, -1.0, 1.0);

    bool in_int  = (df <= 0.0);
    bool in_cap  = (dc > 0.0 && p.x <= ((p.y < 0.0) ? b1 : b2));
    bool in_side = (ds > 0.0 && abs(h * q.y - k * q.x) <= m2 && q.z + k * u1 <= 0.0);
    float d;
    vec3  g;
    if( in_int || in_cap || in_side )
    {
        d = df;
        g = (in_cap || (in_int && dc >= ds))
          ? vec3(0.0, s0.y, 0.0)
          : vec3(h * invM, k * invM, 0.0);
    }
    else
    {
        float u2   = clamp((h * q.y - k * (q.x + q.z)) / (m2 + k * k), -1.0, 1.0);
        vec3  v1   = vec3(q.x + k * u1, q.y - h * u1, max(q.z + k * u1, 0.0));
        vec3  v2   = vec3(q.x + k * u2, q.y - h * u2, q.z + k * u2);
        float d1   = dot(v1, v1), d2 = dot(v2, v2);
        bool  use1 = (q.z + k * u1 <= 0.0 || d1 <= d2);
        vec3  ve   = use1 ? v1 : v2;
        float d2m  = use1 ? d1 : d2;
        float invV = inversesqrt(max(d2m, 1e-30));
        d = d2m * invV;
        g = ve * invV;
    }

    if( sw ) g.xz = g.zx;
    return vec4(d, s0.x * g.x, g.y, s0.z * g.z);
}

#endif // SDF_GEMS_SD_PYRAMID_FRUSTUM_GLSL
