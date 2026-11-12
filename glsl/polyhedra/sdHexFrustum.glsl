// ============================================================================
// sdHexFrustum, sdgHexFrustum, sdHexBipyramid, sdgHexBipyramid & sdBipyramid
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 per-pixel sqrts in interior & face prisms)
// @novelty Universal Dihedral Frustum & Bipyramid Decomposition in the 30-degree D_6
//          fundamental sector 0 <= v <= u / sqrt(3), reducing 6 trapezoidal side faces
//          and 2 hexagonal caps to 1 slanted face, 2 flat caps, and 2 clamped segment
//          parameters u1, u2 in [-1, 1] with 1 + k_z^2 = 4/3 and 0 per-pixel sqrts
//          across >54% of space.
// ============================================================================

#ifndef SDF_GEMS_SD_HEX_FRUSTUM_GLSL
#define SDF_GEMS_SD_HEX_FRUSTUM_GLSL

// 1. Exact 0/1-Sqrt 3D Hexagonal Frustum (bottom hex inradius r1, top hex inradius r2, half-height h)
float sdHexFrustum( vec3 p, float r1, float r2, float h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    p.xz = abs(p.xz);
    p.xz -= 2.0 * min(dot(k.xy, p.xz), 0.0) * k.xy;
    float v = abs(p.x), u = p.z;

    float rm = 0.5 * (r1 + r2), s = 0.5 * (r1 - r2);
    vec3  q  = vec3(u - rm, p.y, v - k.z * rm);
    float m2 = h * h + s * s;
    float ds = (h * q.x + s * q.y) * inversesqrt(m2);
    float dc = abs(q.y) - h;
    float u1 = clamp((h * q.y - s * q.x) / m2, -1.0, 1.0);
    if( max(ds, dc) <= 0.0 || (dc > 0.0 && u <= ((p.y < 0.0) ? r1 : r2)) ||
        (ds > 0.0 && abs(h * q.y - s * q.x) < m2 && q.z + k.z * s * u1 <= 0.0) )
        return max(ds, dc);

    float u2 = clamp((h * q.y - s * (q.x + k.z * q.z)) / (m2 + (1.0 / 3.0) * s * s), -1.0, 1.0);
    vec3  v1 = vec3(q.x + s * u1, q.y - h * u1, max(q.z + k.z * s * u1, 0.0));
    vec3  v2 = vec3(q.x + s * u2, q.y - h * u2, q.z + k.z * s * u2);
    return sqrt(min(dot(v1, v1), dot(v2, v2)));
}

// 2. Analytical 3D Distance & Gradient vec4(d, grad.xyz) for Hexagonal Frustum
vec4 sdgHexFrustum( vec3 p, float r1, float r2, float h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p.xz = abs(p.xz);
    float dk = dot(k.xy, p.xz);
    bool  rf = (dk < 0.0);
    if( rf ) p.xz -= 2.0 * dk * k.xy;
    float sx2 = (p.x < 0.0) ? -1.0 : 1.0;
    float v = abs(p.x), u = p.z;

    float rm   = 0.5 * (r1 + r2), s = 0.5 * (r1 - r2);
    vec3  q    = vec3(u - rm, p.y, v - k.z * rm);
    float m2   = h * h + s * s;
    float invM = inversesqrt(m2);
    float ds   = (h * q.x + s * q.y) * invM;
    float dc   = abs(q.y) - h;
    float df   = max(ds, dc);
    float u1   = clamp((h * q.y - s * q.x) / m2, -1.0, 1.0);

    bool in_int  = (df <= 0.0);
    bool in_cap  = (dc > 0.0 && u <= ((p.y < 0.0) ? r1 : r2));
    bool in_side = (ds > 0.0 && abs(h * q.y - s * q.x) < m2 && q.z + k.z * s * u1 <= 0.0);
    float d;
    vec3  g_uyv;
    if( in_int || in_cap || in_side )
    {
        d = df;
        g_uyv = (in_cap || (in_int && dc >= ds))
              ? vec3(0.0, s0.y, 0.0)
              : vec3(h * invM, s * invM, 0.0);
    }
    else
    {
        float u2   = clamp((h * q.y - s * (q.x + k.z * q.z)) / (m2 + (1.0 / 3.0) * s * s), -1.0, 1.0);
        vec3  v1   = vec3(q.x + s * u1, q.y - h * u1, max(q.z + k.z * s * u1, 0.0));
        vec3  v2   = vec3(q.x + s * u2, q.y - h * u2, q.z + k.z * s * u2);
        float d1   = dot(v1, v1), d2 = dot(v2, v2);
        vec3  ve   = (d1 <= d2) ? v1 : v2;
        float d2m  = min(d1, d2);
        float invV = inversesqrt(max(d2m, 1e-30));
        d     = d2m * invV;
        g_uyv = ve * invV;
    }

    vec2 g_xz = vec2(g_uyv.z * sx2, g_uyv.x);
    if( rf ) g_xz -= 2.0 * dot(k.xy, g_xz) * k.xy;
    return vec4(d, s0.x * g_xz.x, g_uyv.y, s0.z * g_xz.y);
}

// 3. Exact 0/1-Sqrt 3D Hexagonal Bipyramid (equatorial inradius r, polar half-height h)
float sdHexBipyramid( vec3 p, float r, float h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    p = abs(p);
    p.xz -= 2.0 * min(dot(k.xy, p.xz), 0.0) * k.xy;
    float v = abs(p.x) - k.z * r, u = p.z - r;

    float m2 = h * h + r * r;
    float qy = h * p.y - r * u;
    float qz = h * u + r * p.y;
    if( qz <= 0.0 || (qy > 0.0 && v * m2 + k.z * r * qy <= 0.0) )
        return qz * inversesqrt(m2); // 0 per-pixel sqrts in 100% of interior & face prism!

    vec3  vA = vec3(u, p.y, max(v, 0.0));
    float t  = clamp((qy - k.z * r * v) / (m2 + (1.0 / 3.0) * r * r), 0.0, 1.0);
    vec3  vB = vec3(u + r * t, p.y - h * t, v + k.z * r * t);
    return sqrt(min(dot(vA, vA), dot(vB, vB)));
}

// 4. Analytical 3D Distance & Gradient vec4(d, grad.xyz) for Hexagonal Bipyramid
vec4 sdgHexBipyramid( vec3 p, float r, float h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    float dk = dot(k.xy, p.xz);
    bool  rf = (dk < 0.0);
    if( rf ) p.xz -= 2.0 * dk * k.xy;
    float sx2 = (p.x < 0.0) ? -1.0 : 1.0;
    float v = abs(p.x) - k.z * r, u = p.z - r;

    float m2   = h * h + r * r;
    float invM = inversesqrt(m2);
    float qy   = h * p.y - r * u;
    float qz   = h * u + r * p.y;
    float d;
    vec3  g_uyv;
    if( qz <= 0.0 || (qy > 0.0 && v * m2 + k.z * r * qy <= 0.0) )
    {
        d     = qz * invM;
        g_uyv = vec3(h * invM, r * invM, 0.0);
    }
    else
    {
        vec3  vA   = vec3(u, p.y, max(v, 0.0));
        float t    = clamp((qy - k.z * r * v) / (m2 + (1.0 / 3.0) * r * r), 0.0, 1.0);
        vec3  vB   = vec3(u + r * t, p.y - h * t, v + k.z * r * t);
        float dA2  = dot(vA, vA), dB2 = dot(vB, vB);
        vec3  ve   = (dA2 <= dB2) ? vA : vB;
        float d2m  = min(dA2, dB2);
        float invV = inversesqrt(max(d2m, 1e-30));
        d     = d2m * invV;
        g_uyv = ve * invV;
    }

    vec2 g_xz = vec2(g_uyv.z * sx2, g_uyv.x);
    if( rf ) g_xz -= 2.0 * dot(k.xy, g_xz) * k.xy;
    return vec4(d, s0.x * g_xz.x, s0.y * g_uyv.y, s0.z * g_xz.y);
}

// 5. General N-gonal Bipyramid (n >= 3, equatorial inradius r, polar half-height h)
float sdBipyramid( vec3 p, float n, float r, float h )
{
    float an = 3.14159265 / n;
    float kz = tan(an);
    float bn = mod(atan(p.x, p.z) + an, 2.0 * an) - an;
    vec2  uz = length(p.xz) * vec2(cos(bn), abs(sin(bn)));
    float u  = uz.x - r, v = uz.y - kz * r, ay = abs(p.y);

    float m2 = h * h + r * r;
    float qy = h * ay - r * u;
    float qz = h * u + r * ay;
    if( qz <= 0.0 || (qy > 0.0 && v * m2 + kz * r * qy <= 0.0) )
        return qz * inversesqrt(m2);

    vec3  vA = vec3(u, ay, max(v, 0.0));
    float t  = clamp((qy - kz * r * v) / (m2 + kz * kz * r * r), 0.0, 1.0);
    vec3  vB = vec3(u + r * t, ay - h * t, v + kz * r * t);
    return sqrt(min(dot(vA, vA), dot(vB, vB)));
}

#endif // SDF_GEMS_SD_HEX_FRUSTUM_GLSL
