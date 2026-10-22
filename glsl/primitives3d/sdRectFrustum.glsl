// ============================================================================
// sdRectFrustum & sdgRectFrustum — Exact 0/1-Sqrt 3D Rectangular & Hip-Roof Frustum & Gradient
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 per-pixel sqrts in interior, caps, & slanted faces)
// ============================================================================

#ifndef SDF_GEMS_SD_RECT_FRUSTUM_GLSL
#define SDF_GEMS_SD_RECT_FRUSTUM_GLSL

float sdRectFrustum( vec3 p, vec2 b1, vec2 b2, float h )
{
    p.xz = abs(p.xz);
    vec2  bm  = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q   = p - vec3(bm.x, 0.0, bm.y);
    float mx2 = h * h + k.x * k.x, mz2 = h * h + k.y * k.y;
    float dsx = (h * q.x + k.x * q.y) * inversesqrt(mx2);
    float dsz = (h * q.z + k.y * q.y) * inversesqrt(mz2);
    float dc  = abs(q.y) - h;
    float df  = max(max(dsx, dsz), dc);

    float ux  = clamp((h * q.y - k.x * q.x) / mx2, -1.0, 1.0);
    float uz  = clamp((h * q.y - k.y * q.z) / mz2, -1.0, 1.0);
    vec2  bc  = (p.y < 0.0) ? b1 : b2;
    if( df <= 0.0 || (dc > 0.0 && p.x <= bc.x && p.z <= bc.y) ||
        (dsx > 0.0 && abs(h * q.y - k.x * q.x) <= mx2 && q.z + k.y * ux <= 0.0) ||
        (dsz > 0.0 && abs(h * q.y - k.y * q.z) <= mz2 && q.x + k.x * uz <= 0.0) )
        return df;

    vec3  vx = vec3(q.x + k.x * ux, q.y - h * ux, max(q.z + k.y * ux, 0.0));
    vec3  vz = vec3(max(q.x + k.x * uz, 0.0), q.y - h * uz, q.z + k.y * uz);
    float d  = min(dot(vx, vx), dot(vz, vz));
    if( q.z + k.y * ux > 0.0 && q.x + k.x * uz > 0.0 )
    {
        float ue = clamp((h * q.y - dot(k, q.xz)) / (h * h + dot(k, k)), -1.0, 1.0);
        vec3  ve = vec3(q.x + k.x * ue, q.y - h * ue, q.z + k.y * ue);
        d = min(d, dot(ve, ve));
    }
    return sqrt(d);
}

vec4 sdgRectFrustum( vec3 p, vec2 b1, vec2 b2, float h )
{
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p.xz = abs(p.xz);
    vec2  bm    = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q     = p - vec3(bm.x, 0.0, bm.y);
    float mx2   = h * h + k.x * k.x, invMx = inversesqrt(mx2);
    float mz2   = h * h + k.y * k.y, invMz = inversesqrt(mz2);
    float dsx   = (h * q.x + k.x * q.y) * invMx;
    float dsz   = (h * q.z + k.y * q.y) * invMz;
    float dc    = abs(q.y) - h;
    float df    = max(max(dsx, dsz), dc);

    float ux    = clamp((h * q.y - k.x * q.x) / mx2, -1.0, 1.0);
    float uz    = clamp((h * q.y - k.y * q.z) / mz2, -1.0, 1.0);
    vec2  bc    = (p.y < 0.0) ? b1 : b2;
    bool  in_int   = (df <= 0.0);
    bool  in_cap   = (dc > 0.0 && p.x <= bc.x && p.z <= bc.y);
    bool  in_sidex = (dsx > 0.0 && abs(h * q.y - k.x * q.x) <= mx2 && q.z + k.y * ux <= 0.0);
    bool  in_sidez = (dsz > 0.0 && abs(h * q.y - k.y * q.z) <= mz2 && q.x + k.x * uz <= 0.0);
    float d;
    vec3  g;
    if( in_int || in_cap || in_sidex || in_sidez )
    {
        d = df;
        bool use_cap = in_cap || (in_int && dc >= max(dsx, dsz));
        bool use_x   = !use_cap && (in_sidex || (in_int && dsx >= dsz));
        g = use_cap ? vec3(0.0, s0.y, 0.0)
          : (use_x  ? vec3(h * invMx, k.x * invMx, 0.0)
                    : vec3(0.0, k.y * invMz, h * invMz));
    }
    else
    {
        vec3  vx  = vec3(q.x + k.x * ux, q.y - h * ux, max(q.z + k.y * ux, 0.0));
        vec3  vz  = vec3(max(q.x + k.x * uz, 0.0), q.y - h * uz, q.z + k.y * uz);
        float dx2 = dot(vx, vx), dz2 = dot(vz, vz);
        vec3  v   = (dx2 <= dz2) ? vx : vz;
        float d2m = min(dx2, dz2);
        if( q.z + k.y * ux > 0.0 && q.x + k.x * uz > 0.0 )
        {
            float ue  = clamp((h * q.y - dot(k, q.xz)) / (h * h + dot(k, k)), -1.0, 1.0);
            vec3  ve  = vec3(q.x + k.x * ue, q.y - h * ue, q.z + k.y * ue);
            float de2 = dot(ve, ve);
            if( de2 < d2m ) { v = ve; d2m = de2; }
        }
        float invV = inversesqrt(max(d2m, 1e-30));
        d = d2m * invV;
        g = v * invV;
    }
    return vec4(d, s0.x * g.x, g.y, s0.z * g.z);
}

#endif // SDF_GEMS_SD_RECT_FRUSTUM_GLSL
