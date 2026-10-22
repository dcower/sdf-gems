// ============================================================================
// sdgCappedTorus — 1/2-SFU Singularity-Free 3D Capped Torus Distance & Gradient
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (1 inversesqrt in spherical cap sector, 2 in toroidal arc, 0 NaNs)
// ============================================================================

#ifndef SDF_GEMS_SDG_CAPPED_TORUS_GLSL
#define SDF_GEMS_SDG_CAPPED_TORUS_GLSL

// ----------------------------------------------------------------------------
// sdgCappedTorus(p, sc, ra, rb)
//
// Exact 3D signed distance and analytical unit gradient vec4(d, gx, gy, gz)
// for a capped torus of major radius ra, minor tube radius rb, and half-angle
// sc = vec2(sin(theta), cos(theta)).
//
// Key Mathematical Properties:
// 1. In the spherical cap sector (sc.y * |p.x| > sc.x * p.y), the closest point
//    on the spine is the cap endpoint c = vec3(ra * sc.x, ra * sc.y, 0.0).
//    The offset vector w = vec3(|p.x| - ra*sc.x, p.y - ra*sc.y, p.z) gives BOTH
//    the distance and the exact 3D unit gradient in a SINGLE inversesqrt(dot(w,w))!
// 2. In the toroidal arc sector (sc.y * |p.x| <= sc.x * p.y), uses 2 guarded
//    inversesqrt calls with exact analytical fallbacks at the z-axis (r_xy == 0)
//    and the medial spine circle (r_xy == ra && p.z == 0), achieving 0 NaNs
//    everywhere in R^3.
// ----------------------------------------------------------------------------
vec4 sdgCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    float sx  = (p.x < 0.0) ? -1.0 : 1.0;
    p.x       = abs(p.x);

    // Spherical end-cap region: 1 inversesqrt total!
    if( sc.y * p.x > sc.x * p.y )
    {
        vec3  w  = vec3( p.x - ra * sc.x, p.y - ra * sc.y, p.z );
        float w2 = dot(w, w);
        if( w2 > 0.0 )
        {
            float invW = inversesqrt(w2);
            return vec4( w2 * invW - rb, sx * w.x * invW, w.y * invW, w.z * invW );
        }
        return vec4( -rb, sx * sc.y, -sc.x, 0.0 );
    }

    // Toroidal arc region: 2 guarded inversesqrt calls, 0 axis/spine NaNs
    float r2   = dot(p.xy, p.xy);
    float invR = (r2 > 0.0) ? inversesqrt(r2) : 0.0;
    float r_xy = r2 * invR;
    vec2  er   = (r2 > 0.0) ? vec2(sx * p.x, p.y) * invR : vec2(0.0, 1.0);

    float u    = r_xy - ra;
    float l2   = u * u + p.z * p.z;
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4( l2 * invL - rb, u * er.x * invL, u * er.y * invL, p.z * invL );
    }
    return vec4( -rb, er.x, er.y, 0.0 );
}

#endif // SDF_GEMS_SDG_CAPPED_TORUS_GLSL
