// ============================================================================
// sdgDeathStar — Exact 1/2-SFU 3D Gradients for sdDeathStar, sdLens & sdVesicaSegment
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 1 inversesqrt in spherical/cap regions, 2 in rim/body, 0 NaNs)
// ============================================================================

#ifndef SDF_GEMS_SDG_DEATH_STAR_GLSL
#define SDF_GEMS_SDG_DEATH_STAR_GLSL

// ----------------------------------------------------------------------------
// 1. sdgDeathStar(p, ra, rb, d) — 1/2-SFU Sphere Subtracted by Sphere + Gradient
// ----------------------------------------------------------------------------
vec4 sdgDeathStar( vec3 p, float ra, float rb, float d )
{
    float yz2 = dot(p.yz, p.yz), a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2  = max(ra * ra - a * a, 0.0), da = d - a, dx = d - p.x;
    if( b2 * p.x * abs(p.x) > a * abs(a) * yz2 && da * abs(da) * yz2 > b2 * dx * abs(dx) )
    {
        float invY = (yz2 > 0.0) ? inversesqrt(yz2) : 0.0;
        float py   = yz2 * invY;
        vec2  er   = (yz2 > 0.0) ? (p.yz * invY) : vec2(1.0, 0.0);
        float b    = sqrt(b2);
        float wy   = (py + b > 0.0) ? ((yz2 - b2) / (py + b)) : 0.0;
        vec2  w    = vec2( p.x - a, wy );
        float l2   = dot(w, w);
        if( l2 > 0.0 )
        {
            float invL = inversesqrt(l2);
            return vec4( l2 * invL, vec3(w.x, w.y * er.x, w.y * er.y) * invL );
        }
        return vec4( 0.0, 1.0, 0.0, 0.0 );
    }
    float r   = ra + rb, m = clamp(ra * r + d * (p.x - a), 0.0, r * r);
    bool  ca  = (r * r * (p.x * p.x + yz2) >= m * m);
    float sgn = ca ? 1.0 : -1.0;
    vec3  w   = vec3( p.x - (ca ? 0.0 : d), p.yz );
    float r2  = dot(w, w);
    if( r2 > 0.0 )
    {
        float invR = inversesqrt(r2);
        return vec4( r2 * invR * sgn - (ca ? ra : -rb), (sgn * invR) * w );
    }
    return vec4( -(ca ? ra : -rb), sgn, 0.0, 0.0 );
}

// ----------------------------------------------------------------------------
// 2. sdgLens(p, ra, rb, d) — 1/2-SFU Asymmetric Two-Sphere Optical Lens + Gradient
// ----------------------------------------------------------------------------
vec4 sdgLens( vec3 p, float ra, float rb, float d )
{
    float yz2 = dot(p.yz, p.yz), a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2  = max(ra * ra - a * a, 0.0), da = d - a, dx = d - p.x;
    if( a * abs(a) * yz2 > b2 * p.x * abs(p.x) && da * abs(da) * yz2 > b2 * dx * abs(dx) )
    {
        float invY = (yz2 > 0.0) ? inversesqrt(yz2) : 0.0;
        float py   = yz2 * invY;
        vec2  er   = (yz2 > 0.0) ? (p.yz * invY) : vec2(1.0, 0.0);
        float b    = sqrt(b2);
        float wy   = (py + b > 0.0) ? ((yz2 - b2) / (py + b)) : 0.0;
        vec2  w    = vec2( p.x - a, wy );
        float l2   = dot(w, w);
        if( l2 > 0.0 )
        {
            float invL = inversesqrt(l2);
            return vec4( l2 * invL, vec3(w.x, w.y * er.x, w.y * er.y) * invL );
        }
        return vec4( 0.0, 0.0, er.x, er.y );
    }
    float r  = ra - rb, m = ra * r + d * (p.x - a);
    bool  ca = (r * abs(r) * (p.x * p.x + yz2) <= m * abs(m));
    vec3  w  = vec3( p.x - (ca ? 0.0 : d), p.yz );
    float r2 = dot(w, w);
    if( r2 > 0.0 )
    {
        float invR = inversesqrt(r2);
        return vec4( r2 * invR - (ca ? ra : rb), w * invR );
    }
    return vec4( -(ca ? ra : rb), ca ? 1.0 : -1.0, 0.0, 0.0 );
}

// ----------------------------------------------------------------------------
// 3. sdgVesicaSegment(p, a, b, w) — 1/2-SFU 3D Revolved Vesica / Spindle + Gradient
// ----------------------------------------------------------------------------
vec4 sdgVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3  ba = b - a, pc = p - 0.5 * (a + b);
    float l2 = dot(ba, ba), u = dot(pc, ba);
    float d  = (l2 - 4.0 * w * w) / (8.0 * w);
    vec3  up = pc - ba * (u / l2);
    float x2 = dot(up, up), z = d * (2.0 * abs(u) - l2);
    if( z > 0.0 && l2 * l2 * x2 < z * z )
    {
        vec3  vt = pc - ba * ((u >= 0.0) ? 0.5 : -0.5);
        float r2 = dot(vt, vt);
        if( r2 > 0.0 )
        {
            float invR = inversesqrt(r2);
            return vec4( r2 * invR, vt * invR );
        }
        return vec4( 0.0, normalize(ba) * ((u >= 0.0) ? 1.0 : -1.0) );
    }
    float invX = (x2 > 1e-24) ? inversesqrt(x2) : 0.0;
    float qx   = x2 * invX;
    vec3  er   = (x2 > 1e-24) ? (up * invX)
               : normalize((abs(ba.x) < abs(ba.z)) ? vec3(-ba.y, ba.x, 0.0) : vec3(0.0, -ba.z, ba.y));
    vec3  wv   = pc + d * er;
    float w2   = dot(pc, pc) + d * d + 2.0 * d * qx;
    float invW = inversesqrt(max(w2, 1e-30));
    return vec4( w2 * invW - (d + w), wv * invW );
}

#endif // SDF_GEMS_SDG_DEATH_STAR_GLSL
