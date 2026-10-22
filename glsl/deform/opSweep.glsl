// ============================================================================
// opSweep & sdgSweep — Exact Eikonal Planar-Spine, Finite-Segment & Circular-Arc Profile Sweeps + 3D Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere outside spine medial axis; 0 sqrts & 2 muls in sdgSweepPlanar)
// ============================================================================

#ifndef SDF_GEMS_OP_SWEEP_GLSL
#define SDF_GEMS_OP_SWEEP_GLSL

// ----------------------------------------------------------------------------
// 1. Planar-Spine Orthogonal Sweep (opSweepPlanar & sdgSweepPlanar)
//
// By the Planar-Spine Orthogonal Sweep Theorem, for ANY exact 2D spine SDF
// d_spine(p.xy) with ||grad_xy d_spine|| = 1, the mapping
//   p in R^3  |-->  q = vec2(d_spine(p.xy) - o, p.z) in R^2
// has orthonormal Jacobian rows grad_p(u) = (gx_s, gy_s, 0) and grad_p(v) = (0, 0, 1).
// Evaluating ANY exact 2D cross-section profile d_profile(q) produces an exact
// 3D swept solid with unit 3D gradient in 0 sqrts, 0 divisions, and 2 multiplies!
// ----------------------------------------------------------------------------
vec2 opSweepPlanar( float d_spine2d, float z, float o )
{
    return vec2( d_spine2d - o, z );
}

vec2 opSweepPlanar( float d_spine2d, float z )
{
    return vec2( d_spine2d, z );
}

vec4 sdgSweepPlanar( vec3 dg_spine2d, vec3 dg_profile2d )
{
    return vec4( dg_profile2d.x, dg_profile2d.y * dg_spine2d.yz, dg_profile2d.z );
}

// ----------------------------------------------------------------------------
// 2. Revolved-Spine Circular-Arc Sweep (opSweepArcSpine & sdgSweepArcSpine)
//
// Generalizes sdCappedTorus to sweep ANY u-symmetric 2D profile along a circular
// arc of major radius ra and half-angle sc = vec2(sin(theta), cos(theta)) in the
// xy-plane, revolving the 2D profile around the two arc endpoints in the cap zone.
// ----------------------------------------------------------------------------
vec2 opSweepArcSpine( vec3 p, vec2 sc, float ra )
{
    p.x = abs(p.x);
    float u = (sc.y * p.x > sc.x * p.y) ? length(p.xy - ra * sc) : (length(p.xy) - ra);
    return vec2( u, p.z );
}

vec3 sdgArcSpine2D( vec2 p, vec2 sc, float ra )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    if( sc.y * p.x > sc.x * p.y )
    {
        vec2  w  = p - ra * sc;
        float l2 = dot(w, w);
        if( l2 > 0.0 )
        {
            float invL = inversesqrt(l2);
            return vec3( l2 * invL, sx * w.x * invL, w.y * invL );
        }
        return vec3( 0.0, sx * sc.y, -sc.x );
    }
    float r2 = dot(p, p);
    if( r2 > 0.0 )
    {
        float invR = inversesqrt(r2);
        return vec3( r2 * invR - ra, sx * p.x * invR, p.y * invR );
    }
    return vec3( -ra, 0.0, 1.0 );
}

vec4 sdgSweepArcSpine( vec3 p, vec2 sc, float ra, vec3 dg_profile2d )
{
    vec3 dg_s = sdgArcSpine2D( p.xy, sc, ra );
    return sdgSweepPlanar( dg_s, dg_profile2d );
}

// ----------------------------------------------------------------------------
// 3. Finite-Segment Profile Sweep (opSweepSegment & sdgSweepSegment)
//
// Sweeps any 2D profile d2g = sdgProfile2D(p.yz) along a finite x-axis segment
// [-L, L] with exact Flat, Round-Beveled (r), or 45-Degree Chamfered (c) end-caps.
// Uses the 0-SFU median/min check so interior + flat faces use 0 sqrts!
// ----------------------------------------------------------------------------
float opSweepSegment( vec3 p, float L, float d_profile2d )
{
    vec2 w = vec2( d_profile2d, abs(p.x) - L );
    return (min(w.x, w.y) <= 0.0) ? max(w.x, w.y) : length(w);
}

float opSweepSegmentRound( vec3 p, float L, float r, float d_profile2d )
{
    vec2 w = vec2( d_profile2d + r, abs(p.x) - L + r );
    return ((min(w.x, w.y) <= 0.0) ? max(w.x, w.y) : length(w)) - r;
}

float opSweepSegmentChamfer( vec3 p, float L, float c, float d_profile2d )
{
    float u = d_profile2d, v = abs(p.x) - L;
    float a = min(u, v), b = max(u, v);
    const float s2 = 0.70710678;
    return (b > 0.0 && b - a > c) ? length(vec2(max(a + c, 0.0), b))
                                  : max(b, (a + b + c) * s2);
}

vec4 sdgSweepSegment( vec3 p, float L, vec3 dg_profile2d )
{
    float u  = dg_profile2d.x, v = abs(p.x) - L;
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    if( min(u, v) <= 0.0 )
    {
        return (u > v) ? vec4( u, 0.0, dg_profile2d.yz )
                       : vec4( v, sx, 0.0, 0.0 );
    }
    vec2  q    = vec2( u, v );
    float l2   = dot(q, q), invL = inversesqrt(l2);
    vec2  duv  = q * invL;
    return vec4( l2 * invL, duv.y * sx, duv.x * dg_profile2d.yz );
}

vec4 sdgSweepSegmentRound( vec3 p, float L, float r, vec3 dg_profile2d )
{
    vec4 res = sdgSweepSegment( p, L - r, vec3(dg_profile2d.x + r, dg_profile2d.yz) );
    return vec4( res.x - r, res.yzw );
}

vec4 sdgSweepSegmentChamfer( vec3 p, float L, float c, vec3 dg_profile2d )
{
    float u   = dg_profile2d.x, v = abs(p.x) - L;
    float sx  = (p.x < 0.0) ? -1.0 : 1.0;
    float a   = min(u, v), b = max(u, v);
    const float s2 = 0.70710678;
    bool  w   = (b > 0.0 && b - a > c);
    float qa  = max(a + c, 0.0), l = length(vec2(qa, b));
    bool  ch  = (a + b + c) * s2 > b;
    float d   = w ? l : (ch ? (a + b + c) * s2 : b);
    vec2  ab  = w ? vec2(qa, b) / l : (ch ? vec2(s2) : vec2(0.0, 1.0));
    vec2  duv = (u < v) ? ab : ab.yx;
    return vec4( d, duv.y * sx, duv.x * dg_profile2d.yz );
}

// ----------------------------------------------------------------------------
// 4. Orthonormal Cap-Plane Circular-Arc Sweep (opSweepArc & sdgSweepArc)
//
// Sweeps ANY 2D profile (including asymmetric profiles!) along a circular arc
// of radius ra and half-angle sc = vec2(sin(theta), cos(theta)) in the xz-plane:
//   Step 1: vec3 q = opSweepArc(p, ra, sc);
//   Step 2: vec3 d2g = sdgProfile2D(q.xy);
//   Step 3: return sdgSweepArcFlat(p, sc, d2g) / sdgSweepArcRound(p, sc, r, d2g);
// ----------------------------------------------------------------------------
vec3 opSweepArc( vec3 p, float ra, vec2 sc )
{
    vec2  q = vec2(abs(p.x), p.z);
    float u = q.x * sc.x + q.y * sc.y;
    float w = q.x * sc.y - q.y * sc.x;
    float l = length(q);
    return vec3( ((w > 0.0) ? u : l) - ra, p.y, (w <= 0.0 && u < 0.0) ? -l : w );
}

float opSweepArcFlat( float d2d, float w )
{
    return (min(d2d, w) <= 0.0) ? max(d2d, w) : length(vec2(d2d, w));
}

float opSweepArcRound( float d2d, float w, float r )
{
    float u = d2d + r;
    if( w > 0.0 )
        return ((u <= 0.0) ? w : length(vec2(u, w))) - r;
    return max(d2d, w - r);
}

vec4 sdgSweepArcFlat( vec3 p, vec2 sc, vec3 d2g )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  q  = vec2(abs(p.x), p.z);
    float u  = q.x * sc.x + q.y * sc.y;
    float w  = q.x * sc.y - q.y * sc.x;
    vec3  gw = vec3(sx * sc.y, 0.0, -sc.x);
    if( w > 0.0 )
    {
        if( d2g.x <= 0.0 ) return vec4( w, gw );
        float l2   = d2g.x * d2g.x + w * w, invL = inversesqrt(l2);
        vec3  gu   = vec3(sx * d2g.y * sc.x, d2g.z, d2g.y * sc.y);
        return vec4( l2 * invL, (d2g.x * gu + w * gw) * invL );
    }
    if( d2g.x < 0.0 && u > 0.0 && w > d2g.x ) return vec4( w, gw );
    float r2   = dot(q, q);
    float invR = (r2 > 0.0) ? inversesqrt(r2) : 0.0;
    vec2  er   = (r2 > 0.0) ? (q * invR) : vec2(0.0, 1.0);
    return vec4( d2g.x, sx * d2g.y * er.x, d2g.z, d2g.y * er.y );
}

vec4 sdgSweepArcRound( vec3 p, vec2 sc, float r, vec3 d2g )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  q  = vec2(abs(p.x), p.z);
    float u  = q.x * sc.x + q.y * sc.y;
    float w  = q.x * sc.y - q.y * sc.x;
    float di = d2g.x + r;
    vec3  gw = vec3(sx * sc.y, 0.0, -sc.x);
    if( w > 0.0 )
    {
        if( di <= 0.0 ) return vec4( w - r, gw );
        float l2   = di * di + w * w, invL = inversesqrt(l2);
        vec3  gu   = vec3(sx * d2g.y * sc.x, d2g.z, d2g.y * sc.y);
        return vec4( l2 * invL - r, (di * gu + w * gw) * invL );
    }
    if( di < 0.0 && u > 0.0 && w > di ) return vec4( w - r, gw );
    float r2   = dot(q, q);
    float invR = (r2 > 0.0) ? inversesqrt(r2) : 0.0;
    vec2  er   = (r2 > 0.0) ? (q * invR) : vec2(0.0, 1.0);
    return vec4( d2g.x, sx * d2g.y * er.x, d2g.z, d2g.y * er.y );
}

#endif // SDF_GEMS_OP_SWEEP_GLSL
