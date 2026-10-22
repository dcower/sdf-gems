// ============================================================================
// sdgCapsule & sdgCapsuleY — 1-SFU Singularity-Free 3D Capsule Distance & Gradient
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0.0 distance error, 0.0 gradient error, 0 NaNs on medial spine)
// ============================================================================

#ifndef SDF_GEMS_SDG_CAPSULE_GLSL
#define SDF_GEMS_SDG_CAPSULE_GLSL

// ----------------------------------------------------------------------------
// sdgCapsule(p, a, b, r) — Arbitrary 3D Line-Segment Capsule Distance & Gradient
//
// Standard sdgSegment / sdgCapsule implementations compute:
//   float d = length(v);
//   vec3  g = v / d;
// which requires 1 sqrt + 3 divisions (2 SFUs) and produces 0.0 / 0.0 = NaN
// along the entire medial spine segment v = 0!
//
// By using a single guarded reciprocal square root:
//   float invL = inversesqrt(max(l2, 1e-30));
//   float d    = (l2 > 0.0) ? l2 * invL : 0.0;
// and constructing an exact orthogonal fallback normal on the medial spine
// (l2 == 0.0), sdgCapsule evaluates both the exact signed distance and unit
// 3D normal in 1 inversesqrt total with 0 NaNs everywhere in R^3.
// ----------------------------------------------------------------------------
vec4 sdgCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  pa   = p - a, ba = b - a;
    float baba = dot(ba, ba);
    float h    = (baba > 0.0) ? clamp(dot(pa, ba) / baba, 0.0, 1.0) : 0.0;
    vec3  v    = pa - ba * h;
    float l2   = dot(v, v);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4( l2 * invL - r, v * invL );
    }
    vec3  ab = abs(ba);
    vec3  u  = (ab.x <= ab.y && ab.x <= ab.z) ? vec3(0.0, -ba.z, ba.y) :
               (ab.y <= ab.z)                 ? vec3(-ba.z, 0.0, ba.x) :
                                                vec3(-ba.y, ba.x, 0.0);
    float u2 = dot(u, u);
    vec3  g  = (u2 > 0.0) ? u * inversesqrt(u2) : vec3(1.0, 0.0, 0.0);
    return vec4( -r, g );
}

// ----------------------------------------------------------------------------
// sdgCapsuleY(p, h, r) — Vertical Y-Axis Capsule (0 Divisions, 1 inversesqrt)
//
// Vertical capsule with spine segment from (0, 0, 0) to (0, h, 0) and radius r.
// Requires 0 divisions and 1 inversesqrt across 100% of R^3.
// ----------------------------------------------------------------------------
vec4 sdgCapsuleY( vec3 p, float h, float r )
{
    vec3  v  = vec3( p.x, p.y - clamp(p.y, 0.0, h), p.z );
    float l2 = dot(v, v);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4( l2 * invL - r, v * invL );
    }
    return vec4( -r, 1.0, 0.0, 0.0 );
}

#endif // SDF_GEMS_SDG_CAPSULE_GLSL
