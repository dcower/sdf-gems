// ============================================================================
// opRep — 1-Tap Unlimited, Limited, Mirrored, and Clamped-Floor Base Repetition
// @author Inigo Quilez — upgraded with opRepMirror & opRepLimBase using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (symmetric SDFs for 1-tap; arbitrary SDFs with opRepLimBase 2^d loop)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_GLSL
#define SDF_GEMS_OP_REP_GLSL

vec3 opRep( in vec3 p, in vec3 s )
{
    return p - s*round(p/s);
}

vec3 opRepLim( in vec3 p, in float s, in vec3 l )
{
    return p - s*clamp(round(p/s), -l, l);
}

// Pure-floating-point alternating mirror repetition (WebGL1/FP-safe, works for negative cell IDs)
vec3 opRepMirror( in vec3 p, in vec3 s, out vec3 id )
{
    id = round(p / s);
    return (p - s * id) * (1.0 - 2.0 * mod(id, 2.0));
}

// Clamped-floor base index for 2^d-tap asymmetric limited repetition on [lima, limb].
// Unlike clamping inside the loop (which collapses both offsets to the boundary cell
// when p/s > limb or p/s < lima), clamping floor(p/s) to [lima, limb - 1] before the
// loop guarantees 2^d distinct adjacent cells everywhere in R^3!
vec3 opRepLimBase( in vec3 p, in float s, in vec3 lima, in vec3 limb )
{
    return clamp(floor(p / s), lima, limb - 1.0);
}

#endif // SDF_GEMS_OP_REP_GLSL
