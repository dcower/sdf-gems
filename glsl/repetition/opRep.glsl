// ============================================================================
// opRep — 1-Tap Unlimited, Limited, Mirrored, Clamped-Floor Base, and 2-Tap Limited Repetition
// @author Inigo Quilez — upgraded with opRepMirror, opRepLimBase & opRepLim2Tap using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (symmetric SDFs for 1-tap; arbitrary SDFs with opRepLim2Tap or opRepLimBase)
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

// Degenerate-safe clamped-floor base index for 2^d-tap asymmetric limited repetition on [lima, limb].
// Clamping to max(lima, limb - 1.0) prevents GLSL undefined behavior when lima.k == limb.k.
vec3 opRepLimBase( in vec3 p, in float s, in vec3 lima, in vec3 limb )
{
    return clamp(floor(p / s), lima, max(lima, limb - 1.0));
}

// Exact 2-Tap 1D Limited Repetition on [lima, limb] (1st & 2nd closest valid lattice cells)
void opRepLim2Tap( in float p, in float s, in float lima, in float limb, out float id0, out float id1 )
{
    float u     = p / s;
    float b     = clamp(floor(u), lima, max(lima, limb - 1.0));
    float step0 = (u - b >= 0.5 && lima < limb) ? 1.0 : 0.0;
    id0 = b + step0;
    id1 = min(b + (1.0 - step0), limb);
}

// Exact 2-Tap 2D Limited Repetition on [lima, limb] (1st & 2nd closest valid lattice cells)
void opRepLim2Tap( in vec2 p, in float s, in vec2 lima, in vec2 limb, out vec2 id0, out vec2 id1 )
{
    vec2 u     = p / s;
    vec2 b     = clamp(floor(u), lima, max(lima, limb - 1.0));
    vec2 delta = u - (b + 0.5);
    vec2 step0 = vec2((delta.x >= 0.0 && lima.x < limb.x) ? 1.0 : 0.0,
                      (delta.y >= 0.0 && lima.y < limb.y) ? 1.0 : 0.0);
    id0 = b + step0;
    vec2 alt = min(b + (1.0 - step0), limb);
    vec2 pen = vec2((lima.x < limb.x) ? abs(delta.x) : 1e30,
                    (lima.y < limb.y) ? abs(delta.y) : 1e30);
    id1 = (pen.x <= pen.y) ? vec2(alt.x, id0.y) : vec2(id0.x, alt.y);
}

// Exact 2-Tap 3D Limited Repetition on [lima, limb] (1st & 2nd closest valid lattice cells)
void opRepLim2Tap( in vec3 p, in float s, in vec3 lima, in vec3 limb, out vec3 id0, out vec3 id1 )
{
    vec3 u     = p / s;
    vec3 b     = clamp(floor(u), lima, max(lima, limb - 1.0));
    vec3 delta = u - (b + 0.5);
    vec3 step0 = vec3((delta.x >= 0.0 && lima.x < limb.x) ? 1.0 : 0.0,
                      (delta.y >= 0.0 && lima.y < limb.y) ? 1.0 : 0.0,
                      (delta.z >= 0.0 && lima.z < limb.z) ? 1.0 : 0.0);
    id0 = b + step0;
    vec3 alt = min(b + (1.0 - step0), limb);
    vec3 pen = vec3((lima.x < limb.x) ? abs(delta.x) : 1e30,
                    (lima.y < limb.y) ? abs(delta.y) : 1e30,
                    (lima.z < limb.z) ? abs(delta.z) : 1e30);
    id1 = (pen.x <= pen.y && pen.x <= pen.z) ? vec3(alt.x, id0.y, id0.z) :
          (pen.y <= pen.z)                   ? vec3(id0.x, alt.y, id0.z) : vec3(id0.x, id0.y, alt.z);
}

#endif // SDF_GEMS_OP_REP_GLSL
