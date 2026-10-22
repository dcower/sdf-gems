// ============================================================================
// sdSegment2D, sdCapsule2D, sdgSegment2D & sdgCapsule2D — Exact 2D Line Segment / Capsule + Gradients
// @author Inigo Quilez — extended with analytical 1-SFU gradients using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 inversesqrt, 0 NaNs on segment spine)
// ============================================================================

#ifndef SDF_GEMS_SD_SEGMENT2D_GLSL
#define SDF_GEMS_SD_SEGMENT2D_GLSL

float sdSegment2D( vec2 p, vec2 a, vec2 b )
{
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa, ba) / dot(ba, ba), 0.0, 1.0 );
    return length( pa - ba * h );
}

float sdCapsule2D( vec2 p, vec2 a, vec2 b, float r )
{
    return sdSegment2D( p, a, b ) - r;
}

vec3 sdgSegment2D( vec2 p, vec2 a, vec2 b )
{
    vec2  pa = p - a, ba = b - a;
    float h  = clamp( dot(pa, ba) / dot(ba, ba), 0.0, 1.0 );
    vec2  v  = pa - ba * h;
    float l2 = dot(v, v);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec3( l2 * invL, v * invL );
    }
    vec2 n = normalize(vec2(-ba.y, ba.x));
    return vec3( 0.0, n );
}

vec3 sdgCapsule2D( vec2 p, vec2 a, vec2 b, float r )
{
    vec3 res = sdgSegment2D( p, a, b );
    return vec3( res.x - r, res.yz );
}

#endif // SDF_GEMS_SD_SEGMENT2D_GLSL
