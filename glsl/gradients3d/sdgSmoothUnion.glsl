// ============================================================================
// sdgSmoothUnion — Angle-Exact Eikonal Circular Smooth Union, Intersection & Subtraction
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  true (||grad d|| = 1 identically, corner-local for c < 1, 0% overstep)
// ============================================================================

#ifndef SDF_GEMS_SDG_SMOOTH_UNION_GLSL
#define SDF_GEMS_SDG_SMOOTH_UNION_GLSL

vec4 sdgSmoothUnion( in vec4 a, in vec4 b, in float k )
{
    if( a.x > b.x ) { vec4 t = a; a = b; b = t; }
    float c  = max( dot( a.yzw, b.yzw ), -0.95 );
    float c2 = sqrt( 0.5 + 0.5 * c );
    float h  = max( 2.0 * k * ( 1.0 + c2 ) - ( b.x - c * a.x ), 0.0 );
    if( h <= 0.0 ) return a;
    float m  = 0.5 / ( ( 1.0 + c ) * ( 1.0 + c2 ) );
    float w  = m * h * h;
    float u  = k - ( 1.0 - c2 ) * a.x;
    float q  = sqrt( u * u + ( 1.0 - c2 ) * w );
    float b_ = m * h / q;
    return vec4( a.x - w / ( u + q ), ( u / q - c * b_ ) * a.yzw + b_ * b.yzw );
}

// Angle-Exact Eikonal Circular Smooth Intersection (De Morgan dual)
vec4 sdgSmoothIntersection( in vec4 a, in vec4 b, in float k )
{
    return -sdgSmoothUnion( -a, -b, k );
}

// Angle-Exact Eikonal Circular Smooth Subtraction (subtracts a from b)
vec4 sdgSmoothSubtraction( in vec4 a, in vec4 b, in float k )
{
    return -sdgSmoothUnion( a, -b, k );
}

#endif // SDF_GEMS_SDG_SMOOTH_UNION_GLSL
