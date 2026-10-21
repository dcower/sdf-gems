// ============================================================================
// sdgWedgeIntersection — Exact Euclidean Gram-Schmidt Wedge Booleans
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  true (||grad d|| = 1 in exterior normal cone; fixes 10x error of max(a,b))
// ============================================================================

vec4 sdgWedgeIntersection( in vec4 a, in vec4 b )
{
    if( a.x < b.x ) { vec4 t = a; a = b; b = t; }
    float c = dot( a.yzw, b.yzw );
    float h = b.x - c * a.x;
    if( h <= 0.0 ) return a;
    float s2 = max( 1.0 - c * c, 1e-5 );
    float d  = sqrt( a.x * a.x + h * h / s2 );
    return vec4( d, ( a.x * a.yzw + ( h / s2 ) * ( b.yzw - c * a.yzw ) ) / d );
}

vec4 sdgWedgeSubtraction( in vec4 a, in vec4 b ) { return sdgWedgeIntersection( -a, b ); }
vec4 sdgWedgeUnion( in vec4 a, in vec4 b )       { return -sdgWedgeIntersection( -a, -b ); }
