// ============================================================================
// sminLocalOblique — Pure-Scalar Local Oblique Smooth Minimum
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Conservative for theta >= acos(c0) & Strictly Local
// ============================================================================

float sminLocalOblique( float a, float b, float k, float c0 )
{
    float a_s = min( a, b ), b_s = max( a, b );
    float c2  = sqrt( 0.5 + 0.5 * c0 );
    float h   = max( 2.0 * k * ( 1.0 + c2 ) - ( b_s - c0 * a_s ), 0.0 );
    if( h <= 0.0 ) return a_s;
    float w   = 0.5 * h * h / ( ( 1.0 + c0 ) * ( 1.0 + c2 ) );
    float u   = k - ( 1.0 - c2 ) * a_s;
    return a_s - w / ( u + sqrt( u * u + ( 1.0 - c2 ) * w ) );
}
