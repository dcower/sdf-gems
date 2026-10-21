// ============================================================================
// sminCircular — Circular Closest-Distance (CD) Smooth Minimum
// @author Inigo Quilez
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`)
// ============================================================================

float sminCircular( float a, float b, float k )
{
    k *= 1.0/(1.0 - sqrt(0.5));
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - k*0.5*(1.0 + h - sqrt(1.0 - h*(h - 2.0)));
}
