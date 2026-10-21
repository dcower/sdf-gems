// ============================================================================
// sminCubic — Cubic Polynomial Smooth Minimum
// @author Inigo Quilez
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`)
// ============================================================================

float sminCubic( float a, float b, float k )
{
    k *= 6.0;
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - h*h*h*k*(1.0/6.0);
}
