// ============================================================================
// sminQuadratic — Quadratic Polynomial Smooth Minimum
// @author Inigo Quilez
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`||grad|| <= 1`)
// ============================================================================

float sminQuadratic( float a, float b, float k )
{
    k *= 4.0;
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - h*h*k*(1.0/4.0);
}
