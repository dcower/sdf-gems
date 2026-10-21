// ============================================================================
// sminCircularApprox — Polynomial Approximation of Circular CD Smooth Minimum
// @author Inigo Quilez
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative
// ============================================================================

float sminCircularApprox( float a, float b, float k )
{
    k *= 1.0/(1.0 - sqrt(0.5));
    float x = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - k*0.5*x*x*(x*(x*(1.0 - sqrt(0.5)) - (2.0 - sqrt(2.0))) + 1.0);
}
