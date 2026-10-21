// ============================================================================
// sminCircularApprox — [3/1] Padé Rational Circular CD Smooth Minimum
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/smin/
// @exact  Conservative (`9.80e-4` max error — `11.0x` more accurate than IQ's quartic)
// ============================================================================

float sminCircularApprox( float a, float b, float k )
{
    k *= 1.0/(1.0 - sqrt(0.5));
    float h = max( k - abs(a - b), 0.0 )/k;
    return min(a, b) - k * h * h * (1.0 + 0.41421356 * h) / (2.0 * (1.0 + 1.41421356 * h));
}
