// ============================================================================
// sminCircularGeom — Circular Geometrical Smooth Minimum
// @author Inigo Quilez
// @source https://iquilezles.org/articles/smin/
// @exact  Exact at 90 deg (non-conservative for acute angles)
// ============================================================================

float sminCircularGeom( float a, float b, float k )
{
    k *= 1.0/(1.0 - sqrt(0.5));
    return max(k, min(a, b)) - length(max(k - vec2(a, b), 0.0));
}
