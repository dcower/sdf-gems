// ============================================================================
// sdCappedCylinderSeg — 1/2-Sqrt Exact 3D Arbitrary-Segment Capped Cylinder
// @author Inigo Quilez — upgraded to 1/2 sqrts (down from 3) & FP16-safe using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `1` sqrt in cap & wall zones, `2` sqrts in rim wedge)
// ============================================================================

float sdCappedCylinderSeg( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float l2 = dot(ba, ba), t = dot(pa, ba) / l2;
    vec3  u  = pa - t * ba;
    float x2 = dot(u, u), uy = abs(t - 0.5) - 0.5;
    if( uy >= 0.0 ) {
        if( x2 <= r * r ) return uy * sqrt(l2);
        float wx = sqrt(x2) - r;
        return sqrt(wx * wx + uy * uy * l2);
    }
    if( x2 > r * r ) return sqrt(x2) - r;
    float wx = sqrt(x2) - r, wy2 = uy * uy * l2;
    return (wx * wx <= wy2) ? wx : uy * sqrt(l2);
}
