// ============================================================================
// sdSolidAngle — Exact 2-Sqrt 3D Solid Angle (Orthonormal Sector-Dominance Reduction)
// @author Inigo Quilez — upgraded/optimized by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 2 sqrts down from 3, 0.0 error across R^3)
// ============================================================================

float sdSolidAngle( vec3 p, vec2 c, float ra )
{
    vec2  q = vec2( length(p.xz), p.y );
    float u = dot(q, c), v = c.y*q.x - c.x*q.y;
    float l = length(vec2(u - ((v > 0.0) ? clamp(u, 0.0, ra) : 0.0), v));
    return (v > 0.0) ? l : max(l - ra, (u < 0.0) ? -l : v);
}
