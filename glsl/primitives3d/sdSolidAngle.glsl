// ============================================================================
// sdSolidAngle — Exact 3D Solid Angle (Spherical Sector)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdSolidAngle( vec3 p, vec2 c, float ra )
{
    vec2 q = vec2( length(p.xz), p.y );
    float l = length(q) - ra;
    float m = length(q - c*clamp(dot(q,c), 0.0, ra) );
    return max(l, m*sign(c.y*q.x - c.x*q.y));
}
