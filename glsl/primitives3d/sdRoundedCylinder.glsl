// ============================================================================
// sdRoundedCylinder — Exact 3D Rounded Cylinder
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
    vec2 d = vec2( length(p.xz) - 2.0*ra + rb, abs(p.y) - h );
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - rb;
}
