// ============================================================================
// opRevolution — Exact 2D-to-3D Radial Revolution
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

vec2 opRevolution( in vec3 p, float o )
{
    return vec2( length(p.xz) - o, p.y );
}
