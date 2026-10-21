// ============================================================================
// opExtrusion — Exact 2D-to-3D Linear Extrusion
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opExtrusion( in vec3 p, in float sdf2d, in float h )
{
    vec2 w = vec2( sdf2d, abs(p.z) - h );
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}
