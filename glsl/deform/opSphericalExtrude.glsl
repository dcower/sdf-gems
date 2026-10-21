// ============================================================================
// opSphericalExtrude — Polar-to-Cartesian Exact 3D Spherical Shell Extrusion
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opSphericalExtrude( in float r, in vec2 uv, in float R1, in float R2 )
{
    float R0     = 0.5 * (R1 + R2);
    float h      = 0.5 * (R2 - R1);
    float dShell = abs(r - R0) - h;
    float dWall  = length( vec2(uv.x - clamp(uv.x, R1, R2), uv.y) );
    return max( dShell, sign(uv.y) * dWall );
}
