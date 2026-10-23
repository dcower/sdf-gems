// ============================================================================
// sdgSphere — 3D Sphere with Analytical Gradient
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

vec4 sdgSphere( in vec3 p, in float r )
{
    float l = length(p);
    return vec4( l - r, (l > 1e-20) ? (p / l) : vec3(1.0, 0.0, 0.0) );
}
