// ============================================================================
// opTwist — 3D Axial Twist Deformation
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  No (Distorts metric)
// ============================================================================

vec3 opTwist( in vec3 p, in float k )
{
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m*p.xz, p.y);
}
