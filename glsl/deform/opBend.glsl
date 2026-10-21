// ============================================================================
// opBend — 3D Bend Deformation (opCheapBend)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  No (Shear distortion)
// ============================================================================

vec3 opCheapBend( in vec3 p, in float k )
{
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m*p.xy, p.z);
}
