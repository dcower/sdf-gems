// ============================================================================
// opRep3D — 3D Neighbor-Corrected Domain Repetition (8-Tap Cartesian)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes
// ============================================================================

void opRep3D_8Tap( in vec3 p, in float s, out vec3 c[8] )
{
    vec3 q = p / s;
    vec3 id = round(q);
    vec3 o = sign(q - id);
    int idx = 0;
    for( int k = 0; k < 2; k++ )
    for( int j = 0; j < 2; j++ )
    for( int i = 0; i < 2; i++ )
        c[idx++] = id + vec3(i, j, k)*o;
}
