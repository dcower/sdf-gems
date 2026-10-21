// ============================================================================
// opRep3D — 4-Tap 3D BCC Delaunay Tetrahedron Domain Repetition (opRepBCC4)
// @author Inigo Quilez (8-tap Cartesian) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (`0.0` error vs 35-tap brute force)
// ============================================================================

void opRepBCC4( in vec3 p, in float s, out vec3 c0, out vec3 c1, out vec3 c2, out vec3 c3 )
{
    vec3 q  = p / s;
    vec3 b  = round(q);
    vec3 d  = q - b;
    vec3 sg = vec3(greaterThanEqual(d, vec3(0.0))) * 2.0 - 1.0;
    vec3 a  = abs(d);
    vec3 e1 = (a.x >= a.y && a.x >= a.z) ? vec3(sg.x, 0.0, 0.0) :
              (a.y >= a.z)               ? vec3(0.0, sg.y, 0.0) : vec3(0.0, 0.0, sg.z);
    vec3 e3 = (a.x <  a.y && a.x <  a.z) ? vec3(sg.x, 0.0, 0.0) :
              (a.y <  a.z)               ? vec3(0.0, sg.y, 0.0) : vec3(0.0, 0.0, sg.z);
    c0 = b;
    c1 = b + e1;
    c2 = b + 0.5 * sg;
    c3 = b + 0.5 * sg - e3;
}
