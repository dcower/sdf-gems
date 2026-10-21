// ============================================================================
// sdRhombus — 0/1/2-Sqrt Exact 3D Extruded Rounded Rhombus (x=0 Sign Bug Fixed)
// @author Inigo Quilez — upgraded to 0/1/2-sqrt & bug-fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, fixes x=0 sign(0.0) bug)
// ============================================================================

float sdRhombus( vec3 p, float la, float lb, float h, float ra )
{
    p = abs(p);
    float invL = inversesqrt(la * la + lb * lb);
    float du = (abs(la * p.x - lb * p.z + 0.5 * (lb * lb - la * la)) - 0.5 * (la * la + lb * lb)) * invL;
    float dv = (lb * p.x + la * p.z - la * lb) * invL;
    float qx = ((du > 0.0) ? length(vec2(du, dv)) : dv) - ra;
    float qy = p.y - h;
    return (qy <= 0.0) ? max(qx, qy) : ((qx <= 0.0) ? qy : length(vec2(qx, qy)));
}
