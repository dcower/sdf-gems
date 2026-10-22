// ============================================================================
// sdRhombus — 0/1/2-Sqrt Exact 3D Extruded Rounded Rhombus (Direct Vertex-Wedge Form)
// @author Inigo Quilez — upgraded to 0/1/2-sqrt & bug-fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, fixes x=0 sign(0.0) bug)
// ============================================================================

#ifndef SDF_GEMS_SD_RHOMBUS_GLSL
#define SDF_GEMS_SD_RHOMBUS_GLSL

float sdRhombus( vec3 p, float la, float lb, float h, float ra )
{
    p = abs(p);
    float qx = (la * (p.x - la) > lb * p.z) ? length(vec2(p.x - la, p.z)) - ra
             : (lb * (p.z - lb) > la * p.x) ? length(vec2(p.x, p.z - lb)) - ra
             : (lb * p.x + la * p.z - la * lb) * inversesqrt(la * la + lb * lb) - ra;
    float qy = p.y - h;
    return (qy <= 0.0) ? max(qx, qy) : ((qx <= 0.0) ? qy : length(vec2(qx, qy)));
}

#endif // SDF_GEMS_SD_RHOMBUS_GLSL
