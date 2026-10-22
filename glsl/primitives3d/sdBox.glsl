// ============================================================================
// sdBox — Exact 0/1-Sqrt 3D Axis-Aligned Box (0-Sqrt Interior & 6 Face Prisms)
// @author Inigo Quilez — upgraded to 0-sqrt median face-prism check using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0 `sqrt`s in interior + all 6 exterior face prisms, 1 in edges/corners)
// ============================================================================

#ifndef SDF_GEMS_SD_BOX_GLSL
#define SDF_GEMS_SD_BOX_GLSL

float sdBox( vec3 p, vec3 b )
{
    vec3 q = abs(p) - b;
    if( max(min(q.x, q.y), min(max(q.x, q.y), q.z)) <= 0.0 )
        return max(q.x, max(q.y, q.z));
    return length(max(q, 0.0));
}

#endif // SDF_GEMS_SD_BOX_GLSL
