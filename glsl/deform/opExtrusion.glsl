// ============================================================================
// opExtrusion — Exact 2D-to-3D Linear Extrusion & Universal 1-Sqrt Squared Extrusion
// @author Inigo Quilez — upgraded with 1-sqrt opExtrusionSq using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` total in 3D via opExtrusionSq)
// ============================================================================

#ifndef SDF_GEMS_OP_EXTRUSION_GLSL
#define SDF_GEMS_OP_EXTRUSION_GLSL

float opExtrusion( in vec3 p, in float sdf2d, in float h )
{
    vec2 w = vec2( sdf2d, abs(p.z) - h );
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

// Universal 1-Sqrt 3D Extrusion: given any 2D primitive's unsigned squared distance
// d2 >= 0 and outside sign s (s > 0 outside, s <= 0 inside), computes the exact 3D
// extruded prism SDF using only 1 sqrt total (fusing the 2D and 3D sqrts!).
float opExtrusionSq( in float d2, in float s, in float z, in float h )
{
    float qz = abs(z) - h, m = max(qz, 0.0);
    float d  = sqrt((s > 0.0) ? d2 + m * m : d2);
    return (s > 0.0) ? d : max(-d, qz);
}

#endif // SDF_GEMS_OP_EXTRUSION_GLSL
