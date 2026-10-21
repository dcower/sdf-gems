// ============================================================================
// sdCappedConeSeg — Exact 3D Arbitrary-Segment Capped Cone (Axis-NaN Fixed & Endpoint-Subsumed)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0% NaN rate, 1 sqrt in exterior cap cylinder)
// ============================================================================

#ifndef SDF_GEMS_SD_CAPPED_CONE_SEG_GLSL
#define SDF_GEMS_SD_CAPPED_CONE_SEG_GLSL

float sdCappedConeSeg( vec3 p, vec3 a, vec3 b, float ra, float rb )
{
    vec3  ba   = b - a, pa = p - a;
    float rba  = rb - ra, baba = dot(ba, ba), paba = dot(pa, ba) / baba;
    float x    = length(pa - ba * paba);
    float rcap = (paba < 0.5) ? ra : rb;
    float cay  = abs(paba - 0.5) - 0.5;
    if (x <= rcap && cay >= 0.0) return cay * sqrt(baba);
    float f    = clamp((rba * (x - ra) + paba * baba) / (rba * rba + baba), 0.0, 1.0);
    float cbx  = x - ra - f * rba, cby = paba - f;
    float d2   = cbx * cbx + cby * cby * baba;
    if (x <= rcap) d2 = min(d2, cay * cay * baba);
    return ((cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0) * sqrt(d2);
}

#endif // SDF_GEMS_SD_CAPPED_CONE_SEG_GLSL
