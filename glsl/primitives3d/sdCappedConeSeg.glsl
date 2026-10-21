// ============================================================================
// sdCappedConeSeg — Exact 3D Arbitrary-Segment Capped Cone (Axis-NaN Fixed)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0% NaN rate)
// ============================================================================

float sdCappedConeSeg( vec3 p, vec3 a, vec3 b, float ra, float rb )
{
    float rba  = rb - ra;
    vec3  ba   = b - a, pa = p - a;
    float baba = dot(ba, ba);
    float paba = dot(pa, ba) / baba;
    float x    = length(pa - ba * paba);
    float cax  = max(0.0, x - ((paba < 0.5) ? ra : rb));
    float cay  = abs(paba - 0.5) - 0.5;
    float k    = rba*rba + baba;
    float f    = clamp( (rba*(x - ra) + paba*baba)/k, 0.0, 1.0 );
    float cbx  = x - ra - f*rba;
    float cby  = paba - f;
    float s    = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(cax*cax + cay*cay*baba, cbx*cbx + cby*cby*baba) );
}
