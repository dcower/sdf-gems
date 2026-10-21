// ============================================================================
// opRep — 1-Tap Unlimited & Limited Cartesian Repetition
// @author Inigo Quilez
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (symmetric SDFs)
// ============================================================================

vec3 opRep( in vec3 p, in vec3 s )
{
    return p - s*round(p/s);
}

vec3 opRepLim( in vec3 p, in float s, in vec3 l )
{
    return p - s*clamp(round(p/s), -l, l);
}
