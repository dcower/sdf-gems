// ============================================================================
// opRound — Exact Iso-Surface Rounding
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opRound( in float d, in float rad ) { return d - rad; }
