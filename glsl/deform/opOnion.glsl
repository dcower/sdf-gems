// ============================================================================
// opOnion — Exact Hollow Shell / Onion Operator
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opOnion( in float d, in float thickness ) { return abs(d) - thickness; }
