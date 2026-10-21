// ============================================================================
// opSubtraction — Hard Boolean Subtraction
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  No (bound near seam)
// ============================================================================

float opSubtraction( float d1, float d2 ) { return max(-d1, d2); }
