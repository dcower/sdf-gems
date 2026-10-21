// ============================================================================
// opXor — Exact Boolean Symmetric Difference (XOR)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/sdfxor/
// @exact  Yes (up to medial axis)
// ============================================================================

float opXor( float d1, float d2 ) { return max(min(d1, d2), -max(d1, d2)); }
