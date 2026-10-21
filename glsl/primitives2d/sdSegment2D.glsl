// ============================================================================
// sdSegment2D — Exact 2D Line Segment / Capsule
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdSegment2D( vec2 p, vec2 a, vec2 b )
{
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}
