// ============================================================================
// sdVesicaSegment — Exact 3D Revolved Vesica / Spindle
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3 c = (a + b)*0.5;
    float l = length(b - a);
    vec3 v = (b - a)/l;
    float y = dot(p - c, v);
    vec2 q = vec2(length(p - c - y*v), abs(y));
    float r = 0.5*l;
    float d = 0.5*(r*r - w*w)/w;
    vec3 h = (r*q.x < d*(q.y - r)) ? vec3(0.0, r, 0.0) : vec3(-d, 0.0, d + w);
    return length(q - h.xy) - h.z;
}
