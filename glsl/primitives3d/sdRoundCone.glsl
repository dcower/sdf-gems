// ============================================================================
// sdRoundCone & sdRoundConeSeg — Exact 1-Sqrt 3D Round Cone (Signed-Square & O(L^2) Normalized)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` down from 3, O(L^2) scale-safe vs O(L^8) overflow)
// ============================================================================

// 1. Exact 1-Sqrt Vertical Round Cone (down from 3 sqrts in baseline)
float sdRoundCone( vec3 p, float r1, float r2, float h )
{
    float b = (r1 - r2)/h, a2 = 1.0 - b*b;
    float x2 = dot(p.xz, p.xz), y = p.y, z = y - h, s = b*abs(b)*x2;
    return (a2*y*abs(y) < s) ? sqrt(x2 + y*y) - r1 :
           (a2*z*abs(z) > s) ? sqrt(x2 + z*z) - r2 :
                               sqrt(a2*x2) + b*y - r1;
}

// 2. Exact 1-Sqrt O(L^2) Scale-Safe Arbitrary-Segment Round Cone (eliminates O(L^8) FP16/FP32 overflow & 7 muls)
float sdRoundConeSeg( vec3 p, vec3 a, vec3 b, float r1, float r2 )
{
    vec3  ba = b - a, pa = p - a;
    float l2 = dot(ba, ba), il2 = 1.0 / l2;
    float rr = r1 - r2, b2 = rr*rr*il2, a2 = 1.0 - b2;
    float y  = dot(pa, ba)*il2, z = y - 1.0;
    float x2 = dot2(pa - ba*y), y2 = y*y*l2, z2 = z*z*l2;
    float s  = sign(rr)*b2*x2;
    return (sign(y)*a2*y2 < s) ? sqrt(x2 + y2) - r1 :
           (sign(z)*a2*z2 > s) ? sqrt(x2 + z2) - r2 :
                                 sqrt(a2*x2) + y*rr - r1;
}
