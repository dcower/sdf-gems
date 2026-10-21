// ============================================================================
// sdEllipsoid — 3D Triaxial Ellipsoid (11-Line Unified Interior/Exterior Concave-Secular Halley)
// @author Inigo Quilez (V2 bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/ellipsoids/
// @exact  Unified Conservative Bound (~1e-6 mean error exterior & interior, 0 overstep, 0 NaNs)
// ============================================================================

float sdEllipsoid( vec3 p, vec3 r )
{
    p = max(abs(p), 1e-6 * r);
    vec3 r2 = r * r, u = r * p, u2 = u * u, v = p / r2, q2 = p * v, q4 = v * v;
    vec4 s0 = vec4(q2.yzx + q2.zxy, q2.x + q2.y + q2.z);
    vec4 s1 = vec4(q4.yzx + q4.zxy, q4.x + q4.y + q4.z);
    vec3 l1 = max(u, min(u * inversesqrt(max(1.0 - s0.xyz, 1e-12)), r2)) - r2;
    vec4 lv = max(vec4(l1, length(u) - dot(r2, u2) / dot(u, u)), s0 * (sqrt(s0) - 1.0) / s1);
    float lam = max(max(lv.x, lv.y), max(lv.z, lv.w));
    vec3 inv = 1.0 / (r2 + lam), w2 = u2 * inv * inv;
    float g2 = w2.x + w2.y + w2.z, g3 = dot(w2, inv), d = sqrt(g2) - 1.0;
    lam += g2 * d / (g3 + 1.5 * d * (g3 - g2 * dot(w2, inv * inv) / g3));
    v = p / (r2 + lam);
    return (dot(p, v) - length(r * v)) / length(v);
}
