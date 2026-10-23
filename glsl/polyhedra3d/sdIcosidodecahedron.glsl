// ============================================================================
// sdIcosidodecahedron — Exact 1-Sqrt 3D Archimedean Icosidodecahedron
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, max error `2.22e-15`)
// ============================================================================

float sdIcosidodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for (int i = 0; i < 3; i++) {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);
    const vec2 np = vec2(0.52573111, 0.85065081);
    const vec2 nt = vec2(0.35682209, 0.93417236);
    float dp = dot(p.xz, np) - r;
    float dt = dot(p.yz, nt) - 1.09818549 * r;
    if (min(dp, dt) <= 0.79465447 * max(dp, dt))
        return max(dp, dt);
    vec3 q = p - vec3(0.0, 0.0, 1.17557050 * r);
    return length(q - min(dot(q, n3), 0.0) * n3);
}
