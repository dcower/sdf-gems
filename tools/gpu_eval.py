"""
GPU harness for testing the *literal* sdf-gems `.glsl` files.

Compiles GLSL with the platform's real GLSL compiler (OpenGL 4.1 core via `moderngl`)
and evaluates arbitrary expressions on the GPU in FP32 using transform feedback.

Provides:
  - gpu_eval(paths, expr, pts)       : evaluate a vec4 GLSL expression of `vec3 p` at each row of pts (N,3)
  - gpu_fn(paths, expr_fmt)          : wrap gpu_eval into a numpy function usable by tools/render_sdf.py
  - glsl_float(x) / glsl_args(*xs)   : format Python/numpy numbers as GLSL float literals
                                        (never interpolate numpy scalars directly: `np.float64(...)` breaks GLSL)
  - compile_umbrella(new=..., call=...) : single-translation-unit compile of every glsl/**/*.glsl in
                                        forward, reverse and new-file-first order (optionally calling new functions)

Run via UV:
  uv run tools/gpu_eval.py                                   # umbrella compile check (3 orders)
  uv run tools/gpu_eval.py --new glsl/primitives3d/sdFoo.glsl --call "vec4(sdFoo(P, 1.0))"
"""

import argparse
import os
import re

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ctx = None
_progs = {}  # compiled transform-feedback programs, keyed by full source (complex multi-include shaders compile slowly)


def ctx():
    global _ctx
    if _ctx is None:
        import moderngl

        _ctx = moderngl.create_standalone_context(require=410)
    return _ctx


_park = None


def park_framebuffer():
    """Bind a small persistent framebuffer. Call after fragment renders and before releasing their FBOs: binding
    then releasing an FBO in the same context can make later transform-feedback evaluations silently return
    stale buffers."""
    global _park
    if _park is None:
        _park = ctx().simple_framebuffer((4, 4))
    _park.use()


def glsl_float(x):
    """Format a number as a GLSL float literal (works for Python and numpy scalars).
    `repr(float(x))` always yields a valid literal such as '12.0', '0.676', '1e-05'."""
    x = float(x)
    if not np.isfinite(x):
        raise ValueError(f"non-finite GLSL literal: {x}")
    return repr(x)


def glsl_args(*xs):
    """Comma-separated GLSL float literals, e.g. glsl_args(n, rb, ra) -> '12.0, 0.676, 0.84'."""
    return ", ".join(glsl_float(x) for x in xs)


def _strip_includes(src):
    return re.sub(r'^\s*#include\s+"[^"]*"\s*$', "", src, flags=re.M)


def _read(path):
    full = path if os.path.isabs(path) else os.path.join(REPO, path)
    return _strip_includes(open(full).read())


def load_sources(paths, with_common=True):
    srcs = [_read("glsl/common/common.glsl")] if with_common else []
    return "\n".join(srcs + [_read(p) for p in paths])


def gpu_eval(paths, expr, pts, extra="", with_common=True):
    """Evaluate `expr` (a vec4 GLSL expression of `vec3 p`) at each row of pts (N,3). Returns (N,4) float32."""
    import moderngl

    pts = np.ascontiguousarray(np.asarray(pts, dtype=np.float32).reshape(-1, 3))
    vs = (
        "#version 410\n" + load_sources(paths, with_common) + "\n" + extra
        + "\nin vec3 P;\nout vec4 R;\nvoid main(){ vec3 p = P; R = " + expr + "; }\n"
    )
    c = ctx()
    prog = _progs.get(vs)
    if prog is None:
        prog = _progs[vs] = c.program(vertex_shader=vs, varyings=["R"])
    vbo = c.buffer(pts.tobytes())
    out = c.buffer(reserve=max(pts.shape[0], 1) * 16)
    vao = c.vertex_array(prog, [(vbo, "3f", "P")])
    vao.transform(out, mode=moderngl.POINTS, vertices=pts.shape[0])
    res = np.frombuffer(out.read(), dtype=np.float32).reshape(-1, 4)[: pts.shape[0]].copy()
    for o in (vao, vbo, out):
        o.release()
    return res


def gpu_fn(paths, expr, column=0):
    """numpy-callable SDF backed by the literal GLSL: f(P) -> float64 array (for render_sdf / checks).
    `expr` is a vec4 expression of `p`; 2D inputs (N,2) are padded with z = 0."""

    def f(P):
        P = np.asarray(P, dtype=np.float64)
        shp = P.shape[:-1]
        P = P.reshape(-1, P.shape[-1])
        if P.shape[1] == 2:
            P = np.c_[P, np.zeros(len(P))]
        return gpu_eval(paths, expr, P)[:, column].astype(np.float64).reshape(shp)

    return f


def compile_umbrella(new=(), call=None):
    """Compile all headers listed in sdf_gems.glsl (+ `new` files) as ONE translation unit in
    forward, reverse and new-first order. `call` is an optional vec4 expression of `vec3 P`
    that calls the new functions (forces full type-checking of their call paths)."""
    umb = open(os.path.join(REPO, "sdf_gems.glsl")).read()
    incs = re.findall(r'#include\s+"([^"]+)"', umb)
    common, rest = incs[0], [f for f in incs[1:] if f not in new]
    new = list(new)
    body = call if call else "vec4(P, 1.0)"
    orders = {"forward": rest + new, "reverse": (rest + new)[::-1], "new-first": new + rest}
    for name, files in orders.items():
        if name == "new-first" and not new:
            continue
        src = "\n".join(_read(f) for f in [common] + files)
        vs = "#version 410\n" + src + "\nin vec3 P;\nout vec4 R;\nvoid main(){ R = " + body + "; }\n"
        ctx().program(vertex_shader=vs, varyings=["R"]).release()
        print(f"  {name:9s}: OK ({len(files) + 1} files)")
    return len(rest) + len(new) + 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", nargs="*", default=[], help="new .glsl files not yet in sdf_gems.glsl")
    ap.add_argument("--call", default=None, help="vec4 GLSL expression of `vec3 P` calling the new functions")
    a = ap.parse_args()
    n = compile_umbrella(a.new, a.call)
    print(f"Single-translation-unit compile passed ({n} files).")
