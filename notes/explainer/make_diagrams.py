"""Draws the source-layout diagrams of the explainer as inline SVG, from the
coordinates the solvers actually use (layout_data.json, exported from
validation/pec_almond_skeleton.py and pec_cube_lightning3d.py), and builds
index.html from template.html.

    python3 make_diagrams.py
"""
import base64
import io
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "layout_data.json")))


def svg(w, h, body, label):
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">'
            f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0 0L10 5L0 10z" fill="currentColor"/></marker></defs>{body}</svg>')


def text(x, y, s, cls="lab", anchor="start", extra=""):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}" {extra}>{s}</text>'


def dots(pts, cls, r=2.4):
    return "".join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" class="{cls}"/>' for x, y in pts)


def dotpath(pts, cls, width=2.2):
    d = "".join(f"M{x:.1f} {y:.1f}h0" for x, y in pts)
    return f'<path d="{d}" class="{cls}" fill="none" stroke-width="{width}" stroke-linecap="round"/>'


def cross(x, y, cls="singm", s=5):
    return f'<path d="M{x - s:.1f} {y - s:.1f}L{x + s:.1f} {y + s:.1f}M{x - s:.1f} {y + s:.1f}L{x + s:.1f} {y - s:.1f}" class="{cls}" fill="none"/>'


def leader(x1, y1, x2, y2):
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="lead" marker-end="url(#ah)"/>'


# ------------------------------------------------------------------ D1: interior vs scattering

def d_problem():
    b = ""
    # left: thesis
    cx, cy, r = 180, 150, 46
    b += f'<circle cx="{cx}" cy="{cy}" r="{r}" class="field"/>'
    b += f'<circle cx="{cx}" cy="{cy}" r="{r}" class="bnd" fill="none"/>'
    ang = np.arange(16) * 2 * np.pi / 16
    b += dots([(cx + r * np.cos(t), cy + r * np.sin(t)) for t in ang], "col", 2.6)
    R = 2.5 * r
    b += f'<circle cx="{cx}" cy="{cy}" r="{R:.0f}" class="guide" fill="none"/>'
    b += dots([(cx + R * np.cos(t + 0.2), cy + R * np.sin(t + 0.2)) for t in ang], "src", 3)
    b += text(cx, cy + 4, "field wanted here", "labs", "middle")
    b += text(cx, 18, "Thesis: interior problem", "ttl", "middle")
    b += text(cx, 288, "source points outside, R = 2.5", "lab", "middle")
    # right: paper
    cx = 540
    b += f'<rect x="372" y="28" width="336" height="244" class="field" rx="4"/>'
    b += f'<circle cx="{cx}" cy="{cy}" r="{r}" class="body"/>'
    b += f'<circle cx="{cx}" cy="{cy}" r="{r}" class="bnd" fill="none"/>'
    b += dots([(cx + r * np.cos(t), cy + r * np.sin(t)) for t in ang], "col", 2.6)
    R = 0.5 * r
    b += f'<circle cx="{cx}" cy="{cy}" r="{R:.0f}" class="guide" fill="none"/>'
    b += dots([(cx + R * np.cos(t + 0.2), cy + R * np.sin(t + 0.2)) for t in ang[::2]], "src", 3)
    for yy in (70, 100, 130):
        b += f'<line x1="386" y1="{yy}" x2="440" y2="{yy}" class="wave" marker-end="url(#ah)"/>'
    b += text(388, 58, "incident plane wave", "labs")
    b += text(640, 216, "scattered field", "labs", "middle")
    b += text(640, 231, "wanted here", "labs", "middle")
    b += text(cx, 18, "Paper: scattering problem", "ttl", "middle")
    b += text(cx, 288, "source points inside, R = 0.5", "lab", "middle")
    return svg(720, 300, b, "Left: the thesis solved for the field inside the body with source points outside. "
                            "Right: the paper solves for the scattered field outside with source points inside.")


# ------------------------------------------------------------------ D0: the four singular sets

def d_singular_sets():
    b = ""
    w = 180
    names = ("Sphere", "Torus (cut through the axis)", "Almond (cross-section)", "Cube (cross-section)")
    notes = ("one point: any R works", "two rings: R is bounded", "6% from the rim: no copy fits", "reaches the edges: tails")
    for i in range(4):
        cx, cy = 90 + i * w, 96
        b += text(cx, 18, names[i], "ttl", "middle")
        b += text(cx, 186, notes[i], "labs", "middle")
        if i == 0:
            b += f'<circle cx="{cx}" cy="{cy}" r="56" class="body"/><circle cx="{cx}" cy="{cy}" r="56" class="bnd" fill="none"/>'
            b += cross(cx, cy)
        if i == 1:
            for sgn in (-1, 1):
                x0 = cx + sgn * 46
                b += f'<circle cx="{x0}" cy="{cy}" r="23" class="body"/><circle cx="{x0}" cy="{cy}" r="23" class="bnd" fill="none"/>'
                b += cross(x0, cy, s=3.5)
                b += cross(x0 - sgn * 0.134 * 46, cy, s=3.5)
            b += f'<line x1="{cx}" y1="34" x2="{cx}" y2="158" class="axis"/>'
        if i == 2:
            b += f'<ellipse cx="{cx}" cy="{cy}" rx="72" ry="24" class="body"/><ellipse cx="{cx}" cy="{cy}" rx="72" ry="24" class="bnd" fill="none"/>'
            for sgn in (-1, 1):
                b += cross(cx + sgn * 72 * D["eta_f"], cy, s=4)
        if i == 3:
            h = 52
            b += f'<rect x="{cx - h}" y="{cy - h}" width="{2 * h}" height="{2 * h}" class="body"/><rect x="{cx - h}" y="{cy - h}" width="{2 * h}" height="{2 * h}" class="bnd" fill="none"/>'
            b += f'<path d="M{cx - h} {cy - h}L{cx + h} {cy + h}M{cx - h} {cy + h}L{cx + h} {cy - h}" class="singl" fill="none"/>'
    return svg(720, 196, b, "Where the continued scattered field is singular: the centre of the sphere, two rings inside the torus tube, "
                            "the focal points near the rim of the almond, and the diagonal planes that run into the edges of the cube.")


# ------------------------------------------------------------------ D2: torus

def d_torus():
    s, cx, cy = 150.0, 360.0, 150.0
    R, a = 1.0, 0.5
    c = np.sqrt(R * R - a * a)
    b = f'<line x1="{cx}" y1="18" x2="{cx}" y2="282" class="axis"/>'
    b += text(cx + 6, 28, "axis of the torus", "labs")
    for sgn in (-1, 1):
        x0 = cx + sgn * R * s
        b += f'<circle cx="{x0}" cy="{cy}" r="{a * s}" class="body"/><circle cx="{x0}" cy="{cy}" r="{a * s}" class="bnd" fill="none"/>'
        b += f'<circle cx="{x0}" cy="{cy}" r="{0.25 * s}" class="guide" fill="none"/>'
        t = (np.arange(12) + 0.5) * 2 * np.pi / 12
        b += dots([(x0 + 0.25 * s * np.cos(u), cy + 0.25 * s * np.sin(u)) for u in t], "src", 3)
        t = (np.arange(24) + 0.5) * 2 * np.pi / 24
        b += dots([(x0 + a * s * np.cos(u), cy + a * s * np.sin(u)) for u in t], "col", 2.4)
        b += cross(x0, cy, s=4.5)
        b += cross(x0 - sgn * (R - c) * s, cy, s=4.5)
    x0 = cx + R * s
    b += f'<circle cx="{x0}" cy="{cy}" r="{(R - c) * s:.1f}" class="forbid"/>'
    # labels (right-hand disc)
    b += leader(642, 44, x0 + 0.25 * s * np.cos(-1.0) + 3, cy + 0.25 * s * np.sin(-1.0) - 4)
    b += text(646, 40, "source points", "lab", "end")
    b += text(646, 26, "ring of dipoles each", "labs", "end")
    b += leader(668, 262, x0 + a * s * np.cos(0.9) + 3, cy + a * s * np.sin(0.9) + 3)
    b += text(704, 278, "boundary points", "lab", "end")
    b += leader(448, 252, x0 - (R - c) * s - 3, cy + 8)
    b += text(372, 268, "focal ring, \u03c1 = \u221a(R\u00b2\u2212a\u00b2)", "lab")
    b += leader(566, 252, x0 + 2, cy + 8)
    b += text(540, 268, "centre circle", "lab")
    b += text(cx - R * s, 288, "the same rings, cut on the other side", "labs", "middle")
    b += f'<line x1="{cx}" y1="58" x2="{cx + R * s}" y2="58" class="dim"/>'
    b += f'<line x1="{cx + R * s}" y1="53" x2="{cx + R * s}" y2="63" class="dim"/>'
    b += text(cx + 0.5 * R * s, 52, "R", "labs", "middle")
    return svg(720, 300, b, "Cut through a torus along its axis. Source points sit on a smaller circle inside each cut of the tube; "
                            "that circle must contain the centre of the tube and the focal ring, which lies between the centre and the hole.")


# ------------------------------------------------------------------ D3: almond cross-section

def d_almond_section():
    cx, cy, A = 360.0, 122.0, 270.0
    B = A / 3
    b = f'<ellipse cx="{cx}" cy="{cy}" rx="{A}" ry="{B}" class="body"/><ellipse cx="{cx}" cy="{cy}" rx="{A}" ry="{B}" class="bnd" fill="none"/>'
    b += f'<ellipse cx="{cx}" cy="{cy}" rx="{0.65 * A}" ry="{0.65 * B}" class="bad" fill="none"/>'
    b += f'<line x1="{cx - A}" y1="{cy}" x2="{cx + A}" y2="{cy}" class="axis"/>'
    for sgn in (-1, 1):
        b += cross(cx + sgn * D["eta_f"] * A, cy, s=5)
        b += dots([(cx + sgn * e * A, cy) for e in D["etas"]], "src", 3.2)
    b += leader(236, 20, cx - 0.65 * A * np.cos(1.1) - 2, cy - 0.65 * B * np.sin(1.1) - 3)
    b += text(242, 18, "scaled copy of the surface: misses the focal points", "lab")
    b += leader(654, 222, cx + D["eta_f"] * A + 2, cy + 8)
    b += text(700, 238, "focal point, 94% of the half-width", "lab", "end")
    b += leader(140, 222, cx - D["etas"][7] * A, cy + 6)
    b += text(20, 238, "source points on the mid-plane, denser towards the focal point", "lab")
    return svg(720, 250, b, "Cross-section of the almond: an ellipse three times as wide as it is high. Its focal points lie at 94 percent "
                            "of the half-width, outside any scaled-down copy of the surface. The source points sit on the mid-plane between them.")


# ------------------------------------------------------------------ D4: almond mid-plane, real points

def d_almond_plan():
    s = 2650.0
    x0 = D["outline_x"][0]
    X = lambda x: 26 + (np.asarray(x) - x0) * s
    Y = lambda y: 172 - np.asarray(y) * s
    ox, oa = np.array(D["outline_x"]), np.array(D["outline_a"])
    up = " ".join(f"{X(x):.1f},{Y(a):.1f}" for x, a in zip(ox, oa))
    lo = " ".join(f"{X(x):.1f},{Y(-a):.1f}" for x, a in zip(ox[::-1], oa[::-1]))
    b = f'<polygon points="{up} {lo}" class="body"/><polygon points="{up} {lo}" class="bnd" fill="none"/>'
    # focal curve: front part eta_f * a(x); rear part the focal ellipse
    L = D["L"]
    A_, B_, C_ = 0.416667 * L, 0.193333 * L, 0.064444 * L
    psi = np.linspace(0, np.pi / 2, 60)
    rear = [(-np.sqrt(A_ ** 2 - C_ ** 2) * np.sin(p), np.sqrt(B_ ** 2 - C_ ** 2) * np.cos(p)) for p in psi][::-1]
    front = [(x, D["eta_f"] * a) for x, a in zip(ox, oa) if x >= 0]
    fc = rear + front
    for sgn in (1, -1):
        b += '<polyline points="' + " ".join(f"{X(x):.1f},{Y(sgn * y):.1f}" for x, y in fc) + '" class="singl" fill="none"/>'
    src = np.array(D["almond_src"])
    pts = [(X(x), Y(y)) for x, y in src] + [(X(x), Y(-y)) for x, y in src]
    b += dotpath(pts, "srcs", 2.0)
    xt = D["T1"] * L
    b += f'<circle cx="{X(xt):.1f}" cy="172" r="7" class="tailring" fill="none"/>'
    b += leader(650, 62, X(xt) - 3, 163)
    b += text(704, 40, "tip: free source points stop 1 mm away,", "lab", "end")
    b += text(704, 54, "the tip tail continues to 0.0001 mm", "lab", "end")
    b += leader(170, 318, X(0.02), Y(-D["eta_f"] * np.interp(0.02, ox, oa)) + 4)
    b += text(26, 334, "focal curve (dashed): the rows of source points close in on it", "lab")
    b += text(26, 22, "the almond in its mid-plane, seen from above", "lab")
    b += f'<line x1="{X(0.09):.1f}" y1="318" x2="{X(0.14):.1f}" y2="318" class="dim"/>'
    b += text(X(0.115), 334, "5 cm", "labs", "middle")
    return svg(720, 344, b, "The almond seen from above with the 2762 source points used at 1.19 gigahertz. They lie in the mid-plane in ten rows on each side "
                            "that follow the rim and crowd towards the dashed focal curve; spacing along a row is half the local distance to the surface.")


# ------------------------------------------------------------------ D5: cube

def d_cube():
    s, cx, cy = 250.0, 172.0, 176.0
    h = 0.5 * s
    b = f'<rect x="{cx - h}" y="{cy - h}" width="{2 * h}" height="{2 * h}" class="body"/>'
    b += f'<path d="M{cx - h} {cy - h}L{cx + h} {cy + h}M{cx - h} {cy + h}L{cx + h} {cy - h}" class="singl" fill="none"/>'
    b += f'<rect x="{cx - h}" y="{cy - h}" width="{2 * h}" height="{2 * h}" class="bnd" fill="none"/>'
    g = np.array(D["cube_nodes"]) * s
    ticks = ""
    for u in g:
        ticks += f"M{cx + u:.1f} {cy - h - 4}v8M{cx + u:.1f} {cy + h - 4}v8M{cx - h - 4} {cy + u:.1f}h8M{cx + h - 4} {cy + u:.1f}h8"
    b += f'<path d="{ticks}" class="colt" fill="none"/>'
    base = []
    for u in (-0.225, -0.075, 0.075, 0.225):
        for v in (-0.3, 0.3):
            base += [(cx + u * s, cy + v * s), (cx + v * s, cy + u * s)]
    b += "".join(f'<rect x="{x - 2.6:.1f}" y="{y - 2.6:.1f}" width="5.2" height="5.2" class="srco"/>' for x, y in base)
    depths = 0.2 * 0.5 ** np.arange(6)
    for sx in (-1, 1):
        for sy in (-1, 1):
            b += dots([(cx + sx * (0.5 - d / np.sqrt(2)) * s, cy + sy * (0.5 - d / np.sqrt(2)) * s) for d in depths], "src", 2.8)
    b += text(cx, 22, "Cut through the middle of the cube", "ttl", "middle")
    b += text(cx, 330, "dots on the diagonals: edge lines seen end on,", "labs", "middle")
    b += text(cx, 344, "at depths 0.2, 0.1, 0.05, \u2026 from each edge", "labs", "middle")
    b += text(cx, 358, "squares: coarse grid of point dipoles \u00b7 ticks: boundary points", "labs", "middle")
    # oblique view
    ox, oy, t = 545.0, 178.0, 168.0
    P = lambda x, y, z: (ox + t * (x + 0.42 * y), oy - t * (z + 0.30 * y))
    V = [(x, y, z) for x in (-.5, .5) for y in (-.5, .5) for z in (-.5, .5)]
    E = [(i, j) for i in range(8) for j in range(i + 1, 8) if sum(abs(V[i][k] - V[j][k]) > 0.5 for k in range(3)) == 1]
    for i, j in E:
        hidden = any(v[0] < 0 and v[1] > 0 and v[2] < 0 for v in (V[i], V[j]))
        (x1, y1), (x2, y2) = P(*V[i]), P(*V[j])
        b += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="{"hid" if hidden else "bnd"}"/>'
    tri = [P(.5, -.5, -.5), P(.5, -.5, .5), P(0, 0, 0)]
    b += '<polygon points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in tri) + '" class="plane"/>'
    for d in (0.2, 0.1, 0.05, 0.025):
        u = d / np.sqrt(2)
        (x1, y1), (x2, y2) = P(.5 - u, -.5 + u, -(0.5 - d)), P(.5 - u, -.5 + u, 0.5 - d)
        b += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="srcl"/>'
    b += dots([P(*(np.array([1, -1, 1]) * (0.5 - d / np.sqrt(3)))) for d in depths[:4]], "src", 2.6)
    b += text(ox, 22, "One edge and its sources", "ttl", "middle")
    e1 = P(.5, -.5, -0.25)
    b += leader(690, 266, e1[0] + 4, e1[1] + 1)
    b += text(712, 282, "the edge", "lab", "end")
    l1 = P(.5 - 0.2 / np.sqrt(2), -.5 + 0.2 / np.sqrt(2), -0.2)
    b += leader(470, 322, l1[0] - 3, l1[1] + 3)
    b += text(392, 338, "lines of dipoles parallel to the edge,", "lab")
    b += text(392, 352, "in the plane that bisects it (shaded)", "lab")
    v1 = P(*(np.array([1, -1, 1]) * (0.5 - 0.2 / np.sqrt(3))))
    b += leader(470, 58, v1[0] - 4, v1[1] - 3)
    b += text(392, 46, "corner: point dipoles on the diagonal", "lab")
    return svg(720, 368, b, "Left: a cut through the middle of the cube. Source lines along the four vertical edges appear as dots on the diagonals, "
                            "at depths 0.2, 0.1, 0.05 and so on from each edge; a coarse grid of point dipoles sits on an inner cube; boundary points crowd "
                            "towards the edges. Right: one edge with its lines of dipoles lying in the plane that bisects the edge.")


# ------------------------------------------------------------------ D6: the tail

def d_tail():
    J = 5
    free = 0.2 * 0.5 ** np.arange(J + 1)
    tail = free[-1] * 0.5 ** np.arange(1, 10)
    lx = lambda d: 56 + (np.log10(0.28) - np.log10(d)) * 136.0
    yb = 150.0
    b = f'<line x1="40" y1="{yb}" x2="672" y2="{yb}" class="axisl" marker-end="url(#ah)"/>'
    b += text(676, yb + 4, "edge", "lab")
    for d, lab in ((0.1, "0.1"), (0.01, "0.01"), (0.001, "0.001"), (1e-4, "0.0001")):
        b += f'<line x1="{lx(d):.1f}" y1="{yb - 4}" x2="{lx(d):.1f}" y2="{yb + 4}" class="axisl"/>'
        b += text(lx(d), yb + 20, lab, "labs", "middle")
    b += text(40, yb + 38, "distance from the edge (side of the cube = 1), logarithmic scale", "labs")
    for d in free:
        b += f'<circle cx="{lx(d):.1f}" cy="{yb}" r="5.5" class="src"/>'
    e = 2.0 / 3.0
    for i, d in enumerate(tail):
        hgt = 70 * (d / free[-1]) ** e
        b += f'<rect x="{lx(d) - 5:.1f}" y="{yb - 9 - hgt:.1f}" width="10" height="{hgt:.1f}" class="tailbar"/>'
        b += f'<circle cx="{lx(d):.1f}" cy="{yb}" r="4" class="tail"/>'
    b += f'<path d="M{lx(free[0]) - 8:.1f} {yb - 20}v-8H{lx(free[-1]) + 8:.1f}v8" class="brk" fill="none"/>'
    b += text(0.5 * (lx(free[0]) + lx(free[-1])), yb - 36, "free lines: each has its own unknowns α", "lab", "middle")
    b += f'<path d="M{lx(tail[0]) - 8:.1f} {yb - 92}v-8H{lx(tail[-1]) + 8:.1f}v8" class="brk" fill="none"/>'
    b += text(0.5 * (lx(tail[0]) + lx(tail[-1])), yb - 108, "tail: nine more lines that share one set of unknowns", "lab", "middle")
    b += text(lx(tail[2]) + 14, yb - 50, "bar = fixed weight (d / d₅)^(2/3)", "labs")
    return svg(720, 200, b, "Depths of the edge lines on a logarithmic scale. Six free lines each carry their own unknowns. Below the deepest free line, "
                            "nine further lines are tied to one shared set of unknowns with fixed weights that follow a power law of the depth.")


def png_data(path, width=900):
    from PIL import Image
    im = Image.open(path).convert("RGB")
    im = im.resize((width, int(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.quantize(64).save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def main():
    val = os.path.join(HERE, "..", "validation")
    parts = {"D_PROBLEM": d_problem(), "D_SETS": d_singular_sets(), "D_TORUS": d_torus(), "D_ALMOND_SECTION": d_almond_section(),
             "D_ALMOND_PLAN": d_almond_plan(), "D_CUBE": d_cube(), "D_TAIL": d_tail(),
             "IMG_TORUS_COST": png_data(os.path.join(val, "fig_cost_torus.png")),
             "IMG_TORUS_SRC": png_data(os.path.join(val, "fig_torus_sources.png"), 1000),
             "IMG_CUBE_COST": png_data(os.path.join(val, "fig_cost_cube.png"))}
    html = open(os.path.join(HERE, "template.html")).read()
    for k, v in parts.items():
        assert "{{" + k + "}}" in html, k
        html = html.replace("{{" + k + "}}", v)
    assert "{{" not in html
    open(os.path.join(HERE, "index.html"), "w").write(html)
    print("wrote index.html,", len(html) // 1024, "kB")


if __name__ == "__main__":
    main()
