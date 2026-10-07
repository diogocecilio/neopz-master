"""Drawing of the hexahedral finite element models (Hex20-Hex8 on a trilinear geometry) of the examples from the CSV
files of the geometric mesh written by mcc::WriteMeshCSV (<prefix>_nodes.csv, _elements.csv, _faces.csv,
_edges.csv).

Shared by the figure scripts of FLAC3DTriaxial (the single element of the triaxial tests), RS2Triaxial (the same
element in the finite element check) and TerzaghiConsolidation (the column).
The style is that of the model figures of the article (Fig. 8 of the Abaqus benchmark): orthographic projection
with back-face culling (the models are convex boxes), boundary faces coloured by their boundary condition, hidden
edges dashed, vertex nodes (displacement and pore pressure) as squares and mid-edge nodes (displacement only) as
dots.
"""
import os
import sys

import numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon as MPoly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcc_figstyle import C1, INK, INK2, MUTED, read_csv  # noqa: E402

# face colours of the model figures (as fig08_abaqus_model): loaded lateral faces (cell pressure), top faces with a
# prescribed displacement or load, faces with displacement restraints (symmetry planes, base)
FACE_LATERAL = '#c4e5d6'
FACE_TOP = '#b7d3f6'
FACE_RESTRAINED = '#e9e8e3'


class Mesh:
    """Geometric mesh read from the CSV files of mcc::WriteMeshCSV.

    Attributes: X (nodes x 3), faces (list of (4 vertex indices in the order of the element, set of boundary ids));
    coincident boundary elements (e.g. the displacement condition and the drained face) are merged into one face
    with all their ids), edges (n0, n1, nmid) of the volume elements, nelements, vertices (indices of the element
    vertices) and centroid of the body.
    """

    def __init__(self, run, prefix):
        N = read_csv(os.path.join(run, prefix + '_nodes.csv'))
        self.X = np.column_stack([N['x'], N['y'], N['z']])
        F = read_csv(os.path.join(run, prefix + '_faces.csv'))
        merged = {}
        for k in range(len(F['matid'])):
            nodes = tuple(int(F[f'n{j}'][k]) for j in range(4))
            key = tuple(sorted(nodes))
            merged.setdefault(key, [nodes, set()])[1].add(int(F['matid'][k]))
        self.faces = [(np.array(v[0]), v[1]) for v in merged.values()]
        E = read_csv(os.path.join(run, prefix + '_elements.csv'))
        self.nelements = len(E['element'])
        nn = int(E['nnodes'].max())
        self.vertices = sorted({int(E[f'n{j}'][k]) for k in range(self.nelements) for j in range(min(nn, 8))})
        self.centroid = self.X[self.vertices].mean(axis=0)
        G = read_csv(os.path.join(run, prefix + '_edges.csv'))
        self.edges = np.column_stack([G['n0'], G['n1'], G['nmid']]).astype(int)

    def normal(self, nodes):
        """Outward unit normal of a planar boundary face (oriented away from the centroid of the convex body)."""
        P = self.X[nodes]
        n = np.cross(P[2] - P[0], P[3] - P[1])
        n /= np.linalg.norm(n)
        return n if np.dot(n, P.mean(axis=0) - self.centroid) > 0 else -n

    def midpoint(self, edge):
        """Mid-edge node of an edge (the node of a quadratic geometry, or the midpoint of a straight edge)."""
        n0, n1, nm = edge
        return self.X[nm] if nm >= 0 else 0.5 * (self.X[n0] + self.X[n1])

    def faces_with(self, bcid):
        """Faces that carry the boundary id bcid."""
        return [f for f, ids in self.faces if bcid in ids]


class View:
    """Orthographic projection with the camera of matplotlib's view_init(elev, azim), scaled by s and shifted."""

    def __init__(self, elev, azim, scale=(1.0, 1.0, 1.0), shift=(0.0, 0.0)):
        e, a = np.radians(elev), np.radians(azim)
        self.d = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])   # towards the camera
        self.e1 = np.array([-np.sin(a), np.cos(a), 0.0])
        self.e2 = np.cross(self.d, self.e1)
        self.s = np.asarray(scale, dtype=float)
        self.shift = np.asarray(shift, dtype=float)

    def __call__(self, P):
        P = np.atleast_2d(np.asarray(P, dtype=float)) * self.s
        return np.column_stack([P @ self.e1, P @ self.e2]) + self.shift

    def point(self, *x):
        """Screen coordinates of one point."""
        return self(np.array(x, dtype=float))[0]

    def visible(self, n):
        """True for a face whose outward normal n points towards the camera (n is scaled with the coordinates)."""
        return float(np.dot(n / self.s, self.d)) > 1e-9


def draw_model(ax, mesh, view, facecolor, lw=0.6, nodes=True, hidden=True, ms=(2.6, 1.9), open_vertices=()):
    """Draws the visible boundary faces of the mesh coloured by facecolor(ids), the hidden edges dashed and, with
    nodes, the vertex nodes (squares, u and p_w) and the mid-edge nodes (dots, u only); the nodes of the hidden edges
    are drawn lighter. The vertices in open_vertices (e.g. those whose pore pressure is unknown) are drawn as open
    squares. Returns the screen bounding box (xmin, xmax, ymin, ymax)."""
    visible_faces = []
    box = [np.inf, -np.inf, np.inf, -np.inf]
    for nodes_f, ids in mesh.faces:
        if not view.visible(mesh.normal(nodes_f)):
            continue
        P = view(mesh.X[nodes_f])
        ax.add_patch(MPoly(P, closed=True, fc=facecolor(ids), ec=INK2, lw=lw, joinstyle='round', zorder=2))
        visible_faces.append(set(int(v) for v in nodes_f))
        box = [min(box[0], P[:, 0].min()), max(box[1], P[:, 0].max()), min(box[2], P[:, 1].min()),
               max(box[3], P[:, 1].max())]
    # an edge is visible if it belongs to a visible face
    vis_edge = [any(e[0] in f and e[1] in f for f in visible_faces) for e in mesh.edges]
    if hidden:
        for e, v in zip(mesh.edges, vis_edge):
            if not v:
                P = view(mesh.X[[e[0], e[1]]])
                ax.plot(P[:, 0], P[:, 1], color=MUTED, lw=0.55 * lw / 0.6, ls=(0, (3, 2)), zorder=3)
    if nodes:
        vis_vertex = set()
        for f in visible_faces:
            vis_vertex |= f
        for k in mesh.vertices:
            P = view(mesh.X[k])[0]
            seen = k in vis_vertex
            if not (seen or hidden):
                continue
            col = C1 if seen else '#9ec5f4'
            if k in open_vertices:
                ax.plot(*P, 's', ms=ms[0] + 0.4, mfc='white', mec=col, mew=0.8, zorder=6 if seen else 4)
            else:
                ax.plot(*P, 's', ms=ms[0], color=col, mec='none', zorder=6 if seen else 4)
        for e, v in zip(mesh.edges, vis_edge):
            if v or hidden:
                P = view(mesh.midpoint(e))[0]
                ax.plot(*P, 'o', ms=ms[1], color=INK2 if v else '#b9b8b2', mec='none', zorder=6 if v else 4)
    return box


def arrow(ax, p0, p1, color=INK, lw=0.7, ms=5.5, zorder=7):
    """Arrow from the screen point p0 to p1."""
    ax.add_patch(FancyArrowPatch(tuple(p0), tuple(p1), arrowstyle='-|>', mutation_scale=ms, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=zorder))


def triad(ax, view, origin, length, labels=('x', 'y', 'z'), fontsize=7):
    """Axes x, y, z drawn from the 3D point origin (arrows of the given length in model units)."""
    o = view.point(*origin)
    for k in range(3):
        d = np.zeros(3)
        d[k] = length
        p = view.point(*(np.asarray(origin) + d))
        arrow(ax, o, p, color=INK2, lw=0.6, ms=5)
        u = (p - o) / np.linalg.norm(p - o)
        ax.text(*(p + 0.22 * np.linalg.norm(p - o) * u), f'${labels[k]}$', fontsize=fontsize, color=INK2,
                ha='center', va='center')


def node_legend(ax, x, y, dy, fontsize=6.6, vertex='vertex node: $u$ and $p_w$', open_vertex=False):
    """Legend of the node symbols at the screen point (x, y), lines dy apart (vertex: text of the vertex nodes, drawn
    as open squares with open_vertex)."""
    if open_vertex:
        ax.plot(x, y, 's', ms=3.0, mfc='white', mec=C1, mew=0.8)
    else:
        ax.plot(x, y, 's', ms=2.6, color=C1, mec='none')
    ax.text(x + 0.6 * abs(dy), y, vertex, fontsize=fontsize, va='center', color=INK)
    ax.plot(x, y - dy, 'o', ms=1.9, color=INK2, mec='none')
    ax.text(x + 0.6 * abs(dy), y - dy, 'mid-edge node: $u$ only', fontsize=fontsize, va='center', color=INK)
