import os.path
import webbrowser

from ngsolve import *
from netgen.geom2d import SplineGeometry
from ngsolve.webgui import Draw

import matplotlib.pyplot as plt
import numpy as np

"""
[From Gere and Goodno, "Mechanics of Materials"]

Deflection for a horizontal beam with a fixed extremity and a vertical
load applied at the other extremity (gravity neglected):

v(x) = -F*x^2 / (6*E*I) * (3*L - x)

For a rectangular cross-section beam, I = w*h^3/12


"""


def main() -> None:
    geo = SplineGeometry()

    length = 1.0
    height = 0.1
    width = 0.1
    force = -100.0

    p1, p2, p3, p4 = [geo.AppendPoint(x, y) for x, y in [(0, 0), (length, 0), (length, height), (0, height)]]
    geo.Append(["line", p1, p2], bc="bottom")
    geo.Append(["line", p2, p3], bc="right")
    geo.Append(["line", p3, p4], bc="top")
    geo.Append(["line", p4, p1], bc="left")
    mesh = Mesh(geo.GenerateMesh(maxh=height / 5.0))

    E = 70e9  # Young's modulus
    nu = 0.35  # Poisson's ratio

    # Lamé parameters
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))

    def stress(strain):
        return 2.0 * mu * strain + lam * Trace(strain) * Id(2)

    def strain(displacement):
        return Sym(Grad(displacement))

    fes = VectorH1(mesh, order=2, dirichlet="left")
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(InnerProduct(stress(strain(u)), strain(v)) * dx)
    a.Assemble()

    f = LinearForm(CoefficientFunction((0, force)) * v * ds("right"))
    f.Assemble()

    inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data = inv * f.vec

    Draw(gfu, filename="out.html")
    # webbrowser.open("file://" + os.path.abspath("out.html"))

    # Comparison with analytical solution
    x = np.linspace(0.0, 1.0, 1000)
    I_z = width * height**3 / 12.0
    v = force * x * x / (6.0 * E * I_z) * (3.0 * length - x)
    fig, ax = plt.subplots()
    ax.plot(x, v)
    fig.savefig("fig.png")


if __name__ == "__main__":
    main()
