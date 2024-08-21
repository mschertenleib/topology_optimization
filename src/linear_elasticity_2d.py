import os.path
import webbrowser

from ngsolve import *
from netgen.geom2d import SplineGeometry
from ngsolve.webgui import Draw


def main() -> None:
    geo = SplineGeometry()

    p1, p2, p3, p4 = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 0.1), (0, 0.1)]]
    geo.Append(["line", p1, p2], bc="bottom")
    geo.Append(["line", p2, p3], bc="right")
    geo.Append(["line", p3, p4], bc="top")
    geo.Append(["line", p4, p1], bc="left")
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))

    E = 70  # Young's modulus
    nu = 0.35  # Poisson's ratio

    # Lamé parameters
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))

    def stress(strain):
        return 2.0 * mu * strain + lam * Trace(strain) * Id(2)

    def strain(u):
        return Sym(Grad(u))  # NOTE: Sym(A) = 0.5 * (A + A^T)

    fes = VectorH1(mesh, order=2, dirichlet="left")
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(InnerProduct(stress(strain(u)), strain(v)) * dx)
    a.Assemble()

    f = LinearForm(CoefficientFunction((0, -5e-3)) * v * ds("right"))
    f.Assemble()

    inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data = inv * f.vec

    Draw(gfu, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))


if __name__ == "__main__":
    main()
