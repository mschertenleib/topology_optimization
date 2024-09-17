import os.path
import webbrowser

from ngsolve import *
from netgen.occ import *
from ngsolve.webgui import Draw

import matplotlib.pyplot as plt
import numpy as np


def analytical_beam_deflection(width: float, height: float, length: float, E: float, force: float):
    I_z = width * height**3 / 12.0
    return force * length**3 / (3.0 * E * I_z)


def main() -> None:
    length = 0.2
    height = 0.02
    width = 0.02
    force = -100.0
    E = 70e9  # Young's modulus
    nu = 0.35  # Poisson's ratio

    beam = Box(Pnt(0, 0, 0), Pnt(length, height, width))
    beam.faces.Min(X).name = "fix"
    beam.faces.Max(X).name = "force"

    geo = OCCGeometry(beam)
    mesh = Mesh(geo.GenerateMesh(maxh=height / 5.0))

    # Lamé parameters
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))

    def stress(strain):
        return 2.0 * mu * strain + lam * Trace(strain) * Id(3)

    def strain(displacement):
        return Sym(Grad(displacement))

    fes = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(InnerProduct(stress(strain(u)), strain(v)) * dx)
    a.Assemble()

    f = LinearForm(CoefficientFunction((0, force, 0)) * v * ds("force"))
    f.Assemble()

    inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data = inv * f.vec

    Draw(gfu, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))

    deflection = analytical_beam_deflection(width=width, height=height, length=length, E=E, force=force)
    print(f"Analytical deflection: {deflection} m")


if __name__ == "__main__":
    main()
