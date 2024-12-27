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


def stress(strain, mu, lam):
    return 2.0 * mu * strain + lam * Trace(strain) * Id(3)


def strain(displacement):
    return Sym(Grad(displacement))


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

    fes = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(InnerProduct(stress(strain(u), mu=mu, lam=lam), strain(v)) * dx)
    a.Assemble()

    f = LinearForm(CoefficientFunction((0, force / (width * height), 0)) * v * ds("force"))
    f.Assemble()

    inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data = inv * f.vec

    Draw(gfu, filename="out.html")
    # webbrowser.open("file://" + os.path.abspath("out.html"))

    analytical_deflection = analytical_beam_deflection(width=width, height=height, length=length, E=E, force=force)
    print(f"Analytical Y deflection: {analytical_deflection:.9f} m")

    coords = np.asarray([node.point for node in mesh.vertices])
    disp = gfu(mesh(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2]))
    numerical_deflection = np.mean(disp[coords[:, 0] == length, 1])
    print(f"Numerical Y deflection:  {numerical_deflection:.9f} m")


if __name__ == "__main__":
    main()
