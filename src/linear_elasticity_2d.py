import os.path
import webbrowser

import numpy as np
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.webgui import Draw


def analytical_beam_deflection(width: float, height: float, length: float, E: float, force: float):
    I_z = width * height**3 / 12.0
    return force * length**3 / (3.0 * E * I_z)


def stress(strain, mu, lam):
    return 2.0 * mu * strain + lam * Trace(strain) * Id(strain.shape[0])


def strain(displacement):
    return Sym(Grad(displacement))


def main() -> None:
    # NOTE: All values in standard units: m, N, Pa
    length = 0.2
    height = 0.02
    width = 0.03
    force = -100.0
    E = 70e9  # Young's modulus
    nu = 0.35  # Poisson's ratio

    geo = SplineGeometry()
    p1 = geo.AppendPoint(0, 0)
    p2 = geo.AppendPoint(length, 0)
    p3 = geo.AppendPoint(length, height)
    p4 = geo.AppendPoint(0, height)
    geo.Append(["line", p1, p2])
    geo.Append(["line", p2, p3], bc="force")
    geo.Append(["line", p3, p4])
    geo.Append(["line", p4, p1], bc="fix")
    mesh = Mesh(geo.GenerateMesh(maxh=height / 5.0))

    # Lamé parameters
    lam = E * nu / (1.0 - nu**2)  # plane-stress
    mu = E / (2.0 * (1.0 + nu))

    fes = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(InnerProduct(stress(strain(u), mu, lam), strain(v)) * dx)
    a.Assemble()

    f = LinearForm(CoefficientFunction((0, force / (width * height))) * v * ds("force"))
    f.Assemble()

    inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data = inv * f.vec

    # Draw(gfu, filename="out.html")
    # webbrowser.open("file://" + os.path.abspath("out.html"))

    analytical_deflection = analytical_beam_deflection(
        width=width, height=height, length=length, E=E, force=force
    )
    print(f"Analytical Y deflection: {analytical_deflection:.9f} m")

    point = mesh(length, height / 2, VOL_or_BND=BND)
    numerical_deflection = gfu(point)[1]
    print(
        f"Numerical Y deflection:  {numerical_deflection:.9f} m"
        f" ({(numerical_deflection - analytical_deflection) / analytical_deflection * 100.0:+.2f}%)"
    )

    total_force = (
        Integrate(CoefficientFunction(force / (width * height)) * ds("force"), mesh) * width
    )
    print(f"Total integrated force: {total_force:.3f} N")


if __name__ == "__main__":
    main()
