import os.path
import webbrowser

from netgen.geom2d import SplineGeometry
from netgen.meshing import Element1D, Element2D, FaceDescriptor
from netgen.meshing import Mesh as ngMesh
from netgen.meshing import MeshPoint, Pnt
from ngsolve import *
from ngsolve.webgui import Draw


def create_quad_mesh(size_x: float, size_y: float, nx: int, ny: int) -> ngMesh:
    mesh = ngMesh(dim=2)

    # Grid points
    point_ids = []
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            x = ix / nx * size_x
            y = iy / ny * size_y
            point_ids.append(mesh.Add(MeshPoint(Pnt(x, y, 0))))

    # Region
    region = mesh.AddRegion("rectangle", dim=2)

    # Quad elements
    for iy in range(ny):
        for ix in range(nx):
            mesh.Add(
                Element2D(
                    region,
                    [
                        point_ids[iy * (nx + 1) + ix],
                        point_ids[iy * (nx + 1) + ix + 1],
                        point_ids[(iy + 1) * (nx + 1) + ix + 1],
                        point_ids[(iy + 1) * (nx + 1) + ix],
                    ],
                )
            )

    # Left boundary condition (fixed)
    fd_fix = mesh.Add(FaceDescriptor(surfnr=1, domin=region, bc=1))
    mesh.SetBCName(1, "fix")
    for iy in range(ny):
        mesh.Add(
            Element1D([point_ids[iy * (nx + 1)], point_ids[(iy + 1) * (nx + 1)]], index=fd_fix)
        )

    # Right boundary condition (force)
    fd_force = mesh.Add(FaceDescriptor(surfnr=2, domin=region, bc=2))
    mesh.SetBCName(2, "force")
    for iy in range(ny):
        mesh.Add(
            Element1D(
                [point_ids[iy * (nx + 1) + nx], point_ids[(iy + 1) * (nx + 1) + nx]], index=fd_force
            )
        )

    mesh.Compress()

    return mesh


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

    use_quad_mesh = True

    if use_quad_mesh:
        ny = 5
        nx = int(length / height * ny)
        mesh = Mesh(create_quad_mesh(size_x=length, size_y=height, nx=nx, ny=ny))
    else:
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
