import os.path
import webbrowser

from netgen.meshing import Element1D, Element2D, FaceDescriptor
from netgen.meshing import Mesh as ngMesh
from netgen.meshing import MeshPoint, Pnt
from netgen.occ import OCCGeometry, Rectangle, X
from ngsolve import *
from ngsolve.webgui import Draw


def quad_mesh(size_x: float, size_y: float, nx: int, ny: int) -> ngMesh:
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

    # Boundary conditions
    fd_fix = mesh.Add(FaceDescriptor(surfnr=1, domin=region, bc=1))
    fd_force = mesh.Add(FaceDescriptor(surfnr=2, domin=region, bc=2))
    mesh.SetBCName(1, "fix")
    mesh.SetBCName(2, "force")
    for iy in range(ny):
        mesh.Add(
            Element1D([point_ids[iy * (nx + 1)], point_ids[(iy + 1) * (nx + 1)]], index=fd_fix)
        )
        mesh.Add(
            Element1D(
                [point_ids[iy * (nx + 1) + nx], point_ids[(iy + 1) * (nx + 1) + nx]], index=fd_force
            )
        )

    mesh.Compress()

    return mesh


def lame_parameters(E, nu):
    lam = E * nu / (1.0 - nu**2)  # plane-stress
    mu = E / (2.0 * (1.0 + nu))
    return lam, mu


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

    rho_0 = 1.0
    rho_min = 1e-9
    penalty = 3

    use_quad_mesh = False

    if use_quad_mesh:
        ny = 5
        nx = int(length / height * ny)
        mesh = Mesh(quad_mesh(size_x=length, size_y=height, nx=nx, ny=ny))
    else:
        beam = Rectangle(length, height).Face()
        beam.edges.Min(X).name = "fix"
        beam.edges.Max(X).name = "force"
        geo = OCCGeometry(beam, dim=2)
        mesh = Mesh(geo.GenerateMesh(maxh=height / 5.0))

    lam, mu = lame_parameters(E, nu)

    fes_u = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes_u.TrialFunction()
    v = fes_u.TestFunction()
    u_gf = GridFunction(fes_u)

    fes_rho = L2(mesh, order=0)
    rho_gf = GridFunction(fes_rho)

    step_param = Parameter(0)
    density = (
        0.9 / (0.2 * length) * Norm(x - 0.5 * length - 0.3 * length * sin(2 * pi * step_param / 15))
        + 0.1
    )
    density = IfPos(density - 1, 1.0, density)

    rho_factor = rho_min + rho_gf**penalty * (rho_0 - rho_min)

    a = BilinearForm(
        InnerProduct(stress(strain(u), rho_factor * mu, rho_factor * lam), strain(v)) * dx
    )
    f = LinearForm(CoefficientFunction((0, force / (width * height))) * v * ds("force"))
    f.Assemble()
    free_dofs = fes_u.FreeDofs()

    for step in range(100):
        step_param.Set(step)
        rho_gf.Set(density)

        a.Assemble()

        inv = a.mat.Inverse(freedofs=free_dofs, inverse="sparsecholesky")
        u_gf.vec.data = inv * f.vec

        print(f"Step {step}")

        eps_cf = strain(u_gf)
        eps_dev = eps_cf - Trace(eps_cf) / 2 * Id(2)
        eq_strain = sqrt(2 / 3 * InnerProduct(eps_dev, eps_dev))
        Draw(
            eq_strain,
            mesh,
            deformation=2 * u_gf,
            filename="out.html",
            settings={"Colormap": {"ncolors": 32}},
        )

        import time

        time.sleep(0.5)


if __name__ == "__main__":
    main()
