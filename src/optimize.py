from netgen.meshing import Element1D, Element2D, FaceDescriptor
from netgen.meshing import Mesh as ngMesh
from netgen.meshing import MeshPoint, Pnt
from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from tqdm import tqdm


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


def helmholtz_filter(fes, radius):
    u = fes.TrialFunction()
    v = fes.TestFunction()
    A = BilinearForm(u * v * dx + radius**2 * InnerProduct(grad(u), grad(v)) * dx)
    A.Assemble()
    A_inv = A.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")

    M = BilinearForm(u * v * dx)
    M.Assemble()

    return A_inv, M


# FIXME
def optimality_criterion(x, dc, dV, move, volume_fraction, mesh):
    x_lower = x - move
    x_upper = x + move
    sensitivity_frac = -dc / dV
    oc_constant = x * sqrt(IfPos(sensitivity_frac, sensitivity_frac, 0))
    mesh_volume = Integrate(1 * dx, mesh)
    lm_lower = 0
    lm_upper = Integrate(oc_constant * dx) / mesh_volume / volume_fraction
    while (lm_upper - lm_lower) / (lm_upper + lm_lower) > 1e-4:
        lm_middle = 0.5 * (lm_lower + lm_upper)
        x = oc_constant / lm_middle
        x = IfPos(x - 1.0, 1.0, IfPos(-x, 0.0, x))
        x = IfPos(x - x_upper, x_upper, IfPos(x_lower - x, x_lower, x))
        if Integrate(x * dx, mesh) / mesh_volume > volume_fraction:
            lm_lower = lm_middle
        else:
            lm_upper = lm_middle

    return x


def main() -> None:
    # NOTE: All values in standard units: m, N, Pa
    length = 0.2
    height = 0.1
    width = 0.01
    force = -100.0
    E_0 = 70e9  # Young's modulus of solid
    E_min = 1e-9 * E_0  # Young's modulus of void
    nu = 0.35  # Poisson's ratio
    penalty = 3
    filter_radius = 0.005
    num_steps = 100

    use_quad_mesh = False
    use_bracket = False

    # Mesh
    if use_quad_mesh:
        ny = 50
        nx = int(length / height * ny)
        mesh = Mesh(quad_mesh(size_x=length, size_y=height, nx=nx, ny=ny))
    elif use_bracket:
        r1 = 0.045
        r2 = 0.03
        r1_hole = 0.6 * r1
        r2_hole = 0.6 * r2
        length = 0.25
        rod_p1 = (r1 * 0.4, r1 * 0.85)
        rod_p2 = (r1 * 0.9, r1 * 0.1)
        rod_p3 = (length - r2 * 0.9, r2 * 0.1)
        rod_p4 = (length - r2 * 0.4, r2 * 0.85)
        eye_1 = Circle(Pnt(0, 0), r=r1).Face()
        eye_2 = Circle(Pnt(length, 0), r=r2).Face()
        hole_1 = Circle(Pnt(0, 0), r=r1_hole).Face()
        hole_2 = Circle(Pnt(length, 0), r=r2_hole).Face()
        rod_1 = (
            MoveTo(*rod_p1).LineTo(*rod_p2).LineTo(*rod_p3).LineTo(*rod_p4).Close().Face()
            - eye_1
            - eye_2
        )
        rod_2 = (
            (
                MoveTo(rod_p1[0], -rod_p1[1])
                .LineTo(rod_p2[0], -rod_p2[1])
                .LineTo(rod_p3[0], -rod_p3[1])
                .LineTo(rod_p4[0], -rod_p4[1])
                .Close()
                .Reverse()
                .Face()
            )
            - eye_1
            - eye_2
        )
        bracket = Glue([eye_1 - hole_1, eye_2 - hole_2, rod_1, rod_2])
        # bracket.edges.Min(X).name = "fix"
        # bracket.edges.Max(X).name = "force"
        hole_1.edges.name = "fix"
        hole_2.edges.name = "force"
        geo = OCCGeometry(bracket, dim=2)
        mesh = Mesh(geo.GenerateMesh(maxh=0.004))
        # Draw(mesh, filename="out.html")
        # exit()
    else:
        beam = Rectangle(length, height).Face()
        beam.edges.Min(X).name = "fix"
        beam.edges.Max(X).name = "force"
        geo = OCCGeometry(beam, dim=2)
        mesh = Mesh(geo.GenerateMesh(maxh=height / 50.0))

    # Design space (element-wise constant)
    fes_rho = L2(mesh, order=0)
    rho = GridFunction(fes_rho)
    rho_filtered = GridFunction(fes_rho)

    # Filter space (continuous)
    fes_rho_h1 = H1(mesh, order=1)
    rho_h1 = GridFunction(fes_rho_h1)
    rho_filtered_h1 = GridFunction(fes_rho_h1)

    E = E_min + rho_filtered**penalty * (E_0 - E_min)
    lam, mu = lame_parameters(E, nu)

    A_inv, M = helmholtz_filter(fes_rho_h1, filter_radius)

    # Displacement space
    fes_u = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes_u.TrialFunction()
    v = fes_u.TestFunction()
    u_gf = GridFunction(fes_u)

    # Linear elasticity system
    a = BilinearForm(InnerProduct(stress(strain(u), mu, lam), strain(v)).Compile() * dx)
    f = LinearForm(CoefficientFunction((0, force / (width * height))) * v * ds("force"))
    f.Assemble()

    free_dofs = fes_u.FreeDofs()

    # NOTE: temporary explicit density
    step_param = Parameter(0)
    center = 0.5 * length + 0.3 * length * sin(2 * pi * step_param / 15)
    density = IfPos(x - center + 0.03 * length, IfPos(x - center - 0.03 * length, 1.0, 0.1), 1.0)

    with TaskManager():
        for step in tqdm(range(num_steps)):
            step_param.Set(step)
            rho.Set(density)

            rho_h1.Set(rho)
            rho_filtered_h1.vec.data = A_inv * (M.mat * rho_h1.vec)
            rho_filtered.Set(rho_filtered_h1)

            a.Assemble()

            inv = a.mat.Inverse(freedofs=free_dofs, inverse="sparsecholesky")
            u_gf.vec.data = inv * f.vec

            if True:
                Draw(
                    E / E_0,
                    mesh,
                    filename="out.html",
                    settings={"Colormap": {"ncolors": 32}},
                )
            else:
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
