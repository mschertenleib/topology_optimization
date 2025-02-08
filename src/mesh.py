import os.path
import webbrowser

import ngsolve
from netgen.meshing import *
from ngsolve.webgui import Draw


def create_quad_mesh(size_x: float, size_y: float, nx: int, ny: int) -> Mesh:
    mesh = Mesh(dim=2)

    point_ids = []
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            point_ids.append(mesh.Add(MeshPoint(Pnt(ix / nx * size_x, iy / ny * size_y, 0))))

    domain_index = mesh.AddRegion("region", dim=2)
    for iy in range(ny):
        for ix in range(nx):
            mesh.Add(
                Element2D(
                    domain_index,
                    [
                        point_ids[iy * (nx + 1) + ix],
                        point_ids[iy * (nx + 1) + ix + 1],
                        point_ids[(iy + 1) * (nx + 1) + ix + 1],
                        point_ids[(iy + 1) * (nx + 1) + ix],
                    ],
                )
            )

    for ix in range(nx):
        mesh.Add(Element1D([point_ids[ix], point_ids[ix + 1]], index=1))
        mesh.Add(
            Element1D([point_ids[ny * (nx + 1) + ix], point_ids[ny * (nx + 1) + ix + 1]], index=3)
        )
    for iy in range(ny):
        mesh.Add(
            Element1D([point_ids[iy * (nx + 1) + nx], point_ids[(iy + 1) * (nx + 1) + nx]], index=2)
        )
        mesh.Add(Element1D([point_ids[iy * (nx + 1)], point_ids[(iy + 1) * (nx + 1)]], index=4))

    return mesh


def main() -> None:
    mesh = create_quad_mesh(size_x=3.0, size_y=1.5, nx=20, ny=10)
    mesh = ngsolve.Mesh(mesh)
    Draw(mesh, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))


if __name__ == "__main__":
    main()
