import os.path
import webbrowser

import ngsolve
from netgen.meshing import *
from ngsolve.webgui import Draw


def main() -> None:
    ngmesh = Mesh(dim=2)
    N = 5

    point_ids = []
    for i in range(N + 1):
        for j in range(N + 1):
            point_ids.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))

    domain_index = ngmesh.AddRegion("region", dim=2)
    for j in range(N):
        for i in range(N):
            ngmesh.Add(
                Element2D(
                    domain_index,
                    [
                        point_ids[i + j * (N + 1)],
                        point_ids[i + (j + 1) * (N + 1)],
                        point_ids[i + 1 + (j + 1) * (N + 1)],
                        point_ids[i + 1 + j * (N + 1)],
                    ],
                )
            )

    for i in range(N):
        ngmesh.Add(Element1D([point_ids[i], point_ids[i + 1]], index=1))
        ngmesh.Add(
            Element1D([point_ids[N + i * (N + 1)], point_ids[N + (i + 1) * (N + 1)]], index=2)
        )
        ngmesh.Add(Element1D([point_ids[i + N * (N + 1)], point_ids[i + 1 + N * (N + 1)]], index=3))
        ngmesh.Add(
            Element1D([point_ids[0 + i * (N + 1)], point_ids[0 + (i + 1) * (N + 1)]], index=4)
        )

    mesh = ngsolve.Mesh(ngmesh)
    Draw(mesh, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))


if __name__ == "__main__":
    main()
