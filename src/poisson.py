import os.path
import webbrowser

from ngsolve import *
from ngsolve.webgui import Draw


def main() -> None:
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

    fes = H1(mesh, order=2, dirichlet="bottom|right|left")

    u = fes.TrialFunction()
    v = fes.TestFunction()
    gfu = GridFunction(fes)

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += 1 * v * dx
    f.Assemble()

    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    Draw(gfu, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))


if __name__ == "__main__":
    main()
