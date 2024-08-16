import os.path
import webbrowser
from ngsolve import Mesh, unit_square
from ngsolve.webgui import Draw


def main():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
    print(mesh.nv, mesh.ne)
    Draw(mesh, filename="out.html")
    webbrowser.open("file://" + os.path.abspath("out.html"))


if __name__ == "__main__":
    main()
