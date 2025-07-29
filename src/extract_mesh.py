from netgen.meshing import Element2D, Element3D, FaceDescriptor
from netgen.meshing import Mesh as NGMesh
from netgen.meshing import MeshPoint
from netgen.meshing import Pnt as NGPnt
from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw


def extract_solid_mesh(mesh: Mesh, rho: GridFunction, threshold: float = 0.5) -> Mesh:
    old_mesh = mesh.ngmesh
    dim = old_mesh.dim
    new_mesh = NGMesh(dim=dim)

    solid = rho.vec.FV().NumPy() >= threshold

    # Vertices references by solid elements
    used_vertices = set()
    elements = old_mesh.Elements3D() if dim == 3 else old_mesh.Elements2D()
    for idx, el in enumerate(elements):
        if solid[idx]:
            used_vertices.update(el.vertices)

    # Vertices
    # TODO: if vertex IDs are sequential, this map could just be a flat array
    vertex_map = {}  # Vertex ID in source mesh -> vertex ID in new mesh
    points = old_mesh.Points()
    for vertex_id in used_vertices:
        x, y, z = points[vertex_id].p
        vertex_map[vertex_id] = new_mesh.Add(MeshPoint(NGPnt(x, y, z)))

    # TODO: we might want to preserve the volume regions from the original mesh
    region = new_mesh.AddRegion("default", dim=dim)

    # Volume elements
    for idx, el in enumerate(elements):
        if not solid[idx]:
            continue
        vertices = [vertex_map[v] for v in el.vertices]
        if dim == 3:
            new_mesh.Add(Element3D(region, vertices))
        else:
            new_mesh.Add(Element2D(region, vertices))

    # Surface elements for 3D mesh
    # FIXME: This is not sufficient, this only keeps outer faces that were already present in the
    # original mesh
    # TODO: we might want to preserve the face descriptors/boundary conditions from the original
    # mesh
    if dim == 3:
        fd = new_mesh.Add(FaceDescriptor(bc=1, domin=1))
        for face in old_mesh.Elements2D():
            if all(v in used_vertices for v in face.vertices):
                new_mesh.Add(Element2D(fd, [vertex_map[v] for v in face.vertices]))

    return Mesh(new_mesh)


dim = 2
if dim == 3:
    beam = Box(Pnt(0, 0, 0), Pnt(1, 0.5, 0.5))
else:
    beam = Rectangle(1, 0.5).Face()
geo = OCCGeometry(beam, dim=dim)
mesh = Mesh(geo.GenerateMesh(maxh=0.05 if dim == 3 else 0.01))

print(f"Number of elements: {mesh.ne}")

fes = L2(mesh, order=0)
rho = GridFunction(fes)

sigma = 0.05
if dim == 3:
    cf = (
        exp(-((x - 0.5) ** 2) / sigma**2)
        + exp(-((y - 0.25) ** 2) / sigma**2)
        + exp(-((z - 0.25) ** 2) / sigma**2)
    )
else:
    cf = exp(-((x - 0.5) ** 2) / sigma**2) + exp(-((y - 0.25) ** 2) / sigma**2)
rho.Set(IfPos(cf - 1, 1, cf).Compile())

print("Extracting solid mesh...")

new_mesh = extract_solid_mesh(mesh, rho, threshold=0.5)
new_rho = GridFunction(H1(new_mesh, order=1))
new_rho.Set(rho)

Draw(
    new_rho,
    new_mesh,
    filename="out.html",
    settings={
        "Objects": {"Wireframe": False, "Edges": False},
        "Light": {"shininess": 0.0, "specularity": 0.0},
    },
)
