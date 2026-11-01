"""Generate one documentation page per libmat function.

    python docs/gen_api_pages.py

For every function listed in CATEGORIES this writes docs/python/<name>.md,
which renders the numpydoc docstring of python/libmat/<name>.py (via
mkdocstrings) followed by the C++ declaration from include/libmat/<name>.h.
It also writes docs/python/index.md (the API overview) and docs/examples.md
(which embeds every script of examples/). Re-run it after adding
a function; the navigation in mkdocs.yml has to be updated by hand.
"""
import glob
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

CATEGORIES = {
    "Input / output": [
        ("read_tet_mesh", "Read a tet mesh (.msh, .tet, .vtk)"),
        ("read_surface_mesh", "Read a triangle mesh (.obj, .off, .ply, .stl, ...)"),
        ("read_ma", "Read a medial mesh (.ma)"),
        ("write_ma", "Write a medial mesh (.ma)"),
    ],
    "Mesh processing": [
        ("normalize_vertices", "Normalize into [0, 1000]^3, the internal frame"),
        ("extract_surface", "Boundary surface of a tet mesh"),
        ("detect_sharp_features", "Sharp / concave edges, feature lines and corners"),
    ],
    "Medial spheres": [
        ("shrink_spheres", "Medial spheres with the ball-shrinking algorithm"),
        ("regular_triangulation", "Power diagram neighbors of spheres (CGAL)"),
    ],
    "Medial mesh": [
        ("slab_tangent_triangles", "Triangles tangent to medial slabs"),
        ("distance_to_medial_mesh", "Signed distance to a medial mesh (CPU / CUDA)"),
        ("thin_medial_mesh", "Topology-preserving thinning"),
        ("medial_mesh_from_spheres", "Medial mesh from spheres via the RPD (CUDA)"),
    ],
}
CUDA_ONLY = {"medial_mesh_from_spheres"}
CUDA_OPTIONAL = {"distance_to_medial_mesh"}


def cpp_declaration(name):
    """Body of `namespace libmat { ... }` of include/libmat/<name>.h."""
    with open(os.path.join(ROOT, "include", "libmat", name + ".h")) as f:
        text = f.read()
    body = re.search(r"namespace libmat \{\n(.*)\n\}  // namespace libmat",
                     text, re.S).group(1).strip("\n")
    return f"#include <libmat/{name}.h>\n\n{body}\n"


def page(name):
    # the heading (with the signature) is rendered by mkdocstrings
    lines = ["---", f"title: {name}", "---", "",
             f"::: libmat.{name}.{name}",
             "    options:", "      show_root_heading: true",
             "      heading_level: 1", "", "## C++", "", "```cpp",
             cpp_declaration(name).rstrip(), "```", ""]
    note = []
    if name in CUDA_ONLY:
        note = ['!!! warning "Requires CUDA"',
                  "    Only available when libmat is built with "
                  "`LIBMAT_WITH_CUDA=ON`, see [Installation](../installation.md).",
                  ""]
    elif name in CUDA_OPTIONAL:
        note = ['!!! tip "CPU and CUDA"',
                  "    Runs on the CPU in every build, and on the GPU "
                  "(`device=\"cuda\"`) when libmat is built with "
                  "`LIBMAT_WITH_CUDA=ON`.", ""]
    # admonition right before the C++ section
    i = lines.index("## C++")
    return "\n".join(lines[:i] + note + lines[i:])


def index():
    lines = ["---", "title: Python API", "---", "", "# Python API", "",
             "All functions take and return numpy arrays, and are available "
             "directly from the `libmat` module:", "",
             "```python", "import libmat",
             'V, F = libmat.read_surface_mesh("data/mesh_69058.obj")',
             "C, R, types = libmat.shrink_spheres(V, F)", "```", "",
             "Every function has a C++ counterpart with the same name in "
             "namespace `libmat`, declared in `<libmat/<name>.h>`.", ""]
    for category, functions in CATEGORIES.items():
        lines += [f"## {category}", "", "| function | description |",
                  "| --- | --- |"]
        for name, desc in functions:
            badge = " :material-expansion-card:{ title=\"CUDA\" }" if name in CUDA_ONLY else ""
            lines.append(f"| [`{name}`]({name}.md){badge} | {desc} |")
        lines.append("")
    lines += ["## Utilities", "",
              "| function | description |", "| --- | --- |",
              "| `libmat.has_cuda()` | `True` if libmat was built with CUDA |",
              "| `libmat.set_verbose(bool)` | print the log of the C++ modules "
              "(silenced by default) |",
              "| `libmat.SPHERE_TYPES` | names of the sphere types returned by "
              "[`shrink_spheres`](shrink_spheres.md) |",
              "| `libmat.__version__` | version string |", ""]
    return "\n".join(lines)


def examples_page():
    lines = ["---", "title: Examples", "---", "", "# Examples", "",
             "The [`examples/`](https://github.com/ningnawang/libmat/tree/main/examples) "
             "folder contains one Python script per topic. Run them from the "
             "repository root; each one opens a [polyscope](https://polyscope.run/py/) "
             "viewer if it is installed (`pip install polyscope`), or pass `--no-gui`:",
             "", "```bash", "python examples/03_shrink_spheres.py", "```", ""]
    for path in sorted(glob.glob(os.path.join(ROOT, "examples", "[0-9]*.py"))):
        name = os.path.basename(path)
        with open(path) as f:
            title = f.read().split('"""')[1].strip().splitlines()[0].rstrip(".")
        lines += [f"## {title}", "", f"`examples/{name}`", "", "```python",
                  f'--8<-- "examples/{name}"', "```", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    out = os.path.join(ROOT, "docs", "python")
    os.makedirs(out, exist_ok=True)
    for functions in CATEGORIES.values():
        for name, _ in functions:
            with open(os.path.join(out, name + ".md"), "w") as f:
                f.write(page(name))
    with open(os.path.join(out, "index.md"), "w") as f:
        f.write(index())
    with open(os.path.join(ROOT, "docs", "examples.md"), "w") as f:
        f.write(examples_page())
    print("wrote", sum(len(v) for v in CATEGORIES.values()) + 1, "pages to", out)
