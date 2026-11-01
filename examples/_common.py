"""Small helpers shared by the example scripts (not part of libmat)."""
import argparse
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def data(name):
    return os.path.join(DATA, name)


def output(name):
    os.makedirs(OUT, exist_ok=True)
    return os.path.join(OUT, name)


def parse_args(description, **extra):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--no-gui", action="store_true",
                        help="do not open the polyscope viewer")
    for name, kwargs in extra.items():
        parser.add_argument("--" + name.replace("_", "-"), **kwargs)
    return parser.parse_args()


def polyscope_or_none(args):
    """Returns the initialized polyscope module, or None (no GUI)."""
    if args.no_gui:
        return None
    try:
        import polyscope as ps
    except ImportError:
        print("polyscope is not installed (pip install polyscope), skipping GUI")
        return None
    ps.init()
    ps.set_up_dir("z_up")
    return ps
