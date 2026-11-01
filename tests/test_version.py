from importlib.metadata import PackageNotFoundError, version

import pytest

import libmat


def test_version_matches_package():
    # libmat.__version__ comes from CMakeLists.txt, the package version from
    # pyproject.toml: both must be bumped together
    try:
        expected = version("libmat")
    except PackageNotFoundError:
        pytest.skip("libmat is not installed (in-place build)")
    assert libmat.__version__ == expected
