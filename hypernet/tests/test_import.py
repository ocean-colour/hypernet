""" Smoke tests: the package installs and imports cleanly """

def test_import():
    import hypernet


def test_import_wiggles_modules():
    """Every Phase 0 script imports cleanly (no work at import time)."""
    import importlib
    import pkgutil
    import hypernet.wiggles as w
    names = [m.name for m in pkgutil.iter_modules(w.__path__)]
    assert 'phase0a_veit' in names and 'phase0b_index' in names
    for n in names:
        mod = importlib.import_module('hypernet.wiggles.' + n)
        assert hasattr(mod, 'main'), n


def test_import_wiggles_without_os_color():
    """The package imports even when $OS_COLOR is unset (e.g. a fresh install)."""
    import subprocess
    import sys
    import os
    env = {k: v for k, v in os.environ.items() if k != 'OS_COLOR'}
    code = ('import importlib, pkgutil, hypernet.wiggles as w; '
            '[importlib.import_module("hypernet.wiggles." + m.name) '
            'for m in pkgutil.iter_modules(w.__path__)]')
    r = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    assert r.returncode == 0, r.stderr
