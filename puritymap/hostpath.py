# =============================================================================
# hostpath.py   User machine paths ↔ container paths
#
#   The launcher (purity / purity.bat) mounts the user's drives read-only:
#     Linux · macOS   /          → /host
#     Windows         C:\  D:\   → /host/c  /host/d
#   and passes PURITYMAP_HOST_ROOT=/host. Users write their own machine paths in the settings and sheet,
#   and this module maps them to container paths. Without the variable (running outside a container),
#   paths are left unchanged.
# =============================================================================

import os
import re

_DRIVE = re.compile(r"^([A-Za-z]):[\\/]?(.*)$")


def host_root() -> str:
    return os.environ.get("PURITYMAP_HOST_ROOT", "").rstrip("/")


def to_local(path: str) -> str:
    """Path as written by the user → path this process can open."""
    root = host_root()
    if not root:
        return path
    if path.startswith("\\\\"):
        raise ValueError(f"network paths (\\\\server\\share) are not supported → use a mapped drive letter: {path}")
    m = _DRIVE.match(path)
    if m:                                           # Windows drive path
        rest = m.group(2).replace("\\", "/").lstrip("/")
        return f"{root}/{m.group(1).lower()}/{rest}".rstrip("/")
    if path.startswith("/"):                        # Linux · macOS absolute path
        return path if path == root or path.startswith(root + "/") else root + path
    return path                                     # relative path: the caller joins it to its base folder


def display(text: str) -> str:
    """Turn container paths inside a message back into the user's machine paths."""
    root = host_root()
    if not root:
        return text
    if os.environ.get("PURITYMAP_HOST_OS") == "windows":
        return re.sub(re.escape(root) + r"/([a-z])(?=/|\b)", lambda m: m.group(1).upper() + ":", text)
    return re.sub(re.escape(root) + r"(?=/)", "", text)
