__version__ = "0.1.0"

import os as _os

# Earlier names of the tool (PurityAgent, then ngsPurityMap): settings made under them keep working.
#   PURITYAGENT_* / NGSPURITYMAP_* environment variables → PURITYMAP_* (the new name wins when both are set)
OLD_NAMES = ("purityagent", "ngspuritymap")
for _old in OLD_NAMES:
    for _k, _v in list(_os.environ.items()):
        if _k.startswith(_old.upper() + "_"):
            _os.environ.setdefault("PURITYMAP_" + _k[len(_old) + 1:], _v)
