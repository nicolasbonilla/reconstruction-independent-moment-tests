# -*- coding: utf-8 -*-
"""Check the present manuscript (main.tex + supplement.tex of this folder) against the deposited records.
Part 1 runs the repository's own check, ../src/verify.py, unchanged: its physics part (moments, thresholds, the d=97
catch) is recomputed from source, and its record check, which normally reads ../paper/main.tex, is given the text of
this folder instead, so every anchored sentence is searched in the main text and its Supplemental Material together.
Anchors written for sentences that this text words differently are reported as ABSENT. Runs in seconds.
Part 2 runs closing_records.py on the same text: numbers that this text prints in a form ../src/verify.py does not
anchor, each checked against its JSON key with the direction of rounding that its sentence requires (see that file).
Exit code 1 if either part fails."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, '..', 'src'))
PRB = HERE
sys.path.insert(0, SRC)
os.chdir(SRC)
import verify  # noqa: E402

_orig = verify._load


def _load(rel):
    if rel.replace('\\', '/') == 'paper/main.tex':
        return ''.join(open(os.path.join(PRB, f), encoding='utf-8').read() + '\n' for f in ('main.tex', 'supplement.tex'))
    return _orig(rel)


verify._load = _load
sys.path.insert(0, HERE)
import closing_records  # noqa: E402

code = 0
try:
    verify.main()
except SystemExit as e:
    code = e.code if isinstance(e.code, int) else 1
code = closing_records.report(_load('paper/main.tex')) or code
sys.exit(code)
