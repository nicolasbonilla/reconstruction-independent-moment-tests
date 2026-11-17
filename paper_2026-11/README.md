# The present manuscript (November 2026)

*What low-order sum rules can and cannot detect in sample-based spectral functions*, N. Bonilla Vargas. This folder
holds the LaTeX sources of the manuscript as of deposit 2.0.3 and three scripts that check printed numbers of that
text against the deposited records. The earlier manuscript (October 2026, *Reconstruction-independent moment tests
for sample-based spectral functions*) stays in [`../paper/`](../paper/).

## Files

- `main.tex`: the manuscript in the two-column REVTeX layout; `main_jpcm.tex`: the same manuscript in a
  single-column layout at 12 pt (it inputs `main.tex`). Besides the layout, the switch between the two changes the
  name given to the Supplemental Material, the way it is pointed to and the label of the e-mail address; in the
  single-column layout the note on the Supplemental Material cites no references, so its reference list is shorter.
  `supplement.tex`: the Supplemental Material; `refs.bib`; `figstyle.tex` (shared figure style).
- The figures are the files of [`../paper/figs/`](../paper/figs/): the sources read `figs/...`, so to compile,
  copy that folder here as `figs/` (see Build below).
- `verify_prb.py`: runs [`../src/verify.py`](../src/verify.py) unchanged (its recomputation of the `L=6` moments,
  the thresholds and the `d=97` catch, then its record check, here on `main.tex` and `supplement.tex` of this
  folder) and then `closing_records.py`. The anchors of `../src/verify.py` were written for the earlier manuscript:
  on this text 74 printed numbers match their records and five anchors are reported as `ABSENT` (sentences this
  text words differently); the numbers of four of them (eight numbers) are checked by `closing_records.py`, and the
  fifth anchor's numbers are checked by another anchor of the same record.
- `closing_records.py`: makes 37 checks under 15 anchors of its own (passages of this text that `../src/verify.py`
  does not anchor in their present form), against five deposited data files. 35 of the checks compare a printed
  value with its record, with the direction of rounding that its sentence requires (an upper bound is rounded up,
  the lower end of an enclosure down, a point value to nearest). The other two read no text: one checks that the
  two maxima stored in the summary of the inverse-moment record (the relative shifts of `m_-1` and `m_1`) equal the
  values at its last sweep point, the other that the Krylov record holds the default-rule thresholds of orders
  1 to 10 in each of its 16 rows (160 entries; the four rows of one U/t carry the same values), the set over which
  the threshold range is stated.
- `check_krylov_numbers.py`: checks numbers that the text prints from the two records of 2026-11-08. From
  `../data/2026-11-08_krylov_detection_order.json`: numbers of the Krylov detection-order analysis (Sec. V B of
  the manuscript, its Table III, and Table SII and the threshold range of the Supplemental Material). From
  `../data/2026-11-08_operator_form_tie_scan.json`: three strings of the operator-form example (Sec. II C), which
  carry the residual 0.81 and the first-moment error 0.48, the 11% deficit of the deterministic cut, and the
  8.7-12.0% range of that deficit over the draws of the tie-order scan; the script matches no other number of that
  example in the text. Of its 43 checks, 42 rebuild a printed string from one of the two records and require it
  verbatim in the text; the remaining one reads no text: it checks that the tie-scan record counts the same two
  verdicts (the `O_1` form rejects, the commutator form does not) in all of its 300 draws. The script does not
  read the plotted points of Fig. 5: those are the `.dat` files that `../src/export_krylov_order_dat.py` writes
  from the Krylov record.

## Run

```bash
python verify_prb.py              # from this folder; exits 0 when every check passes (seconds)
python check_krylov_numbers.py    # exits 0 when every check passes (seconds)
```

The three scripts check the numbers listed above, not every number of the text.

The two programs that write the data files of Fig. 5 and the rows of Table SII from the Krylov record are
`../src/export_krylov_order_dat.py` and `../src/export_krylov_table_tex.py` (`make figures` runs them).

## Build

```bash
cp -r ../paper/figs figs          # from this folder: the sources read figs/...
for d in supplement main main_jpcm; do pdflatex $d; done
for d in supplement main main_jpcm; do bibtex $d; done
for d in supplement main main_jpcm; do pdflatex $d; done
for d in supplement main main_jpcm; do pdflatex $d; done
bibtex main
for d in supplement main main_jpcm; do pdflatex $d; done
for d in supplement main main_jpcm; do pdflatex $d; done
```

The three documents cross-reference each other, so they are compiled in turn. `main` needs the second `bibtex` run:
its note on the Supplemental Material lists references that reach BibTeX only after that note has been typeset.
Run from a clean copy of this folder on 2026-11-16 (pdfTeX of MiKTeX 25.12), these commands give `main.pdf` (16
pages), `main_jpcm.pdf` (34 pages) and `supplement.pdf` (26 pages) with no LaTeX error, no overfull box and no
undefined citation or reference; BibTeX prints one warning per document, `jnrlst (dependency: not reversed) set 1`,
which comes from the REVTeX bibliography style.
