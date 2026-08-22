# arXiv submission

The upload-ready source bundle is [`../paper/arxiv-submission.tar.gz`](../paper/arxiv-submission.tar.gz).
It is the **only** file you upload.

## Verified (clean-room)
Extracted to a fresh directory and compiled with `pdflatex` alone (arXiv's workflow):
- **18 pages, 0 undefined references, 0 missing files, 0 fatal errors.**
- Source is **100 % ASCII**; `\pdfoutput=1` forces pdfLaTeX (pgfplots + PNG raster fields).
- `main.bbl` is **included**, so arXiv uses it and does **not** run BibTeX (`refs.bib` is omitted).
- Bundle root contains `main.tex`, `figstyle.tex`, `main.bbl`, `figs/` — no `main.pdf`, no aux.

## Metadata (copy-paste)
- **Title:** Self-falsifying quantum spectroscopy: a transportable necessary-condition screen for quantum-computed dynamical spectra
- **Authors:** Nicolás Bonilla Vargas
- **Primary category:** `quant-ph`  ·  **Cross-list:** `cond-mat.str-el`
- **License:** CC BY 4.0
- **Comments:** `18 pages, 10 figures. Companion to arXiv:2608.16436.`
- **Abstract:** the condensed, ASCII, `<1920`-char version is in the repo's release notes; the full
  abstract lives in the PDF.

## Steps
1. arxiv.org → **Submit** → license **CC BY 4.0**.
2. Upload `arxiv-submission.tar.gz`; let AutoTeX build (it uses `main.bbl`; no BibTeX — expected).
3. Paste the metadata above.
4. Preview the arXiv-generated PDF end to end.
5. Submit; note the assigned `arXiv:XXXX.XXXXX` and update the badge / `CITATION.cff` / `README.md` here
   and the companion cross-links in the two companion repositories.

Endorsement is **not** needed for `quant-ph` once the companion `arXiv:2608.16436` has announced.
