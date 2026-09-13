# Wet-lab validated peptide measurements

`df_wetlab_validated.tsv` is a strict, long-form table with the requested schema:

```text
sequence  wetlab_experiment  mic  paper_title  paper_doi
```

It contains 1,551 complete peptide–assay rows, representing 280 distinct peptide
sequences from 10 primary papers in `df_cesar.tsv`. Every DOI and paper title was
matched exactly to that 32-paper catalogue. Each source heatmap peptide × organism
cell is represented independently, so a sequence can occur in several rows when it
was tested against several strains. Rows with missing sequence/result fields,
predicted MICs, no-inhibition/censored results, reviews and preprints are excluded.

`df_cesar_wetlab_round1.tsv` covers all 32 papers (22 primary studies and 10
reviews/perspectives) and records the source-read cohort counts. A blank count means
the abstract/intro or inspected source did not support a defensible number; it does
not mean zero. Complete sequence-to-MIC extraction has not yet been achieved for all
22 primary papers, so this release must not be described as exhaustive across the
entire catalogue.

## AMP Challenge relationship

The AMP Challenge submission rules apply to generated output. The 50,000-sequence
library may not contain an exact sequence from the organizers' `antibacterial.fasta`,
and no top-100 sequence may have a Python `Levenshtein.ratio` above 0.80 to that
reference. Training data must be disclosed.

The September 12, 2026 audit used the official 39,448-sequence reference. Among the
280 unique sequences here:

- 15 are exact reference matches;
- 33 have a maximum Levenshtein ratio above 0.80 (including exact matches);
- 265 pass the full-library no-exact-match rule;
- 247 pass the stricter top-100 similarity rule.

`df_wetlab_challenge_similarity.tsv` records the result for every unique sequence.
This wet-lab table is evidence/training data, not a generated submission library.
Generated candidates must be screened independently with the official Challenge
validator; do not copy the 15 exact matches into the 50,000-sequence output or any
of the 33 >80%-similar sequences into the top 100.

Sources: the paper DOI in each measurement row, the source tables named in
`wetlab_experiment`, and the official AMP Challenge repository and validator:
https://github.com/szczurek-lab/amp-challenge-2027
