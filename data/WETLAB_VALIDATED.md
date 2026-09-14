# Wet-lab validated antimicrobial-peptide release

## Deliverable

`df_wetlab_validated.tsv` is the mentor-facing table with the requested schema:

```text
sequence  wetlab_experiment  mic  paper_title  paper_doi
```

The release contains 2,449 peptide–assay rows, 596 unique amino-acid sequences,
and all 22 primary-paper DOIs in `df_cesar.tsv`. Reviews and perspectives are
represented only in the 32-paper Round-1 audit. Preprints and predicted MICs are
excluded from the validated table.

Of the 2,449 rows, 2,342 contain an exact source-reported MIC. The other 107
retain a blank MIC rather than an inferred number: 55 human-microbiome SEPs, 36
CRPS wasp derivatives, ten AMPSphere mouse candidates, five dF-AndD1 peptides,
and Mast-MO. Their sequence-specific antimicrobial evidence and source figure or
table are recorded in `wetlab_experiment`.

## What was extracted

- `df_wetlab_validated.tsv` — positive antimicrobial evidence only, in the five
  requested columns.
- `df_wetlab_heatmap_results.tsv` — 6,448 source assay cells, including censored
  and no-inhibition results useful for negative/control analyses.
- `df_cesar_wetlab_round1.tsv` — all 32 final journal papers, classified as 22
  primary studies and ten reviews/perspectives, with source-read cohort counts.
- `df_cesar_murine_validated.tsv` — 38 source-linked murine-efficacy records;
  every row now has a resolved sequence.

Exact tables or publisher source workbooks were used for human-proteome mining,
molecular de-extinction, Venomics, Archaeasins, non-immune proteins, ApoB,
MMP-19, prionins, ApexGO, key-cutting, the scorpion hybrids, and the nanopore
designs. Figure-supported activity was retained without fabricated MIC values
for the remaining cohorts. Terminal amidation and all-D stereochemistry are
encoded in `wetlab_experiment` because the requested `sequence` field is the
canonical one-letter amino-acid string.

## Scope and limitations

This is complete at the paper-coverage level, not a claim that publishers expose
every reported active identity in machine-readable form. Two source limitations
remain explicit:

1. AMPSphere reports 79 active candidates among 100 tested. The public material
   resolves accession/sequence mappings for the ten candidates advanced to the
   mouse experiment, but not for the other 69 active display names. Only the ten
   defensibly mapped sequences are included.
2. AMP-Diffusion reports 35 active candidates among 46 synthesized. Twenty-six
   sequence identities have unambiguous positive cells in the published Figure
   2a heatmap and are included; nine reported-active identities cannot be
   assigned without guessing.

These gaps are documented rather than silently filled. A publisher-provided
mapping table or author response can be appended later without changing the
five-column schema.

## Quality controls

- All sequence strings contain only the 20 standard one-letter amino-acid codes.
- Every DOI and title matches the 32-paper César catalogue.
- All 22 primary-paper DOIs occur in the validated file.
- Censored/no-inhibition cells are excluded from the positive-only release and
  retained in the companion heatmap-results table.
- No unresolved sequence remains in the murine-efficacy table.
- Rows are deterministically sorted and duplicate rows are removed.

Use this as curated evidence/training data, not as a generated AMP Challenge
submission. Cite the original DOI associated with every reused measurement and
conduct a final domain-expert spot check before publication or model benchmarking.
