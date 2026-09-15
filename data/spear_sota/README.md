# SPEAR SOTA wet-lab and segmentation release

This directory separates two questions that must not be conflated:

1. **Did a synthesized compound show antimicrobial activity?** Use the assay
   tables. These include synthetic and chemically modified compounds.
2. **Is this exact native sequence a defensible cleavage span inside a longer
   biological parent?** Use the strict segmentation gold table.

## What to use

- `SPEAR_596_wetlab_validated.tsv` is the broad one-row-per-sequence view. All
  596 sequences have positive wet-lab evidence. Modified and synthetic peptides
  are intentionally retained; repeated assays and papers are aggregated.
- `SPEAR_596_activity_peptides.tsv` is the model index for the activity task.
  It assigns each sequence to a 60%-identity peptide cluster and keeps every
  cluster within one train/validation/test split.
- `SPEAR_2449_activity_measurements.tsv` is the linked target-conditioned MIC
  table. Each row is one peptide-assay-target observation, with normalized
  target text, numeric MIC fields, a regression-eligibility flag, and a sample
  weight that prevents heavily tested peptides from dominating each unit task.
- `SPEAR_270_strict_segmentation.tsv` is the strict one-row-per-label view for
  training or evaluating parent-to-fragment segmentation. It includes the
  parent protein and exact peptide coordinates.
- `df_wetlab_validated_sota.tsv` is the requested Vinay-facing five-column
  table: `sequence, wetlab_experiment, mic, paper_title, paper_doi`. It has
  2,449 assay rows, 596 unique canonical strings, and all 22 primary papers.
- `df_wetlab_assays_structured.tsv` is the authoritative analysis table. It
  separates assay target, MIC value/unit, chemistry, source locator, murine
  evidence, and AMP Challenge similarity.
- `df_wetlab_compound_decisions.tsv` explains, for every one of the 601
  sequence-paper records, why it is included in or excluded from segmentation.
- `df_spear_segmentation_gold.tsv` contains 270 strict labels over 217 parent
  sequences. Every peptide is an exact substring of a unique longer natural
  parent and every coordinate is verified.
- `spear_segmentation_gold.jsonl` is the grouped model-ready version of that
  strict gold set.
- `segmentation_splits/*.jsonl` contains the full 12,278-parent training corpus.
  MMseqs2 clusters are kept intact across train/validation/test/holdout.
- `df_spear_homology_splits.tsv` is the auditable split assignment.
- `df_amp_challenge_similarity.tsv` is the sequence-level challenge audit.
- `SPEAR_SOTA_release.xlsx` is a review-friendly workbook; TSV/JSONL remain the
  canonical machine-readable files.

## False-positive policy

The strict segmentation target requires all of the following: direct positive
antimicrobial evidence, canonical unmodified sequence evidence, an exact span,
a unique resolved parent sequence, and a parent longer than the peptide.
Therefore 331 sequence-paper records remain valid AMP assay evidence but are
not cleavage labels: 261 lack a resolved biological parent (mostly synthetic
designs), 35 have ambiguous parent sequences, 23 are supported only as modified
chemistry/stereochemistry, and 12 were tested as engineered constructs.

Modification is not the main reason the broad set has 596 sequences while the
strict set has 270. Allowing modified and engineered compounds while retaining
the exact-parent requirements raises the usable mapping set only to 280
sequence-paper mappings (277 unique sequences). Most of the remaining peptides
cannot define a segmentation target because no unique longer parent is known.

The five-column table intentionally preserves the requested schema, but a plain
amino-acid string cannot encode C-terminal amidation, all-D stereochemistry, or
engineered termini. Always use `chemical_form` and `tested_sequence_note` from
the structured table when compound identity matters.

## Leakage-safe split

Parent sequences were deduplicated, then clustered with MMseqs2 at 30% sequence
identity and 80% coverage (`cov-mode 1`). Entire clusters—not rows—were assigned
to a split. Any cluster containing an existing wet-lab label or a strict-gold
parent was forced into `wetlab_holdout`. The remaining clusters were assigned
deterministically (seed 0) to train/validation/test within proteome strata.

This conservative rule leaves 5,481 parents in train, 650 in validation, 673 in
test, and 5,474 in the wet-lab holdout. It is intentionally stricter than an
accession-only split because homologous parents otherwise leak across folds.

The activity tables use a separate peptide-level split because many of the 596
active compounds have no biological parent. Peptides are clustered with
MMseqs2 at 60% identity and 80% coverage; all assay rows inherit the split of
their peptide cluster. Do not randomly split the 2,449 assay rows. The same or
closely related peptide would otherwise appear in both training and evaluation.

The activity measurements intentionally remain long rather than creating one
MIC column per strain. A wide strain matrix would be mostly missing and would
discard assay context. Train a target-conditioned MIC regressor on rows marked
`regression_eligible=True`; use `sequence_task_sample_weight` when each peptide
should contribute equal total weight. Use `MIC_uM` as the primary comparable
regression task; keep `MIC_ug_per_mL` as a separate task rather than mixing raw
values from different units. The 97 positive rows without an exact MIC remain
useful as evidence but are not numeric regression labels. These tables
contain only positive compounds and therefore do not by themselves support a
binary AMP-versus-non-AMP classifier.

## AMP Challenge status

Against the official `antibacterial.fasta` at commit
`5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e` (SHA-256
`cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9`):

- 23 of 596 sequences are exact reference matches.
- 93 of 596 have maximum Levenshtein similarity greater than 0.80, including
  the exact matches.
- 573 pass the full-library no-exact-match rule.
- 503 pass the top-100 similarity rule.

These checks describe eligibility of generated challenge submissions. They do
not by themselves decide whether a published sequence may be used as training
data; challenge terms and provenance must be followed separately.

## Validation and rebuild

Run:

```bash
python scripts/validate_sota_release.py
```

The deterministic builder records the challenge reference hash and exact split
parameters in `release_manifest.json`. It requires Python, MMseqs2, and
python-Levenshtein:

```bash
python scripts/build_sota_release.py \
  --repo . \
  --challenge-fasta /path/to/amp-challenge/data/antibacterial.fasta \
  --challenge-commit 5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e
```

Source-reported exact MICs are never invented. Blank MICs mean the source
supports activity but an exact numerical value was not defensibly transcribed.
