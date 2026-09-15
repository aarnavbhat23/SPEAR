# SPEAR SOTA dataset report

## Executive conclusion

The release is ready to use as a defensible SPEAR baseline. It contains a broad
model-development corpus, a chemistry-aware assay table, a strict wet-lab gold
set, explicit exclusion decisions, and homology-safe splits. The key design
choice is to keep **AMP activity evidence** separate from **segmentation labels**.
A compound can be experimentally active without proving that its written
sequence is a native fragment of a longer parent protein.

For simple handoff, use `SPEAR_596_wetlab_validated.tsv` for the complete
wet-lab-positive sequence set and `SPEAR_270_strict_segmentation.tsv` for exact
parent-to-fragment labels.

For activity modeling, use `SPEAR_596_activity_peptides.tsv` as the peptide and
split index and join it to `SPEAR_2449_activity_measurements.tsv` by
`peptide_id`. The measurement table contains 1,957 exact µM labels, 395 exact
µg/mL labels, and 97 positive observations without a numerical MIC. Unit tasks
are kept separate, and 60%-identity peptide clusters never cross splits.

## Release statistics

| Layer | Count | Intended use |
|---|---:|---|
| Positive peptide-assay rows | 2,449 | Activity evidence and MIC analyses |
| Unique canonical sequence strings | 596 | Sequence-level AMP analysis |
| Sequence-paper compound records | 601 | Inclusion/exclusion audit |
| Primary journal papers | 22 | Literature coverage |
| Rows with exact MIC | 2,352 | Quantitative potency analyses |
| Rows with source-supported activity but no exact MIC | 97 | Qualitative positive evidence |
| Strict segmentation labels | 270 | Final wet-lab gold evaluation |
| Strict-gold parent sequences | 217 | Model-level holdout |
| All deduplicated parent sequences | 12,278 | Segmentation development corpus |
| Parent homology clusters | 3,811 | Leakage-safe partitioning |

The apparent asymmetry—2,449 assay rows for 596 sequences—is expected. One
peptide may be tested against several organisms, strains, conditions, or in
multiple models. Each peptide-assay measurement is one row; it is not a new
peptide.

## What changed from the first five-column release

1. The ten AMPSphere mouse-candidate MICs were resolved against *A. baumannii*
   and added as exact values: 1, 2, 8, 16, or 64 µM depending on candidate.
2. The 36 CRPS derivatives were manually confirmed in the supplementary
   heatmaps as active against at least one tested strain. Their MIC field stays
   blank because no exact value is invented from a qualitative inclusion call.
3. Chemical form is now explicit. C-terminal amidation, all-D forms, and
   engineered constructs are preserved as assay evidence but cannot silently
   act as unmodified sequence labels.
4. Every sequence-paper record has a segmentation inclusion decision and human-
   readable reason.
5. Parent proteins are clustered before splitting, which prevents homologous
   parents from crossing model partitions.
6. AMP Challenge exact-match and Levenshtein checks cover all 596 sequences.

## Strict segmentation criteria

A label enters `df_spear_segmentation_gold.tsv` only when all five conditions
hold:

1. Direct positive antimicrobial evidence is linked to the sequence.
2. Evidence supports the canonical, unmodified peptide rather than only an
   amidated, all-D, retro-inverso, or added-residue construct.
3. The peptide occurs as an exact contiguous substring of a recovered parent.
4. The parent sequence is unique for the curated mapping.
5. The parent is longer than the peptide, so a real boundary-prediction target
   exists.

The 331 excluded sequence-paper records remain in the assay dataset:

| Reason | Records |
|---|---:|
| No resolved biological parent, usually synthetic design | 261 |
| Ambiguous parent sequence | 35 |
| Activity supported only for modified chemistry/stereochemistry | 23 |
| Evidence is for an engineered construct | 12 |

This exclusion is not a claim that those compounds are inactive. It means they
do not provide a clean native cleavage label for a segmentation model.

## Recommended model protocol

1. Train only on `segmentation_splits/train.jsonl`.
2. Select checkpoints and thresholds on `val.jsonl`.
3. Use `test.jsonl` for the mining-corpus test result.
4. Treat `wetlab_holdout.jsonl` as untouched until model and threshold choices
   are frozen.
5. Report a second, most stringent result on
   `spear_segmentation_gold.jsonl`. Never merge these 217 strict-gold parents
   back into training for the headline evaluation.
6. Measure exact-span precision/recall/F1, boundary tolerance F1 (for example
   ±1 and ±3 residues), peptide intersection-over-union, and parent-level
   precision-recall. Residue accuracy alone will be misleading because most
   residues are negative.
7. Report results both with and without the 93 challenge-similarity-restricted
   sequences when comparing to the AMP Challenge.

## Leakage control

The original accession-grouped split prevents the same accession from appearing
twice, but it does not prevent related proteins from appearing on opposite sides
of the evaluation. This release deduplicates exact parent sequences and clusters
them with MMseqs2 at 30% identity and 80% coverage. All members of a cluster
share a split. Clusters containing any wet-lab-labelled or strict-gold parent are
forced into the holdout.

The resulting source counts are 5,481 train, 650 validation, 673 test, and 5,474
wet-lab holdout. The large holdout is a consequence of conservative homology
blocking in repetitive and closely related protein families, not a row-counting
error.

## AMP Challenge result

Against the pinned official reference:

- 23 sequences are exact matches and therefore fail the no-exact-match rule for
  the full generated library.
- 93 sequences have maximum Levenshtein similarity greater than 0.80, including
  those exact matches, and therefore fail the stricter top-100 similarity rule.
- 573 sequences pass the exact-match screen.
- 503 sequences pass the top-100 similarity screen.

This does **not** mean the dataset is banned as training data. It means those
specific published sequences should not be submitted as generated challenge
entries under the corresponding screens. The challenge's current rules and
data-use terms remain authoritative.

## Remaining limitations

- A canonical one-letter sequence does not fully specify an experimental
  compound. The structured table must be used for chemical identity.
- Ninety-seven positive rows lack a defensible exact MIC; they are retained
  with blank values rather than guessed numbers.
- The segmentation development corpus includes literature-mined candidates, not
  239,610 independently wet-lab-confirmed cleavage events. The strict 270-label
  gold set is the honest experimental benchmark.
- Public supplementary materials do not expose every sequence identity reported
  as active in every study. Unresolvable identities are omitted rather than
  hallucinated.
- A domain expert should spot-check a stratified sample before a formal paper or
  leaderboard submission.

## Suggested message to Vinay

> I rebuilt the release as two separate layers so AMP activity cannot be
> confused with a cleavage label. The full evidence table has 2,449 positive
> peptide-assay rows for 596 unique canonical sequences across all 22 primary
> papers, with the exact five columns you requested. I also added a structured
> chemistry-aware version and an explicit decision table. For SPEAR itself,
> 270 peptides pass a strict gold standard: direct antimicrobial evidence, an
> unmodified canonical sequence, an exact unique span, and a longer biological
> parent. Synthetic, modified, engineered, whole-source, or parent-ambiguous
> records stay in the AMP evidence layer but are excluded from segmentation.
> The 12,278 parent sequences are split by MMseqs2 homology cluster, with every
> wet-lab cluster held out, so related proteins cannot leak across train/test.
> I also audited all 596 sequences against the pinned official AMP Challenge
> reference: 23 are exact matches and 93 exceed 0.80 similarity, so the release
> includes explicit eligibility flags rather than claiming every peptide is
> challenge-submittable.
