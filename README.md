# Project SPEAR: <img src="data/logo.jpeg" height="100">

## Segmentation of Peptides for Encrypted Antimicrobial Revelation 

Goal: To train a model to automatically cleave peptides strategically to reveal antimicrobial peptides encrypted inside them (cryptomics mining for antimicrobial peptides).

# AMP encryptome dataset

Curated 8 September 2026 from final journal publications associated with [César de la Fuente-Nunez's lab](https://delafuentelab.engineering.upenn.edu/publications-2/). 
## Disclaimers:
- This is a source-audited release, not a claim that every paper's reported candidate pool or historical parent proteome is publicly recoverable. No new antimicrobial predictions were generated.
- The main co-author was GPT 6 Astra 

## Data Contents and counts

1. `df_amp_encryptome.tsv`: **239,610 rows**, representing **98,790 peptide–paper pairs** and **98,771 distinct peptide sequences**, across 10 core papers. Exactly the requested seven columns.
2. `df_cesar.tsv`: **32 final journal papers**, with the requested four columns. Includes 10 core mining papers, 2 adjacent smORF papers, 10 synthetic-design papers and 10 reviews/perspectives. No preprints.
3. `df_amp_provenance.tsv`: one row for every main-table row, including source locations, reported and recovered accessions, mapping evidence, sequence hashes, positions and annotation conflicts.
4. `df_amp_coverage.tsv`: paper-by-paper extraction coverage and limitations. `df_cesar_audit.tsv` preserves the 37-paper discovery audit and explicitly excludes its 5 preprints.
5. `df_amp_smorf.tsv`: **863,821** standalone smORF records, comprising all **863,498 AMPSphere** sequences and **323 human-microbiome SEPs**. These do not establish a longer host protein.
6. `df_amp_engineered.tsv`: **14 constructs** from the core mining papers, including added-Pro and D/retro-inverso variants. This is not a complete catalogue of de novo or engineered peptides from the 10 synthetic-design context papers.
7. `df_amp_nonfragment_candidates.tsv`: **22 candidates** whose only resolved source equals the whole deposited sequence. They remain available without being represented as verified fragments of a longer protein.
8. `source_manifest.json`, `dataset_summary.json`, `validation_report.json` and `reproducibility/`: provenance, machine-readable counts, validation and reconstruction code.
9. `df_wetlab_validated.tsv`: **1,551 complete measured peptide–assay rows** for **280 distinct sequences** from the journal-only César catalogue. See [`data/WETLAB_VALIDATED.md`](data/WETLAB_VALIDATED.md) for extraction coverage and the separate AMP Challenge similarity audit.

## Main-table interpretation

A row represents `(encrypted_peptide, distinct source_peptide sequence, doi_paper)`. Identical full parent sequences are collapsed even if several accessions share them. Different matching parent sequences produce different rows; they must not be interpreted as independent discoveries or independent experimental observations. The provenance `main_data_row` counts data rows starting at 1, excluding the header. `record_id` is a stable SHA-256-derived identifier for that three-field key.

Every nonempty `source_peptide` is a literal, case-sensitive contiguous match containing the encrypted sequence and is longer than it. **274 candidate–paper records have an empty parent sequence**. Unresolved candidates are retained provisionally to avoid discarding released sequences. The presence of an exact reference match does not, by itself, establish which historical protein the authors actually used. The provenance distinguishes reported accessions, author-released parents, domain-reference matches, broader reference searches and explicitly documented alternative mappings.

`source_peptide` contains the entire available reference record. Some records are partial proteins or contain ambiguous residues; these are flagged. Sequence versions and available isoforms are preserved in the parent headers. No parent protein was fabricated by extending a peptide or replacing residues. `source_proteome` is the study's domain, not independent proof of the taxonomy or extinction status of every matched parent.

`model_used` refers specifically to MIC regression. `apex_1.0` normalizes the original APEX used in the 2024 extinctome and 2025 Venomics papers. `apex_1.1` denotes APEX-Pathogen in the archaeasin and prionin papers. Classifiers such as panCleave, Macrel and AmPEP, and physicochemical scores, are not mislabeled as MIC regressors; these rows use `not_applicable`. Unverified models remain explicit in the broader paper catalogue.

`mic_threshold` preserves the paper's selection statistic, inequality and units. Here `uM` means micromolar, equal to micromoles per litre. It is a **prediction-selection rule**, not a measured per-peptide MIC or a universal definition of activity. The core APEX cutoffs are median MIC <=80 uM for the extinctome, median MIC <=32 uM for Venomics, mean MIC <=100 uM for archaea, and median MIC <=64 uM for prionins, across the reported 11-strain panels. Separately scored controls do not inherit an APEX threshold.

`wet_lab` describes sequence-specific experimental evidence. `not_reported` means no assay could be confidently assigned to that exact native sequence in this extraction; it does not mean the peptide is inactive or was never tested. In-vitro testing is distinct from mouse testing, and testing does not imply success. Assays on an added-Pro or D analogue are assigned to the engineered table. Terminal chemical modifications and experimental formulation are not exhaustively encoded in the seven-column amino-acid dataset; use the original methods for synthesis-level identity.

## Coverage and material limitations

| Study | Recovered peptide–paper records | Main candidates | Parent unresolved | In-vitro labels | Mouse labels |
|---|---:|---:|---:|---:|---:|
| human2021 | 43,053 | 43,053 | 7 | 56 | 2 |
| extinctome | 37,228 | 37,228 | 78 | 118 | 5 |
| venomics | 4,664 | 4,645 | 40 | 55 | 3 |
| archaea | 12,623 | 12,620 | 0 | 80 | 3 |
| prion | 1,179 | 1,179 | 147 | 75 | 0 |
| nonimmune | 39 | 39 | 0 | 39 | 16 |
| pancleave | 14 | 14 | 0 | 14 | 6 |
| mmp19 | 3 | 3 | 0 | 0 | 0 |
| fibrinogen | 6 | 6 | 0 | 6 | 0 |
| plasma | 3 | 3 | 2 | 0 | 0 |

Counts describe released sequence records, not a uniform selection stage. The main-table counts include both initial released prediction/ranking pools and published experimental candidates; `candidate_stages` permits filtering them.

1. **Human proteome, 2021:** the released CSV contains 43,000 ranked peptides (1,000 per length from 8 to 50), not an explicitly marked set of the final 2,603 mentioned in the abstract. The dataset retains that pool and 53 additional tested sequences. All 56 tested sequences are represented, consistent with the 2022 correction (10.1038/s41551-022-00967-2). The exact 2020 human reference FASTA bundled with the paper's code was used.
2. **Extinctome, 2024:** all 37,176 released APEX predictions, three additional APEX-tested sequences, and 49 separately scored controls are retained. Several common-name/accession/source annotations are inconsistent. Mouse labels use Figure 6's accession-fragment identifiers, not the conflicting common-name column. Current reference matching recovers some sequences from the extant taxon *Equus quagga boehmi*; those rows retain the paper-level `extinctome` label with a provenance warning. They do not establish an extinct origin.
3. **Venomics, 2025:** Supplementary Data 1 contains only the 4,618-sequence UniProt initial pool. The initial ConoServer, ArachnoServer and ISOB pools underlying the paper's 7,379 count were not recovered. The released 386-member diverse pool and all 58 tested instances are included, giving 4,664 unique recovered sequences before separating 19 whole-source matches. The 58 tested instances correspond to 55 unique sequences. Only UniprotKB-7, Conoserver-14 and Arachnoserver-5 receive the mouse skin-infection label.
4. **Archaeasins:** all 12,623 released predictions and 80 tested sequences are captured. The caption's 12,632 is inconsistent with the actual table and main text. Three complete-source matches are separated.
5. **Prionins:** all 1,179 released predictions and 75 tested sequences are captured. The tested-sequence table lacks the prionin number/name key, so the mouse evidence for prionin-7 and prionin-38 cannot safely be assigned to sequence rows. There are 147 unresolved parent sequences.
6. **panCleave:** the 14 active hits have recoverable published sequences and parents. The larger predicted/selected pool was not recovered as a sequence catalogue. The companion repository's wrapped fragment sequence column is misaligned across several rows; the final journal Supplementary Tables S3/S5 resolve those errors.
7. **Non-immune proteins:** all 39 sequences are captured. Table S1's parent assignments for Eyesin-3 and Eyesin-5 conflict with exact matches to the other BRINP family member; both the reported IDs and inferred alternatives are preserved. CD8B isoforms are retained.
8. **MMP19/plasma:** native candidates do not inherit added-Pro or retro-inverso construct assays. The full parent of the ApoB Pro900 variant was not reconstructed; only the canonical Ala parent is filled.
9. **Standalone smORFs:** all released sequences are included, but peptide-specific experimental linkage is incomplete. AMPSphere has 31 confidently linked tested records out of 100 reported tested, including all 10 mouse lead accessions. SEPs has 55 linked active records out of 78 reported tested, including prevotellin-2. Do not interpret the other `not_reported` rows as known negatives.

## Use and reproducibility

TSVs are UTF-8, tab-delimited and quote fields when needed. The longest parent is 34,350 residues, beyond Excel's per-cell text limit; use a dataframe or text reader to preserve such sequences. Read them as strings, preserving empty fields:

```python
import pandas as pd
amp = pd.read_csv('df_amp_encryptome.tsv', sep='\t', dtype=str, keep_default_na=False)
provenance = pd.read_csv('df_amp_provenance.tsv', sep='\t', dtype=str, keep_default_na=False)
```

Count peptide discoveries by sequence/DOI, not by expanded parent rows. For a stricter analysis, use the provenance to restrict to reported parent accessions or author-released parent sequences, exclude unresolved/partial parents and separate initial pools from downstream selection stages. Do not combine MIC thresholds from different models/panels as if they were measured potencies.

The bundle contains frozen candidate annotations, parent sequences/matches and the published AMPSphere FASTA. `python3 reproducibility/rebuild.py --output rebuilt` reconstructs the tables from those frozen inputs using the standard Python library. The output ordering and values are deterministic. Source-extraction scripts and the source manifest document how the curation was assembled, including manual transcriptions from final tables. Re-fetching historical source databases is not guaranteed to reproduce their current contents; source hashes and frozen records preserve this release. No model weights or inference are required.

Primary article and supplementary links are in the paper audit, coverage table and source manifest. Cite the original papers when reusing their sequence data. This compilation makes no claim to relicense the source datasets.


