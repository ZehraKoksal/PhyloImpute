# Changelog

All notable changes to this project will be documented here.

## [v1.3] - 2025-12-09
### Added
- ability to input also multisample vcf files
- uncertainty warning for user: Asterisk * before predicted haplogroup (if Confidence value < Penalty value 1)
- automatically generated log file
- possibility to output imputed file in multisample vcf format (0:999 indicates imputed Reference allele, 1:999 indicates imputed Alternative 1 allele, 2:999 etc.)
- In "_haplogroups.csv" file add "downstream_ancestral" and "downstream_unknown" columns that contains the downstream haplogroup names with the specific SNPs that are ancestral or not known to allow user to verify if the predicted haplogroup is correct or simply lacking resolution due to sample quality
### Fixed
- bug for SNPs without any derived alleles
- replace asterisk (*) in NAMQY haplogroup names by " to avoid confusion in nomenclature meaning for asterisks
---


## [v1.2] – 2025-07-24
### Added
- Functionality for reference genome T2T-aligned data
- Generate allele frequency maps
- Option to receive phyloimputed nucleotides (A, C, G, T, N) instead of default ancestral states (A, D, X).
### Fixed
- Bug fix for ISOGG 2020 tree, when using vcf file. The database contains ancestral/derived allele switch-ups that leads to wrong haplogroup predictions for haplogroups afar from reference haplogroup R, e.g., A and B. >90 switch up positions were found and changed (see "switched_alleles_isogg.csv")
- Functioning for vcf files that contain non-polymorphic variants, i.e. alleles corresponding to reference allele.
- erroneous entries of "nan" SNPs in ISOGG tree were removed
- recurrent mutation-aware imputation
---

## [v1.1] – 2024-10-23
### Added
- New phylogenetic tree added for Y chromosome: ISOGG 2019/2020 
### Fixed
- Bug fix for errors caused by big phylogenetic trees: removing missing values using additionally np.notnull instead of only np.nan
- Bug fix for Hg Penalty value 2 calculation: exclude variants that are derived in the data, but not present in the tree from categories "derived variants in parallel branches"
- Bug fix for Hg Penalty value 2 calculation: reduce missing values from counting total SNPs in tree (denominator of Hg Penalty value 2)

---

## [v1.0] – 2024-06-12
### Added
- Initial release with core functionality.
- Input file formats:
-- csv: character-by-taxon data matrix with ancestral allelic state ("A"), derived allelic state ("D"), or missing information ("X")
-- vcf: traditional vcf files

- Phylogenetic tree:
-- Minimal Y tree
-- custom trees

- Output files:
-- phyloimputed.csv: Imputed missing data
-- haplogroups.csv: haplogroup prediction, confidence values, penalty value1, penalty value2
-- conflicting_SNPs.csv: list of variants causing contradictions (counted in penalty values 1 and 2)
