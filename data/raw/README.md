# Data files

Large genomic files are gitignored and must be downloaded before running any scripts.

## 1. VCF (chr22, ~425 MB)

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

Source: 1000 Genomes 2022 high-coverage (30x) release, GRCh38.
Contains 1,066,557 variants across 3,202 samples (2,504 unrelated + 698 trio relatives).

## 2. VCF index

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
```

## 3. Population panel (already tracked in repo)

`integrated_call_samples_v3.20130502.ALL.panel` — maps each of the 2,504 unrelated Phase 3 samples to population (e.g. YRI) and superpopulation (e.g. AFR).

Download URL if needed:
```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

## 4. IGSR population metadata (already tracked in repo)

`igsr_populations.tsv` — official sampling coordinates and population descriptions from the International Genome Sample Resource. Used by `scripts/07_africa_map.R`.

Download from: https://www.internationalgenome.org/data-portal/population (click "Download the list").
