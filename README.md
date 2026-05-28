# Epigenetic diploidisation and conserved subgenome asymmetry in post-WGD cyprinids


## Repository layout

Each top-level folder corresponds to a section of the Methods.

```
.
├── 01_assembly/                 Hifiasm + purge_dups + HapHiC + BUSCO
├── 02_repeat_annotation/        Earl Grey vertebrata pipeline
├── 03_gene_annotation/          BRAKER3 + PASA + InterProScan + eggNOG
├── 04_subgenome_partitioning/   SubPhaser + WGDI Ks validation
├── 05_phylogenomics/            OrthoFinder + IQ-TREE + treePL
├── 06_wgd_ks_analysis/          Paralog Ks distributions
├── 07_synteny/                  JCVI MCscan macro-/micro-synteny
├── 08_fractionation/            Sliding-window retention vs Danio
├── 09_homoeolog_expression/     HEB across tissues, TE-flanking density
├── 10_te_age_analysis/          Kimura age of subgenome-flanking TEs
├── 11_dnds_analysis/            PAML branch-model ω per homoeolog
├── 12_structural_variation/     MUMmer + SyRI inter-subgenome SVs
├── 13_homoeologous_exchanges/   Read-depth + synteny HE detection
├── 14_wgbs_methylation/         Trim Galore + Bismark + metaplots
├── 15_epigenetic_diploidisation/ Silencing-direction classification + Fisher
├── 16_3d_genome/                BWA + HiCExplorer + CALDER (A/B + TADs)
└── utils/                       Common helpers (GFF→BED, pair table, plots)
```

### Mapping of Methods sections to scripts

| Methods subsection | Script(s) |
|---|---|
| Genome assembly | `01_assembly/assembly_pipeline.sh` |
| Repeat annotation | `02_repeat_annotation/earl_grey_repeat_annotation.sh` |
| Gene structure & function annotation | `03_gene_annotation/gene_annotation_pipeline.sh`, `rnaseq_quantification.sh`, `aggregate_tpm_matrix.py` |
| Subgenome partitioning | `04_subgenome_partitioning/subphaser_run.sh`, `ks_based_subgenome_assignment.sh` |
| Phylogenomics & divergence dating | `05_phylogenomics/phylogenomics_pipeline.sh`, `extract_cds_for_og.py` |
| Ks paralog-pair analysis | `06_wgd_ks_analysis/ks_distribution_analysis.py` |
| Synteny | `07_synteny/mcscan_jcvi_synteny.sh` |
| Fractionation | `08_fractionation/fractionation_analysis.py` |
| Homoeolog expression bias | `09_homoeolog_expression/homoeolog_expression_bias.py`, `te_flanking_density.py` |
| TE-age comparison | `10_te_age_analysis/te_age_comparison.py` |
| dN/dS / ω analysis | `11_dnds_analysis/branch_specific_omega.py` |
| Structural variation | `12_structural_variation/syri_subgenome_sv.sh` |
| Homoeologous exchanges | `13_homoeologous_exchanges/he_detection_read_depth.py` |
| WGBS bisulfite analysis | `14_wgbs_methylation/wgbs_bismark_pipeline.sh`, `methylation_metaplot.py`, `gene_promoter_methylation.py` |
| Epigenetic diploidisation | `15_epigenetic_diploidisation/epigenetic_diploidisation.py`, `run_cross_species_fisher.py` |
| 3D genome (Hi-C compartments / TADs) | `16_3d_genome/hic_3d_genome_pipeline.sh`, `calder_compartments.R`, `ab_compartment_gene_analysis.py` |

---

## Dependencies

The pipelines combine command-line bioinformatics tools and Python/R
analysis code. Versions match the manuscript.

**Command-line tools**

- Hifiasm v0.25.0, purge_dups v1.2.5, HapHiC v1.0.6, BUSCO (Actinopterygii_odb10)
- Earl Grey v6.1.1 (RepeatModeler, LTRfinder, RepeatCraft)
- BRAKER3 v3.0.8 (AUGUSTUS, GeneMark-ETP, TSEBRA), PASA, InterProScan, eggNOG-mapper
- HISAT2, StringTie, samtools, BWA-MEM, samblaster, bedtools v2.27.1
- SubPhaser, WGDI v0.74, DIAMOND, OrthoFinder v3.0.1, MAFFT v7.453, pal2nal,
  AMAS, IQ-TREE v2.2.6, treePL, PAML v4.9 (codeml)
- JCVI MCscan, MUMmer v4.0, SyRI v1.7.1
- Trim Galore v0.6.10, Bismark v0.24.2, HiCExplorer v3.7.6, CALDER

**Python**

See `environment.yml` / `requirements.txt`. Core packages: numpy, scipy,
pandas, scikit-learn, biopython, pysam, matplotlib, seaborn, statsmodels.

**R**

R ≥ 4.2, with CALDER and its declared dependencies.

---

## License

MIT — see `LICENSE`.
