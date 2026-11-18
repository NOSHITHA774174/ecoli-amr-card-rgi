# E. coli antibiotic-resistance genes: reference vs clinical genomes (CARD/RGI)

## Question
Do clinical (patient) E. coli isolates carry acquired resistance genes that common reference strains lack, and which drug classes do those genes cover?

## Data
Complete genomes from NCBI RefSeq.

| Group | Strain | Accession |
|---|---|---|
| Reference | K-12 MG1655 | GCF_000005845.2 |
| Reference | EDL933 | GCF_000006665.1 |
| Reference | CFT073 | GCF_000007445.1 |
| Reference | Sakai | GCF_000008865.2 |
| Reference | UTI89 | GCF_000013265.1 |
| Clinical | 755 (urine, Taiwan) | GCF_048403405.1 |
| Clinical | 405 (urine, Houston) | GCF_054047845.1 |
| Clinical | AH19B041 (blood, Melbourne) | GCF_051619495.1 |
| Clinical | p11B (blood, Houston, ST131) | GCF_024918915.1 |

The four clinical genomes were picked by hand from NCBI complete E. coli genomes with clinical source information. Genome files and the CARD database are not included here; download them with the accessions above.

## Methods
- NCBI `datasets` CLI 18.41.1 to download genomes
- RGI 6.0.8 with CARD 4.0.2: `rgi main --input_sequence <genome>.fna --input_type contig --local --clean -n 4`
- Only Perfect and Strict hits were kept
- R 4.6.2 (tidyverse, pheatmap) for counting, comparison and plots: see `analysis.R`

## Results
Resistance genes detected per genome: references 54-55; AH19B041 57; 405 63; p11B 75; 755 88.

33 genes were found in clinical isolates but in none of the five references (`results/clinical_only_genes.csv`; which isolate carries which: `results/gene_by_isolate.csv`).

Genes per drug class (aminoglycoside / penicillin beta-lactam / tetracycline):
- references: 6-7 / 19-20 / 16-17
- 405: 7 / 20 / 19
- AH19B041: 7 / 20 / 17
- p11B: 12 / 24 / 18
- 755: 15 / 23 / 19

Main observations:
- Clinical-only genes include acquired beta-lactamases (TEM-1, OXA-1, CTX-M-14, CTX-M-15, CMY-181), the carbapenemase NDM-5, aminoglycoside-modifying enzymes (e.g. APH(3'')-Ib, APH(6)-Id, aadA2, aadA5), tet(A), tet(B), sul1, sul2, dfrA12 and dfrA17.
- NDM-5 was found in both urine isolates (755, 405) and in neither blood isolate.
- 755 carries 27 of the 33 genes; p11B 16; 405 9; AH19B041 3.
- sul1, qacEdelta1, aadA5 and dfrA17 occur together in 755 and p11B, which is consistent with a class 1 integron (not verified here).
- Gene counts alone can mislead: 405 has few extra genes but one is a carbapenemase.

Plots: `plots/genes_per_genome_v2.png` and `plots/variable_genes_heatmap.png`.

## Caveats
- "Absent" means no Perfect/Strict RGI hit, not proof the gene is missing from the genome.
- A resistance gene does not prove the isolate is resistant in a lab test; no susceptibility data were used.
- Four hand-picked isolates is a small sample. This analysis is descriptive, with no statistics.
- gyrA, parC, UhpT, PtsI and cyaA hits are mutation models in core genes, not acquired genes.
- RGI prints a non-fatal "Requested rname ... does not exist" message on some genomes; output files are still produced.
- Gene location (plasmid vs chromosome) was not analysed.

## Reproduce
1. Download the genomes with the accessions above and place them in `clinical/` and `genomes/`.
2. Load CARD 4.0.2 locally and run RGI as shown in Methods, writing outputs to `results/`.
3. Run `analysis.R` in R.
