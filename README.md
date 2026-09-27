# Cancer Genomics → CRISPR Targets

Going from a public tumour mutation dataset to sequence-level CRISPR guide candidates: which gene does the data point at, and what would a guide against it have to look like?

> **How this was made:** a guided learning project built and run with Claude Code (an AI assistant), which wrote the code and the explanations and executed each step. It is a learning exercise, not independent research. The *My notes* sections are left for me to fill in as I work through the code.

**Full write-up: [report.md](report.md)**

## The question

TCGA-LUAD (lung adenocarcinoma, 517 patients) — start from real somatic mutation calls and work towards a guide design, without skipping the step where you check whether the obvious answer is actually right.

| Step | Script | What it does |
|---|---|---|
| 1 | `steps/01_mutation_landscape.R` | Load TCGA-LUAD MC3 calls; oncoplot and top mutated genes |
| 2 | `steps/02_normalise_and_choose_target.py` | Correct for gene length; fetch CDS lengths from NCBI |
| 3 | `steps/03_positional_clustering.R` | OncodriveCLUST: are mutations clustered or scattered? |
| 4 | `steps/04_design_allele_specific_guides.py` | Design guides against the KRAS codon-12 hotspot |
| 5 | `steps/05_offtarget_and_report.py` | Off-target search on chr12; build the report |

## What the data showed

**The obvious answer is wrong.** Ranking by how many patients carry a mutation puts TP53 first — but then TTN, MUC16, CSMD3, RYR2 and USH2A. Those are not lung cancer drivers. They are very large genes that collect passenger mutations in proportion to their length ([Lawrence et al. 2013](https://doi.org/10.1038/nature12213)).

Dividing by coding sequence length reorders things:

| Gene | CDS (kb) | Rank by frequency | Rank by density |
|---|---|---|---|
| **KRAS** | 0.6 | 9 | **1** |
| **TP53** | 1.2 | 1 | **2** |
| **KEAP1** | 1.9 | 22 | **3** |
| TTN | 108.0 | 2 | 25 |
| MUC16 | 46.5 | 3 | 23 |

TTN encodes titin at 35,992 amino acids — the largest human protein. It falls from 2nd to 25th. KEAP1, a genuine and well-described lung adenocarcinoma driver, rises from 22nd to 3rd, which is a useful sign the correction is doing something real rather than just shuffling.

![Frequency versus density](results/frequency_vs_density.png)

**Clustering separates the two kinds of driver.** OncodriveCLUST puts KRAS top (FDR 0.034) with 152 of 161 mutations in a single cluster. The lollipop plots show why that matters:

![KRAS lollipop](results/kras_lollipop.png)

KRAS: one enormous spike at codon 12, nearly all missense — an oncogene hotspot. TP53 (`results/tp53_lollipop.png`): mutations scattered along the protein, many of them truncating — a tumour suppressor being broken.

## Why KRAS needs a different guide design

That difference changes the design problem completely.

TP53 is already broken in the tumour, so disrupting it further is a sensible knockout — that is what the [crispr-guide-design](https://github.com/mhommii/crispr-guide-design) project does. KRAS is the opposite: it is driven by **one** mutant copy, the other copy is normal, and healthy cells need it. A guide matching both copies cuts the normal one too.

So the guide has to tell a single base apart. Whether that is possible depends on where the change sits relative to the PAM, because Cas9 tolerates mismatches poorly in the ~12 bases next to it.

| Mutation | Guide (mutant allele) | PAM | Distance from PAM | Off-targets on chr12 (≤4 mm) |
|---|---|---|---|---|
| G12A | `CTTGTGGTAGTTGGAGCTGC` | TGG | 0 | 8 |
| G12C | `CTTGTGGTAGTTGGAGCTTG` | TGG | 1 | 10 |
| G12D | `CTTGTGGTAGTTGGAGCTGA` | TGG | 0 | 11 |
| G12R | `CTTGTGGTAGTTGGAGCTCG` | TGG | 1 | 3 |
| G12S | `CTTGTGGTAGTTGGAGCTAG` | TGG | 1 | 7 |
| G12V | `CTTGTGGTAGTTGGAGCTGT` | TGG | 0 | 9 |

Two findings worth stating plainly:

- **Only two PAMs place a protospacer over codon 12.** The region is PAM-poor, which is a practical constraint, not an oversight.
- **No codon-12 substitution creates a new PAM.** That would have been the strongest form of discrimination — Cas9 simply cannot cut an allele with no PAM — and none of these mutations provides it.

### The wild-type allele is the off-target

Every one of these guides has **exactly one single-mismatch site on chromosome 12, and it is the KRAS locus itself** — the normal copy, in every healthy cell. Verified by coordinate for all 12 guides.

That is not a flaw to be filtered out. It is the design problem. The question is whether one mismatch, in that position, is enough for Cas9 to cut the mutant and spare the normal copy — which is why the distance-from-PAM column is the one to read.

## Checks built in

Each step stops rather than producing quietly wrong numbers:

- **Step 2** requires TTN to have the longest CDS in the table. It catches a real failure: taking the first RefSeq hit per gene returned 444 bp for TTN instead of ~108,000, which inverted the whole length argument. The fix takes the longest annotated CDS.
- **Step 4** requires KRAS codon 12 to read `GGT` (glycine) after mapping transcript coordinates back onto the minus-strand genome. A mapping error here would design guides against the wrong codon.
- **Step 5** checks that no mutant-allele guide matches the reference exactly (the reference is wild-type), and that each guide's single-mismatch site really is the KRAS locus.

## Running it

```powershell
# R side (maftools, TCGAmutations)
Rscript steps\01_mutation_landscape.R
Rscript steps\03_positional_clustering.R     # several minutes

# Python side
py -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
# put your own email in config.toml (NCBI requires one)
.venv\Scripts\python.exe steps\02_normalise_and_choose_target.py
.venv\Scripts\python.exe steps\04_design_allele_specific_guides.py
# step 5 needs chr12 (~38 MB) in data/chr12.fa.gz, from Ensembl:
#   https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/
.venv\Scripts\python.exe steps\05_offtarget_and_report.py
```

R packages: `BiocManager::install("maftools")` and `BiocManager::install("PoisonAlien/TCGAmutations")`.

## My notes

*(to be written by me)*

### On the gene-length artefact

### On oncogenes versus tumour suppressors

### On allele-specific targeting

## Limitations

- **Nothing here is validated.** These are sequence analyses. No guide has been tested, and none is a recommendation for any use.
- **Recurrent mutation is not dependency.** A gene being mutated often does not show the tumour needs it. Showing that requires functional work, such as a CRISPR screen.
- **Length correction is crude.** Mutations per kb ignores sequence context, replication timing and expression — all of which MutSigCV models properly. It is used here to demonstrate the effect, not as a substitute.
- **One chromosome.** Off-target search covers chr12 only (~2.6% of the genome), and searches mismatches but not bulges.
- **Discrimination is predicted from position alone.** Whether one mismatch is actually enough depends on the guide, the chromatin context and the Cas9 variant.
- **Reference sources differ between projects.** chr12 here comes from Ensembl (UCSC was unreachable); the companion project used UCSC hg38. Both are GRCh38 primary assembly, so coordinates agree.

## Sources

- Ellrott, K. et al. (2018). Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. *Cell Systems*. https://doi.org/10.1016/j.cels.2018.03.002
- Mayakonda, A. et al. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*. https://doi.org/10.1101/gr.239244.118
- Lawrence, M. S. et al. (2013). Mutational heterogeneity in cancer and the search for new cancer-associated genes. *Nature*. https://doi.org/10.1038/nature12213
- Tamborero, D., Gonzalez-Perez, A. & Lopez-Bigas, N. (2013). OncodriveCLUST: exploiting the positional clustering of somatic mutations to identify cancer genes. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btt395
- Bae, S., Park, J. & Kim, J.-S. (2014). Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btu048
- Cock, P. J. A. et al. (2009). Biopython. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btp163
- Sequence data: NCBI Gene and RefSeq (https://www.ncbi.nlm.nih.gov/gene/); chromosome 12 from Ensembl GRCh38 (https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/)
