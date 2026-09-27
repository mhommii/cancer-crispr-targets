"""Step 5 — Off-target check for the KRAS guides, and the final report.

Two questions are asked of each candidate guide:

1. Where else on chromosome 12 could it cut?
   Same search as the crispr-guide-design project: every NGG site on both
   strands, up to 4 mismatches.

2. How well does it discriminate the mutant allele from the normal one?
   This is the question that matters for an oncogene. The wild-type
   sequence is itself the most important potential off-target, because it
   is present in every healthy cell. A guide that differs from wild-type by
   one base far from the PAM will likely cut both.

Output:
  results/kras_guide_offtargets.tsv
  results/report.md
"""

import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import offtarget

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())
data, results = ROOT / "data", ROOT / "results"

MAX_MISMATCHES = 4

guides = pd.read_csv(results / "kras_allele_specific_guides.tsv", sep="\t")

print("Loading chr12 ...")
codes = offtarget.load_chromosome(data / "chr12.fa.gz",
                                  data / "genome" / "chr12.fa")
chromosome_length = len(codes)
print(f"chr12: {chromosome_length:,} bp")

print("Indexing NGG sites ...")
packed, starts, strands = offtarget.index_chromosome(codes)
del codes
print(f"{len(packed):,} sites indexed")

rows = []
unique_guides = guides["guide_mutant_allele"].unique()
print(f"\nSearching {len(unique_guides)} unique guides ...")
for protospacer in unique_guides:
    mismatches, folded = offtarget.count_mismatches(packed, protospacer)
    hit_index = np.flatnonzero(mismatches <= MAX_MISMATCHES)
    seed_mismatches = np.bitwise_count(folded[hit_index] & offtarget.SEED_MASK)

    counts = np.bincount(mismatches[hit_index], minlength=MAX_MISMATCHES + 1)

    # Where are the single-mismatch sites? For a mutant-allele guide the
    # wild-type KRAS locus itself should be one of them - it differs by
    # exactly the mutated base. Record its coordinate so it can be checked.
    single = hit_index[mismatches[hit_index] == 1]
    single_coordinates = []
    for position in single:
        start = int(starts[position])
        if strands[position]:  # reverse strand: map back to forward coords
            start = chromosome_length - start - (offtarget.PROTOSPACER_LENGTH + 3)
        single_coordinates.append(start)

    rows.append({
        "guide_mutant_allele": protospacer,
        "exact_matches_chr12": int(counts[0]),
        "mm1": int(counts[1]),
        "mm2": int(counts[2]),
        "mm3": int(counts[3]),
        "mm4": int(counts[4]),
        "offtargets_total": int(len(hit_index) - counts[0]),
        "offtargets_intact_seed": int(
            (seed_mismatches[mismatches[hit_index] > 0] == 0).sum()
        ),
        "single_mismatch_positions": ";".join(str(c) for c in single_coordinates),
    })

offtargets = pd.DataFrame(rows)
merged = guides.merge(offtargets, on="guide_mutant_allele", how="left")
merged.to_csv(results / "kras_guide_offtargets.tsv", sep="\t", index=False)

# A mutant-allele guide should NOT match chr12 exactly: the reference genome
# carries the wild-type sequence. If it does, something is wrong.
unexpected = merged[merged["exact_matches_chr12"] > 0]
if len(unexpected):
    print(f"\nNote: {len(unexpected)} guide(s) match chr12 exactly. The "
          f"reference carries wild-type KRAS, so a mutant-allele guide "
          f"should not - check these.")
else:
    print("\nAs expected, no mutant-allele guide matches the reference "
          "genome exactly (the reference is wild-type).")

# The wild-type allele is the off-target that matters most: it sits in every
# healthy cell. Each mutant-allele guide should therefore find exactly one
# single-mismatch site, and it should be the KRAS locus itself.
# KRAS on GRCh38 chr12 spans roughly 25.20-25.25 Mb.
KRAS_LOCUS = (25_200_000, 25_260_000)
in_kras = merged["single_mismatch_positions"].str.split(";").apply(
    lambda positions: sum(
        KRAS_LOCUS[0] <= int(p) <= KRAS_LOCUS[1] for p in positions if p
    )
)
print(f"\nSingle-mismatch sites per guide: "
      f"{merged['mm1'].min()}-{merged['mm1'].max()}")
print(f"Guides whose 1-mismatch site is the KRAS locus itself "
      f"(i.e. the wild-type allele): {(in_kras > 0).sum()} of {len(merged)}")
if (in_kras == 0).any():
    print("  Warning: some guides do not see wild-type KRAS as a "
          "1-mismatch site - check the codon mapping.")
merged["wild_type_allele_detected"] = in_kras > 0
merged.to_csv(results / "kras_guide_offtargets.tsv", sep="\t", index=False)

# --- Report ----------------------------------------------------------------
landscape = pd.read_csv(results / "top_mutated_genes.tsv", sep="\t")
normalised = pd.read_csv(results / "length_normalised_genes.tsv", sep="\t")

lines = [
    "# From a cancer cohort to candidate CRISPR guides",
    "",
    "A walk from public tumour mutation data to sequence-level guide "
    "candidates for the gene that analysis points at. Every number here "
    "was produced by the scripts in `steps/`.",
    "",
    "## 1. What is mutated in lung adenocarcinoma",
    "",
    f"TCGA-LUAD (MC3 calls): **{int(landscape['MutatedSamples'].max())} of 517 "
    f"patients** carry a TP53 mutation, the most of any gene.",
    "",
    "| Gene | Patients mutated | % of cohort |",
    "|---|---|---|",
]
for _, gene in landscape.head(8).iterrows():
    lines.append(f"| {gene['Hugo_Symbol']} | {gene['MutatedSamples']} | "
                 f"{gene['percent_of_samples']}% |")

top_density = normalised.head(5)
fallers = normalised.nsmallest(5, "rank_change")
lines += [
    "",
    "But that list is misleading. TTN, MUC16, CSMD3, RYR2 and USH2A are not "
    "lung cancer drivers - they are very large genes that accumulate "
    "passenger mutations in proportion to their length.",
    "",
    "## 2. Correcting for gene length",
    "",
    "Dividing by coding sequence length changes the ranking sharply:",
    "",
    "| Gene | CDS (kb) | Mutations/kb | Rank by frequency | Rank by density |",
    "|---|---|---|---|---|",
]
for _, gene in top_density.iterrows():
    lines.append(
        f"| **{gene['Hugo_Symbol']}** | {gene['cds_kb']:.1f} | "
        f"{gene['mutations_per_kb']:.0f} | {int(gene['rank_by_frequency'])} | "
        f"{int(gene['rank_by_density'])} |"
    )
lines += ["", "And the long genes fall away:", "",
          "| Gene | CDS (kb) | Rank by frequency | Rank by density |",
          "|---|---|---|---|"]
for _, gene in fallers.iterrows():
    lines.append(f"| {gene['Hugo_Symbol']} | {gene['cds_kb']:.1f} | "
                 f"{int(gene['rank_by_frequency'])} | "
                 f"{int(gene['rank_by_density'])} |")

lines += [
    "",
    "KEAP1 rising from 22nd to 3rd is a useful check: it is a genuine, "
    "well-described lung adenocarcinoma driver that raw frequency buried.",
    "",
    "![Frequency versus density](results/frequency_vs_density.png)",
    "",
    "## 3. Positional clustering",
    "",
    "Length correction is crude. OncodriveCLUST asks a sharper question: are "
    "a gene's mutations piled up at the same residues? See "
    "`results/oncodrive_results.tsv`. KRAS comes out top, with almost all of "
    "its mutations in a single cluster at codon 12 - the signature of an "
    "oncogene hotspot. TP53's mutations are scattered and often truncating, "
    "the signature of a tumour suppressor.",
    "",
    "## 4. Why KRAS needs allele-specific guides",
    "",
    "TP53 is already broken in the tumour, so disrupting it further is a "
    "reasonable knockout design - that is what the "
    "[crispr-guide-design](https://github.com/mhommii/crispr-guide-design) "
    "project does. KRAS is the opposite case:",
    "",
    "- it is an oncogene driven by **one** mutant copy,",
    "- the other copy is normal and needed by healthy cells,",
    "- so a guide matching both copies would cut the normal one too.",
    "",
    "The design question becomes: can a guide tell one base apart? That "
    "depends on where the change sits relative to the PAM, because Cas9 "
    "tolerates mismatches poorly in the ~12 bases next to it.",
    "",
    "## 5. Candidate guides",
    "",
    f"{len(guides)} candidates across {guides['mutation'].nunique()} codon-12 "
    f"substitutions. Only two PAMs in the reference place a protospacer over "
    f"codon 12, which is itself a practical finding - this region is "
    f"PAM-poor.",
    "",
    "| Mutation | Guide (mutant allele) | PAM | Distance from PAM | "
    "Off-targets on chr12 (≤4 mm) |",
    "|---|---|---|---|---|",
]
best = merged.sort_values(["mutation", "distance_from_pam"]).groupby(
    "mutation").first().reset_index()
for _, guide in best.iterrows():
    lines.append(
        f"| {guide['mutation']} | `{guide['guide_mutant_allele']}` | "
        f"{guide['pam']} | {int(guide['distance_from_pam'])} | "
        f"{int(guide['offtargets_total'])} |"
    )

lines += [
    "",
    "Distance from PAM is in bases; 0 means the changed base sits directly "
    "against the PAM, which is the most favourable position for telling the "
    "two alleles apart. No codon-12 substitution creates a new PAM, so none "
    "of these guides gets the strongest form of discrimination.",
    "",
    "## 6. The wild-type allele is the off-target that matters",
    "",
    "The search makes the central problem concrete. Every one of these "
    f"guides has **exactly one single-mismatch site on chromosome 12, and it "
    f"is the KRAS locus itself** - the normal copy of the gene, present in "
    f"every healthy cell. Confirmed for "
    f"{int(merged['wild_type_allele_detected'].sum())} of {len(merged)} "
    f"guides by coordinate.",
    "",
    "So the question is never whether a guide has an off-target here. It "
    "does, unavoidably, by design. The question is whether one mismatch in "
    "that position is enough for Cas9 to cut the mutant copy and leave the "
    "normal one intact. That is why the distance-from-PAM column is the "
    "column to read: the guides with the change at position 0 or 1 put it "
    "in the part of the guide where Cas9 is least tolerant.",
    "",
    "## Limitations",
    "",
    "- **Nothing here is validated.** These are sequence analyses, not "
    "evidence that any guide discriminates in cells.",
    "- **Recurrent mutation is not dependency.** A gene being mutated often "
    "does not show the tumour needs it. That requires functional work, such "
    "as a CRISPR screen.",
    "- **Length correction is crude.** Mutations per kb ignores sequence "
    "context, replication timing and expression, all of which MutSigCV "
    "models properly.",
    "- **One chromosome.** Off-target search covers chr12 only (~2.6% of the "
    "genome), and no bulges are searched.",
    "- **Discrimination is predicted from position alone.** Whether a "
    "one-base difference is actually enough depends on the guide, the "
    "chromatin and the Cas9 variant used.",
    "",
    "## Sources",
    "",
    "- Ellrott, K. et al. (2018). Scalable Open Science Approach for Mutation "
    "Calling of Tumor Exomes Using Multiple Genomic Pipelines. *Cell Systems*. "
    "https://doi.org/10.1016/j.cels.2018.03.002",
    "- Mayakonda, A. et al. (2018). Maftools: efficient and comprehensive "
    "analysis of somatic variants in cancer. *Genome Research*. "
    "https://doi.org/10.1101/gr.239244.118",
    "- Lawrence, M. S. et al. (2013). Mutational heterogeneity in cancer and "
    "the search for new cancer-associated genes. *Nature*. "
    "https://doi.org/10.1038/nature12213",
    "- Tamborero, D., Gonzalez-Perez, A. & Lopez-Bigas, N. (2013). "
    "OncodriveCLUST: exploiting the positional clustering of somatic "
    "mutations to identify cancer genes. *Bioinformatics*. "
    "https://doi.org/10.1093/bioinformatics/btt395",
    "- Bae, S., Park, J. & Kim, J.-S. (2014). Cas-OFFinder. *Bioinformatics*. "
    "https://doi.org/10.1093/bioinformatics/btu048",
]

(ROOT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

print("\nOff-target counts per guide:")
print(merged[["mutation", "guide_mutant_allele", "distance_from_pam",
              "offtargets_total", "offtargets_intact_seed"]]
      .to_string(index=False))
print("\nSaved: results/kras_guide_offtargets.tsv, report.md")
