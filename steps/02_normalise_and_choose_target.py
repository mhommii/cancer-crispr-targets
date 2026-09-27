"""Step 2 — Separate likely drivers from long genes, then choose a target.

The problem with the step 1 ranking
-----------------------------------
Step 1 ranked genes by how many patients carry a mutation. The top of that
list includes TTN, MUC16, CSMD3, RYR2 and USH2A. These are not lung cancer
drivers. They are simply enormous genes: TTN encodes titin, the largest
human protein at over 30,000 amino acids. A gene with ten times the coding
sequence collects roughly ten times the passenger mutations, so raw
frequency puts big genes at the top regardless of whether they matter.

This is a well-described problem. Lawrence et al. (2013) showed that
naive frequency analysis produces long, highly-expressed-tissue-specific
false positives, and that correcting for gene length and background
mutation rate removes most of them. Proper tools for this (MutSigCV,
OncodriveCLUST, dNdScv) model the background rate carefully.

What this step does is the simplest version of that correction: divide the
mutation count by the length of the coding sequence, giving mutations per
kilobase. That is much cruder than MutSigCV - it ignores sequence context,
replication timing and expression - but it is enough to show the effect and
to see the length artefacts fall away.

Output:
  results/length_normalised_genes.tsv
  results/frequency_vs_density.png
  results/cds_lengths.tsv           cached, so NCBI is queried only once
"""

import os
import sys
import time
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Entrez, SeqIO

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())

# NCBI asks every script that downloads data to identify itself with a
# contact address. Read it from the environment first so that no real email
# has to be committed to this public repository.
email = os.environ.get("NCBI_EMAIL") or config.get("email", "")
if "@" not in email:
    sys.exit("Set NCBI_EMAIL in your environment (or edit config.toml) - "
             "NCBI requires a contact address for downloads.")
Entrez.email = email

results = ROOT / "results"
genes = pd.read_csv(results / "top_mutated_genes.tsv", sep="\t")


TRANSCRIPTS_PER_GENE = 20


def fetch_cds_length(symbol: str) -> int | None:
    """Length of the longest annotated RefSeq coding sequence, in bases.

    Most genes have several RefSeq transcripts, and simply taking the first
    one returns whichever record NCBI happens to rank first - sometimes a
    short isoform. For TTN that gave 444 bp instead of roughly 100,000,
    which reverses the whole point of this step. Taking the longest CDS
    across the returned transcripts is a stable, stated rule.
    """
    query = (f'{symbol}[Gene Name] AND "Homo sapiens"[Organism] '
             f'AND srcdb_refseq_known[PROP] AND biomol_mrna[PROP]')
    with Entrez.esearch(db="nuccore", term=query,
                        retmax=TRANSCRIPTS_PER_GENE) as handle:
        ids = Entrez.read(handle)["IdList"]
    if not ids:
        return None
    with Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb",
                       retmode="text") as handle:
        lengths = [
            len(feature.location)
            for record in SeqIO.parse(handle, "genbank")
            for feature in record.features
            if feature.type == "CDS"
        ]
    return max(lengths) if lengths else None


cache_path = results / "cds_lengths.tsv"
if cache_path.exists():
    lengths = pd.read_csv(cache_path, sep="\t").set_index("Hugo_Symbol")["cds_bp"]
    print(f"Using cached CDS lengths for {len(lengths)} genes.")
else:
    print(f"Fetching CDS lengths from NCBI for {len(genes)} genes ...")
    collected = {}
    for symbol in genes["Hugo_Symbol"]:
        try:
            collected[symbol] = fetch_cds_length(symbol)
            print(f"  {symbol:<10} {collected[symbol]}")
        except Exception as error:                      # network hiccups
            print(f"  {symbol:<10} failed ({error})")
            collected[symbol] = None
        time.sleep(0.4)                                 # respect NCBI rate limits
    lengths = pd.Series(collected, name="cds_bp")
    lengths.index.name = "Hugo_Symbol"
    lengths.to_csv(cache_path, sep="\t")

genes["cds_bp"] = genes["Hugo_Symbol"].map(lengths)
genes = genes.dropna(subset=["cds_bp"]).copy()

# Sanity check on the lengths themselves. TTN encodes titin, by a wide
# margin the largest human protein, so if TTN is not the longest coding
# sequence in the table then the wrong transcripts were fetched and every
# number below would be wrong.
if "TTN" in set(genes["Hugo_Symbol"]):
    longest = genes.loc[genes["cds_bp"].idxmax(), "Hugo_Symbol"]
    if longest != "TTN":
        sys.exit(
            f"CDS lengths look wrong: the longest is {longest}, not TTN.\n"
            f"Delete {cache_path.name} and re-run to fetch them again."
        )
    # cds_bp / 3 counts codons, which includes the stop codon, so this is
    # one more than the protein length in residues.
    titin_codons = genes.loc[genes["Hugo_Symbol"] == "TTN",
                             "cds_bp"].iloc[0] / 3
    print(f"Sanity check: TTN is the longest CDS "
          f"({titin_codons:,.0f} codons).")
genes["cds_kb"] = genes["cds_bp"] / 1000
genes["mutations_per_kb"] = genes["total"] / genes["cds_kb"]

by_density = genes.sort_values("mutations_per_kb", ascending=False).copy()
by_density["rank_by_frequency"] = (
    genes["MutatedSamples"].rank(ascending=False).astype(int)
)
by_density["rank_by_density"] = range(1, len(by_density) + 1)
by_density["rank_change"] = (
    by_density["rank_by_frequency"] - by_density["rank_by_density"]
)
by_density.to_csv(results / "length_normalised_genes.tsv", sep="\t", index=False)

print("\nRanked by mutations per kb of coding sequence:")
print(by_density[["Hugo_Symbol", "MutatedSamples", "cds_kb",
                  "mutations_per_kb", "rank_by_frequency",
                  "rank_by_density"]].head(12).to_string(index=False))

dropped = by_density.nsmallest(5, "rank_change")
print("\nGenes that fall furthest once length is accounted for:")
print(dropped[["Hugo_Symbol", "cds_kb", "rank_by_frequency",
               "rank_by_density"]].to_string(index=False))

# --- Figure: frequency versus density --------------------------------------
figure, axis = plt.subplots(figsize=(9, 7), constrained_layout=True)
axis.scatter(genes["cds_kb"], genes["mutations_per_kb"],
             s=genes["MutatedSamples"] / 3, alpha=0.65, color="#3d7ea6",
             edgecolor="white", linewidth=0.5)
# Label only the genes worth reading: the densest few, and the ones that
# fall furthest once length is accounted for. Labelling all 25 collides.
worth_labelling = set(by_density.head(5)["Hugo_Symbol"]) | set(
    by_density.nsmallest(5, "rank_change")["Hugo_Symbol"]
)
for _, gene in genes.iterrows():
    if gene["Hugo_Symbol"] in worth_labelling:
        axis.annotate(gene["Hugo_Symbol"],
                      (gene["cds_kb"], gene["mutations_per_kb"]),
                      fontsize=9, fontweight="bold",
                      xytext=(5, 4), textcoords="offset points")
axis.set_xscale("log")
axis.set_xlabel("coding sequence length (kb, log scale)")
axis.set_ylabel("mutations per kb of coding sequence")
axis.set_title("TCGA-LUAD: long genes collect mutations without driving cancer\n"
               "(bubble size = number of patients mutated)")
figure.savefig(results / "frequency_vs_density.png", dpi=150)

print("\nSaved: results/length_normalised_genes.tsv, "
      "results/frequency_vs_density.png")
