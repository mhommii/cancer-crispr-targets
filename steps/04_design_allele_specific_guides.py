"""Step 4 — Design guides against the KRAS codon 12 hotspot.

Why KRAS needs a different approach from TP53
---------------------------------------------
Steps 2 and 3 both put KRAS at the top: it has the highest mutation density
per kilobase, and OncodriveCLUST finds nearly all of its mutations packed
into a single cluster. That cluster is codon 12.

That changes what a sensible guide looks like. TP53 is a tumour suppressor
already broken in the tumour, so a guide that disrupts the coding sequence
anywhere early is a reasonable knockout design. KRAS is different:

  - It is an oncogene. The tumour is driven by one mutant copy.
  - The other copy is normal, and KRAS is needed by healthy cells.
  - A guide matching both copies would cut the normal allele too.

So the useful design question is not "where can I cut KRAS" but "can I cut
only the mutant copy". That is allele-specific editing, and whether it
works depends on where the single changed base falls inside the guide.

Two things make discrimination easier:

  1. The mutation sits in the seed region, the ~12 bases next to the PAM,
     where Cas9 tolerates mismatches poorly. A guide matching the mutant
     seed should cut the mutant and tolerate the wild-type badly.
  2. The mutation creates a PAM that does not exist in the wild-type
     sequence. This is the strongest case: with no NGG, Cas9 cannot cut the
     normal allele at that position at all.

This step finds guides overlapping codon 12 for the common G12 substitutions
and classifies them on exactly those two criteria.

Nothing here is validated. It is an analysis of sequence, not evidence that
any of these guides discriminate in cells.

Output:
  results/kras_allele_specific_guides.tsv
  results/kras_codon12_context.txt
"""

import os
import sys
import tomllib
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

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

GENE = "KRAS"
CODON = 12
PROTOSPACER_LENGTH = 20
SEED_LENGTH = 12

# The substitutions seen most often at KRAS codon 12. The wild-type codon is
# GGT (glycine); each entry is the mutant codon and the resulting residue.
G12_MUTATIONS = {
    "G12C": "TGT",
    "G12D": "GAT",
    "G12V": "GTT",
    "G12A": "GCT",
    "G12S": "AGT",
    "G12R": "CGT",
}

data, results = ROOT / "data", ROOT / "results"
data.mkdir(exist_ok=True)
results.mkdir(exist_ok=True)

# --- Fetch KRAS ------------------------------------------------------------
record_path = data / f"{GENE}.gb"
if not record_path.exists():
    print(f"Fetching {GENE} from NCBI ...")
    query = f'{GENE}[Gene Name] AND "Homo sapiens"[Organism] AND alive[property]'
    with Entrez.esearch(db="gene", term=query) as handle:
        gene_id = Entrez.read(handle)["IdList"][0]
    with Entrez.esummary(db="gene", id=gene_id) as handle:
        info = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]
    location = info["GenomicInfo"][0]
    start, stop = sorted((int(location["ChrStart"]), int(location["ChrStop"])))
    with Entrez.efetch(db="nuccore", id=location["ChrAccVer"],
                       rettype="gbwithparts", retmode="text",
                       seq_start=start + 1, seq_stop=stop + 1) as handle:
        SeqIO.write(SeqIO.read(handle, "genbank"), record_path, "genbank")

record = SeqIO.read(record_path, "genbank")
sequence = str(record.seq).upper()
print(f"{GENE}: {len(sequence):,} bp downloaded "
      f"({record.annotations.get('accessions', ['?'])[0]})")

# --- Locate codon 12 in genomic coordinates --------------------------------
# Build the genomic position of every coding base, in transcript order, so a
# codon number can be mapped back onto the chromosome. Done explicitly
# because KRAS is on the minus strand and getting this wrong silently
# targets the wrong codon.
coding_features = [f for f in record.features if f.type == "CDS"]
if not coding_features:
    sys.exit("No CDS feature found in the downloaded record.")

# Prefer the longest CDS (the principal isoform).
coding = max(coding_features, key=lambda f: len(f.location))
strand = coding.location.strand

positions = []
for part in coding.location.parts:
    if strand == -1:
        positions.extend(range(int(part.end) - 1, int(part.start) - 1, -1))
    else:
        positions.extend(range(int(part.start), int(part.end)))

codon_positions = positions[(CODON - 1) * 3:(CODON - 1) * 3 + 3]
codon_bases = "".join(sequence[p] for p in codon_positions)
if strand == -1:
    codon_bases = str(Seq(codon_bases).complement())

# Self-check: codon 12 of KRAS must be GGT, coding for glycine.
if codon_bases != "GGT":
    sys.exit(
        f"Codon {CODON} reads {codon_bases}, expected GGT (glycine).\n"
        f"The coordinate mapping is wrong, so every guide below would be "
        f"designed against the wrong position."
    )
print(f"Codon {CODON} = {codon_bases} (glycine), "
      f"strand {'-' if strand == -1 else '+'}  [self-check passed]")

codon_span = sorted(codon_positions)
codon_start, codon_end = codon_span[0], codon_span[-1] + 1


def make_mutant(mutant_codon: str) -> str:
    """The genomic sequence carrying one codon-12 substitution."""
    genomic_codon = (str(Seq(mutant_codon).complement()) if strand == -1
                     else mutant_codon)
    # Genomic positions run the other way on the minus strand.
    ordered = codon_positions if strand == 1 else codon_positions[::-1]
    mutated = list(sequence)
    for position, base in zip(sorted(ordered), genomic_codon if strand == 1
                              else genomic_codon[::-1]):
        mutated[position] = base
    return "".join(mutated)


def find_guides_over_codon(target_sequence: str):
    """Guides whose protospacer or PAM covers codon 12."""
    found = []
    window_start = max(0, codon_start - 30)
    window_end = min(len(target_sequence), codon_end + 30)
    for position in range(window_start, window_end):
        # Forward strand: PAM = NGG at position..position+2
        if (target_sequence[position + 1:position + 3] == "GG"
                and position >= PROTOSPACER_LENGTH):
            protospacer = target_sequence[position - PROTOSPACER_LENGTH:position]
            span = (position - PROTOSPACER_LENGTH, position + 3)
            if span[0] < codon_end and span[1] > codon_start:
                found.append(("+", protospacer,
                              target_sequence[position:position + 3], span))
        # Reverse strand: CCN here reads NGG on the other strand
        if target_sequence[position:position + 2] == "CC":
            end = position + 3 + PROTOSPACER_LENGTH
            if end <= len(target_sequence):
                protospacer = str(Seq(
                    target_sequence[position + 3:end]).reverse_complement())
                pam = str(Seq(
                    target_sequence[position:position + 3]).reverse_complement())
                span = (position, end)
                if span[0] < codon_end and span[1] > codon_start:
                    found.append(("-", protospacer, pam, span))
    return found


wild_type_guides = {
    (strand_sign, protospacer)
    for strand_sign, protospacer, _, _ in find_guides_over_codon(sequence)
}
wild_type_pam_sites = {
    (strand_sign, span)
    for strand_sign, _, _, span in find_guides_over_codon(sequence)
}

rows = []
for name, mutant_codon in G12_MUTATIONS.items():
    mutant_sequence = make_mutant(mutant_codon)
    for strand_sign, protospacer, pam, span in find_guides_over_codon(mutant_sequence):
        # The equivalent wild-type protospacer at the same coordinates.
        if strand_sign == "+":
            wild_type = sequence[span[0]:span[0] + PROTOSPACER_LENGTH]
        else:
            wild_type = str(Seq(
                sequence[span[1] - PROTOSPACER_LENGTH:span[1]]
            ).reverse_complement())

        differences = [i for i, (a, b) in enumerate(zip(protospacer, wild_type))
                       if a != b]
        # Does this PAM exist in the wild-type sequence at all? If not, the
        # mutation created it, and Cas9 simply cannot cut the normal allele
        # here - the strongest form of discrimination there is.
        creates_new_pam = (strand_sign, span) not in wild_type_pam_sites

        if not differences and not creates_new_pam:
            # Same protospacer, same PAM as wild-type: no way to tell the
            # two alleles apart.
            continue

        if differences:
            # Distance from the PAM: index 19 sits immediately next to it.
            distance_from_pam = min(PROTOSPACER_LENGTH - 1 - i
                                    for i in differences)
            in_seed = distance_from_pam < SEED_LENGTH
        else:
            # The change is in the PAM itself, not the protospacer.
            distance_from_pam = -1
            in_seed = False

        rows.append({
            "mutation": name,
            "mutant_codon": mutant_codon,
            "strand": strand_sign,
            "guide_mutant_allele": protospacer,
            "wild_type_sequence": wild_type,
            "pam": pam,
            "mismatches_vs_wild_type": len(differences),
            "distance_from_pam": distance_from_pam,
            "mutation_in_seed": in_seed,
            "creates_new_pam": creates_new_pam,
            # Best case first: a PAM the wild-type does not have, then a
            # mutation sitting closest to the PAM.
            "discrimination": ("new PAM" if creates_new_pam
                               else "seed mismatch" if in_seed
                               else "distal mismatch"),
        })

guides = pd.DataFrame(rows)
if guides.empty:
    sys.exit("No allele-specific guides found - check the codon mapping.")

priority = {"new PAM": 0, "seed mismatch": 1, "distal mismatch": 2}
guides["priority"] = guides["discrimination"].map(priority)
guides.sort_values(["mutation", "priority", "distance_from_pam"], inplace=True)
guides.drop(columns="priority", inplace=True)
guides.to_csv(results / "kras_allele_specific_guides.tsv", sep="\t", index=False)

# --- Readable context ------------------------------------------------------
context_lines = [
    f"{GENE} codon {CODON} in genomic context",
    "=" * 44,
    f"Gene on the {'minus' if strand == -1 else 'plus'} strand; "
    f"codon {CODON} at positions {codon_start}-{codon_end - 1} "
    f"of the downloaded region.",
    "",
    f"wild-type codon: {codon_bases} (glycine)",
    "",
    "Guides found per substitution:",
]
for name in G12_MUTATIONS:
    subset = guides[guides["mutation"] == name]
    new_pam = (subset["discrimination"] == "new PAM").sum()
    seed = (subset["discrimination"] == "seed mismatch").sum()
    context_lines.append(
        f"  {name}: {len(subset):>2} guides  "
        f"({new_pam} create a new PAM, {seed} put the change in the seed)"
    )
(results / "kras_codon12_context.txt").write_text(
    "\n".join(context_lines) + "\n", encoding="utf-8")

print(f"\n{len(guides)} candidate allele-specific guides across "
      f"{len(G12_MUTATIONS)} substitutions.")
print("\nBy discrimination mechanism:")
print(guides["discrimination"].value_counts().to_string())
print("\nBest candidate per substitution:")
best = guides.groupby("mutation").first().reset_index()
print(best[["mutation", "strand", "guide_mutant_allele", "pam",
            "distance_from_pam", "discrimination"]].to_string(index=False))
print("\nSaved: results/kras_allele_specific_guides.tsv, "
      "results/kras_codon12_context.txt")
