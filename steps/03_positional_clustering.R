# Step 3 - Which genes have mutations clustered in one spot?
#
# Step 2 corrected for gene length, which is crude. This step uses a
# published method that asks a sharper question: are a gene's mutations
# scattered along the protein, or piled up at the same few residues?
#
# The distinction matters biologically:
#
#   Oncogene          mutations cluster at specific residues, because only
#                     particular substitutions switch the protein on.
#                     KRAS codon 12 is the textbook example.
#
#   Tumour suppressor mutations are scattered and often truncating
#                     (nonsense, frameshift), because any change that
#                     breaks the protein will do. TP53 behaves this way.
#
# OncodriveCLUST (Tamborero et al. 2013), available in maftools as
# oncodrive(), scores that clustering. Passenger mutations in long genes
# are spread out, so genes like TTN and MUC16 score poorly here even though
# they topped the raw frequency list in step 1.
#
# This takes several minutes to run. The result is written to disk so the
# later steps do not need to repeat it.
#
# Output:
#   results/oncodrive_results.tsv
#   results/kras_lollipop.png   where KRAS mutations sit along the protein
#   results/tp53_lollipop.png   the contrasting tumour suppressor pattern

suppressMessages({
  library(TCGAmutations)
  library(maftools)
})

cohort <- "LUAD"
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

cat("Loading TCGA-", cohort, " (MC3) ...\n", sep = "")
maf <- tcga_load(study = cohort, source = "MC3")

cat("Scoring positional clustering (this takes a few minutes) ...\n")
clustering <- oncodrive(maf = maf, minMut = 5, pvalMethod = "zscore")
clustering <- clustering[order(fdr)]

write.table(clustering, file = "results/oncodrive_results.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

significant <- clustering[fdr < 0.1]
cat("\nGenes with significantly clustered mutations (FDR < 0.1):\n")
print(as.data.frame(
  significant[, c("Hugo_Symbol", "clusters", "muts_in_clusters", "total", "fdr")]
))

# Did the long genes from step 1 survive this test?
length_artefacts <- c("TTN", "MUC16", "CSMD3", "RYR2", "USH2A")
cat("\nHow the long genes from step 1 score here:\n")
print(as.data.frame(
  clustering[Hugo_Symbol %in% length_artefacts,
             c("Hugo_Symbol", "muts_in_clusters", "total", "fdr")]
))

# --- Lollipop plots: mutation position along the protein -------------------
# The shape of these two plots is the whole oncogene/tumour-suppressor point.
png("results/kras_lollipop.png", width = 1300, height = 600, res = 130)
lollipopPlot(maf = maf, gene = "KRAS", showMutationRate = TRUE,
             labelPos = "all")
invisible(dev.off())

png("results/tp53_lollipop.png", width = 1300, height = 600, res = 130)
lollipopPlot(maf = maf, gene = "TP53", showMutationRate = TRUE)
invisible(dev.off())

cat("\nSaved: results/oncodrive_results.tsv, results/kras_lollipop.png,",
    "results/tp53_lollipop.png\n")
