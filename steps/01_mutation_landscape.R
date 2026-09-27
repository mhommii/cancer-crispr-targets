# Step 1 - What is recurrently mutated in this cancer type?
#
# Before designing a CRISPR guide you need a reason to target a particular
# gene. This step takes a real tumour cohort and asks which genes are
# mutated in many patients, which is the usual starting point for picking
# a candidate.
#
# Data: TCGA MC3, the harmonised somatic mutation calls across TCGA
# (Ellrott et al. 2018), loaded through the TCGAmutations package.
# Analysis: maftools (Mayakonda et al. 2018).
#
# A caution that matters for the whole project: recurrent mutation means a
# gene is altered often. It does not prove the tumour depends on that gene,
# and it does not make the gene a validated drug or editing target. Showing
# dependency needs functional work - CRISPR screens, for instance - not
# mutation frequency alone.
#
# Output:
#   results/oncoplot.png          the classic mutated-genes-by-patient view
#   results/top_mutated_genes.tsv counts per gene
#   results/mutation_summary.png  variant types and burden across the cohort

suppressMessages({
  library(TCGAmutations)
  library(maftools)
})

cohort <- "LUAD"   # lung adenocarcinoma
top_n  <- 25

# Run this from the repository root, so results/ lands in the right place.
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

cat("Loading TCGA-", cohort, " (MC3) ...\n", sep = "")
maf <- tcga_load(study = cohort, source = "MC3")

samples <- as.numeric(maf@summary[ID == "Samples", summary])
cat("Samples in cohort:", samples, "\n")
cat("Mutations recorded:", nrow(maf@data), "\n\n")

# --- Which genes are mutated in the most patients? -------------------------
gene_summary <- getGeneSummary(maf)
top_genes <- head(gene_summary, top_n)
top_genes$percent_of_samples <- round(100 * top_genes$MutatedSamples / samples, 1)

write.table(
  top_genes[, .(Hugo_Symbol, MutatedSamples, percent_of_samples,
                total, Missense_Mutation, Nonsense_Mutation,
                Frame_Shift_Del, Frame_Shift_Ins)],
  file = "results/top_mutated_genes.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Top 10 genes by number of patients mutated:\n")
print(top_genes[1:10, .(Hugo_Symbol, MutatedSamples, percent_of_samples)])

# --- Oncoplot: genes down the side, patients across the bottom -------------
png("results/oncoplot.png", width = 1400, height = 900, res = 130)
oncoplot(maf = maf, top = 20,
         titleText = paste0("TCGA-", cohort, " (MC3): 20 most frequently mutated genes"))
invisible(dev.off())

# --- Cohort-level summary --------------------------------------------------
png("results/mutation_summary.png", width = 1400, height = 900, res = 130)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
invisible(dev.off())

cat("\nSaved: results/top_mutated_genes.tsv, results/oncoplot.png,",
    "results/mutation_summary.png\n")
