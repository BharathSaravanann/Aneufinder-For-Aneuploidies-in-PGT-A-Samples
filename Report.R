#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 2) {
  cat("Usage: Rscript analyze_aneufinder.R <input.RData> <output_directory>\n")
  cat("Example: Rscript analyze_aneufinder.R sample.RData /path/to/output/\n")
  quit(status=1)
}

input_file <- args[1]
output_dir <- args[2]

if(!file.exists(input_file)) {
  cat("Error: Input file does not exist:", input_file, "\n")
  quit(status=1)
}

if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
}

cat("Loading data from:", input_file, "\n")
results <- get(load(input_file))

library(dplyr)

segments_df <- as.data.frame(results$segments)

report <- segments_df %>%
  group_by(seqnames) %>%
  summarise(
    Copy_Number = round(sum(copy.number * width) / sum(width), 2),
    total_length = sum(width),
    length_0somy = sum(width[state == "0-somy"]),
    length_1somy = sum(width[state == "1-somy"]),
    length_2somy = sum(width[state == "2-somy"]),
    length_3somy = sum(width[state == "3-somy"])
  ) %>%
  mutate(
    percent_0somy = round(length_0somy / total_length * 100, 1),
    percent_1somy = round(length_1somy / total_length * 100, 1),
    percent_2somy = round(length_2somy / total_length * 100, 1),
    percent_3somy = round(length_3somy / total_length * 100, 1),
    dominant_state = case_when(
      percent_3somy == pmax(percent_0somy, percent_1somy, percent_2somy, percent_3somy) ~ "3-somy",
      percent_2somy == pmax(percent_0somy, percent_1somy, percent_2somy, percent_3somy) ~ "2-somy",
      percent_1somy == pmax(percent_0somy, percent_1somy, percent_2somy, percent_3somy) ~ "1-somy",
      TRUE ~ "0-somy"
    ),
    Call = case_when(
      dominant_state == "3-somy" ~ "TRISOMY",
      dominant_state == "2-somy" ~ "NORMAL",
      dominant_state == "1-somy" ~ "MONOSOMY",
      dominant_state == "0-somy" ~ "NULLISOMY"
    )
  ) %>%
  select(seqnames, Copy_Number, percent_0somy, percent_1somy, percent_2somy, percent_3somy, dominant_state, Call)

output_csv <- file.path(output_dir, "AneuFinder_Report.csv")
write.csv(report, output_csv, row.names=FALSE)
cat("Report saved to:", output_csv, "\n\n")

print(report)

cat("\n=== SUMMARY ===\n")
cat("Sample:", basename(input_file), "\n")
cat("Total reads:", results$qualityInfo$total.read.count, "\n")
cat("Normal:", sum(report$Call == "NORMAL"), "\n")
cat("Trisomy:", sum(report$Call == "TRISOMY"), "\n")
cat("Monosomy:", sum(report$Call == "MONOSOMY"), "\n")
cat("Nullisomy:", sum(report$Call == "NULLISOMY"), "\n")

abnormal <- report %>% filter(Call != "NORMAL")
if(nrow(abnormal) > 0) {
  cat("\n=== ABNORMAL CHROMOSOMES ===\n")
  print(abnormal)

  output_abnormal <- file.path(output_dir, "Abnormal_Chromosomes.csv")
  write.csv(abnormal, output_abnormal, row.names=FALSE)
  cat("\nAbnormal chromosomes saved to:", output_abnormal, "\n")

  cat("\n=== CLINICAL INTERPRETATION ===\n")
  cat("RESULT: ABNORMAL EMBRYO\n")
  cat("RECOMMENDATION: NOT SUITABLE FOR TRANSFER\n")
} else {
  cat("\n=== CLINICAL INTERPRETATION ===\n")
  cat("RESULT: EUPLOID EMBRYO\n")
  cat("RECOMMENDATION: SUITABLE FOR TRANSFER\n")
}

output_summary <- file.path(output_dir, "Summary.txt")
sink(output_summary)
cat("AneuFinder Analysis Summary\n")
cat("============================\n")
cat("Sample:", basename(input_file), "\n")
cat("Total reads:", results$qualityInfo$total.read.count, "\n")
cat("Normal chromosomes:", sum(report$Call == "NORMAL"), "\n")
cat("Trisomies:", sum(report$Call == "TRISOMY"), "\n")
cat("Monosomies:", sum(report$Call == "MONOSOMY"), "\n")
cat("Nullisomies:", sum(report$Call == "NULLISOMY"), "\n")
if(nrow(abnormal) > 0) {
  cat("\nAbnormal Chromosomes:\n")
  for(i in 1:nrow(abnormal)) {
    cat(sprintf("%s: %s (CN=%.2f)\n", abnormal$seqnames[i], abnormal$Call[i], abnormal$Copy_Number[i]))
  }
}
sink()

cat("\nSummary saved to:", output_summary, "\n")
cat("\nAnalysis complete!\n")
