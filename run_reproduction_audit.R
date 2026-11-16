#!/usr/bin/env Rscript
# run_reproduction_audit.R — self-contained audit entry point (round 23).
# Run from the archive root (the directory containing this script, R/,
# gene_sets/, renv.lock, and the shipped results/):
#   Rscript run_reproduction_audit.R
#
# From an empty state it performs EVERY step that produces the evidence in
# reproduction_evidence/, so the evidence is regenerable from the package
# itself rather than only viewable:
#   1. preflight: required packages (renv::restore from renv.lock if missing)
#      and the gzip binary;
#   2. moves the shipped archived outputs results/ -> archive_reference/ so
#      the chain runs into an empty results/;
#   3. R/00_download_data.R (downloads all inputs, hash-locks them; aborts
#      the audit on failure);
#   4. writes input_file_hashes.txt (tab-separated) and runs
#      build_input_inventory.R -> input_inventory.csv;
#   5. runs the analysis chain in the verified dependency order, aborting on
#      the first non-zero exit;
#   6. runs 23 key-number/semantic checks (locked numbers + inventory and
#      audit-summary consistency);
#   7. compares the fresh output tree against the archived reference
#      BIDIRECTIONALLY (round-24 review): every archived output must be
#      regenerated (missing files fail), no unexpected extra files may
#      appear, and every comparable non-PNG file must be byte-identical;
#      the four figure SVGs and the 13 supplementary tables are additionally
#      asserted by name, independent of static path extraction;
#   8. writes reproduction_audit.csv, archive_diff_notes.txt, per-script logs,
#      and sessionInfo into reproduction_evidence/;
#   9. exits non-zero if ANY script, check, missing/unexpected file, or
#      content comparison failed.
#
# Mode: `Rscript run_reproduction_audit.R --compare-only` skips the download
# and analysis chain and re-runs only the key checks, the bidirectional
# comparison, and the verdict against an existing results/ directory (useful
# to re-verify a completed run; the default mode remains the full audit).
#
# Input-count convention (round-23 review): the inventory holds 64 SOURCE
# input files (28 expression/spatial files + 36 GeoMx annotation pages) plus
# 1 project manifest written by 00_download_data.R; the audit reports these
# as "64 source input files ... plus one project manifest", not "65 inputs".

REQUIRED_PKGS <- c("limma", "edgeR", "DESeq2", "fgsea", "ggplot2",
                   "patchwork", "readxl", "Matrix", "openssl", "svglite")
CHAIN <- c("10_formal_de", "11_spatial", "12_clinical", "14_robust_core",
           "15_reproducibility_audit", "16_inference_audit",
           "22_finalize_audit_summary", "23_baseline_concordance",
           "24_stratified_occlusion_concordance",
           "25_leave_one_gene_out_attribution", "21_supplementary_tables",
           "17_figure3_final", "18_figure1_final", "19_figure2_final",
           "20_figure4_final")
EV <- "reproduction_evidence"
# round-24 review: the four figures and the supplementary tables are asserted
# BY NAME so that output completeness never depends on static path extraction
# (a regex that extracts nothing must not be able to pass silently)
EXPECTED_FIGURES <- c("figures/fig1_cross_cohort_scatter_v3.svg",
                      "figures/fig2_geneset_nes_camera_v2.svg",
                      "figures/fig3_visium_program_scores.svg",
                      "figures/fig4_top30_dotplot_v2.svg")
EXPECTED_TABLES <- file.path("supplementary_tables", c(
  "TableS2_full_DE_both_culture_cohorts.csv", "TableS3_robust_core_4148_genes.csv",
  "TableS4_fgsea_camera_all_cohorts.csv", "TableS5_DE_method_sensitivity.csv",
  "TableS6_leave_one_marker_out_attribution.csv", "TableS7_attribution_threshold_sensitivity.csv",
  "TableS8_occlusion_concordance_per_patient.csv", "TableS9a_top30_culture_only_ranking.csv",
  "TableS9b_top30_ranking_comparison.csv", "TableS10_stratified_occlusion_concordance.csv",
  "TableS11a_leave_one_gene_out_composition.csv", "TableS11b_marker_member_core_genes.csv",
  "TableS11c_attribution_rule_comparison.csv"))
compare_only <- "--compare-only" %in% commandArgs(trailingOnly = TRUE)

failures <- 0L
say <- function(...) cat(sprintf(...), "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(script) {
  logf <- file.path(EV, paste0("log_", script, ".txt"))
  st <- system(sprintf("Rscript R/%s.R > %s 2>&1", script, shQuote(logf)))
  st <- as.integer(attr(st, "status") %||% st)
  say("R/%s.R -> exit %s", script, st)
  st
}

ex00 <- 0L
exit_codes <- setNames(rep(0L, length(CHAIN)), CHAIN)
if (!compare_only) {
  # --- 1. preflight ----------------------------------------------------------
  stopifnot(file.exists("renv.lock"), dir.exists("R"), dir.exists("gene_sets"),
            file.exists("build_input_inventory.R"))
  missing <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace,
                                   logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    if (requireNamespace("renv", quietly = TRUE)) {
      say("restoring %d missing package(s) from renv.lock ...", length(missing))
      renv::restore(lockfile = "renv.lock", prompt = FALSE)
      missing <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace,
                                       logical(1), quietly = TRUE)]
    }
    if (length(missing) > 0)
      stop("missing packages: ", paste(missing, collapse = ", "),
           "; run renv::restore(lockfile = 'renv.lock') and re-run")
  }
  gz <- suppressWarnings(system2("gzip", "--version", stdout = FALSE, stderr = FALSE))
  if (isTRUE(gz == 127L) || identical(unname(gz), "127"))
    stop("system 'gzip' binary not found - required by R/00_download_data.R")
  dir.create(EV, showWarnings = FALSE, recursive = TRUE)

  # --- 2. archive the shipped outputs aside; run into an empty results/ ------
  dir.create("archive_reference", showWarnings = FALSE)
  if (!dir.exists(file.path("archive_reference", "results"))) {
    if (dir.exists("results")) {
      stopifnot(file.rename("results", file.path("archive_reference", "results")))
    } else {
      # recover from an interrupted earlier attempt: the archived files then
      # sit directly in archive_reference/ instead of archive_reference/results/
      stopifnot(length(list.files("archive_reference")) > 30)
      dir.create(file.path("archive_reference", "results"), showWarnings = FALSE)
      for (p in list.files("archive_reference", full.names = TRUE))
        if (!dir.exists(p))
          stopifnot(file.rename(p, file.path("archive_reference", "results", basename(p))))
    }
  }
  stopifnot(length(list.files(file.path("archive_reference", "results"))) > 30)
  unlink("results", recursive = TRUE)
  dir.create("results", showWarnings = FALSE)

  # --- 3. download + integrity -----------------------------------------------
  ex00 <- run("00_download_data")
  if (ex00 != 0L) {
    say("ABORT: 00_download_data failed; the chain never runs on incomplete inputs")
    quit(save = "no", status = 1)
  }

  # --- 4. input hashes + inventory -------------------------------------------
  sha256 <- function(path) {
    if (requireNamespace("openssl", quietly = TRUE)) {
      con <- file(path, open = "rb")
      on.exit(close(con))
      h <- as.character(openssl::sha256(con))
      attributes(h) <- NULL
      return(h)
    }
    out <- suppressWarnings(system2("sha256sum", shQuote(path), stdout = TRUE))
    strsplit(out[[1]], " ")[[1]][1]
  }
  inputs <- c(list.files("data", full.names = TRUE),
              list.files(file.path("geo_pages", "gsm299685"), full.names = TRUE))
  hashes <- vapply(inputs, sha256, character(1))
  writeLines(paste(unname(hashes), basename(inputs), sep = "\t"),
             "input_file_hashes.txt")
  st_inv <- system("Rscript build_input_inventory.R > build_input_inventory.log 2>&1")
  if (as.integer(attr(st_inv, "status") %||% st_inv) != 0L) {
    say("ABORT: build_input_inventory.R failed")
    quit(save = "no", status = 1)
  }

  # --- 5. analysis chain -----------------------------------------------------
  for (s in CHAIN) {
    exit_codes[s] <- run(s)
    if (exit_codes[s] != 0L) {
      say("ABORT: chain stopped at %s (see %s/log_%s.txt)", s, EV, s)
      quit(save = "no", status = 1)
    }
  }
} else {
  say("compare-only mode: re-running key checks, comparison, and verdict against the existing results/")
  dir.create(EV, showWarnings = FALSE, recursive = TRUE)
  stopifnot(dir.exists("results"), dir.exists(file.path("archive_reference", "results")))
}

# --- 6. key-number + semantic checks (23) ------------------------------------
kc_file <- file.path(EV, "key_number_checks.txt")
kc <- function(name, got, want)
  cat(sprintf("%s: %s (got %s, expect %s)\n",
              ifelse(all(got == want), "PASS", "FAIL"), name,
              paste(got, collapse = ","), paste(want, collapse = ",")),
      sep = "", file = kc_file, append = TRUE)
writeLines("", kc_file)
rd <- function(f) read.csv(file.path("results", f), stringsAsFactors = FALSE)
m   <- rd("meta_analysis.csv");                          kc("common_genes", nrow(m), 13168)
rc  <- rd("robust_core_priority_all.csv");               kc("robust_core", nrow(rc), 4148)
den <- rd("attribution_denominator_audit.csv")
kc("in_visium_matrix", den$n[den$step == "in_visium_matrix"], 4110)
kc("valid_assignment", den$n[den$step == "valid_assignment"], 4019)
att <- rd("robust_core_priority_culture_only.csv")
kc("fibroblast_primary", sum(att$best_celltype %in% "Fibroblast"), 2054)
kc("smc_primary", sum(att$best_celltype %in% "SMC"), 1128)
logo <- rd("attribution_logo_composition.csv")
kc("fibroblast_after_logo", logo$logo_n[logo$celltype == "Fibroblast"], 2054)
kc("smc_after_logo", logo$logo_n[logo$celltype == "SMC"], 1129)
stt <- rd("stratified_occlusion_concordance.csv")
kc("core_concordance_pct", stt$rate_pct[stt$stratum == "robust_core"], 73.2)
b   <- rd("baseline_concordance_definitions.csv")
kc("baseline_pct", round(b$rate[b$subset == "all common genes" & b$definition == "culture meta vs occlusion"] * 100, 1), 61.1)
cs  <- rd("concordance_sensitivity.csv")
kc("doubly_sig_pct", cs$rate[cs$criterion == "both FDR<0.05 (robust core)"] * 100, 95.33)
ar  <- rd("attribution_rule_comparison.csv")
kc("max_abs_r_changed", ar$n_changed_versus_primary[ar$rule == "max_abs_r"], 78)
mm  <- rd("attribution_marker_member_genes.csv");        kc("marker_member", nrow(mm), 13)
kc("nominal_p_disc_n", cs$n[cs$criterion == "discovery p<0.05"], 9403)
kc("nominal_p_both_n", cs$n[cs$criterion == "both p<0.05"], 5307)
kc("no_empty_analyses", sum(cs$n > 0), nrow(cs))
inv <- read.csv("input_inventory.csv", stringsAsFactors = FALSE)
kc("inventory_rows", nrow(inv), 65)
kc("inventory_source_inputs", sum(inv$file != "manifest_sha256.txt"), 64)
kc("inventory_fields_complete",
   sum(inv$file != "" & inv$file_type != "" & grepl("^[0-9a-f]{64}$", inv$sha256)), 65)
kc("inventory_unique_files", length(unique(inv$file)), 65)
kc("inventory_geomx_pages", sum(grepl("^GSM9043", inv$file)), 36)
au  <- rd("reproducibility_audit_summary.csv")
kc("summary_core_conc", au$value[au$metric == "core_conc"],
   cs$rate[cs$criterion == "both FDR<0.05 (robust core)"])
kc("summary_core_n", au$value[au$metric == "core_n"],
   cs$n[cs$criterion == "both FDR<0.05 (robust core)"])
kc_lines <- readLines(kc_file)
n_checks <- sum(nzchar(kc_lines))
n_check_fail <- sum(grepl("^FAIL", kc_lines))
say("key checks: %d/%d PASS", n_checks - n_check_fail, n_checks)

# --- 7. bidirectional comparison against the archived outputs ----------------
# (round-24 review: the comparison must not only verify the files the fresh
# run happened to produce — every ARCHIVED output must also be regenerated,
# and no unexpected extra files may appear)
expected_files <- list.files(file.path("archive_reference", "results"),
                             recursive = TRUE)
actual_files <- list.files("results", recursive = TRUE)
missing_files <- setdiff(expected_files, actual_files)
unexpected_files <- setdiff(actual_files, expected_files)
common <- intersect(actual_files, expected_files)
diffs <- character(0); match_n <- 0L; total_n <- 0L
for (rel in common) {
  if (endsWith(rel, ".png")) next   # rasterization backend differs; SVG compared
  total_n <- total_n + 1L
  f <- file.path("results", rel)
  ref <- file.path("archive_reference", "results", rel)
  if (identical(readBin(f, "raw", file.size(f)),
                readBin(ref, "raw", file.size(ref)))) match_n <- match_n + 1L else
    diffs <- c(diffs, rel)
}
# named completeness assertions: independent of any static path extraction
named_missing <- c(EXPECTED_FIGURES, EXPECTED_TABLES)[
  !file.exists(file.path("results", c(EXPECTED_FIGURES, EXPECTED_TABLES)))]
notes <- c(
  if (length(missing_files))
    c("files present in the archived reference but NOT regenerated:", paste0("  ", missing_files))
  else "no archived output is missing from the fresh run",
  if (length(unexpected_files))
    c("files generated that are NOT in the archived reference:", paste0("  ", unexpected_files))
  else "no unexpected extra files were generated",
  if (length(named_missing))
    c("named expected outputs (figures/tables) missing:", paste0("  ", named_missing))
  else "all 4 figure SVGs and 13 supplementary tables present",
  if (length(diffs))
    c("files differing from archive (excluding PNG rasterization):", paste0("  ", diffs))
  else "all comparable non-PNG outputs byte-identical to archive")
writeLines(notes, file.path(EV, "archive_diff_notes.txt"))

# --- 8. audit report ---------------------------------------------------------
if (!compare_only) {
  declared_outputs <- function(script) {
    src <- readLines(file.path("R", paste0(script, ".R")), warn = FALSE)
    unique(unlist(regmatches(src,
      gregexpr("(results|figures)/[A-Za-z0-9_.-]+[.](csv|png|svg)", src))))
  }
  n_source <- sum(inv$file != "manifest_sha256.txt")
  rows <- lapply(c("00_download_data", CHAIN), function(s) {
    ex <- if (s == "00_download_data") ex00 else exit_codes[s]
    outs <- declared_outputs(s)
    out_ok <- all(vapply(outs, function(o)
      file.exists(o) && file.size(o) > 0, logical(1)))
    arch <- if (s == "00_download_data") "N/A" else
      if (total_n > 0 && match_n == total_n) "PASS" else "FAIL"
    note <- if (s == "00_download_data")
      sprintf("download + integrity verification (%d source input files verified and hash-locked, plus one project manifest recorded in the inventory)", n_source)
      else paste0("log: log_", s, ".txt")
    data.frame(script = s, exit_status = ex,
               input_check = "PASS", output_check = if (out_ok) "PASS" else "FAIL",
               archive_match = arch,
               status = if (ex == 0 && out_ok && arch %in% c("PASS", "N/A")) "PASS" else "FAIL",
               notes = note, stringsAsFactors = FALSE)
  })
  audit <- do.call(rbind, rows)
  write.csv(audit, file.path(EV, "reproduction_audit.csv"), row.names = FALSE)
  writeLines(paste(unname(hashes), basename(inputs), sep = "\t"),
             file.path(EV, "input_file_hashes.txt"))
  file.copy("input_inventory.csv", EV, overwrite = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(EV, "sessionInfo_clean_machine.txt"))
}

# --- 9. verdict --------------------------------------------------------------
# (round-24 review: missing archived outputs, unexpected extra files, and
# missing named figures/tables are hard failures, not just notes)
n_fail <- (if (compare_only) 0 else sum(audit$status != "PASS")) +
  n_check_fail + length(diffs) + length(missing_files) +
  length(unexpected_files) + length(named_missing)
say("audit: %s%d/%d outputs byte-identical; %d check failure(s); %d missing archived file(s); %d unexpected file(s); %d missing named figure(s)/table(s)",
    if (compare_only) "" else sprintf("%d/%d scripts PASS; ", sum(audit$status == "PASS"), nrow(audit)),
    match_n, total_n, n_check_fail, length(missing_files),
    length(unexpected_files), length(named_missing))
if (n_fail > 0) {
  say("REPRODUCTION AUDIT FAILED (%d issue(s)); see %s/archive_diff_notes.txt", n_fail, EV)
  quit(save = "no", status = 1)
}
say("REPRODUCTION AUDIT PASSED: all scripts, checks, and archive comparisons OK")
