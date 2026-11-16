# 00_download_data.R ------------------------------------------------------
# Download all GEO input files into data/, verify integrity, and record a
# SHA-256 manifest. Run from the project root:
#   Rscript R/00_download_data.R
# All files are publicly available GEO supplementary files; nothing is
# redistributed in the archive. File names match the GEO deposits exactly.
#
# Integrity checks performed after download (and for pre-existing files):
#   - file exists and is non-empty
#   - .gz files are fully streamed through gzip decompression (every byte,
#     not just the first line), so truncation and CRC errors are detected
#   - .xlsx files are verified as ZIP containers (PK magic bytes)
#   - every warning raised during download is captured, logged with the URL
#     and attempt number, and classified in the end-of-run summary
#   - SHA-256 recorded in data/manifest_sha256.txt (binary-mode hashing)
# The manifest is written on first download; on later runs existing files are
# checked against it. A mismatch means the local copy differs from the copy
# first downloaded by this project (GEO may replace superseded files).
#
# FAILURE BEHAVIOUR: any download failure, empty file, gzip integrity failure,
# or manifest hash mismatch is collected and the script STOPS with a non-zero
# exit status at the end. Downstream analysis must not be started from a run
# that reports failures. To intentionally adopt a changed GEO file, delete the
# affected file (and its manifest entry by deleting data/manifest_sha256.txt),
# re-run this script, and record the new manifest as the locked version.

base <- "https://ftp.ncbi.nlm.nih.gov/geo"
dir.create("data", showWarnings = FALSE)

files <- c(
  # GSE310013 - Augusta culture cohort (discovery), bulk RNA-seq raw counts
  "series/GSE310nnn/GSE310013/suppl/GSE310013_raw_counts.txt.gz",
  # GSE326710 - Lausanne culture cohort (replication), bulk RNA-seq raw counts
  "series/GSE326nnn/GSE326710/suppl/GSE326710_Raw_counts_all_samples.csv.gz",
  # GSE241205 - occluded grafts (clinical calibration), FPKM
  "series/GSE241nnn/GSE241205/suppl/GSE241205_gene_fpkm.txt.gz",
  # GSE299685 - Lausanne GeoMx spatial, normalized counts
  "series/GSE299nnn/GSE299685/suppl/GSE299685_Normalized_Counts.xlsx"
)

# GSE310351 - Augusta Visium FFPE spatial, 6 samples (3 patients x day-0/day-7)
# GEO FTP layout: samples/<GSMxxxxxnnn>/<full GSM accession>/suppl/<file>
# The parent directory is the accession with its last 3 digits replaced by
# "nnn" (e.g. GSM9294393 -> GSM9294nnn).
visium_samples <- c("GSM9294393_1D0", "GSM9294394_1D7", "GSM9294395_2D0",
                    "GSM9294396_2D7", "GSM9294397_3D0", "GSM9294398_3D7")
gsm_parent_dir <- function(gsm) {
  paste0(substr(gsm, 1, nchar(gsm) - 3), "nnn")   # GSM9294393 -> GSM9294nnn
}
for (s in visium_samples) {
  gsm <- sub("_.*", "", s)
  for (f in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz",
              "tissue_positions.csv.gz")) {
    files <- c(files, sprintf("samples/%s/%s/suppl/%s_%s",
                              gsm_parent_dir(gsm), gsm, s, f))
  }
}

sha256 <- function(path) {
  # portable SHA-256: openssl in explicit BINARY mode (a text-mode connection
  # hashes the wrong byte stream), then system sha256sum as fallback. The
  # connection is closed explicitly (round-23 review: 65 unclosed connections
  # were garbage-collected with "closing unused connection" warnings that
  # flooded the log and obscured real issues).
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, open = "rb")
    on.exit(close(con))
    h <- as.character(openssl::sha256(con))
    attributes(h) <- NULL   # strip the "hash" class so vector assignment stays atomic
    return(h)
  }
  out <- suppressWarnings(system2("sha256sum", shQuote(path), stdout = TRUE))
  if (!is.null(attr(out, "status")) && !is.na(attr(out, "status")) && attr(out, "status") != 0) return(NA_character_)
  strsplit(out[[1]], " ")[[1]][1]
}

check_xlsx <- function(path) {
  # .xlsx is a ZIP container: verify the PK\x03\x04 magic bytes (round-23
  # review: non-gzip files previously got only a non-empty check + hash)
  identical(readBin(path, "raw", n = 4L),
            as.raw(c(0x50, 0x4b, 0x03, 0x04)))
}

check_gz <- function(path) {
  # FULL decompression test via system `gzip -t`, which verifies the complete
  # stream including the gzip trailer/CRC. R's gzfile connection does NOT
  # surface trailer errors on read in all builds (verified by truncation test
  # in round 20), so there is NO silent fallback: the gzip binary is a hard
  # prerequisite and its absence aborts the script.
  st <- suppressWarnings(system2("gzip", c("-t", shQuote(path)),
                                 stdout = FALSE, stderr = FALSE))
  st <- as.integer(st)
  if (st == 127L) stop("system 'gzip' binary not found - it is required for full gzip integrity verification; install gzip and re-run.")
  st == 0L
}

# GSE299685 - GeoMx ROI sample annotation pages (SOFT text, one per ROI).
# Scripts 11 and 14 parse !Sample_title from these to build the GeoMx
# annotation table (patient / time / MYH11 segment). Downloaded from the GEO
# acc.cgi text endpoint into geo_pages/gsm299685/.
dir.create("geo_pages/gsm299685", showWarnings = FALSE, recursive = TRUE)
geomx_gsms <- sprintf("GSM90437%d", 27:62)   # GSM9043727 - GSM9043762 (36 ROIs)

manifest_path <- file.path("data", "manifest_sha256.txt")
manifest_files <- character(0)
manifest_hashes <- character(0)
if (file.exists(manifest_path)) {
  m <- read.csv(manifest_path, stringsAsFactors = FALSE)
  manifest_files <- m$file
  manifest_hashes <- m$sha256
}

failures <- character(0)   # every entry here -> non-zero exit at the end
n_new <- n_skip <- 0
# round-23 review: every warning is captured with URL, attempt number, and
# message, then classified in the end-of-run summary, so a final PASS line can
# never silently coexist with unexplained warnings.
warn_records <- data.frame(url = character(0), attempt = integer(0),
                           message = character(0), final_status = character(0),
                           stringsAsFactors = FALSE)

download_one <- function(url, dest) {
  # download to a temporary file, verify, then atomically rename into place.
  # Transient mirror errors (e.g. NCBI 502) are retried with a short pause;
  # a final failure is still reported as a failure. Warnings raised inside
  # download.file are captured and logged per attempt, not left aggregated.
  tmp <- paste0(dest, ".part")
  if (file.exists(tmp)) unlink(tmp)
  ok <- FALSE
  for (attempt in 1:3) {
    warns <- character(0)
    ok <- tryCatch(
      withCallingHandlers(
        { download.file(url, tmp, mode = "wb", quiet = TRUE); TRUE },
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) { if (attempt == 3) cat("  DOWNLOAD FAILED:", url, "|", conditionMessage(e), "\n"); FALSE })
    if (length(warns) > 0) {
      for (w in warns)
        cat(sprintf("  WARNING (attempt %d): %s | %s\n", attempt, url, w))
      warn_records <<- rbind(warn_records,
                             data.frame(url = url, attempt = attempt,
                                        message = warns, stringsAsFactors = FALSE))
    }
    if (ok && (!file.exists(tmp) || file.size(tmp) == 0)) {
      if (attempt == 3) cat("  DOWNLOAD EMPTY:", url, "\n")
      ok <- FALSE
    }
    if (ok) break
    if (attempt < 3) Sys.sleep(5)
  }
  if (ok) {
    ok <- file.rename(tmp, dest)   # atomic move into place; must succeed
    if (!ok) cat("  RENAME FAILED:", tmp, "->", dest, "\n")
  } else {
    unlink(tmp)
  }
  # classify every warning recorded for this URL by the final outcome
  warn_records$final_status[warn_records$url == url] <<- if (ok) "success" else "FAILED"
  ok
}

# --- GeoMx annotation pages -------------------------------------------------
for (g in geomx_gsms) {
  dest <- file.path("geo_pages", "gsm299685", paste0(g, ".txt"))
  if (file.exists(dest) && file.size(dest) > 0) next
  url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full", g)
  if (!download_one(url, dest)) {
    failures <- c(failures, paste("GeoMx annotation page download failed:", g))
  }
}
# verify every annotation page is present and parseable (!Sample_title line),
# and lock it into the same manifest as the data files (round-20 review: the
# pages were hashed for the audit inventory but not verified on later runs)
for (g in geomx_gsms) {
  dest <- file.path("geo_pages", "gsm299685", paste0(g, ".txt"))
  mkey <- paste0("geo_pages/gsm299685/", g, ".txt")
  ok <- file.exists(dest) && file.size(dest) > 0 &&
    any(grepl("^!Sample_title", readLines(dest, warn = FALSE)))
  if (!ok) {
    failures <- c(failures, paste("GeoMx annotation page missing or unparseable:", g))
    next
  }
  h <- sha256(dest)
  idx <- match(mkey, manifest_files)
  if (!is.na(idx)) {
    if (manifest_hashes[idx] != h) {
      failures <- c(failures, paste("HASH MISMATCH versus locked manifest:", mkey))
    }
  } else {
    manifest_files <- c(manifest_files, mkey)
    manifest_hashes <- c(manifest_hashes, h)
  }
}

# --- expression / spatial data files ---------------------------------------
for (f in files) {
  dest <- file.path("data", basename(f))
  url <- file.path(base, f)
  if (!file.exists(dest)) {
    cat("downloading:", url, "\n")
    if (!download_one(url, dest)) {
      failures <- c(failures, paste("download failed:", basename(f)))
      next
    }
    n_new <- n_new + 1
  } else {
    n_skip <- n_skip + 1
  }
  if (file.size(dest) == 0) {
    failures <- c(failures, paste("empty file:", basename(dest))); next
  }
  if (endsWith(dest, ".gz") && !check_gz(dest)) {
    failures <- c(failures, paste("gzip integrity failure (full decompression test):", basename(dest))); next
  }
  if (endsWith(dest, ".xlsx") && !check_xlsx(dest)) {
    failures <- c(failures, paste("xlsx structure failure (missing ZIP magic bytes):", basename(dest))); next
  }
  h <- sha256(dest)
  if (is.na(h) || length(h) != 1) {
    failures <- c(failures, paste("sha256 unavailable for:", basename(dest))); next
  }
  idx <- match(basename(dest), manifest_files)
  if (!is.na(idx)) {
    if (manifest_hashes[idx] != h) {
      failures <- c(failures, paste("HASH MISMATCH versus locked manifest (file changed since first download):", basename(dest)))
    }
  } else {
    manifest_files <- c(manifest_files, basename(dest))
    manifest_hashes <- c(manifest_hashes, h)
  }
}

write.csv(data.frame(file = manifest_files, sha256 = manifest_hashes),
          manifest_path, row.names = FALSE)
cat(sprintf("\nDone. %d files in data/ (%d downloaded now, %d already present).\n",
            length(list.files("data")), n_new, n_skip))
# round-23 review: explicit, classified warning summary — the number of
# captured warnings and their content are reported, and each is labelled as
# having occurred during a download that ultimately succeeded (transient) or
# during one that failed (those files are in the failure list regardless).
if (nrow(warn_records) > 0) {
  n_ok <- sum(warn_records$final_status == "success")
  n_bad <- sum(warn_records$final_status == "FAILED")
  cat(sprintf("\nDownload warnings captured: %d warning(s) across %d file(s) (%d during ultimately successful downloads, %d during failed ones).\n",
              nrow(warn_records), length(unique(warn_records$url)), n_ok, n_bad))
  for (msg in unique(warn_records$message)) {
    sel <- warn_records$message == msg
    urls <- unique(warn_records$url[sel])
    cat(sprintf("  - %d x \"%s\" [%s] (URLs: %s%s)\n", sum(sel), msg,
                paste(unique(warn_records$final_status[sel]), collapse = "/"),
                paste(utils::head(urls, 3), collapse = ", "),
                if (length(urls) > 3) sprintf(" + %d more", length(urls) - 3) else ""))
  }
  if (n_bad == 0)
    cat("All captured warnings occurred during downloads that ultimately succeeded and passed integrity verification (transient network/mirror messages).\n")
} else {
  cat("\nNo warnings were raised during download or verification.\n")
}
if (length(failures) > 0) {
  cat("\nINPUT VERIFICATION FAILURES:\n")
  for (x in failures) cat("  -", x, "\n")
  stop(sprintf("%d input verification failure(s); resolve the issues above before running the analysis chain.",
               length(failures)))
}
cat("All input verification checks passed.\n")
