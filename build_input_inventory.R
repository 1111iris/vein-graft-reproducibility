#!/usr/bin/env Rscript
# build_input_inventory.R — enriched input inventory: GEO accession, file role,
# SHA-256, and download date for every input the download step produced.
# Run from the project root after 00_download_data.R and the runner's hashing
# step. Round-22: parses sha256sum output robustly (whitespace-separated, no
# naive splitting) and hard-asserts field completeness, filename uniqueness,
# hash format, and the expected row count before writing.
h <- read.table("input_file_hashes.txt", sep = "\t",
                col.names = c("sha256", "path"), stringsAsFactors = FALSE,
                quote = "", comment.char = "")
# the runner writes "hash<TAB>path"; if a tab is missing, fail loudly instead
# of producing a broken table (round-22 review: the shipped inventory had NA
# fields because a space-separated parse silently merged hash and filename)
stopifnot(nrow(h) > 0, !any(is.na(h$sha256)), !any(is.na(h$path)),
          !any(grepl(" ", h$sha256)))
h$file <- basename(h$path)

role <- c(
  "GSE310013_raw_counts.txt.gz"             = "bulk RNA-seq raw counts | Augusta culture cohort (discovery), 6 patients, paired day-0/day-7",
  "GSE326710_Raw_counts_all_samples.csv.gz" = "bulk RNA-seq raw counts | Lausanne culture cohort (replication), 4 patients, paired day-0/day-7",
  "GSE241205_gene_fpkm.txt.gz"              = "gene FPKM | occluded grafts (clinical calibration), 3 patients, occluded vs spare vein",
  "GSE299685_Normalized_Counts.xlsx"        = "normalized counts | Lausanne GeoMx spatial, 4 patients, media vs intima/adventitia ROIs")
h$geo_accession <- ifelse(grepl("^GSM9294", h$file), "GSE310351",
                   ifelse(grepl("^GSM9043", h$file), "GSE299685",
                    sub("_.*", "", h$file)))
h$file_type <- ifelse(h$file == "manifest_sha256.txt",
                      "project manifest (data/manifest_sha256.txt; written by 00_download_data.R, not a GEO input)",
               ifelse(grepl("^GSM9294", h$file),
                      paste0("Visium FFPE spatial | ", sub("^GSM929439[0-9]_", "", h$file)),
               ifelse(grepl("^GSM9043", h$file),
                      "GeoMx ROI SOFT annotation page (!Sample_title parsing; 36 ROIs)",
                      unname(role[h$file]))))
h$download_date <- format(Sys.Date(), "%Y-%m-%d")
inv <- h[, c("file", "geo_accession", "file_type", "sha256", "download_date")]

# hard assertions: complete, unique, well-formed, expected size
stopifnot(nrow(inv) == 65,                                   # 28 data + 36 pages + 1 manifest
          all(!is.na(inv$file) & inv$file != ""),
          all(!is.na(inv$geo_accession) | inv$file == "manifest_sha256.txt"),
          all(!is.na(inv$file_type) & inv$file_type != ""),
          all(grepl("^[0-9a-f]{64}$", inv$sha256)),          # 64-hex SHA-256
          !any(duplicated(inv$file)),
          sum(grepl("^GSM9043", inv$file)) == 36,            # all 36 GeoMx pages
          sum(grepl("^GSM9294", inv$file)) == 24)            # all 24 Visium files
write.csv(inv, "input_inventory.csv", row.names = FALSE)
cat("input_inventory.csv written:", nrow(inv), "rows; all assertions passed\n")
