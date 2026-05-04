#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 IBDmix report: sample-level summaries by archaic reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic:
# Rscript ibdmix_report.R --in /mnt/d/analysis/gu/ibdmix/report/len1000.lod5.segments.tsv.gz --outdir /mnt/d/analysis/gu/ibdmix/report
#
# Optional sample-info report:
# Rscript ibdmix_report.R --in ...segments.tsv.gz --outdir .../report --sample_info 1kg.v3.sample.txt --summary_by super_pop,gender
#
# Back-compatible old phenotype report:
# Rscript ibdmix_report.R --in ...segments.tsv.gz --outdir .../report --pheno phenotype.tsv --id_col ID --group_col pop

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NA_character_) {
	i <- match(key, args)
	if (is.na(i) || i == length(args)) return(default)
	args[i + 1]
}

infile <- get_arg("--in")
outdir <- get_arg("--outdir", ".")
genome_length <- as.numeric(get_arg("--genome_length", "3000000000"))
ref_divisor_arg <- get_arg("--nref", NA_character_)
pheno_file <- get_arg("--pheno", NA_character_)
id_col <- get_arg("--id_col", "ID")
group_col <- get_arg("--group_col", "pop")
sample_info_file <- get_arg("--sample_info", NA_character_)
sample_id_col <- get_arg("--sample_id_col", "sample")
summary_by_arg <- get_arg("--summary_by", NA_character_)

if (is.na(infile)) stop("Please provide --in input.tsv.gz")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (is.na(sample_info_file) && !is.na(pheno_file)) {
	sample_info_file <- pheno_file
	sample_id_col <- id_col
	if (is.na(summary_by_arg)) summary_by_arg <- group_col
}

summary_by <- character(0)
if (!is.na(summary_by_arg) && nzchar(summary_by_arg)) {
	summary_by <- trimws(unlist(strsplit(summary_by_arg, ",", fixed = TRUE)))
	summary_by <- summary_by[nzchar(summary_by)]
}

find_input <- function(f, outdir) {
	bn <- basename(f)
	bn2 <- sub("\\.tsv\\.gz$", ".segments.tsv.gz", bn)
	cand <- unique(c(
		f,
		file.path(outdir, bn),
		sub("\\.tsv\\.gz$", ".segments.tsv.gz", f),
		file.path(outdir, bn2),
		file.path(dirname(f), "report", bn),
		file.path(dirname(f), "report", bn2)
	))
	cand[file.exists(cand)][1]
}

in0 <- find_input(infile, outdir)
if (is.na(in0)) stop("Input file not found: ", infile)
infile <- in0

prefix <- basename(infile)
prefix <- sub("\\.segments\\.tsv\\.gz$", "", prefix)
prefix <- sub("\\.tsv\\.gz$", "", prefix)
prefix <- sub("\\.gz$", "", prefix)

message("Input:  ", infile)
message("Outdir: ", outdir)
message("Prefix: ", prefix)
message("Genome length: ", genome_length)
if (!is.na(sample_info_file)) message("Sample info: ", sample_info_file)
if (length(summary_by) > 0) message("Summary by: ", paste(summary_by, collapse = ","))

count_do <- function(f) {
	cmd <- paste("gzip -dc", shQuote(f), "| awk 'BEGIN{n=0} $1==\"DO\"{n++} END{print n}'")
	as.integer(system(cmd, intern = TRUE))
}
n_do <- tryCatch(count_do(infile), error = function(e) NA_integer_)

cmd <- paste("gzip -dc", shQuote(infile), "| awk 'BEGIN{FS=OFS=\"\\t\"} $1!=\"DO\"'")
dat <- fread(cmd = cmd, sep = "\t", header = TRUE, fill = TRUE, showProgress = TRUE)

need <- c("ID", "chrom", "start", "end", "length", "slod", "sites", "positive_lods", "negative_lods", "anc")
miss <- setdiff(need, names(dat))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

if (!("sample_set" %in% names(dat))) {
	if ("pop" %in% names(dat)) dat[, sample_set := as.character(pop)] else dat[, sample_set := "sample_keep"]
}

dat <- dat[!is.na(ID) & ID != "" & ID != "ID" & ID != "DO"]
dat[, chrom := as.character(chrom)]
dat[, start := as.numeric(start)]
dat[, end := as.numeric(end)]
dat[, seg_len := as.numeric(length)]
dat[, slod := as.numeric(slod)]
dat[, sites := as.numeric(sites)]
dat[, positive_lods := as.numeric(positive_lods)]
dat[, negative_lods := as.numeric(negative_lods)]
dat[, anc := as.character(anc)]
dat[, sample_set := as.character(sample_set)]
dat <- dat[!is.na(start) & !is.na(end) & !is.na(seg_len) & !is.na(anc) & anc != ""]

n0 <- nrow(dat)
dat <- unique(dat, by = c("ID", "anc", "chrom", "start", "end"))
n_dup_removed <- n0 - nrow(dat)

ancs <- sort(unique(dat$anc))
if (!is.na(ref_divisor_arg) && nzchar(ref_divisor_arg)) {
	ref_divisor <- as.numeric(ref_divisor_arg)
	ref_divisor_source <- "--nref"
} else {
	ref_divisor <- length(ancs)
	ref_divisor_source <- "auto_from_anc"
}
if (is.na(ref_divisor) || ref_divisor <= 0) stop("Bad --nref value or no archaic references found in anc column")

message("Rows after cleanup: ", nrow(dat))
message("Archaic references: ", paste(ancs, collapse = ", "))
message("Total panel divisor: ", ref_divisor, " (", ref_divisor_source, ")")
if (length(ancs) != ref_divisor) {
	message("Note: observed ", length(ancs), " archaic references, but total panel divisor is ", ref_divisor, ".")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Input check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input.check <- data.table(
	input = infile,
	n_rows = nrow(dat),
	n_samples = uniqueN(dat$ID),
	n_archaic_references = length(ancs),
	archaic_references = paste(ancs, collapse = ","),
	n_DO_lines_skipped = n_do,
	n_exact_duplicate_segments_removed = n_dup_removed,
	genome_length = genome_length,
	ref_divisor_for_total = ref_divisor,
	ref_divisor_source = ref_divisor_source,
	sample_info = ifelse(is.na(sample_info_file), "", sample_info_file),
	sample_id_col = sample_id_col,
	summary_by = paste(summary_by, collapse = ",")
)
fwrite(input.check, file.path(outdir, paste0(prefix, ".input.check.tsv")), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Sample-level summaries by archaic reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample.by_anc <- dat[, .(
	n_segments = .N,
	total_length = sum(seg_len, na.rm = TRUE),
	percentage = sum(seg_len, na.rm = TRUE) / genome_length * 100,
	mean_slod = mean(slod, na.rm = TRUE),
	median_slod = median(slod, na.rm = TRUE),
	mean_segment_length = mean(seg_len, na.rm = TRUE),
	median_segment_length = median(seg_len, na.rm = TRUE)
), by = .(ID, anc)]

setorder(sample.by_anc, ID, anc)
fwrite(sample.by_anc, file.path(outdir, paste0(prefix, ".sample.by_anc.sum.tsv")), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Sample-level total summary across archaic references, divided by nref
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample.total <- dat[, .(
	n_ref_observed = uniqueN(anc),
	n_segments_all_refs = .N,
	total_length_all_refs = sum(seg_len, na.rm = TRUE)
), by = ID]

sample.total[, total_length_div_ref := total_length_all_refs / ref_divisor]
sample.total[, percentage_div_ref := total_length_div_ref / genome_length * 100]
setorder(sample.total, ID)
fwrite(sample.total, file.path(outdir, paste0(prefix, ".sample.total_div", ref_divisor, ".sum.tsv")), sep = "\t")

sample.mean_ref <- sample.by_anc[, .(
	n_ref_observed = .N,
	total_length_mean_ref = mean(total_length),
	percentage_mean_ref = mean(percentage),
	percentage_median_ref = median(percentage),
	percentage_min_ref = min(percentage),
	percentage_max_ref = max(percentage)
), by = ID]
setorder(sample.mean_ref, ID)
fwrite(sample.mean_ref, file.path(outdir, paste0(prefix, ".sample.mean_by_anc.sum.tsv")), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: Overall mean percentage table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean.by_anc <- sample.by_anc[, .(
	n_samples = .N,
	mean_percentage = mean(percentage),
	median_percentage = median(percentage),
	sd_percentage = sd(percentage),
	min_percentage = min(percentage),
	q25_percentage = as.numeric(quantile(percentage, 0.25)),
	q75_percentage = as.numeric(quantile(percentage, 0.75)),
	max_percentage = max(percentage),
	mean_total_length = mean(total_length),
	mean_n_segments = mean(n_segments)
), by = anc]
setnames(mean.by_anc, "anc", "group")

mean.total <- sample.total[, .(
	group = paste0("TOTAL_div", ref_divisor),
	n_samples = .N,
	mean_percentage = mean(percentage_div_ref),
	median_percentage = median(percentage_div_ref),
	sd_percentage = sd(percentage_div_ref),
	min_percentage = min(percentage_div_ref),
	q25_percentage = as.numeric(quantile(percentage_div_ref, 0.25)),
	q75_percentage = as.numeric(quantile(percentage_div_ref, 0.75)),
	max_percentage = max(percentage_div_ref),
	mean_total_length = mean(total_length_div_ref),
	mean_n_segments = mean(n_segments_all_refs / ref_divisor)
)]

mean.percentage <- rbindlist(list(mean.by_anc, mean.total), fill = TRUE)
setorder(mean.percentage, group)
fwrite(mean.percentage, file.path(outdir, paste0(prefix, ".mean_percentage.tsv")), sep = "\t")
print(mean.percentage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s5: Overall PNG with 5 archaic-reference histograms + total/nref histogram
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.file <- file.path(outdir, paste0(prefix, ".percentage.histograms.3x2.png"))

vals <- list()
for (a in ancs) vals[[a]] <- sample.by_anc[anc == a, percentage]
vals[[paste0("TOTAL_div", ref_divisor)]] <- sample.total$percentage_div_ref

plot.names <- c(ancs, paste0("TOTAL_div", ref_divisor))
if (length(plot.names) < 6) plot.names <- c(plot.names, rep(NA_character_, 6 - length(plot.names)))
if (length(plot.names) > 6) plot.names <- plot.names[1:6]

all.x <- unlist(vals[!is.na(names(vals))], use.names = FALSE)
xlim <- range(all.x, na.rm = TRUE)
if (!is.finite(xlim[1]) || !is.finite(xlim[2]) || xlim[1] == xlim[2]) xlim <- c(0, max(1, xlim[2], na.rm = TRUE))
breaks <- seq(xlim[1], xlim[2], length.out = 41)

png(plot.file, width = 2400, height = 1800, res = 220)
op <- par(mfrow = c(3, 2), mar = c(4.2, 4.2, 3.6, 1.2), oma = c(0, 0, 2, 0))
for (nm in plot.names) {
	if (is.na(nm) || !(nm %in% names(vals))) {
		plot.new()
		next
	}
	x <- vals[[nm]]
	main <- if (grepl("^TOTAL_div", nm)) paste0("All references / ", ref_divisor) else nm
	hist(x, breaks = breaks, xlim = xlim, main = main, xlab = "Inherited segment percentage (%)", ylab = "Samples")
	abline(v = mean(x, na.rm = TRUE), lwd = 2, lty = 2)
	mtext(paste0("mean = ", round(mean(x, na.rm = TRUE), 3), "%"), side = 3, line = 0.2, cex = 0.85)
}
mtext("IBDmix inherited segment percentage by archaic reference", outer = TRUE, cex = 1.1, font = 2)
par(op)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s6: Optional sample-info/group-specific report
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!is.na(sample_info_file) && file.exists(sample_info_file) && length(summary_by) > 0) {
	info <- fread(sample_info_file)

	if (!(sample_id_col %in% names(info))) stop("sample_id_col not found in sample_info file: ", sample_id_col)
	miss <- setdiff(summary_by, names(info))
	if (length(miss) > 0) stop("summary_by column(s) not found in sample_info file: ", paste(miss, collapse = ", "))

	ph <- unique(info[, c(sample_id_col, summary_by), with = FALSE])
	setnames(ph, sample_id_col, "ID")
	ph[, ID := as.character(ID)]
	for (x in summary_by) ph[, (x) := as.character(get(x))]
	ph <- ph[!is.na(ID) & ID != ""]
	for (x in summary_by) ph <- ph[!is.na(get(x)) & get(x) != ""]

	sba <- merge(sample.by_anc, ph, by = "ID", all.x = FALSE)
	st <- merge(sample.total, ph, by = "ID", all.x = FALSE)
	if (nrow(sba) == 0 || nrow(st) == 0) stop("No matching sample IDs between --in and --sample_info")

	group_tag <- paste(summary_by, collapse = "_")
	by_anc_cols <- c(summary_by, "anc")

	group.by_anc <- sba[, .(
		n_samples = uniqueN(ID),
		mean_percentage = mean(percentage),
		median_percentage = median(percentage),
		sd_percentage = sd(percentage),
		min_percentage = min(percentage),
		q25_percentage = as.numeric(quantile(percentage, 0.25)),
		q75_percentage = as.numeric(quantile(percentage, 0.75)),
		max_percentage = max(percentage),
		mean_total_length = mean(total_length),
		mean_n_segments = mean(n_segments)
	), by = by_anc_cols]
	setnames(group.by_anc, "anc", "reference")

	group.total <- st[, .(
		reference = paste0("TOTAL_div", ref_divisor),
		n_samples = uniqueN(ID),
		mean_percentage = mean(percentage_div_ref),
		median_percentage = median(percentage_div_ref),
		sd_percentage = sd(percentage_div_ref),
		min_percentage = min(percentage_div_ref),
		q25_percentage = as.numeric(quantile(percentage_div_ref, 0.25)),
		q75_percentage = as.numeric(quantile(percentage_div_ref, 0.75)),
		max_percentage = max(percentage_div_ref),
		mean_total_length = mean(total_length_div_ref),
		mean_n_segments = mean(n_segments_all_refs / ref_divisor)
	), by = summary_by]

	group.report <- rbindlist(list(group.by_anc, group.total), fill = TRUE)
	setcolorder(group.report, c(summary_by, "reference", setdiff(names(group.report), c(summary_by, "reference"))))
	setorderv(group.report, c(summary_by, "reference"))
	group.report.file <- file.path(outdir, paste0(prefix, ".by_", group_tag, ".mean_percentage.tsv"))
	fwrite(group.report, group.report.file, sep = "\t")

	plot.group.file <- file.path(outdir, paste0(prefix, ".by_", group_tag, ".percentage.boxplots.3x2.png"))
	make_group_label <- function(d) {
		if (length(summary_by) == 1) {
			d[, group_label := get(summary_by)]
		} else {
			d[, group_label := do.call(paste, c(.SD, sep = " | ")), .SDcols = summary_by]
		}
		d
	}

	sba.plot <- make_group_label(copy(sba))
	st.plot <- make_group_label(copy(st))

	png(plot.group.file, width = 2400, height = 1800, res = 220)
	op <- par(mfrow = c(3, 2), mar = c(6.5, 4.2, 3.6, 1.2), oma = c(0, 0, 2, 0))
	for (nm in plot.names) {
		if (is.na(nm)) {
			plot.new()
			next
		}
		if (grepl("^TOTAL_div", nm)) {
			d <- st.plot[, .(percentage = percentage_div_ref, group_label)]
			main <- paste0("All references / ", ref_divisor)
		} else {
			d <- sba.plot[anc == nm, .(percentage, group_label)]
			main <- nm
		}
		if (nrow(d) == 0) {
			plot.new()
			next
		}
		boxplot(percentage ~ group_label, data = d, las = 2, main = main, xlab = "", ylab = "Inherited segment percentage (%)")
	}
	mtext(paste0("IBDmix percentage by ", paste(summary_by, collapse = ", ")), outer = TRUE, cex = 1.1, font = 2)
	par(op)
	dev.off()

	message("Sample-info-specific report:")
	message("  ", group.report.file)
	message("  ", plot.group.file)
} else if (!is.na(sample_info_file) && !file.exists(sample_info_file)) {
	message("Sample info file not found, skip group-specific report: ", sample_info_file)
} else if (!is.na(sample_info_file) && length(summary_by) == 0) {
	message("--sample_info was provided, but --summary_by was not provided. Skip group-specific report.")
}

message("Done:")
message("  ", file.path(outdir, paste0(prefix, ".input.check.tsv")))
message("  ", file.path(outdir, paste0(prefix, ".sample.by_anc.sum.tsv")))
message("  ", file.path(outdir, paste0(prefix, ".sample.total_div", ref_divisor, ".sum.tsv")))
message("  ", file.path(outdir, paste0(prefix, ".sample.mean_by_anc.sum.tsv")))
message("  ", file.path(outdir, paste0(prefix, ".mean_percentage.tsv")))
message("  ", plot.file)
