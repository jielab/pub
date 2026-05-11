#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 IBDmix report: genome-wide burden + optional locus-level summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic genome-wide window burden:
# Rscript ibdmix_report.R --in len1000.lod5.segments.tsv.gz --outdir report
#
# With sample groups:
# Rscript ibdmix_report.R --in len1000.lod5.segments.tsv.gz --outdir report \
#   --sample_info /mnt/d/files/1kg.v3.sample.txt --stats_by_group super_pop,gender
#
# With focal loci:
# Rscript ibdmix_report.R --in len1000.lod5.segments.tsv.gz --outdir report \
#   --sample_info /mnt/d/files/1kg.v3.sample.txt --stats_by_group super_pop,gender \
#   --stats_for_locus '3:45859651-45909024;9:136130562-136150630'

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NA_character_) {
	i <- match(key, args)
	if (is.na(i) || i == length(args)) return(default)
	args[i + 1]
}

deprecated_args <- intersect(args, c("--summary_by", "--stats_by"))
if (length(deprecated_args) > 0) {
	stop("Deprecated option(s) not supported: ", paste(deprecated_args, collapse = ", "),
		". Please use only --stats_by_group and --stats_for_locus.")
}
valid_args <- c("--in", "--outdir", "--sample_info", "--sample_keep", "--sample_id_col", "--stats_by_group",
	"--stats_for_locus", "--window_size", "--locus_min_cover", "--top_n_windows")
unknown_args <- args[grepl("^--", args) & !(args %in% valid_args)]
if (length(unknown_args) > 0) stop("Unknown option(s): ", paste(unknown_args, collapse = ", "))

infile <- get_arg("--in")
outdir <- get_arg("--outdir", ".")
sample_info_file <- get_arg("--sample_info", NA_character_)
sample_keep_file <- get_arg("--sample_keep", NA_character_)
sample_id_col <- get_arg("--sample_id_col", "sample")
stats_by_group_arg <- get_arg("--stats_by_group", NA_character_)
stats_for_locus_arg <- get_arg("--stats_for_locus", NA_character_)
window_size <- as.integer(get_arg("--window_size", "1000000"))
locus_min_cover <- as.numeric(get_arg("--locus_min_cover", "0.5"))
top_n_windows <- as.integer(get_arg("--top_n_windows", "100"))

if (is.na(infile)) stop("Please provide --in input.segments.tsv.gz")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (is.na(window_size) || window_size <= 0) stop("Bad --window_size: ", window_size)
if (is.na(locus_min_cover) || locus_min_cover <= 0 || locus_min_cover > 1) stop("Bad --locus_min_cover: ", locus_min_cover)

stats_by_group <- character(0)
if (!is.na(stats_by_group_arg) && nzchar(stats_by_group_arg)) {
	stats_by_group <- trimws(unlist(strsplit(stats_by_group_arg, ",", fixed = TRUE)))
	stats_by_group <- stats_by_group[nzchar(stats_by_group)]
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

input_prefix <- basename(infile)
input_prefix <- sub("\\.segments\\.tsv\\.gz$", "", input_prefix)
input_prefix <- sub("\\.tsv\\.gz$", "", input_prefix)
input_prefix <- sub("\\.gz$", "", input_prefix)
out_file <- function(suffix) file.path(outdir, suffix)
out_dir <- function(dirname) file.path(outdir, dirname)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
write_xlsx_sheets <- function(sheets, path) {
	if (!requireNamespace("writexl", quietly = TRUE)) {
		stop("R package 'writexl' is required to write Excel files: ", path)
	}
	writexl::write_xlsx(sheets, path)
}
unlink(file.path(outdir, c(
	"check.tsv",
	"check.csv",
	"burden.tsv.gz",
	"burden.csv.gz",
	"burden.csv",
	"burden_by_group.tsv.gz",
	"burden_by_group.csv.gz",
	"burden_by_group.csv",
	"top_windows.tsv",
	"top_windows.csv",
	"top.csv",
	"burden.png",
	"all_gw.csv",
	"all_gw.png",
	"all.xlsx",
	"group_burden.tsv",
	"group_burden.csv",
	"group_burden.png",
	"all_percent.csv",
	"all_percent.png",
	"sample_burden.tsv",
	"sample_burden.csv",
	"sample_boxplot.png",
	"all_length.csv",
	"all_length.png",
	"locus_matches.tsv",
	"locus_matches.csv",
	"locus_summary.tsv",
	"locus_summary.csv",
	"locus_summary.png",
	"selected_overall.csv",
	"selected_overall.png",
	"selected.png",
	"selected.xlsx",
	"locus_group.tsv",
	"locus_group.csv",
	"locus_group.png",
	"selected_percent.csv",
	"selected_percent.png",
	"input.check.tsv",
	"group.by_anc.segment_percent.tsv",
	"sample.by_anc.segment_burden.tsv",
	"group.by_anc.segment_percent.png",
	"sample.by_anc.segment_mb.boxplot.png",
	"locus.sample_matches.tsv",
	"locus.by_anc.summary.tsv",
	"locus.by_anc.summary.png",
	"locus.by_super_pop_gender.summary.tsv",
	"locus.by_group_anc.summary.png",
	"locus.by_group_anc.summary.tsv",
	"genomewide.window1000000.burden.tsv.gz",
	"genomewide.window1000000.by_super_pop_gender.burden.tsv.gz",
	"genomewide.window1000000.top_windows.tsv",
	"genomewide.window1000000.burden.3x2.png"
)), force = TRUE)
unlink(c(out_dir("sample.by_anc.segment_mb.boxplot"), out_dir("locus.by_group_anc.summary"), out_dir(paste0("bedGraph.window", window_size))), recursive = TRUE, force = TRUE)

message("Input:  ", infile)
message("Outdir: ", outdir)
message("Input prefix ignored for output names: ", input_prefix)
message("Window size: ", window_size)
if (!is.na(sample_info_file)) message("Sample info: ", sample_info_file)
if (!is.na(sample_keep_file)) message("Sample keep: ", sample_keep_file)
if (length(stats_by_group) > 0) message("Stats by group: ", paste(stats_by_group, collapse = ","))
if (!is.na(stats_for_locus_arg) && nzchar(stats_for_locus_arg)) message("Stats for locus: ", stats_for_locus_arg)

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

if (!("sample_set" %in% names(dat))) dat[, sample_set := "sample_keep"]

dat <- dat[!is.na(ID) & ID != "" & ID != "ID" & ID != "DO"]
dat[, ID := as.character(ID)]
dat[, chrom := as.character(chrom)]
dat[, chrom := sub("^chr", "", chrom, ignore.case = TRUE)]
dat[, start := as.integer(start)]
dat[, end := as.integer(end)]
dat[, seg_len := as.integer(length)]
dat[, slod := as.numeric(slod)]
dat[, sites := as.numeric(sites)]
dat[, positive_lods := as.numeric(positive_lods)]
dat[, negative_lods := as.numeric(negative_lods)]
dat[, anc := as.character(anc)]
dat[, sample_set := as.character(sample_set)]
dat <- dat[!is.na(start) & !is.na(end) & !is.na(seg_len) & !is.na(anc) & anc != ""]
dat <- dat[end >= start]
if (nrow(dat) == 0) stop("No usable segment rows after cleanup")

n0 <- nrow(dat)
dat <- unique(dat, by = c("ID", "anc", "chrom", "start", "end"))
n_dup_removed <- n0 - nrow(dat)

seg_ids <- sort(unique(dat$ID))
ids <- seg_ids
if (!is.na(sample_keep_file) && file.exists(sample_keep_file)) {
	ids0 <- fread(sample_keep_file, header = FALSE, showProgress = FALSE)[[1]]
	ids0 <- sort(unique(as.character(ids0[!is.na(ids0) & nzchar(ids0)])))
	if (length(ids0) == 0) stop("sample_keep has no sample IDs: ", sample_keep_file)
	dat <- dat[ID %in% ids0]
	ids <- ids0
	seg_ids <- sort(unique(dat$ID))
} else if (!is.na(sample_keep_file)) {
	stop("sample_keep file not found: ", sample_keep_file)
}
if (nrow(dat) == 0) stop("No segment rows remain after applying sample_keep: ", sample_keep_file)
ancs <- sort(unique(dat$anc))
chrom_order <- c(as.character(1:22), "X", "Y", "MT", "M")
chrom_rank <- function(x) {
	r <- match(x, chrom_order)
	other <- match(x, sort(unique(x)))
	ifelse(is.na(r), length(chrom_order) + other, r)
}

message("Rows after cleanup: ", nrow(dat))
message("Samples in segment file: ", length(seg_ids))
message("Samples in denominator: ", length(ids))
message("Archaic references: ", paste(ancs, collapse = ", "))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Optional sample-info table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ph <- NULL
if (!is.na(sample_info_file) && file.exists(sample_info_file)) {
	info <- fread(sample_info_file)
	if (!(sample_id_col %in% names(info))) stop("sample_id_col not found in sample_info file: ", sample_id_col)
	if (length(stats_by_group) > 0) {
		miss <- setdiff(stats_by_group, names(info))
		if (length(miss) > 0) stop("stats_by_group column(s) not found in sample_info file: ", paste(miss, collapse = ", "))
	}
		sex_col <- intersect(c("sex", "gender", "Sex", "Gender", "SEX", "GENDER"), names(info))[1]
		keep_cols <- unique(c(sample_id_col, stats_by_group, sex_col))
		ph <- unique(info[, keep_cols, with = FALSE])
		setnames(ph, sample_id_col, "ID")
		if (!is.na(sex_col)) setnames(ph, sex_col, ".sex_for_denom")
		ph[, ID := as.character(ID)]
		ph <- ph[ID %in% ids]
		for (x in stats_by_group) ph[, (x) := as.character(get(x))]
		if (!(".sex_for_denom" %in% names(ph))) ph[, .sex_for_denom := NA_character_]
		ph[, .sex_for_denom := tolower(trimws(as.character(.sex_for_denom)))]
		for (x in stats_by_group) ph <- ph[!is.na(get(x)) & get(x) != ""]
	if (nrow(ph) == 0) stop("No matching sample IDs between --in and --sample_info")
	message("Matched sample_info IDs: ", uniqueN(ph$ID))
} else if (!is.na(sample_info_file)) {
	message("Sample info file not found, skip group-specific outputs: ", sample_info_file)
	stats_by_group <- character(0)
}

if (is.null(ph)) {
	ph <- data.table(ID = ids)
	ph[, .sex_for_denom := NA_character_]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Input check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input.check <- data.table(
	input = infile,
	n_rows = nrow(dat),
	n_samples = length(ids),
	n_samples_with_segments = length(seg_ids),
	n_archaic_references = length(ancs),
	archaic_references = paste(ancs, collapse = ","),
	n_DO_lines_skipped = n_do,
	n_exact_duplicate_segments_removed = n_dup_removed,
	window_size = window_size,
	sample_info = ifelse(is.na(sample_info_file), "", sample_info_file),
	sample_keep = ifelse(is.na(sample_keep_file), "", sample_keep_file),
	sample_id_col = sample_id_col,
	stats_by_group = paste(stats_by_group, collapse = ","),
	stats_for_locus = ifelse(is.na(stats_for_locus_arg), "", stats_for_locus_arg),
	locus_min_cover = locus_min_cover
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Genome-wide window burden table + bedGraph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_windows <- function(d, w) {
	chr.max <- d[, .(chr_end = max(end, na.rm = TRUE)), by = chrom]
	chr.max[, nwin := ceiling(chr_end / w)]
	wins <- chr.max[, .(win_id = seq_len(nwin)), by = .(chrom, chr_end)]
	wins[, start := (win_id - 1L) * w + 1L]
	wins[, end := pmin(win_id * w, chr_end)]
	wins[, window := paste0(chrom, ":", start, "-", end)]
	wins[, chr_rank := chrom_rank(chrom)]
	setorder(wins, chr_rank, start)
	wins[, genome_pos := cumsum(as.numeric(c(1L, head(end - start + 1L, -1L))))]
	wins[]
}

windows <- make_windows(dat, window_size)
genome_bp <- windows[, sum(as.numeric(end - start + 1L), na.rm = TRUE)]
chr_lengths_grch37 <- c(
	"1" = 249250621, "2" = 243199373, "3" = 198022430, "4" = 191154276,
	"5" = 180915260, "6" = 171115067, "7" = 159138663, "8" = 146364022,
	"9" = 141213431, "10" = 135534747, "11" = 135006516, "12" = 133851895,
	"13" = 115169878, "14" = 107349540, "15" = 102531392, "16" = 90354753,
	"17" = 81195210, "18" = 78077248, "19" = 59128983, "20" = 63025520,
	"21" = 48129895, "22" = 51304566, "X" = 155270560, "23" = 155270560,
	"Y" = 59373566, "MT" = 16569, "M" = 16569
)
denom_chroms <- sort(unique(windows$chrom))
known_denom_chroms <- denom_chroms[denom_chroms %in% names(chr_lengths_grch37)]
unknown_denom_chroms <- setdiff(denom_chroms, known_denom_chroms)
haploid_denom_bp <- sum(as.numeric(chr_lengths_grch37[known_denom_chroms]), na.rm = TRUE)
if (length(unknown_denom_chroms) > 0) {
	haploid_denom_bp <- haploid_denom_bp + windows[chrom %in% unknown_denom_chroms, sum(as.numeric(end - start + 1L), na.rm = TRUE)]
}
if (haploid_denom_bp <= 0) haploid_denom_bp <- genome_bp
diploid_genome_bp <- 6e9
sample.denom <- unique(ph[, .(ID, .sex_for_denom)])
sample.denom <- merge(data.table(ID = ids), sample.denom, by = "ID", all.x = TRUE)
sample.denom[is.na(.sex_for_denom), .sex_for_denom := ""]
sample.denom[, denominator_bp := diploid_genome_bp]
sample.denom[, denominator_gb := denominator_bp / 1e9]
wins.map <- windows[, .(chrom, win_id, window, win_start = start, win_end = end, chr_rank, genome_pos)]
setkey(wins.map, chrom, win_id)

# Expand each inherited segment to the window(s) it overlaps. With the default
# 1 Mb window, most segments contribute to only one row here.
dseg <- dat[, .(
	seg_id = .I, ID, anc, chrom,
	start_win = ((start - 1L) %/% window_size) + 1L,
	end_win = ((end - 1L) %/% window_size) + 1L
)]
ov <- dseg[, .(win_id = seq.int(start_win[1], end_win[1])), by = .(seg_id, ID, anc, chrom)]
ov[, seg_id := NULL]
ov <- unique(ov)
ov <- wins.map[ov, on = .(chrom, win_id), nomatch = 0L]
ov <- ov[, .(ID, anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)]

den_all <- length(ids)
win.by_anc <- ov[, .(n_samples_with_segment = uniqueN(ID)), by = .(anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)]
win.by_anc[, n_samples_total := den_all]
win.by_anc[, burden_fraction := n_samples_with_segment / n_samples_total]
win.by_anc[, burden_percent := burden_fraction * 100]

win.any <- unique(ov[, .(ID, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)])
win.total <- win.any[, .(n_samples_with_segment = uniqueN(ID)), by = .(chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)]
win.total[, anc := "ANY_REF"]
win.total[, n_samples_total := den_all]
win.total[, burden_fraction := n_samples_with_segment / n_samples_total]
win.total[, burden_percent := burden_fraction * 100]

win.burden <- rbindlist(list(win.by_anc, win.total), fill = TRUE)
setcolorder(win.burden, c("anc", "chrom", "win_id", "window", "win_start", "win_end", "n_samples_with_segment", "n_samples_total", "burden_fraction", "burden_percent", "chr_rank", "genome_pos"))
setorder(win.burden, anc, chr_rank, win_start)

win.file <- out_file("all_gw.csv")
fwrite(win.burden, win.file)

bedgraph.dir <- out_dir("bedgraph_tmp")
dir.create(bedgraph.dir, recursive = TRUE, showWarnings = FALSE)
for (a in sort(unique(win.burden$anc))) {
	bg <- win.burden[anc == a, .(
		chrom = ifelse(chrom %in% c("X", "Y", "M", "MT"), paste0("chr", chrom), paste0("chr", chrom)),
		start0 = win_start - 1L,
		end = win_end,
		burden_percent,
		chr_rank
	)]
	setorder(bg, chr_rank, start0)
	bg[, chr_rank := NULL]
	fwrite(bg, file.path(bedgraph.dir, paste0(a, ".window", window_size, ".burden_percent.bedGraph")), sep = "\t", col.names = FALSE)
}

top.file <- out_file("top.csv")
top.windows <- win.burden[, head(.SD[order(-burden_percent, -n_samples_with_segment)], top_n_windows), by = anc]
fwrite(top.windows, top.file)

# Optional group-specific window burden table, when --sample_info and --stats_by_group are provided.
if (length(stats_by_group) > 0) {
	ov.g <- merge(ov[, .(ID, anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)], ph, by = "ID", all.x = FALSE)
	if (nrow(ov.g) > 0) {
		den.group <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
		win.g <- ov.g[, .(n_samples_with_segment = uniqueN(ID)), by = c(stats_by_group, "anc", "chrom", "win_id", "window", "win_start", "win_end", "chr_rank", "genome_pos")]
		win.g <- merge(win.g, den.group, by = stats_by_group, all.x = TRUE)
		win.g[, burden_fraction := n_samples_with_segment / n_samples_total]
		win.g[, burden_percent := burden_fraction * 100]
		setorderv(win.g, c(stats_by_group, "anc", "chr_rank", "win_start"))
		group_tag <- paste(stats_by_group, collapse = "_")
	}
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seg.sample <- dat[, .(
	n_segments = .N,
	total_bp = sum(seg_len, na.rm = TRUE),
	total_mb = sum(seg_len, na.rm = TRUE) / 1e6,
	mean_length = mean(seg_len, na.rm = TRUE),
	mean_slod = mean(slod, na.rm = TRUE),
	median_slod = median(slod, na.rm = TRUE)
), by = .(ID, anc)]

all.sample.anc <- CJ(ID = ids, anc = ancs, unique = TRUE)
seg.sample <- merge(all.sample.anc, seg.sample, by = c("ID", "anc"), all.x = TRUE)
for (cc in c("n_segments", "total_bp", "total_mb")) seg.sample[is.na(get(cc)), (cc) := 0]
seg.sample <- merge(seg.sample, ph, by = "ID", all.x = TRUE)
seg.sample <- merge(seg.sample, sample.denom[, .(ID, denominator_bp, denominator_gb)], by = "ID", all.x = TRUE)
seg.sample[, segment_percent := 100 * total_bp / denominator_bp]
fwrite(seg.sample, out_file("all_length.csv"))

if (length(stats_by_group) > 0) {
	group.cols <- stats_by_group
	seg.sample.group <- copy(seg.sample)
	for (g in group.cols) seg.sample.group <- seg.sample.group[!is.na(get(g)) & get(g) != ""]
	if (nrow(seg.sample.group) == 0) stop("No samples with complete stats_by_group columns: ", paste(group.cols, collapse = ","))
	group.percent <- seg.sample.group[, .(
		n_samples_total = uniqueN(ID),
		n_samples_with_segment = uniqueN(ID[n_segments > 0]),
		mean_total_bp_per_sample = mean(total_bp, na.rm = TRUE),
		mean_segments_per_sample = as.numeric(mean(n_segments, na.rm = TRUE)),
		median_segments_per_sample = as.numeric(median(n_segments, na.rm = TRUE)),
		mean_total_mb_per_sample = mean(total_mb, na.rm = TRUE),
		median_total_mb_per_sample = median(total_mb, na.rm = TRUE),
		mean_denominator_gb = mean(denominator_gb, na.rm = TRUE),
		segment_percent = mean(segment_percent, na.rm = TRUE)
	), by = c(group.cols, "anc")]
	group.percent[, sample_prevalence_percent := 100 * n_samples_with_segment / n_samples_total]
	setorderv(group.percent, c(group.cols, "anc"))
	fwrite(group.percent, out_file("all_percent.csv"))

	png(out_file("all_percent.png"), width = if (length(group.cols) == 1) 1700 else 1700 * length(group.cols), height = 1250, res = 220)
	op <- par(mfrow = c(1, length(group.cols)), mar = c(5, 5, 3.2, 1.2), oma = c(0, 0, 4, 0))
	for (g in group.cols) {
		dg <- seg.sample.group[, .(
			n_samples_total = uniqueN(ID),
			n_samples_with_segment = uniqueN(ID[n_segments > 0]),
			mean_total_bp_per_sample = mean(total_bp, na.rm = TRUE),
			segment_percent = mean(segment_percent, na.rm = TRUE)
		), by = c(g, "anc")]
		dg[, sample_prevalence_percent := 100 * n_samples_with_segment / n_samples_total]
		dg[, group_value := get(g)]
		tab <- dcast(dg, group_value ~ anc, value.var = "segment_percent", fill = 0)
		mat <- as.matrix(tab[, setdiff(names(tab), "group_value"), with = FALSE])
		rownames(mat) <- tab$group_value
		barplot(mat, beside = TRUE, las = 1, ylab = "Archaic segment burden (% of diploid genome)",
			xlab = "Archaic reference", main = paste("Grouped by", g),
			legend.text = rownames(mat), args.legend = list(x = "topright", cex = 0.75, title = g))
	}
	mtext("Genome-wide archaic segment burden", outer = TRUE, line = 2.2, font = 2)
	mtext(paste0("Percentage = total IBDmix segment length / diploid genome size (",
		round(diploid_genome_bp / 1e9, 2), " Gb)."), outer = TRUE, line = 0.9, cex = 0.75)
	par(op)
	dev.off()

	}

all_sheets <- list(all_gw = win.burden, all_length = seg.sample)
if (exists("group.percent")) all_sheets$all_percent <- group.percent
write_xlsx_sheets(all_sheets, out_file("all.xlsx"))
unlink(out_file(c("all_gw.csv", "all_length.csv", "all_percent.csv")), force = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: Genome-wide burden plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.file <- out_file("all_gw.png")
plot.names <- c(ancs, "ANY_REF")
if (length(plot.names) < 6) plot.names <- c(plot.names, rep(NA_character_, 6 - length(plot.names)))
if (length(plot.names) > 6) plot.names <- plot.names[1:6]

chr.mid <- windows[, .(mid = mean(genome_pos), chr_rank = unique(chr_rank)), by = chrom]
setorder(chr.mid, chr_rank)

png(plot.file, width = 2600, height = 1800, res = 220)
op <- par(mfrow = c(3, 2), mar = c(4.2, 4.4, 3.4, 1.2), oma = c(0, 0, 2, 0))
for (a in plot.names) {
	if (is.na(a)) {
		plot.new()
		next
	}
	d <- win.burden[anc == a]
	if (nrow(d) == 0) {
		plot.new()
		next
	}
	main <- ifelse(a == "ANY_REF", "Any archaic reference", a)
	plot(d$genome_pos, d$burden_percent, pch = 20, cex = 0.35, xaxt = "n",
		xlab = "Chromosome", ylab = "Samples with segment (%)", main = main)
	axis(1, at = chr.mid$mid, labels = chr.mid$chrom, cex.axis = 0.55)
	q99 <- as.numeric(quantile(d$burden_percent, 0.99, na.rm = TRUE))
	abline(h = q99, lty = 2, lwd = 1)
	mtext(paste0("99th pct = ", round(q99, 2), "%"), side = 3, line = 0.15, cex = 0.75)
}
mtext(paste0("Genome-wide IBDmix window burden, window = ", format(window_size, big.mark = ","), " bp"), outer = TRUE, cex = 1.1, font = 2)
par(op)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s5: Optional focal-locus sample table and summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parse_loci <- function(x) {
	if (is.na(x) || !nzchar(x)) return(data.table())
	x <- gsub("–|—|−", "-", x)
	parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
	parts <- parts[nzchar(parts)]
	out <- rbindlist(lapply(seq_along(parts), function(i) {
		p <- parts[i]
		m <- regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", p)
		r <- regmatches(p, m)[[1]]
		if (length(r) != 4) stop("Bad locus format: ", p, ". Expected chr:start-end;chr:start-end")
		data.table(
			locus_id = paste0("locus", i),
			locus = p,
			chrom = sub("^chr", "", r[2], ignore.case = TRUE),
			locus_start = as.integer(gsub(",", "", r[3])),
			locus_end = as.integer(gsub(",", "", r[4]))
		)
	}), fill = TRUE)
	out[locus_end < locus_start, c("locus_start", "locus_end") := .(locus_end, locus_start)]
	out[, locus_len := locus_end - locus_start + 1L]
	out
}

loci <- parse_loci(stats_for_locus_arg)
if (nrow(loci) > 0) {
	loci[, chrom := as.character(chrom)]
	lo <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
		l <- loci[i]
		z <- dat[chrom == l$chrom & end >= l$locus_start & start <= l$locus_end,
			.(ID, anc, chrom, seg_start = start, seg_end = end, seg_len, slod, sites)]
		if (nrow(z) == 0) return(data.table())
		z[, `:=`(
			locus_id = l$locus_id,
			locus = l$locus,
			locus_start = l$locus_start,
			locus_end = l$locus_end,
			locus_len = l$locus_len
		)]
		z
	}), fill = TRUE)
	if (nrow(lo) > 0) {
		lo[, overlap_start := pmax(seg_start, locus_start)]
		lo[, overlap_end := pmin(seg_end, locus_end)]
		lo[, overlap_bp := pmax(0L, overlap_end - overlap_start + 1L)]
		lo[, locus_cover := overlap_bp / locus_len]
		lo <- lo[locus_cover >= locus_min_cover]
	}

	if (nrow(lo) == 0) {
		message("No locus matches at locus_min_cover >= ", locus_min_cover)
		matches2 <- data.table(
			locus_id = character(), locus = character(), chrom = character(),
			locus_start = integer(), locus_end = integer(), locus_len = integer(),
			ID = character(), anc = character(), n_matching_segments = integer(),
			max_locus_cover = numeric(), total_overlap_bp = integer(),
			best_segment_start = integer(), best_segment_end = integer(),
			best_segment_length = integer(), best_slod = numeric(), best_sites = numeric()
		)
	} else {
		matches <- lo[, .(
			n_matching_segments = .N,
			max_locus_cover = max(locus_cover, na.rm = TRUE),
			total_overlap_bp = sum(overlap_bp, na.rm = TRUE),
			best_segment_start = seg_start[which.max(locus_cover)],
			best_segment_end = seg_end[which.max(locus_cover)],
			best_segment_length = seg_len[which.max(locus_cover)],
			best_slod = slod[which.max(locus_cover)],
			best_sites = sites[which.max(locus_cover)]
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, ID, anc)]
		setorder(matches, locus_id, anc, ID)
		matches2 <- matches

		den.locus <- uniqueN(ph$ID)
		locus.by_anc <- matches[, .(n_samples_match = uniqueN(ID)), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc)]
		locus.by_anc[, n_samples_total := den.locus]
		locus.by_anc[, match_fraction := n_samples_match / n_samples_total]
		locus.by_anc[, match_percent := match_fraction * 100]
		setorder(locus.by_anc, locus_id, anc)
		fwrite(locus.by_anc, out_file("selected_overall.csv"))

		if (length(stats_by_group) > 0) {
			mg <- merge(matches, ph, by = "ID", all.x = FALSE)
			den.g <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
			locus.g <- mg[, .(n_samples_match = uniqueN(ID)), by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
			locus.g <- merge(locus.g, den.g, by = stats_by_group, all.x = TRUE)
			locus.g[, match_fraction := n_samples_match / n_samples_total]
			locus.g[, match_percent := match_fraction * 100]
			group_tag <- paste(stats_by_group, collapse = "_")
			setorderv(locus.g, c(stats_by_group, "locus_id", "anc"))
			fwrite(locus.g, out_file("selected_percent.csv"))
		}
	}


	if (nrow(matches2) > 0) matches2[, chrom := as.character(chrom)]
	locus.grid <- merge(copy(loci)[, tmp := 1], data.table(anc = ancs, tmp = 1), by = "tmp", allow.cartesian = TRUE)[, tmp := NULL]
	if (nrow(matches2) > 0) {
		locus.by_anc <- matches2[, .(
			n_samples_match = uniqueN(ID),
			n_matching_segments = as.integer(sum(n_matching_segments, na.rm = TRUE)),
			median_max_locus_cover = as.numeric(median(max_locus_cover, na.rm = TRUE)),
			median_best_slod = as.numeric(median(best_slod, na.rm = TRUE)),
			median_best_segment_length = as.numeric(median(best_segment_length, na.rm = TRUE))
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc)]
	} else {
		locus.by_anc <- locus.grid[0, .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc,
			n_samples_match = integer(), n_matching_segments = integer(), median_max_locus_cover = numeric(),
			median_best_slod = numeric(), median_best_segment_length = numeric())]
	}
	locus.by_anc <- merge(locus.grid, locus.by_anc,
		by = c("locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc"), all.x = TRUE)
	for (cc in c("n_samples_match", "n_matching_segments")) locus.by_anc[is.na(get(cc)), (cc) := 0]
	locus.by_anc[, n_samples_total := uniqueN(ph$ID)]
	locus.by_anc[, match_fraction := n_samples_match / n_samples_total]
	locus.by_anc[, match_percent := 100 * match_fraction]
	setorder(locus.by_anc, locus_id, anc)
	fwrite(locus.by_anc, out_file("selected_overall.csv"))

	if (length(stats_by_group) > 0) {
		group_tag <- paste(stats_by_group, collapse = "_")
		den.g <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
		if (nrow(matches2) > 0) {
			mg <- merge(matches2, ph, by = "ID", all.x = FALSE)
			locus.g <- mg[, .(n_samples_match = uniqueN(ID)), by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
		} else {
			locus.g <- data.table()
		}
		group.levels <- unique(ph[, stats_by_group, with = FALSE])
		group.levels[, tmp := 1]
		locus.grid.g <- merge(merge(copy(loci)[, tmp := 1], data.table(anc = ancs, tmp = 1), by = "tmp", allow.cartesian = TRUE),
			group.levels, by = "tmp", allow.cartesian = TRUE)[, tmp := NULL]
		locus.g <- merge(locus.grid.g, locus.g,
			by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc"), all.x = TRUE)
		locus.g[is.na(n_samples_match), n_samples_match := 0]
		locus.g <- merge(locus.g, den.g, by = stats_by_group, all.x = TRUE)
		locus.g[, match_fraction := n_samples_match / n_samples_total]
		locus.g[, match_percent := 100 * match_fraction]
		setorderv(locus.g, c(stats_by_group, "locus_id", "anc"))
		fwrite(locus.g, out_file("selected_percent.csv"))
		}

	selected_sheets <- list(selected_overall = locus.by_anc)
	if (exists("locus.g")) selected_sheets$selected_percent <- locus.g
	selected_sheets$locus_matches <- matches2
	write_xlsx_sheets(selected_sheets, out_file("selected.xlsx"))

	plot_selected_panel <- function(d, main_title, show_legend = FALSE) {
		tab <- dcast(d, anc ~ locus, value.var = "match_percent", fill = 0)
		mat <- as.matrix(tab[, setdiff(names(tab), "anc"), with = FALSE])
		rownames(mat) <- tab$anc
		ymax <- max(mat, na.rm = TRUE) * 1.2
		if (!is.finite(ymax) || ymax <= 0) ymax <- 1
		barplot(mat, beside = TRUE, las = 1, cex.names = 0.82,
			ylab = "Samples carrying selected-locus segment (%)",
			xlab = "Selected locus",
			main = main_title,
			ylim = c(0, ymax),
			legend.text = if (show_legend) rownames(mat) else NULL,
			args.legend = list(x = "topright", cex = 0.68, title = "Archaic reference"))
	}

	panel_data <- list(list(title = "Overall", data = locus.by_anc))
	if (exists("locus.g")) {
		locus.g[, group_label := do.call(paste, c(.SD, sep = " / ")), .SDcols = stats_by_group]
		group_title <- paste(stats_by_group, collapse = " / ")
		for (gv in sort(unique(locus.g$group_label))) {
			panel_data[[length(panel_data) + 1L]] <- list(
				title = paste0(group_title, " = ", gv),
				data = locus.g[group_label == gv]
			)
		}
	}
	png(out_file("selected.png"), width = 1850, height = 780 * length(panel_data), res = 220)
	op <- par(mfrow = c(length(panel_data), 1), mar = c(5, 7, 3.5, 1.2), oma = c(0, 0, 3.8, 0))
	for (i in seq_along(panel_data)) {
		plot_selected_panel(panel_data[[i]]$data, panel_data[[i]]$title, show_legend = i == 1L)
	}
	mtext("Selected loci: carrier percentage", outer = TRUE, line = 2.1, font = 2)
	mtext(paste0("A sample is counted once per archaic reference when any IBDmix segment covers at least ",
		round(100 * locus_min_cover), "% of the selected locus."), outer = TRUE, line = 0.7, cex = 0.75)
	par(op)
	dev.off()
	unlink(out_file(c("selected_overall.csv", "selected_percent.csv")), force = TRUE)
	}

message("Done:")
unlink(bedgraph.dir, recursive = TRUE, force = TRUE)
message("  ", out_file("all.xlsx"))
message("  ", top.file)
message("  ", plot.file)
if (nrow(loci) > 0) {
	message("  ", out_file("selected.xlsx"))
	message("  ", out_file("selected.png"))
}
