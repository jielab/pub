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
valid_args <- c("--in", "--outdir", "--sample_info", "--sample_id_col", "--stats_by_group",
	"--stats_for_locus", "--window_size", "--locus_min_cover", "--top_n_windows")
unknown_args <- args[grepl("^--", args) & !(args %in% valid_args)]
if (length(unknown_args) > 0) stop("Unknown option(s): ", paste(unknown_args, collapse = ", "))

infile <- get_arg("--in")
outdir <- get_arg("--outdir", ".")
sample_info_file <- get_arg("--sample_info", NA_character_)
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

message("Input:  ", infile)
message("Outdir: ", outdir)
message("Input prefix ignored for output names: ", input_prefix)
message("Window size: ", window_size)
if (!is.na(sample_info_file)) message("Sample info: ", sample_info_file)
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

n0 <- nrow(dat)
dat <- unique(dat, by = c("ID", "anc", "chrom", "start", "end"))
n_dup_removed <- n0 - nrow(dat)

ancs <- sort(unique(dat$anc))
ids <- sort(unique(dat$ID))
chrom_order <- c(as.character(1:22), "X", "Y", "MT", "M")
chrom_rank <- function(x) {
	r <- match(x, chrom_order)
	other <- match(x, sort(unique(x)))
	ifelse(is.na(r), length(chrom_order) + other, r)
}

message("Rows after cleanup: ", nrow(dat))
message("Samples in segment file: ", length(ids))
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
	keep_cols <- unique(c(sample_id_col, stats_by_group))
	ph <- unique(info[, keep_cols, with = FALSE])
	setnames(ph, sample_id_col, "ID")
	ph[, ID := as.character(ID)]
	ph <- ph[ID %in% ids]
	for (x in stats_by_group) ph[, (x) := as.character(get(x))]
	for (x in stats_by_group) ph <- ph[!is.na(get(x)) & get(x) != ""]
	if (nrow(ph) == 0) stop("No matching sample IDs between --in and --sample_info")
	message("Matched sample_info IDs: ", uniqueN(ph$ID))
} else if (!is.na(sample_info_file)) {
	message("Sample info file not found, skip group-specific outputs: ", sample_info_file)
}

if (is.null(ph)) {
	ph <- data.table(ID = ids)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Input check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input.check <- data.table(
	input = infile,
	n_rows = nrow(dat),
	n_samples = length(ids),
	n_archaic_references = length(ancs),
	archaic_references = paste(ancs, collapse = ","),
	n_DO_lines_skipped = n_do,
	n_exact_duplicate_segments_removed = n_dup_removed,
	window_size = window_size,
	sample_info = ifelse(is.na(sample_info_file), "", sample_info_file),
	sample_id_col = sample_id_col,
	stats_by_group = paste(stats_by_group, collapse = ","),
	stats_for_locus = ifelse(is.na(stats_for_locus_arg), "", stats_for_locus_arg),
	locus_min_cover = locus_min_cover
)
fwrite(input.check, out_file("input.check.tsv"), sep = "\t")

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
	wins[, genome_pos := cumsum(c(1L, head(end - start + 1L, -1L)))]
	wins[]
}

windows <- make_windows(dat, window_size)
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

den_all <- uniqueN(ph$ID)
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

win.file <- out_file(paste0("genomewide.window", window_size, ".burden.tsv.gz"))
fwrite(win.burden, win.file, sep = "\t")

bedgraph.dir <- out_dir(paste0("bedGraph.window", window_size))
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

top.file <- out_file(paste0("genomewide.window", window_size, ".top_windows.tsv"))
top.windows <- win.burden[, head(.SD[order(-burden_percent, -n_samples_with_segment)], top_n_windows), by = anc]
fwrite(top.windows, top.file, sep = "\t")

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
		fwrite(win.g, out_file(paste0("genomewide.window", window_size, ".by_", group_tag, ".burden.tsv.gz")), sep = "\t")
	}
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: Genome-wide burden plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.file <- out_file(paste0("genomewide.window", window_size, ".burden.3x2.png"))
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
		empty <- loci[, .(locus_id, locus, chrom, locus_start, locus_end, locus_len)]
		fwrite(empty[0], out_file("locus.sample_matches.tsv"), sep = "\t")
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
		fwrite(matches, out_file("locus.sample_matches.tsv"), sep = "\t")

		den.locus <- uniqueN(ph$ID)
		locus.by_anc <- matches[, .(n_samples_match = uniqueN(ID)), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc)]
		locus.by_anc[, n_samples_total := den.locus]
		locus.by_anc[, match_fraction := n_samples_match / n_samples_total]
		locus.by_anc[, match_percent := match_fraction * 100]
		setorder(locus.by_anc, locus_id, anc)
		fwrite(locus.by_anc, out_file("locus.by_anc.summary.tsv"), sep = "\t")

		if (length(stats_by_group) > 0) {
			mg <- merge(matches, ph, by = "ID", all.x = FALSE)
			den.g <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
			locus.g <- mg[, .(n_samples_match = uniqueN(ID)), by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
			locus.g <- merge(locus.g, den.g, by = stats_by_group, all.x = TRUE)
			locus.g[, match_fraction := n_samples_match / n_samples_total]
			locus.g[, match_percent := match_fraction * 100]
			group_tag <- paste(stats_by_group, collapse = "_")
			setorderv(locus.g, c(stats_by_group, "locus_id", "anc"))
			fwrite(locus.g, out_file(paste0("locus.by_", group_tag, ".summary.tsv")), sep = "\t")
		}
	}
}

message("Done:")
message("  ", out_file("input.check.tsv"))
message("  ", win.file)
message("  ", top.file)
message("  ", bedgraph.dir)
message("  ", plot.file)
if (nrow(loci) > 0) {
	message("  ", out_file("locus.sample_matches.tsv"))
	message("  ", out_file("locus.by_anc.summary.tsv"))
	if (length(stats_by_group) > 0) message("  ", out_file(paste0("locus.by_", paste(stats_by_group, collapse = "_"), ".summary.tsv")))
}
