library(data.table)
args <- commandArgs(TRUE)
root <- normalizePath(if(length(args)) args[1] else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
hapf <- file.path(root, "hap_match.tsv"); mat0 <- file.path(root, "mat"); phy0 <- file.path(root, "phy")
dir.create(phy0, recursive = TRUE, showWarnings = FALSE)
if(!file.exists(hapf) || file.info(hapf)$size == 0) quit(save = "no", status = 1)
base_set <- c("A", "C", "G", "T")
hap <- fread(hapf)
oldf <- list.files(phy0, pattern = "\\.(phy|meta\\.tsv|phy_phyml_tree\\.txt|phy_phyml_stats\\.txt)$", recursive = TRUE, full.names = TRUE)
if(length(oldf)) unlink(oldf, force = TRUE)

arch_label <- function(f){ x <- tools::file_path_sans_ext(basename(f)); paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))) }
arch_files <- function(d0){ fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE); fs <- fs[basename(fs) != "kg.tsv"]; setNames(fs, arch_label(fs)) }
gt_hap <- function(x){ x <- gsub("\\|", "/", x); sp <- tstrsplit(x, "/", fixed = TRUE); list(a = sp[[1]], b = sp[[2]]) }
to_base <- function(h, ref, alt){ out <- rep("N", length(h)); out[h == "0"] <- ref[h == "0"]; out[h == "1"] <- alt[h == "1"]; out }
arch_base <- function(gt, ref, alt){
	gt <- gsub("\\|", "/", gt); sp <- strsplit(gt, "/", fixed = TRUE)
	vapply(seq_along(sp), function(i){
		z <- sp[[i]]; z <- z[z != "."]
		if(length(z) != 2L || z[1] != z[2]) return(NA_character_)
		a <- z[1]; alts <- if(alt[i] %in% c(".", "")) character(0) else strsplit(alt[i], ",", fixed = TRUE)[[1]]
		if(a == "0") return(ref[i]); ai <- suppressWarnings(as.integer(a))
		if(is.na(ai) || ai < 1L || ai > length(alts)) return(NA_character_)
		alts[ai]
	}, character(1))
}
read_mat <- function(f, type){
	x <- fread(f, header = FALSE)
	if(!nrow(x)) return(NULL)
	if(type == "kg") setnames(x, c("chr", "pos", "ref", "alt", "aa", paste0("s", seq_len(ncol(x) - 5)))) else setnames(x, c("chr", "pos", "ref", "alt", "gt"))
	x[, chr := sub("^chr", "", as.character(chr), ignore.case = TRUE)]
	x[, `:=`(ref = toupper(trimws(ref)), alt = toupper(trimws(alt)))]
	if("aa" %in% names(x)) x[, aa := toupper(trimws(aa))]
	x
}
read_arch <- function(f, lab){
	x <- read_mat(f, "arch")
	if(is.null(x)) return(NULL)
	x[, allele := arch_base(gt, ref, alt)]
	x <- x[, .(chr, pos, allele)]
	setnames(x, "allele", lab)
	x
}
merge_arch <- function(kg, fs, keep_arch){
	x <- kg; keep <- character(0)
	for(lab in keep_arch){
		if(!lab %in% names(fs)) next
		a <- read_arch(fs[[lab]], lab)
		if(is.null(a)) next
		x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
		keep <- c(keep, lab)
	}
	setorder(x, pos)
	list(x = x, arch = keep)
}
clean_phy <- function(tr0, id0){
	od <- file.path(phy0, tr0); dir.create(od, TRUE, FALSE)
	unlink(file.path(od, paste0(id0, c(".full.phy", ".full.meta.tsv", ".full.phy_phyml_tree.txt", ".full.phy_phyml_stats.txt", ".main.phy", ".main.meta.tsv", ".main.phy_phyml_tree.txt", ".main.phy_phyml_stats.txt"))), force = TRUE)
}

keys <- unique(hap[, .(trait, id, matched_archaics)])
for(i in seq_len(nrow(keys))){
	tr0 <- keys$trait[i]; id0 <- keys$id[i]; clean_phy(tr0, id0)
	keep_arch <- unlist(strsplit(keys$matched_archaics[i], ";", fixed = TRUE)); keep_arch <- keep_arch[nzchar(keep_arch)]
	if(!length(keep_arch)) next
	h <- hap[trait == tr0 & id == id0]
	if(!nrow(h)) next
	d0 <- file.path(mat0, tr0, id0); kgf <- file.path(d0, "kg.tsv"); fs <- arch_files(d0)
	if(!file.exists(kgf) || !length(fs)) next
	kg <- read_mat(kgf, "kg")
	if(is.null(kg)) next
	ma <- merge_arch(kg, fs, keep_arch); x <- ma$x; keep_arch <- ma$arch
	if(!length(keep_arch)) next

	samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", keep_arch))
	H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
	keep_modern <- apply(H, 1, function(z){ z <- z[!is.na(z) & z %chin% base_set]; if(!length(z)) return(FALSE); tab <- table(z); length(tab) >= 2L && min(tab) >= 2L })
	keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set))
	x2 <- x[keep_modern & keep_arch_called]
	if(!nrow(x2)) next
	anc <- x2$aa; anc[!anc %in% base_set] <- "N"
	od <- file.path(phy0, tr0); dir.create(od, TRUE, FALSE)

	make_one <- function(sub, tag){
		if(!nrow(sub)) return(invisible(NULL))
		seqs <- setNames(sub$seq, sub$hap_id)
		for(an in keep_arch){ z <- x2[[an]]; z[!z %in% base_set] <- "N"; seqs <- c(seqs, setNames(paste0(z, collapse = ""), an)) }
		seqs <- c(seqs, Ancestral = paste0(anc, collapse = ""))
		len <- unique(nchar(seqs)); if(length(len) != 1L) return(invisible(NULL))
		phyf <- file.path(od, paste0(id0, ".", tag, ".phy")); con <- file(phyf, "w")
		writeLines(sprintf("%d %d", length(seqs), len), con)
		for(j in seq_along(seqs)) writeLines(sprintf("%-12s%s", substr(names(seqs)[j], 1, 12), seqs[j]), con)
		close(con)
		meta <- data.table(label = names(seqs), type = fifelse(names(seqs) %in% keep_arch, "archaic", fifelse(names(seqs) == "Ancestral", "ancestral", "modern")))
		meta <- merge(meta, sub[, .(label = hap_id, n, best_lineage, best_arch, best_match)], by = "label", all.x = TRUE)
		fwrite(meta, file.path(od, paste0(id0, ".", tag, ".meta.tsv")), sep = "\t")
	}
	make_one(h, "full")
	make_one(h[n > 10], "main")
}
