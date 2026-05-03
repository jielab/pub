library(data.table)
args0 <- commandArgs(TRUE)
has_root <- length(args0) && grepl("[/\\\\]", args0[1])
root <- normalizePath(if(has_root) args0[1] else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
args <- if(has_root) args0[-1] else args0
mat0 <- file.path(root, "mat"); hap0 <- file.path(root, "hap"); ld0 <- file.path(root, "ld")
riskf <- file.path(root, "lead", "lead.assoc")
dir.create(hap0, recursive = TRUE, showWarnings = FALSE)
base_set <- c("A", "C", "G", "T")

pick1 <- function(nm, cand){ x <- cand[cand %chin% nm]; if(length(x)) x[1] else NA_character_ }
arch_label <- function(f){ x <- tools::file_path_sans_ext(basename(f)); paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))) }
arch_files <- function(d0){ fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE); fs <- fs[basename(fs) != "kg.tsv"]; setNames(fs, arch_label(fs)) }
lineage_group <- function(x){ y <- tolower(x); fifelse(grepl("vindija|altai|chagyr|neand", y), "Neanderthal", fifelse(grepl("denisova|denisovan", y), "Denisovan", "Archaic")) }

read_risk <- function(f){
	x <- fread(f); nm <- names(x)
	tr_col <- pick1(nm, c("trait","Trait")); sn_col <- pick1(nm, c("lead_snp","SNP","snp","rsid","RSID","rsID","ID","id"))
	rk_col <- pick1(nm, c("risk_allele","RiskAllele","riskAllele")); ea_col <- pick1(nm, c("effect_allele","EA","ea","A1","a1","ALT","alt","Allele1","allele1","tested_allele","TestedAllele"))
	oa_col <- pick1(nm, c("other_allele","OA","oa","NEA","nea","A2","a2","REF","ref","Allele2","allele2","non_effect_allele","Non_Effect_Allele"))
	bt_col <- pick1(nm, c("beta","BETA","Beta","effect","Effect","Estimate","estimate","B","b")); or_col <- pick1(nm, c("OR","or","OddsRatio","odds_ratio","oddsratio"))
	if(is.na(tr_col) || is.na(sn_col)) stop("lead.assoc must contain trait and lead_snp")
	keep <- unique(na.omit(c(tr_col, sn_col, rk_col, ea_col, oa_col, bt_col, or_col))); x <- x[, ..keep]
	map <- c(trait = tr_col, lead_snp = sn_col, risk_allele = rk_col, effect_allele = ea_col, other_allele = oa_col, beta = bt_col, OR = or_col)
	for(nm0 in names(map)) if(!is.na(map[nm0])) setnames(x, map[nm0], nm0)
	if(!"risk_allele" %in% names(x)){
		if(!all(c("effect_allele", "other_allele") %in% names(x))) stop("lead.assoc missing risk_allele and effect_allele/other_allele")
		if("beta" %in% names(x)) x[, risk_allele := fifelse(as.numeric(beta) >= 0, effect_allele, other_allele)]
		else if("OR" %in% names(x)) x[, risk_allele := fifelse(as.numeric(OR) >= 1, effect_allele, other_allele)]
		else stop("lead.assoc missing risk_allele and beta/OR")
	}
	x[, `:=`(trait = as.character(trait), lead_snp = as.character(lead_snp), risk_allele = toupper(as.character(risk_allele)))]
	unique(x[, .(trait, lead_snp, risk_allele)])
}
risk <- read_risk(riskf)

parse_id <- function(id) list(chr = suppressWarnings(as.integer(sub("\\..*$", "", id))), snp = sub("\\.[0-9]+$", "", sub("^[0-9]+\\.", "", id)), bp = suppressWarnings(as.integer(sub("^.*\\.", "", id))))
gt_hap <- function(x){ x <- gsub("\\|", "/", x); sp <- tstrsplit(x, "/", fixed = TRUE); list(a = sp[[1]], b = sp[[2]]) }
roman_n <- function(x) vapply(x, function(i) as.character(as.roman(i)), character(1))
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
score_seq <- function(x, y){ ok <- x != "N" & y != "N" & !is.na(x) & !is.na(y); c(match = sum(x[ok] == y[ok]), n = sum(ok)) }
ils_p <- function(size_bp, recomb_cM_Mb = 0.53, split_years = 550000, archaic_age_years = 50000, gen_years = 29){ r <- recomb_cM_Mb * 1e-8; L <- 1 / (r * ((split_years / gen_years) + ((split_years - archaic_age_years) / gen_years))); 1 - pgamma(size_bp, shape = 2, rate = 1 / L) }

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
merge_arch <- function(kg, fs){
	x <- kg; keep <- character(0)
	for(lab in names(fs)){
		a <- read_arch(fs[[lab]], lab)
		if(is.null(a)) next
		x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
		keep <- c(keep, lab)
	}
	setorder(x, pos)
	list(x = x, arch = keep)
}
comp_allele <- function(x){ x <- toupper(x); ifelse(nchar(x) == 1L, chartr("ACGT", "TGCA", x), x) }
harmonize_lead <- function(a, ref, alt){ a <- toupper(trimws(a)); if(!nzchar(a)) return(NA_character_); if(a == ref || a == alt) return(a); b <- comp_allele(a); if(b == ref || b == alt) b else NA_character_ }
pick_lineage <- function(arch_stat, p, p_cut = 0.1, min_snp = 2L, min_prop = 0.5){
	if(!is.finite(p) || p >= p_cut || !nrow(arch_stat)) return(list(best_lineage = NA_character_, keep_arch = character(0)))
	z <- copy(arch_stat); z[, group := lineage_group(archaic)]; z[, ok := n_compared_risk >= min_snp & n_match_risk >= min_snp & prop_match_risk >= min_prop]
	z <- z[ok == TRUE]
	if(!nrow(z)) return(list(best_lineage = NA_character_, keep_arch = character(0)))
	grp <- z[, .(prop = max(prop_match_risk, na.rm = TRUE), match = max(n_match_risk, na.rm = TRUE)), by = group][order(-prop, -match)][1, group]
	list(best_lineage = grp, keep_arch = z[group == grp][order(-prop_match_risk, -n_match_risk), archaic])
}
write_one <- function(x, f){ if(nrow(x)) fwrite(x, f, sep = "\t") }
clean_one <- function(tr0, id0){ od <- file.path(hap0, tr0); dir.create(od, TRUE, FALSE); unlink(file.path(od, paste0(id0, c(".hap.tsv", ".site.tsv", ".arch.tsv", ".region.tsv", ".core.tsv"))), force = TRUE) }

process_one <- function(tr0, id0){
	od <- file.path(hap0, tr0); dir.create(od, TRUE, FALSE)
	p0 <- parse_id(id0); chr0 <- p0$chr; snp0 <- p0$snp; bp0 <- p0$bp
	ldf <- file.path(ld0, tr0, "ld.tsv"); bdf <- file.path(ld0, tr0, "block.tsv")
	if(!file.exists(ldf) || !file.exists(bdf)) return(invisible(NULL))
	ld <- fread(ldf); ld[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), pos = as.integer(pos), R2 = as.numeric(R2))]; ld <- ld[!is.na(lead_chr) & !is.na(lead_bp) & !is.na(pos)]
	blk <- fread(bdf); blk[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), start = as.integer(start), end = as.integer(end), n = as.integer(n), size_bp = as.integer(size_bp))]
	blk <- blk[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0]
	if(!nrow(blk)) return(invisible(NULL)); blk <- blk[1]
	rk_raw <- risk[trait == tr0 & lead_snp == snp0, risk_allele]; rk_raw <- if(length(rk_raw)) rk_raw[1] else NA_character_
	p <- ils_p(as.numeric(blk$size_bp))
	region <- data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, core_start = blk$start, core_end = blk$end, n_ld_snp = blk$n, core_size_bp = blk$size_bp, p_ils = p, best_lineage = NA_character_, matched_archaics = "")
	site <- data.table(trait = tr0, id = id0, n_site_raw = 0L, n_site_called = 0L, n_site_keep = 0L, n_hap_raw = 0L, n_hap_keep = 0L, best_lineage = NA_character_, matched_archaics = "")
	d0 <- file.path(mat0, tr0, id0); kgf <- file.path(d0, "kg.tsv"); afs <- arch_files(d0)
	if(!file.exists(kgf) || !length(afs)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	kg <- read_mat(kgf, "kg"); if(is.null(kg)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	ma <- merge_arch(kg, afs); x <- ma$x; arch_names <- ma$arch
	if(!length(arch_names)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	site[, `:=`(n_site_raw = nrow(x), n_site_called = nrow(x))]

	core_pos <- ld[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0, sort(unique(pos))]
	xcore <- x[pos %in% core_pos]; setorder(xcore, pos)
	if(!nrow(xcore)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	samp_core <- setdiff(names(xcore), c("chr", "pos", "ref", "alt", "aa", arch_names))
	Hcore <- do.call(cbind, c(lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$a, xcore$ref, xcore$alt)), lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$b, xcore$ref, xcore$alt))))
	colnames(Hcore) <- c(paste0(samp_core, "_1"), paste0(samp_core, "_2"))
	lead_i <- which(xcore$pos == bp0); rk <- if(length(lead_i) == 1L) harmonize_lead(rk_raw, xcore$ref[lead_i], xcore$alt[lead_i]) else NA_character_
	carry <- if(length(lead_i) == 1L && !is.na(rk)) Hcore[lead_i, ] == rk else rep(FALSE, ncol(Hcore))
	core <- copy(xcore)[, .(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, pos, ref, alt)]
	core[, `:=`(lead_risk_raw = rk_raw, lead_risk_harmonized = rk, n_lead_risk_haps = sum(carry, na.rm = TRUE), risk_core_allele = NA_character_, carry_freq = NA_real_, noncarry_freq = NA_real_)]
	if(sum(carry, na.rm = TRUE) > 0L){
		rka <- vapply(seq_len(nrow(Hcore)), function(i){ z <- Hcore[i, carry, drop = TRUE]; z <- z[!is.na(z) & z != "N"]; if(!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] }, character(1))
		core[, risk_core_allele := rka]
		core[, carry_freq := vapply(seq_len(nrow(Hcore)), function(i){ z <- Hcore[i, carry, drop = TRUE]; z <- z[!is.na(z) & z != "N"]; if(!length(z) || is.na(rka[i])) NA_real_ else mean(z == rka[i]) }, numeric(1))]
		core[, noncarry_freq := vapply(seq_len(nrow(Hcore)), function(i){ z <- Hcore[i, !carry, drop = TRUE]; z <- z[!is.na(z) & z != "N"]; if(!length(z) || is.na(rka[i])) NA_real_ else mean(z == rka[i]) }, numeric(1))]
	}
	arch_stat <- rbindlist(lapply(arch_names, function(an){
		ca <- xcore[[an]]; ok1 <- ca %chin% base_set; ok2 <- ok1 & !is.na(core$risk_core_allele); z1 <- core[ok1]; z2 <- core[ok2]; a2 <- ca[ok2]
		data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, archaic = an, lineage_group = lineage_group(an), n_core_ld_snp = length(core_pos), n_risk_defined = sum(!is.na(core$risk_core_allele)), n_callable = nrow(z1), n_compared_risk = nrow(z2), n_match_risk = sum(a2 == z2$risk_core_allele, na.rm = TRUE), n_match_ref = sum(ca[ok1] == z1$ref, na.rm = TRUE), n_match_alt = sum(ca[ok1] == z1$alt, na.rm = TRUE), prop_match_risk = fifelse(nrow(z2) > 0L, sum(a2 == z2$risk_core_allele, na.rm = TRUE) / nrow(z2), NA_real_))
	}), fill = TRUE)
	cls <- pick_lineage(arch_stat, p); keep_arch <- cls$keep_arch; best_lineage_val <- cls$best_lineage; matched_archaics_val <- paste(keep_arch, collapse = ";")
	region[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]; site[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
	write_one(core, file.path(od, paste0(id0, ".core.tsv"))); write_one(arch_stat, file.path(od, paste0(id0, ".arch.tsv")))
	if(!length(keep_arch)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }

	x <- x[nchar(ref) == 1L & nchar(alt) == 1L & ref %chin% base_set & alt %chin% base_set]
	samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", arch_names))
	H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
	colnames(H) <- c(paste0(samp, "_1"), paste0(samp, "_2"))
	keep_modern <- apply(H, 1, function(z){ z <- z[!is.na(z) & z %chin% base_set]; if(!length(z)) return(FALSE); tab <- table(z); length(tab) >= 2L && min(tab) >= 2L })
	keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set)); keep <- keep_modern & keep_arch_called
	x2 <- x[keep]; H2 <- H[keep, , drop = FALSE]; site[, n_site_keep := nrow(x2)]
	if(!nrow(x2)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	hap <- data.table(copy = colnames(H2), seq = apply(H2, 2, paste0, collapse = ""))[, .(n = .N, copies = paste(copy, collapse = ";")), by = seq][order(-n, seq)]
	site[, n_hap_raw := nrow(hap)]; hap <- hap[n > 1L][order(-n, seq)]
	if(!nrow(hap)){ write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv"))); return(invisible(NULL)) }
	hap[, hap_id := roman_n(.I)]; site[, n_hap_keep := nrow(hap)]
	for(an in keep_arch){ aseq <- paste0(x2[[an]], collapse = ""); hap[, (paste0(an, "_match")) := vapply(seq, function(sq) score_seq(strsplit(sq, "", fixed = TRUE)[[1]], strsplit(aseq, "", fixed = TRUE)[[1]])["match"], numeric(1))] }
	cols <- paste0(keep_arch, "_match"); hap[, best_arch := keep_arch[max.col(.SD, ties.method = "first")], .SDcols = cols]; hap[, best_match := apply(.SD, 1, max), .SDcols = cols]
	risk_i <- which(x2$pos == bp0)
	if(length(risk_i) == 1L){ risk_a <- harmonize_lead(rk_raw, x2$ref[risk_i], x2$alt[risk_i]); hap[, carry_risk := if(!is.na(risk_a)) substring(seq, risk_i, risk_i) == risk_a else NA] } else hap[, carry_risk := NA]
	hap[, `:=`(trait = tr0, id = id0, best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
	setcolorder(hap, c("trait", "id", "hap_id", "n", "best_lineage", "matched_archaics", "best_arch", "best_match", "carry_risk", "seq", "copies", cols))
	write_one(hap, file.path(od, paste0(id0, ".hap.tsv"))); write_one(site, file.path(od, paste0(id0, ".site.tsv"))); write_one(region, file.path(od, paste0(id0, ".region.tsv")))
}

merge_all <- function(){
	gather <- function(pat){ fs <- list.files(hap0, pattern = pat, recursive = TRUE, full.names = TRUE); fs <- fs[file.info(fs)$size > 0]; if(!length(fs)) return(data.table()); rbindlist(lapply(fs, fread), fill = TRUE) }
	hap_dt <- gather("\\.hap\\.tsv$"); site_dt <- gather("\\.site\\.tsv$"); arch_dt <- gather("\\.arch\\.tsv$"); region_dt <- gather("\\.region\\.tsv$"); core_dt <- gather("\\.core\\.tsv$")
	unlink(file.path(root, c("hap_match.tsv", "hap_site_count.tsv", "core_archaic_match.tsv", "region_summary.tsv", "core_risk.tsv")), force = TRUE)
	if(nrow(hap_dt)) fwrite(hap_dt, file.path(root, "hap_match.tsv"), sep = "\t"); if(nrow(site_dt)) fwrite(site_dt, file.path(root, "hap_site_count.tsv"), sep = "\t")
	if(nrow(arch_dt)) fwrite(arch_dt, file.path(root, "core_archaic_match.tsv"), sep = "\t"); if(nrow(region_dt)) fwrite(region_dt, file.path(root, "region_summary.tsv"), sep = "\t"); if(nrow(core_dt)) fwrite(core_dt, file.path(root, "core_risk.tsv"), sep = "\t")
}
run_trait <- function(tr0){ td <- file.path(mat0, tr0); if(!dir.exists(td)) stop("trait not found: ", tr0); for(id0 in list.dirs(td, recursive = FALSE, full.names = FALSE)){ clean_one(tr0, id0); process_one(tr0, id0) } }
run_all <- function(){ for(tr0 in list.dirs(mat0, recursive = FALSE, full.names = FALSE)) run_trait(tr0); merge_all() }
if(!length(args)) run_all() else if(length(args) == 1L && args[1] == "merge") merge_all() else if(length(args) == 1L) run_trait(args[1]) else if(length(args) == 2L){ clean_one(args[1], args[2]); process_one(args[1], args[2]) } else stop("usage: Rscript make_hap.R <res_dir> [trait | trait id | merge]")
