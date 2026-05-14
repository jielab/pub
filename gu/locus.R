#!/usr/bin/env Rscript

# Unified standalone R entry point for the locus workflow.
# Commands: prep_input, prep_archaic, make_hap, filter_hap, make_phy, make_tree.

usage <- function() {
    stop(paste("Usage:", "  Rscript locus.R prep_input --dirgwas DIR --dirout DIR --dirmod DIR --traits trait ...", 
        "  Rscript locus.R prep_archaic kg.vcf.gz archaic.raw.vcf.gz out.vcf.gz sample_name", 
        "  Rscript locus.R make_hap RES_DIR [trait | trait id | merge]", "  Rscript locus.R filter_hap RES_DIR sample.tsv POP MAX_COUNT", 
        "  Rscript locus.R make_phy RES_DIR [trait | trait id]", "  Rscript locus.R make_tree RES_DIR", 
        sep = "\n"), call. = FALSE)
}

run_prep_input <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    fphe <- "/mnt/d/scripts/f/phe.f.R"
    if (file.exists(fphe)) 
        source(fphe)
    args <- commandArgs(TRUE)
    arg <- function(k, multi = FALSE) {
        i <- which(args %in% k)[1]
        if (is.na(i) || i == length(args)) 
            return(NULL)
        j <- i + 1
        if (!multi) 
            return(args[j])
        z <- args[j:length(args)]
        z[seq_len(which(c(grepl("^-", z[-1]), TRUE))[1])]
    }
    pick1 <- function(nm, z) {
        x <- z[z %chin% nm]
        if (length(x)) 
            x[1]
        else NA_character_
    }
    chr_int <- function(x) {
        x <- toupper(sub("^CHR", "", as.character(x)))
        x[x == "X"] <- "23"
        x[x == "Y"] <- "24"
        suppressWarnings(as.integer(x))
    }
    id_mode <- function(x) {
        x <- unique(na.omit(as.character(x)))
        x <- x[nzchar(x) & x != "."]
        x <- head(x, 5000)
        if (!length(x)) 
            return("missing")
        y <- toupper(x)
        rs <- mean(grepl("^RS[0-9]+$", y))
        cp <- mean(grepl("^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$", y))
        if (rs > 0.8) 
            "rsid"
        else if (cp > 0.8) 
            "chrpos"
        else "other"
    }
    usable_id <- function(x) {
        !is.na(x) & nzchar(x) & x != "."
    }
    allele_score <- function(ea, nea, ref, alt) {
        ea <- toupper(as.character(ea))
        nea <- toupper(as.character(nea))
        ref <- toupper(as.character(ref))
        alt <- toupper(as.character(alt))
        pair <- ((ea == ref & nea == alt) | (ea == alt & nea == ref))
        one <- (ea == ref | ea == alt | nea == ref | nea == alt)
        fifelse(pair, 2L, fifelse(one, 1L, 0L))
    }
    rd_gwas <- function(f, keep_snps = NULL) {
        h <- names(fread(f, nrows = 0))
        cc <- c(SNP = pick1(h, c("SNP", "rsid", "RSID", "ID")), EA = pick1(h, c("EA", "A1", "effect_allele", 
            "ALT", "alt")), NEA = pick1(h, c("NEA", "A2", "other_allele", "REF", "ref")), BETA = pick1(h, 
            c("BETA", "beta", "b")))
        if (anyNA(cc)) 
            stop("Missing GWAS columns in ", f, ": ", paste(h, collapse = ","))
        if (!is.null(keep_snps)) {
            keep_snps <- unique(as.character(keep_snps[!is.na(keep_snps) & nzchar(keep_snps)]))
            if (!length(keep_snps)) 
                return(data.table(SNP = character(), EA = character(), NEA = character(), BETA = numeric()))
            tmp <- tempfile()
            fwrite(data.table(SNP = keep_snps), tmp, col.names = FALSE, sep = "\t")
            idx <- match(unname(cc), h)
            cmd <- sprintf("gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{k[$1]=1; next} FNR==1{next} ($%d in k){print $%d,$%d,$%d,$%d}' %s -", 
                shQuote(f), idx[1], idx[1], idx[2], idx[3], idx[4], shQuote(tmp))
            x <- fread(cmd = cmd, col.names = names(cc), showProgress = FALSE)
            unlink(tmp)
        }
        else {
            x <- fread(f, select = unname(cc), showProgress = FALSE)
            setnames(x, unname(cc), names(cc))
        }
        unique(x[, .(SNP = as.character(SNP), EA = toupper(EA), NEA = toupper(NEA), BETA = as.numeric(BETA))], 
            by = "SNP")
    }
    pvar_files <- function(dirmod, chr) {
        if (chr == 23L) 
            unique(c(paste0(dirmod, "/EUR.male.chrX.par.pvar"), paste0(dirmod, "/EUR.male.chrX.nonPar.pvar"), 
                paste0(dirmod, "/EUR.chrX.pvar")))
        else paste0(dirmod, "/EUR.chr", chr, ".pvar")
    }
    pvar_ids <- function(f, n = 2000) {
        if (!file.exists(f)) 
            return(character())
        cmd <- sprintf("awk '/^#/{next} $3!=\".\" && $3!=\"\"{print $3; if(++n==%d) exit}' %s", 
            as.integer(n), shQuote(f))
        system(cmd, intern = TRUE)
    }
    pvar_mode <- function(fs) {
        id_mode(unlist(lapply(fs, pvar_ids), use.names = FALSE))
    }
    rd_pvar <- function(dirmod, lead) {
        out <- list()
        k <- 0L
        for (chr in sort(unique(lead$lead_chr))) {
            fs <- pvar_files(dirmod, chr)
            fs <- fs[file.exists(fs)]
            if (!length(fs)) 
                stop("No 1000G .pvar found for chr ", chr, " under ", dirmod)
            key <- unique(lead[lead_chr == chr, .(lead_bp, lead_snp)])
            tmp <- tempfile()
            fwrite(key, tmp, col.names = FALSE, sep = "\t")
            for (f in fs) {
                awk_script <- paste0(
                  "BEGIN{FS=OFS=\"\\t\"} ",
                  "NR==FNR{p[$1]=1; id[$2]=1; next} ",
                  "/^##/{next} ",
                  "/^#CHROM/{print \"CHR\",$2,$3,$4,$5; next} ",
                  "FNR==1 && ($1==\"CHR\" || $1==\"CHROM\")",
                  "{print \"CHR\",$2,$3,$4,$5; next} ",
                  "(($2 in p) || ($3 in id))",
                  "{gsub(/^chr/,\"\",$1); print $1,$2,$3,$4,$5}"
                )
                cmd <- sprintf("awk %s %s %s", shQuote(awk_script),
                  shQuote(tmp), shQuote(f))
                x <- fread(cmd = cmd, header = TRUE, showProgress = FALSE)
                if (nrow(x) <= 0) 
                  next
                cc <- c(CHR = pick1(names(x), c("CHR", "CHROM")), POS = pick1(names(x), c("POS", 
                  "BP")), ID = pick1(names(x), c("ID", "SNP")), REF = pick1(names(x), c("REF")), 
                  ALT = pick1(names(x), c("ALT")))
                if (anyNA(cc)) 
                  stop("Bad pvar columns after skipping ## lines: ", f, "; names=", paste(names(x), 
                    collapse = ","))
                k <- k + 1L
                out[[k]] <- x[, .(CHR = chr_int(get(cc["CHR"])), POS = as.integer(get(cc["POS"])), 
                  ID = as.character(get(cc["ID"])), REF = toupper(get(cc["REF"])), ALT = toupper(get(cc["ALT"])), 
                  pvar_file = f)]
            }
            unlink(tmp)
        }
        if (!length(out)) 
            return(data.table(CHR = integer(), POS = integer(), ID = character(), REF = character(), 
                ALT = character(), pvar_file = character()))
        unique(rbindlist(out, fill = TRUE)[!is.na(CHR) & !is.na(POS)])
    }
    align_1kg <- function(z, dirmod, lead0, tr) {
        z[, `:=`(rid, .I)]
        if (is.null(dirmod)) {
            z[, `:=`(lead_snp0 = lead_snp, pvar_id = lead_snp, match_type_1kg = "not_checked", 
                pvar_file = NA_character_, pvar_bp = lead_bp)]
            return(z)
        }
        fs <- unique(unlist(lapply(sort(unique(z$lead_chr)), function(chr) pvar_files(dirmod, 
            chr))))
        fs <- fs[file.exists(fs)]
        mode <- pvar_mode(fs)
        cojo_mode <- id_mode(z$lead_snp)
        p <- rd_pvar(dirmod, z)
        fwrite(data.table(source = c("COJO_lead", "1000G_pvar"), id_mode = c(cojo_mode, mode), 
            example = c(paste(head(z$lead_snp, 8), collapse = ","), paste(head(p$ID, 8), collapse = ","))), 
            paste0(lead0, "/", tr, ".id_sanity.tsv"), sep = "\t")
        if (!mode %chin% c("rsid", "chrpos")) {
            stop(
              "ERROR: 1000G pvar ID mode is ", mode, " for ", tr,
              ". Expected rsID or CHR:POS:REF:ALT. Fix the 1kg VCF/pfile IDs first; ",
              "this script will not run --set-all-var-ids automatically."
            )
        }
        z[, `:=`(lead_snp0 = lead_snp, lead_bp0 = lead_bp, pvar_id = NA_character_, pvar_bp = NA_integer_, 
            pvar_file = NA_character_, match_type_1kg = NA_character_)]
        p[, `:=`(pvar_rank, fifelse(grepl("male\\.chrX\\.(par|nonPar)\\.pvar$", pvar_file), 1L, 
            2L))]
        m <- data.table(rid = integer(), pvar_id = character(), pvar_bp = integer(), pvar_file = character(), 
            match_type_1kg = character(), priority = integer())
        idm <- merge(z[usable_id(lead_snp), .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)], p[usable_id(ID)], 
            by.x = "lead_snp", by.y = "ID", allow.cartesian = TRUE)
        idm <- idm[lead_chr == CHR]
        if (nrow(idm)) {
            idm[, `:=`(allele_match, allele_score(EA, NEA, REF, ALT))]
            idm <- idm[allele_match > 0]
            idm[, `:=`(pos_same, as.integer(lead_bp == POS))]
            if (nrow(idm)) {
                setorder(idm, rid, -allele_match, -pos_same, pvar_rank)
                idm <- idm[, .SD[1], by = rid]
                m <- rbind(m, idm[, .(rid, pvar_id = lead_snp, pvar_bp = POS, pvar_file, match_type_1kg = fifelse(pos_same == 
                  1L, "ID_allele_POS_same", "ID_allele_POS_from_pvar"), priority = 0L)], fill = TRUE)
            }
        }
        unmatched <- z[!(rid %in% m$rid)]
        if (nrow(unmatched)) {
            pm <- merge(unmatched[, .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)], p, by.x = c("lead_chr", 
                "lead_bp"), by.y = c("CHR", "POS"), allow.cartesian = TRUE)
            pm <- pm[(!usable_id(lead_snp)) | (!usable_id(ID))]
            if (nrow(pm)) {
                pm[, `:=`(allele_match, allele_score(EA, NEA, REF, ALT))]
                pm <- pm[allele_match > 0]
                if (nrow(pm)) {
                  setorder(pm, rid, -allele_match, pvar_rank)
                  pm <- pm[, .SD[1], by = rid]
                  pm[, `:=`(out_id, fifelse(usable_id(ID), ID, lead_snp))]
                  m <- rbind(m, pm[, .(rid, pvar_id = out_id, pvar_bp = lead_bp, pvar_file, match_type_1kg = "POS_allele_missing_ID", 
                    priority = 1L)], fill = TRUE)
                }
            }
        }
        h <- m
        h <- h[!is.na(pvar_id) & nzchar(pvar_id) & pvar_id != "."]
        if (nrow(h)) {
            setorder(h, rid, priority)
            h <- h[, .SD[1], by = rid]
            z[h, `:=`(pvar_id = i.pvar_id, pvar_bp = i.pvar_bp, pvar_file = i.pvar_file, match_type_1kg = i.match_type_1kg), 
                on = .(rid)]
        }
        z[, `:=`(lead_snp, fifelse(!is.na(pvar_id), pvar_id, lead_snp))]
        z[, `:=`(lead_bp, fifelse(!is.na(pvar_bp), pvar_bp, lead_bp))]
        z[]
    }
    dirgwas <- arg(c("--dirgwas", "-dirgwas"))
    dirout <- arg(c("--dirout", "-dirout"))
    dirmod <- arg(c("--dirmod", "-dirmod"))
    traits <- arg(c("--traits", "-traits"), TRUE)
    if (is.null(dirgwas) || is.null(dirout) || is.null(traits)) 
        stop("Usage: Rscript locus.R prep_input --dirgwas DIR --dirout DIR --dirmod 1KG_DIR --traits bald bald12 ...")
    traits <- unlist(strsplit(paste(traits, collapse = " "), "[, ]+"))
    traits <- traits[nzchar(traits)]
    lead0 <- paste0(normalizePath(dirout, winslash = "/", mustWork = FALSE), "/lead")
    dir.create(lead0, recursive = TRUE, showWarnings = FALSE)
    all_match <- list()
    prep_version <- "prep_input_v3_id_first_pos_if_missing_id"
    for (tr in traits) {
        cj_f <- paste0(dirgwas, "/cojo/", tr, "/", tr, ".jma.cojo")
        gw_f <- paste0(dirgwas, "/clean/", tr, "/", tr, ".gz")
        cj <- fread(cj_f)
        cc <- c(chr = pick1(names(cj), c("Chr", "CHR", "chr", "#CHROM", "CHROM")), bp = pick1(names(cj), 
            c("bp", "BP", "POS", "pos")))
        if (anyNA(cc) || !"SNP" %chin% names(cj)) 
            stop("COJO needs Chr/SNP/bp: ", cj_f)
        lead <- unique(cj[, .(lead_chr = chr_int(get(cc["chr"])), lead_snp = as.character(SNP), 
            lead_bp = as.integer(get(cc["bp"])))])[lead_chr %between% c(1L, 23L) & !is.na(lead_bp)]
        g <- rd_gwas(gw_f, lead$lead_snp)
        z <- merge(lead, g, by.x = "lead_snp", by.y = "SNP", all.x = TRUE, sort = FALSE)
        z[, `:=`(fail, is.na(EA) | !nzchar(EA) | is.na(NEA) | !nzchar(NEA) | is.na(BETA))]
        fail <- z[fail == TRUE, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele = EA, 
            other_allele = NEA, beta = BETA, fail_reason = "not_in_gwas_or_missing_allele_beta")]
        z <- z[fail == FALSE]
        if (nrow(z)) 
            z <- align_1kg(z, dirmod, lead0, tr)
        else z[, `:=`(lead_snp0 = lead_snp, lead_bp0 = lead_bp, pvar_id = NA_character_, pvar_bp = NA_integer_, 
            pvar_file = NA_character_, match_type_1kg = NA_character_)]
        all_match[[tr]] <- z[, .(trait = tr, lead_chr, lead_snp0, lead_bp0, lead_snp, lead_bp, 
            EA, NEA, BETA, pvar_id, pvar_bp, pvar_file, match_type_1kg)]
        fail <- rbind(fail, z[is.na(pvar_id), .(trait = tr, lead_chr, lead_snp = lead_snp0, lead_bp = lead_bp0, 
            effect_allele = EA, other_allele = NEA, beta = BETA, fail_reason = "not_in_1000G_pvar_or_allele_mismatch")], 
            fill = TRUE)
        z <- z[!is.na(pvar_id)][, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele = EA, 
            other_allele = NEA, beta = BETA, risk_allele = fifelse(BETA >= 0, EA, NEA), lead_snp0, 
            match_type_1kg)]
        fwrite(z[, .(lead_chr, lead_snp, lead_bp)], paste0(lead0, "/", tr, ".lead.3col"), sep = "\t")
        if (nrow(fail)) 
            fwrite(fail, paste0(lead0, "/", tr, ".lead.fail.tsv"), sep = "\t")
        else unlink(paste0(lead0, "/", tr, ".lead.fail.tsv"))
        fwrite(z, paste0(lead0, "/", tr, ".lead.assoc"), sep = "\t")
    }
    fwrite(rbindlist(lapply(traits, function(tr) fread(paste0(lead0, "/", tr, ".lead.assoc"))), 
        fill = TRUE), paste0(lead0, "/lead.assoc"), sep = "\t")
    fwrite(rbindlist(all_match, fill = TRUE), paste0(lead0, "/lead.match.tsv"), sep = "\t")
    writeLines(prep_version, paste0(lead0, "/prep.version"))
}

run_prep_archaic <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    args <- commandArgs(TRUE)
    if (length(args) < 4) {
        stop("Usage: Rscript locus.R prep_archaic kg.vcf.gz archaic.raw.vcf.gz out.vcf.gz sample_name")
    }
    kg_vcf <- args[1]
    arc_vcf <- args[2]
    out_vcf <- args[3]
    sample <- args[4]
    if (!file.exists(kg_vcf)) 
        stop("Missing kg_vcf: ", kg_vcf)
    if (!file.exists(arc_vcf)) 
        stop("Missing arc_vcf: ", arc_vcf)
    bcftools <- Sys.which("bcftools")
    bgzip <- Sys.which("bgzip")
    tabix <- Sys.which("tabix")
    if (bcftools == "") 
        stop("bcftools not found in PATH")
    if (bgzip == "") 
        stop("bgzip not found in PATH")
    if (tabix == "") 
        stop("tabix not found in PATH")
    message("kg_vcf  = ", kg_vcf)
    message("arc_vcf = ", arc_vcf)
    message("out_vcf = ", out_vcf)
    message("sample  = ", sample)
    read_query <- function(vcf, fmt, cols, label) {
        tmp <- tempfile()
        rc <- system2(bcftools, c("query", "-f", shQuote(fmt), shQuote(vcf)), stdout = tmp)
        if (rc != 0) 
            stop("bcftools query failed for ", label, ": ", vcf)
        if (!file.exists(tmp) || file.info(tmp)$size == 0) {
            unlink(tmp)
            return(as.data.table(setNames(rep(list(character()), length(cols)), cols)))
        }
        x <- fread(tmp, header = FALSE, col.names = cols, showProgress = FALSE)
        unlink(tmp)
        x
    }
    kg <- read_query(kg_vcf, "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n", c("CHR", "POS", "ID", "REF_kg", 
        "ALT_kg"), "1000G template")
    arc <- read_query(arc_vcf, "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", c("CHR", "POS", "REF_arc", 
        "ALT_arc", "GT_arc"), "archaic raw")
    if (nrow(kg) == 0) 
        stop("kg_vcf has 0 variants: ", kg_vcf)
    if (nrow(arc) == 0) 
        warning("arc_vcf has 0 records: ", arc_vcf)
    kg[, `:=`(idx = .I, CHR = as.character(CHR), POS = as.integer(POS), ID = fifelse(is.na(ID) | 
        ID == "", ".", as.character(ID)), REF_kg = toupper(REF_kg), ALT_kg = toupper(ALT_kg))]
    arc[, `:=`(CHR = as.character(CHR), POS = as.integer(POS), REF_arc = toupper(REF_arc), ALT_arc = toupper(ALT_arc), 
        GT_arc = gsub("\\|", "/", as.character(GT_arc)))]
    x <- merge(kg, arc, by = c("CHR", "POS"), all.x = TRUE, sort = FALSE)
    setorder(x, idx)
    x[, `:=`(GT, "./.")]
    x[, `:=`(project_match, fifelse(is.na(REF_arc), "missing_archaic_record", "unmatched_allele"))]
    x[!is.na(REF_arc) & REF_arc == REF_kg & ALT_arc == ".", `:=`(GT = "0/0", project_match = "REF_arc_matches_REF_kg")]
    x[!is.na(REF_arc) & REF_arc == ALT_kg & ALT_arc == ".", `:=`(GT = "1/1", project_match = "REF_arc_matches_ALT_kg")]
    flip_gt <- function(gt) {
        gt <- gsub("\\|", "/", as.character(gt))
        out <- rep("./.", length(gt))
        out[gt == "0/0"] <- "1/1"
        out[gt == "1/1"] <- "0/0"
        out[gt %chin% c("0/1", "1/0")] <- "0/1"
        out[gt %chin% c(".", "./.") | is.na(gt)] <- "./."
        out
    }
    x[!is.na(REF_arc) & REF_arc == REF_kg & ALT_arc == ALT_kg, `:=`(GT = GT_arc, project_match = "REF_ALT_exact")]
    x[!is.na(REF_arc) & REF_arc == ALT_kg & ALT_arc == REF_kg, `:=`(GT = flip_gt(GT_arc), project_match = "REF_ALT_swapped")]
    x[, `:=`(project_priority, fcase(project_match == "REF_ALT_exact", 1L, project_match == "REF_ALT_swapped", 
        2L, project_match == "REF_arc_matches_REF_kg", 3L, project_match == "REF_arc_matches_ALT_kg", 
        4L, default = 9L))]
    setorder(x, idx, project_priority)
    x <- x[, .SD[1], by = idx]
    setorder(x, idx)
    x[GT == "0", `:=`(GT, "0/0")]
    x[GT == "1", `:=`(GT, "1/1")]
    x[GT == ".", `:=`(GT, "./.")]
    x[is.na(GT) | !(GT %chin% c("0/0", "0/1", "1/0", "1/1", "./.")), `:=`(GT, "./.")]
    x <- x[REF_kg %chin% c("A", "C", "G", "T") & ALT_kg %chin% c("A", "C", "G", "T")]
    vcf <- x[, .(CHR, POS, ID, REF = REF_kg, ALT = ALT_kg, QUAL = ".", FILTER = ".", INFO = ".", 
        FORMAT = "GT", GT)]
    dir.create(dirname(out_vcf), recursive = TRUE, showWarnings = FALSE)
    tmp_vcf <- tempfile(fileext = ".vcf")
    cat("##fileformat=VCFv4.2\n", file = tmp_vcf)
    cat("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", file = tmp_vcf, append = TRUE)
    cat(paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample, 
        sep = "\t"), "\n", file = tmp_vcf, append = TRUE)
    fwrite(vcf, tmp_vcf, sep = "\t", append = TRUE, col.names = FALSE)
    if (file.exists(out_vcf)) 
        file.remove(out_vcf)
    if (file.exists(paste0(out_vcf, ".tbi"))) 
        file.remove(paste0(out_vcf, ".tbi"))
    tmp_gz <- tempfile(fileext = ".vcf.gz")
    rc1 <- system2(bgzip, c("-c", shQuote(tmp_vcf)), stdout = tmp_gz)
    if (rc1 != 0) 
        stop("bgzip failed for: ", out_vcf)
    if (!file.copy(tmp_gz, out_vcf, overwrite = TRUE)) 
        stop("failed to copy bgzip output to: ", out_vcf)
    unlink(tmp_gz)
    if (!file.exists(out_vcf) || file.info(out_vcf)$size == 0) 
        stop("bgzip did not create output: ", out_vcf)
    rc2 <- system2(tabix, c("-f", "-p", "vcf", shQuote(out_vcf)))
    if (rc2 != 0) 
        stop("tabix failed for: ", out_vcf)
    if (!file.exists(paste0(out_vcf, ".tbi"))) 
        stop("tabix did not create index: ", paste0(out_vcf, ".tbi"))
    message("Wrote: ", out_vcf)
    message("n_template=", nrow(kg))
    message("n_archaic_raw=", nrow(arc))
    message("n_output=", nrow(vcf))
    message("GT_0_0=", sum(vcf$GT == "0/0"))
    message("GT_0_1=", sum(vcf$GT %chin% c("0/1", "1/0")))
    message("GT_1_1=", sum(vcf$GT == "1/1"))
    message("GT_missing=", sum(vcf$GT == "./."))
    message("projection_counts=", paste(names(table(x$project_match)), as.integer(table(x$project_match)), 
        sep = ":", collapse = ","))
}

run_make_hap <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    args0 <- commandArgs(TRUE)
    has_root <- length(args0) && grepl("[/\\\\]", args0[1])
    root <- normalizePath(if (has_root) 
        args0[1]
    else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
    args <- if (has_root) 
        args0[-1]
    else args0
    mat0 <- file.path(root, "mat")
    hap0 <- file.path(root, "hap")
    ld0 <- file.path(root, "ld")
    report0 <- file.path(root, "report")
    riskf <- file.path(root, "lead", "lead.assoc")
    dir.create(hap0, recursive = TRUE, showWarnings = FALSE)
    dir.create(report0, recursive = TRUE, showWarnings = FALSE)
    base_set <- c("A", "C", "G", "T")
    pick1 <- function(nm, cand) {
        x <- cand[cand %chin% nm]
        if (length(x)) 
            x[1]
        else NA_character_
    }
    arch_label <- function(f) {
        x <- tools::file_path_sans_ext(basename(f))
        paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    }
    arch_files <- function(d0) {
        fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE)
        fs <- fs[!basename(fs) %in% c("kg.tsv", "kg.samples.tsv")]
        setNames(fs, arch_label(fs))
    }
    lineage_group <- function(x) {
        y <- tolower(x)
        fifelse(grepl("vindija|altai|chagyr|neand", y), "Neanderthal", fifelse(grepl("denisova|denisovan", 
            y), "Denisovan", "Archaic"))
    }
    read_risk <- function(f) {
        x <- fread(f)
        nm <- names(x)
        tr_col <- pick1(nm, c("trait", "Trait"))
        sn_col <- pick1(nm, c("lead_snp", "SNP", "snp", "rsid", "RSID", "rsID", "ID", "id"))
        rk_col <- pick1(nm, c("risk_allele", "RiskAllele", "riskAllele"))
        ea_col <- pick1(nm, c("effect_allele", "EA", "ea", "A1", "a1", "ALT", "alt", "Allele1", 
            "allele1", "tested_allele", "TestedAllele"))
        oa_col <- pick1(nm, c("other_allele", "OA", "oa", "NEA", "nea", "A2", "a2", "REF", "ref", 
            "Allele2", "allele2", "non_effect_allele", "Non_Effect_Allele"))
        bt_col <- pick1(nm, c("beta", "BETA", "Beta", "effect", "Effect", "Estimate", "estimate", 
            "B", "b"))
        or_col <- pick1(nm, c("OR", "or", "OddsRatio", "odds_ratio", "oddsratio"))
        if (is.na(tr_col) || is.na(sn_col)) 
            stop("lead.assoc must contain trait and lead_snp")
        keep <- unique(na.omit(c(tr_col, sn_col, rk_col, ea_col, oa_col, bt_col, or_col)))
        x <- x[, ..keep]
        map <- c(trait = tr_col, lead_snp = sn_col, risk_allele = rk_col, effect_allele = ea_col, 
            other_allele = oa_col, beta = bt_col, OR = or_col)
        for (nm0 in names(map)) if (!is.na(map[nm0])) 
            setnames(x, map[nm0], nm0)
        if (!"risk_allele" %in% names(x)) {
            if (!all(c("effect_allele", "other_allele") %in% names(x))) 
                stop("lead.assoc missing risk_allele and effect_allele/other_allele")
            if ("beta" %in% names(x)) 
                x[, `:=`(risk_allele, fifelse(as.numeric(beta) >= 0, effect_allele, other_allele))]
            else if ("OR" %in% names(x)) 
                x[, `:=`(risk_allele, fifelse(as.numeric(OR) >= 1, effect_allele, other_allele))]
            else stop("lead.assoc missing risk_allele and beta/OR")
        }
        x[, `:=`(trait = as.character(trait), lead_snp = as.character(lead_snp), risk_allele = toupper(as.character(risk_allele)))]
        unique(x[, .(trait, lead_snp, risk_allele)])
    }
    risk <- read_risk(riskf)
    chr_id <- function(x) {
        x <- toupper(sub("^CHR", "", as.character(x)))
        x[x == "X"] <- "23"
        suppressWarnings(as.integer(x))
    }
    parse_id <- function(id) {
        z <- strsplit(id, "\\.")[[1]]
        list(chr = chr_id(z[1]), snp = paste(z[2:(length(z) - 1)], collapse = "."), bp = suppressWarnings(as.integer(z[length(z)])))
    }
    gt_hap <- function(x) {
        x <- gsub("\\|", "/", x)
        sp <- tstrsplit(x, "/", fixed = TRUE)
        list(a = sp[[1]], b = sp[[2]])
    }
    roman_n <- function(x) vapply(x, function(i) as.character(as.roman(i)), character(1))
    to_base <- function(h, ref, alt) {
        out <- rep("N", length(h))
        out[h == "0"] <- ref[h == "0"]
        out[h == "1"] <- alt[h == "1"]
        out
    }
    arch_base <- function(gt, ref, alt) {
        gt <- gsub("\\|", "/", gt)
        sp <- strsplit(gt, "/", fixed = TRUE)
        vapply(seq_along(sp), function(i) {
            z <- sp[[i]]
            z <- z[z != "."]
            if (length(z) != 2L || z[1] != z[2]) 
                return(NA_character_)
            a <- z[1]
            alts <- if (alt[i] %in% c(".", "")) 
                character(0)
            else strsplit(alt[i], ",", fixed = TRUE)[[1]]
            if (a == "0") 
                return(ref[i])
            ai <- suppressWarnings(as.integer(a))
            if (is.na(ai) || ai < 1L || ai > length(alts)) 
                return(NA_character_)
            alts[ai]
        }, character(1))
    }
    score_seq <- function(x, y) {
        ok <- x != "N" & y != "N" & !is.na(x) & !is.na(y)
        c(match = sum(x[ok] == y[ok]), n = sum(ok))
    }
    ils_p <- function(size_bp, recomb_cM_Mb = 0.53, split_years = 550000, archaic_age_years = 50000, 
        gen_years = 29) {
        r <- recomb_cM_Mb * 1e-08
        L <- 1/(r * ((split_years/gen_years) + ((split_years - archaic_age_years)/gen_years)))
        1 - pgamma(size_bp, shape = 2, rate = 1/L)
    }
    read_mat <- function(f, type) {
        x <- fread(f, header = FALSE)
        if (!nrow(x)) 
            return(NULL)
        if (type == "kg") 
            setnames(x, c("chr", "pos", "ref", "alt", "aa", paste0("s", seq_len(ncol(x) - 5))))
        else setnames(x, c("chr", "pos", "ref", "alt", "gt"))
        x[, `:=`(chr, sub("^chr", "", as.character(chr), ignore.case = TRUE))]
        x[, `:=`(ref = toupper(trimws(ref)), alt = toupper(trimws(alt)))]
        if ("aa" %in% names(x)) 
            x[, `:=`(aa, toupper(trimws(aa)))]
        x
    }
    read_arch <- function(f, lab) {
        x <- read_mat(f, "arch")
        if (is.null(x)) 
            return(NULL)
        x[, `:=`(allele, arch_base(gt, ref, alt))]
        x <- x[, .(chr, pos, allele)]
        setnames(x, "allele", lab)
        x
    }
    merge_arch <- function(kg, fs) {
        x <- kg
        keep <- character(0)
        for (lab in names(fs)) {
            a <- read_arch(fs[[lab]], lab)
            if (is.null(a)) 
                next
            x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
            keep <- c(keep, lab)
        }
        setorder(x, pos)
        list(x = x, arch = keep)
    }
    comp_allele <- function(x) {
        x <- toupper(x)
        ifelse(nchar(x) == 1L, chartr("ACGT", "TGCA", x), x)
    }
    harmonize_lead <- function(a, ref, alt) {
        a <- toupper(trimws(a))
        if (!nzchar(a)) 
            return(NA_character_)
        if (a == ref || a == alt) 
            return(a)
        b <- comp_allele(a)
        if (b == ref || b == alt) 
            b
        else NA_character_
    }
    pick_lineage <- function(arch_stat, p, p_cut = 0.1, min_snp = 2L, min_prop = 0.5) {
        if (!is.finite(p) || p >= p_cut || !nrow(arch_stat)) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        z <- copy(arch_stat)
        z[, `:=`(group, lineage_group(archaic))]
        z[, `:=`(ok, n_compared_risk >= min_snp & n_match_risk >= min_snp & prop_match_risk >= 
            min_prop)]
        z <- z[ok == TRUE]
        if (!nrow(z)) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        grp <- z[, .(prop = max(prop_match_risk, na.rm = TRUE), match = max(n_match_risk, na.rm = TRUE)), 
            by = group][order(-prop, -match)][1, group]
        list(best_lineage = grp, keep_arch = z[group == grp][order(-prop_match_risk, -n_match_risk), 
            archaic])
    }
    write_one <- function(x, f) {
        if (nrow(x)) 
            fwrite(x, f, sep = "\t")
    }
    clean_one <- function(tr0, id0) {
        od <- file.path(hap0, tr0)
        dir.create(od, TRUE, FALSE)
        unlink(file.path(od, paste0(id0, c(".hap.tsv", ".site.tsv", ".arch.tsv", ".region.tsv", 
            ".core.tsv"))), force = TRUE)
    }
    process_one <- function(tr0, id0) {
        od <- file.path(hap0, tr0)
        dir.create(od, TRUE, FALSE)
        p0 <- parse_id(id0)
        chr0 <- p0$chr
        snp0 <- p0$snp
        bp0 <- p0$bp
        ldf <- file.path(ld0, tr0, "ld.tsv")
        bdf <- file.path(ld0, tr0, "block.tsv")
        if (!file.exists(ldf) || !file.exists(bdf)) 
            return(invisible(NULL))
        ld <- fread(ldf)
        ld[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), pos = as.integer(pos), 
            R2 = as.numeric(R2))]
        ld <- ld[!is.na(lead_chr) & !is.na(lead_bp) & !is.na(pos)]
        blk <- fread(bdf)
        blk[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), start = as.integer(start), 
            end = as.integer(end), n = as.integer(n), size_bp = as.integer(size_bp))]
        blk <- blk[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0]
        if (!nrow(blk)) 
            return(invisible(NULL))
        blk <- blk[1]
        rk_raw <- risk[trait == tr0 & lead_snp == snp0, risk_allele]
        rk_raw <- if (length(rk_raw)) 
            rk_raw[1]
        else NA_character_
        p <- ils_p(as.numeric(blk$size_bp))
        region <- data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
            core_start = blk$start, core_end = blk$end, n_ld_snp = blk$n, core_size_bp = blk$size_bp, 
            p_ils = p, best_lineage = NA_character_, matched_archaics = "")
        site <- data.table(trait = tr0, id = id0, n_site_raw = 0L, n_site_called = 0L, n_site_keep = 0L, 
            n_hap_raw = 0L, n_hap_keep = 0L, best_lineage = NA_character_, matched_archaics = "")
        d0 <- file.path(mat0, tr0, id0)
        kgf <- file.path(d0, "kg.tsv")
        afs <- arch_files(d0)
        if (!file.exists(kgf) || !length(afs)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        kg <- read_mat(kgf, "kg")
        if (is.null(kg)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        ma <- merge_arch(kg, afs)
        x <- ma$x
        arch_names <- ma$arch
        if (!length(arch_names)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        site[, `:=`(n_site_raw = nrow(x), n_site_called = nrow(x))]
        core_pos <- ld[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0, sort(unique(pos))]
        xcore <- x[pos %in% core_pos]
        setorder(xcore, pos)
        if (!nrow(xcore)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        samp_core <- setdiff(names(xcore), c("chr", "pos", "ref", "alt", "aa", arch_names))
        Hcore <- do.call(cbind, c(lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$a, 
            xcore$ref, xcore$alt)), lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$b, 
            xcore$ref, xcore$alt))))
        colnames(Hcore) <- c(paste0(samp_core, "_1"), paste0(samp_core, "_2"))
        lead_i <- which(xcore$pos == bp0)
        rk <- if (length(lead_i) == 1L) 
            harmonize_lead(rk_raw, xcore$ref[lead_i], xcore$alt[lead_i])
        else NA_character_
        carry <- if (length(lead_i) == 1L && !is.na(rk)) 
            Hcore[lead_i, ] == rk
        else rep(FALSE, ncol(Hcore))
        core <- copy(xcore)[, .(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
            pos, ref, alt)]
        core[, `:=`(lead_risk_raw = rk_raw, lead_risk_harmonized = rk, n_lead_risk_haps = sum(carry, 
            na.rm = TRUE), risk_core_allele = NA_character_, carry_freq = NA_real_, noncarry_freq = NA_real_)]
        if (sum(carry, na.rm = TRUE) > 0L) {
            rka <- vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z)) 
                  NA_character_
                else names(sort(table(z), decreasing = TRUE))[1]
            }, character(1))
            core[, `:=`(risk_core_allele, rka)]
            core[, `:=`(carry_freq, vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z) || is.na(rka[i])) 
                  NA_real_
                else mean(z == rka[i])
            }, numeric(1)))]
            core[, `:=`(noncarry_freq, vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, !carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z) || is.na(rka[i])) 
                  NA_real_
                else mean(z == rka[i])
            }, numeric(1)))]
        }
        arch_stat <- rbindlist(lapply(arch_names, function(an) {
            ca <- xcore[[an]]
            ok1 <- ca %chin% base_set
            ok2 <- ok1 & !is.na(core$risk_core_allele)
            z1 <- core[ok1]
            z2 <- core[ok2]
            a2 <- ca[ok2]
            data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
                archaic = an, lineage_group = lineage_group(an), n_core_ld_snp = length(core_pos), 
                n_risk_defined = sum(!is.na(core$risk_core_allele)), n_callable = nrow(z1), n_compared_risk = nrow(z2), 
                n_match_risk = sum(a2 == z2$risk_core_allele, na.rm = TRUE), n_match_ref = sum(ca[ok1] == 
                  z1$ref, na.rm = TRUE), n_match_alt = sum(ca[ok1] == z1$alt, na.rm = TRUE), 
                prop_match_risk = fifelse(nrow(z2) > 0L, sum(a2 == z2$risk_core_allele, na.rm = TRUE)/nrow(z2), 
                  NA_real_))
        }), fill = TRUE)
        cls <- pick_lineage(arch_stat, p)
        keep_arch <- cls$keep_arch
        best_lineage_val <- cls$best_lineage
        matched_archaics_val <- paste(keep_arch, collapse = ";")
        region[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        site[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        write_one(core, file.path(od, paste0(id0, ".core.tsv")))
        write_one(arch_stat, file.path(od, paste0(id0, ".arch.tsv")))
        if (!length(keep_arch)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        x <- x[nchar(ref) == 1L & nchar(alt) == 1L & ref %chin% base_set & alt %chin% base_set]
        samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", arch_names))
        H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), 
            lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
        colnames(H) <- c(paste0(samp, "_1"), paste0(samp, "_2"))
        keep_modern <- apply(H, 1, function(z) {
            z <- z[!is.na(z) & z %chin% base_set]
            if (!length(z)) 
                return(FALSE)
            tab <- table(z)
            length(tab) >= 2L && min(tab) >= 2L
        })
        keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set))
        keep <- keep_modern & keep_arch_called
        x2 <- x[keep]
        H2 <- H[keep, , drop = FALSE]
        site[, `:=`(n_site_keep, nrow(x2))]
        if (!nrow(x2)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        hap <- data.table(copy = colnames(H2), seq = apply(H2, 2, paste0, collapse = ""))[, .(n = .N, 
            copies = paste(copy, collapse = ";")), by = seq][order(-n, seq)]
        site[, `:=`(n_hap_raw, nrow(hap))]
        hap <- hap[n >= 1L][order(-n, seq)]
        if (!nrow(hap)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        hap[, `:=`(hap_id, roman_n(.I))]
        site[, `:=`(n_hap_keep, nrow(hap))]
        for (an in keep_arch) {
            aseq <- paste0(x2[[an]], collapse = "")
            hap[, `:=`((paste0(an, "_match")), vapply(seq, function(sq) score_seq(strsplit(sq, 
                "", fixed = TRUE)[[1]], strsplit(aseq, "", fixed = TRUE)[[1]])["match"], numeric(1)))]
        }
        cols <- paste0(keep_arch, "_match")
        hap[, `:=`(best_arch, keep_arch[max.col(.SD, ties.method = "first")]), .SDcols = cols]
        hap[, `:=`(best_match, apply(.SD, 1, max)), .SDcols = cols]
        risk_i <- which(x2$pos == bp0)
        if (length(risk_i) == 1L) {
            risk_a <- harmonize_lead(rk_raw, x2$ref[risk_i], x2$alt[risk_i])
            hap[, `:=`(carry_risk, if (!is.na(risk_a)) 
                substring(seq, risk_i, risk_i) == risk_a
            else NA)]
        }
        else hap[, `:=`(carry_risk, NA)]
        hap[, `:=`(trait = tr0, id = id0, best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        setcolorder(hap, c("trait", "id", "hap_id", "n", "best_lineage", "matched_archaics", 
            "best_arch", "best_match", "carry_risk", "seq", "copies", cols))
        write_one(hap, file.path(od, paste0(id0, ".hap.tsv")))
        write_one(site, file.path(od, paste0(id0, ".site.tsv")))
        write_one(region, file.path(od, paste0(id0, ".region.tsv")))
    }
    merge_all <- function() {
        gather <- function(pat) {
            fs <- list.files(hap0, pattern = pat, recursive = TRUE, full.names = TRUE)
            fs <- fs[file.info(fs)$size > 0]
            if (!length(fs)) 
                return(data.table())
            rbindlist(lapply(fs, fread), fill = TRUE)
        }
        hap_dt <- gather("\\.hap\\.tsv$")
        site_dt <- gather("\\.site\\.tsv$")
        arch_dt <- gather("\\.arch\\.tsv$")
        region_dt <- gather("\\.region\\.tsv$")
        core_dt <- gather("\\.core\\.tsv$")
        out_files <- c("hap_match.csv", "site_count.csv", "archaic_match.csv", "region_summary.csv", 
            "risk_core.csv", "inherited_segments.csv", "hap_match.tsv", "site_count.tsv", "archaic_match.tsv", 
            "region_summary.tsv", "risk_core.tsv", "inherited_segments.tsv", "hap_site_count.tsv", 
            "core_archaic_match.tsv", "core_risk.tsv")
        unlink(c(file.path(root, out_files), file.path(report0, out_files)), force = TRUE)
        write_report <- function(dt, stem) {
            if (!nrow(dt)) 
                return(invisible(NULL))
            fwrite(dt, file.path(report0, paste0(stem, ".tsv")), sep = "\t")
        }
        write_report(hap_dt, "hap_match")
        write_report(site_dt, "site_count")
        write_report(arch_dt, "archaic_match")
        write_report(region_dt, "region_summary")
        write_report(core_dt, "risk_core")
        if (nrow(region_dt)) {
            inherited <- region_dt[!is.na(best_lineage) & nzchar(best_lineage) & nzchar(matched_archaics)]
            write_report(inherited, "inherited_segments")
        }
    }
    run_trait <- function(tr0) {
        td <- file.path(mat0, tr0)
        if (!dir.exists(td)) 
            stop("trait not found: ", tr0)
        for (id0 in list.dirs(td, recursive = FALSE, full.names = FALSE)) {
            clean_one(tr0, id0)
            process_one(tr0, id0)
        }
    }
    run_all <- function() {
        for (tr0 in list.dirs(mat0, recursive = FALSE, full.names = FALSE)) run_trait(tr0)
        merge_all()
    }
    if (!length(args)) 
        run_all()
    else if (length(args) == 1L && args[1] == "merge") 
        merge_all()
    else if (length(args) == 1L) 
        run_trait(args[1])
    else if (length(args) == 2L) {
        clean_one(args[1], args[2])
        process_one(args[1], args[2])
    }
    else stop("usage: Rscript locus.R make_hap <res_dir> [trait | trait id | merge]")
}

run_filter_hap <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table, writexl)
    args <- commandArgs(TRUE)
    root <- normalizePath(if (length(args) >= 1L) 
        args[1]
    else "/mnt/d/analysis/gu/locus", winslash = "/", mustWork = FALSE)
    sample_file <- normalizePath(if (length(args) >= 2L) 
        args[2]
    else "/mnt/d/files/1kg.v3.sample.txt", winslash = "/", mustWork = FALSE)
    target_pop <- if (length(args) >= 3L) 
        args[3]
    else Sys.getenv("FILTER_POP", "YRI")
    max_count <- as.integer(if (length(args) >= 4L) args[4] else Sys.getenv("FILTER_MAX_COUNT", 
        "1"))
    if (is.na(max_count)) 
        stop("FILTER_MAX_COUNT must be an integer")
    report0 <- file.path(root, "report")
    mat0 <- file.path(root, "mat")
    lead0 <- file.path(root, "lead")
    dir.create(report0, recursive = TRUE, showWarnings = FALSE)
    hapf <- file.path(report0, "hap_match.tsv")
    segf <- file.path(report0, "inherited_segments.tsv")
    if (!file.exists(hapf)) 
        stop("missing hap_match.tsv: ", hapf)
    if (!file.exists(segf)) 
        stop("missing inherited_segments.tsv: ", segf)
    if (!file.exists(sample_file)) 
        stop("missing sample file: ", sample_file)
    hap <- fread(hapf)
    seg <- fread(segf)
    sample_info <- fread(sample_file)
    if (!all(c("sample", "pop") %in% names(sample_info))) 
        stop("sample file must contain columns: sample,pop")
    if (!"super_pop" %in% names(sample_info)) 
        sample_info[, `:=`(super_pop, NA_character_)]
    sample_info[, `:=`(pop = as.character(pop), super_pop = as.character(super_pop))]
    pop_levels <- sort(unique(na.omit(sample_info$pop)))
    super_levels <- sort(unique(na.omit(sample_info$super_pop)))
    pop_cols <- paste0("pop_", pop_levels)
    super_cols <- paste0("super_", super_levels)
    main_hap <- hap[n > 10]
    main_hap <- main_hap[trait %chin% seg$trait & id %chin% seg$id]
    filter_scope <- "n_gt_10"
    if ("carry_risk" %in% names(main_hap)) {
        main_hap[, `:=`(carry_risk, as.logical(carry_risk))]
        filter_hap <- main_hap[carry_risk == TRUE]
        filter_scope <- "carry_risk_TRUE_and_n_gt_10"
    }
    else {
        filter_hap <- main_hap
    }
    sample_cache <- new.env(parent = emptyenv())
    get_sample_map <- function(trait, id) {
        key <- paste(trait, id, sep = "\r")
        if (exists(key, envir = sample_cache, inherits = FALSE)) 
            return(get(key, envir = sample_cache))
        f <- file.path(mat0, trait, id, "kg.samples.tsv")
        if (!file.exists(f)) 
            stop("missing kg sample list: ", f)
        x <- fread(f, header = FALSE, col.names = "sample")
        x[, `:=`(sidx, paste0("s", .I))]
        x <- merge(x, sample_info[, .(sample, pop, super_pop)], by = "sample", all.x = TRUE, 
            sort = FALSE)
        assign(key, x, envir = sample_cache)
        x
    }
    copy_table <- function(trait, id, hap_id, n, copies) {
        cp <- unlist(strsplit(copies, ";", fixed = TRUE), use.names = FALSE)
        cp <- cp[nzchar(cp)]
        if (!length(cp)) {
            return(data.table(trait = trait, id = id, hap_id = hap_id, hap_n = n, copy = character(), 
                sidx = character(), sample = character(), pop = character(), super_pop = character()))
        }
        x <- data.table(copy = cp, sidx = sub("_[12]$", "", cp))
        x <- merge(x, get_sample_map(trait, id)[, .(sidx, sample, pop, super_pop)], by = "sidx", 
            all.x = TRUE, sort = FALSE)
        x[, `:=`(trait = trait, id = id, hap_id = hap_id, hap_n = n)]
        setcolorder(x, c("trait", "id", "hap_id", "hap_n", "copy", "sidx", "sample", "pop", "super_pop"))
        x
    }
    copies_long <- if (nrow(filter_hap)) {
        rbindlist(lapply(seq_len(nrow(filter_hap)), function(i) {
            copy_table(filter_hap$trait[i], filter_hap$id[i], filter_hap$hap_id[i], filter_hap$n[i], 
                filter_hap$copies[i])
        }), fill = TRUE)
    }
    else data.table(trait = character(), id = character(), hap_id = character(), hap_n = integer(), 
        copy = character(), sidx = character(), sample = character(), pop = character(), super_pop = character())
    wide_counts <- function(x, by_cols) {
        base <- unique(seg[, .(trait, id)])
        if (!nrow(x)) {
            for (cc in c(pop_cols, super_cols)) base[, `:=`((cc), 0L)]
            return(base)
        }
        p <- dcast(x[!is.na(pop), .N, by = c(by_cols, "pop")], paste(paste(by_cols, collapse = " + "), 
            "~ pop"), value.var = "N", fill = 0)
        s <- dcast(x[!is.na(super_pop), .N, by = c(by_cols, "super_pop")], paste(paste(by_cols, 
            collapse = " + "), "~ super_pop"), value.var = "N", fill = 0)
        for (z in list(p, s)) {
            if (!nrow(z)) 
                next
        }
        if (nrow(p)) 
            setnames(p, setdiff(names(p), by_cols), paste0("pop_", setdiff(names(p), by_cols)))
        if (nrow(s)) 
            setnames(s, setdiff(names(s), by_cols), paste0("super_", setdiff(names(s), by_cols)))
        out <- Reduce(function(a, b) merge(a, b, by = by_cols, all = TRUE, sort = FALSE), Filter(nrow, 
            list(p, s)))
        if (!length(out)) 
            out <- unique(x[, ..by_cols])
        for (cc in c(pop_cols, super_cols)) if (!cc %in% names(out)) 
            out[, `:=`((cc), 0L)]
        for (cc in c(pop_cols, super_cols)) set(out, which(is.na(out[[cc]])), cc, 0L)
        out
    }
    locus_counts <- wide_counts(copies_long, c("trait", "id"))
    hap_counts <- if (nrow(copies_long)) 
        wide_counts(copies_long, c("trait", "id", "hap_id"))
    else data.table()
    locus_meta <- filter_hap[, .(n_filter_hap = .N, n_filter_copies = sum(n)), by = .(trait, 
        id)]
    locus_summary <- merge(seg, locus_meta, by = c("trait", "id"), all.x = TRUE, sort = FALSE)
    locus_summary <- merge(locus_summary, locus_counts, by = c("trait", "id"), all.x = TRUE, 
        sort = FALSE)
    locus_summary[is.na(n_filter_hap), `:=`(n_filter_hap, 0L)]
    locus_summary[is.na(n_filter_copies), `:=`(n_filter_copies, 0L)]
    for (cc in c(pop_cols, super_cols)) {
        if (!cc %in% names(locus_summary)) 
            locus_summary[, `:=`((cc), 0L)]
        set(locus_summary, which(is.na(locus_summary[[cc]])), cc, 0L)
    }
    target_col <- paste0("pop_", target_pop)
    if (!target_col %in% names(locus_summary)) 
        locus_summary[, `:=`((target_col), 0L)]
    locus_summary[, `:=`(keep, n_filter_hap > 0L & get(target_col) <= max_count)]
    front <- c("trait", "id", "lead_chr", "lead_snp", "lead_bp", "core_start", "core_end", "n_ld_snp", 
        "core_size_bp", "p_ils", "best_lineage", "matched_archaics", "n_filter_hap", "n_filter_copies", 
        target_col, "keep")
    setcolorder(locus_summary, c(front[front %in% names(locus_summary)], setdiff(names(locus_summary), 
        front)))
    if (nrow(hap_counts)) {
        hap_counts <- merge(filter_hap[, setdiff(names(filter_hap), "copies"), with = FALSE], 
            hap_counts, by = c("trait", "id", "hap_id"), all.x = TRUE, sort = FALSE)
        for (cc in c(pop_cols, super_cols)) {
            if (!cc %in% names(hap_counts)) 
                hap_counts[, `:=`((cc), 0L)]
            set(hap_counts, which(is.na(hap_counts[[cc]])), cc, 0L)
        }
    }
    filtered_loci <- locus_summary[keep == TRUE]
    filtered_keys <- filtered_loci[, .(trait, id)]
    filtered_hap <- merge(filter_hap, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    filtered_seg <- merge(seg, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    criteria <- data.table(filter_pop = target_pop, max_count = max_count, filter_scope = filter_scope, 
        input_inherited_loci = nrow(seg), input_filter_hap = nrow(filter_hap), kept_hap = nrow(filtered_hap), 
        kept_loci = nrow(filtered_seg))
    read_optional <- function(stem) {
        f <- file.path(report0, paste0(stem, ".tsv"))
        if (file.exists(f) && file.info(f)$size > 0) 
            fread(f)
        else data.table()
    }
    filter_dt <- function(x) {
        if (!nrow(x) || !all(c("trait", "id") %in% names(x)) || !nrow(filtered_keys)) 
            return(x[0])
        merge(x, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    }
    write_book <- function(sheets, path) {
        sheets <- sheets[vapply(sheets, function(x) is.data.frame(x) && nrow(x) >= 0, logical(1))]
        write_xlsx(sheets, path)
    }
    all_sheets <- list(inherited_loci_counts = locus_summary, inherited_haplotype_counts = hap_counts, 
        inherited_segments = seg, archaic_match = read_optional("archaic_match"), risk_core = read_optional("risk_core"), 
        site_count = read_optional("site_count"), region_summary = read_optional("region_summary"), 
        filter_criteria = criteria)
    filtered_sheets <- list(filtered_loci = filtered_loci, filtered_haplotypes = hap_counts[filtered_keys, 
        on = .(trait, id), nomatch = 0], filtered_inherited_segments = filtered_seg, filtered_archaic_match = filter_dt(read_optional("archaic_match")), 
        filtered_risk_core = filter_dt(read_optional("risk_core")), filtered_site_count = filter_dt(read_optional("site_count")), 
        filter_criteria = criteria)
    selected_sheets <- list(selected_loci = if (file.exists(file.path(lead0, "pick.tsv"))) fread(file.path(lead0, 
        "pick.tsv")) else data.table(), single_snp_loci = if (file.exists(file.path(lead0, "pick.single.tsv"))) fread(file.path(lead0, 
        "pick.single.tsv")) else data.table(), lead_assoc = if (file.exists(file.path(lead0, 
        "lead.assoc"))) fread(file.path(lead0, "lead.assoc")) else data.table())
    fwrite(locus_summary, file.path(report0, "hap_pop_counts.tsv"), sep = "\t")
    fwrite(filtered_hap, file.path(report0, "filtered_hap_match.tsv"), sep = "\t")
    fwrite(filtered_seg, file.path(report0, "filtered_inherited_segments.tsv"), sep = "\t")
    write_book(all_sheets, file.path(report0, "all.xlsx"))
    write_book(filtered_sheets, file.path(report0, "filtered.xlsx"))
    write_book(selected_sheets, file.path(report0, "selected.xlsx"))
    cat(sprintf("filter_pop=%s max_count=%d input_inherited_loci=%d input_filter_hap=%d kept_hap=%d kept_loci=%d\n", 
        target_pop, max_count, nrow(seg), nrow(filter_hap), nrow(filtered_hap), nrow(filtered_seg)))
}

run_make_phy <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    library(data.table)
    args0 <- commandArgs(TRUE)
    has_root <- length(args0) && grepl("[/\\\\]", args0[1])
    root <- normalizePath(if (has_root) 
        args0[1]
    else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
    args <- if (has_root) 
        args0[-1]
    else args0
    target_trait <- if (length(args) >= 1L) 
        args[1]
    else NA_character_
    target_id <- if (length(args) >= 2L) 
        args[2]
    else NA_character_
    report0 <- file.path(root, "report")
    hapf <- file.path(report0, "filtered_hap_match.tsv")
    if (!file.exists(hapf)) 
        hapf <- file.path(report0, "hap_match.tsv")
    if (!file.exists(hapf)) 
        hapf <- file.path(root, "hap_match.tsv")
    mat0 <- file.path(root, "mat")
    phy0 <- file.path(root, "phy")
    dir.create(phy0, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(hapf) || file.info(hapf)$size == 0) 
        quit(save = "no", status = 1)
    base_set <- c("A", "C", "G", "T")
    hap <- fread(hapf)
    if (!is.na(target_trait)) 
        hap <- hap[trait == target_trait]
    if (!is.na(target_id)) 
        hap <- hap[id == target_id]
    if (!nrow(hap)) 
        quit(save = "no", status = 0)
    old_root <- if (is.na(target_trait)) 
        phy0
    else file.path(phy0, target_trait)
    oldf <- list.files(old_root, pattern = "\\.(phy|meta\\.tsv|phy_phyml_tree\\.txt|phy_phyml_stats\\.txt)$", 
        recursive = TRUE, full.names = TRUE)
    if (is.na(target_id) && length(oldf)) 
        unlink(oldf, force = TRUE)
    arch_label <- function(f) {
        x <- tools::file_path_sans_ext(basename(f))
        paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    }
    arch_files <- function(d0) {
        fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE)
        fs <- fs[!basename(fs) %in% c("kg.tsv", "kg.samples.tsv")]
        setNames(fs, arch_label(fs))
    }
    gt_hap <- function(x) {
        x <- gsub("\\|", "/", x)
        sp <- tstrsplit(x, "/", fixed = TRUE)
        list(a = sp[[1]], b = sp[[2]])
    }
    to_base <- function(h, ref, alt) {
        out <- rep("N", length(h))
        out[h == "0"] <- ref[h == "0"]
        out[h == "1"] <- alt[h == "1"]
        out
    }
    arch_base <- function(gt, ref, alt) {
        gt <- gsub("\\|", "/", gt)
        sp <- strsplit(gt, "/", fixed = TRUE)
        vapply(seq_along(sp), function(i) {
            z <- sp[[i]]
            z <- z[z != "."]
            if (length(z) != 2L || z[1] != z[2]) 
                return(NA_character_)
            a <- z[1]
            alts <- if (alt[i] %in% c(".", "")) 
                character(0)
            else strsplit(alt[i], ",", fixed = TRUE)[[1]]
            if (a == "0") 
                return(ref[i])
            ai <- suppressWarnings(as.integer(a))
            if (is.na(ai) || ai < 1L || ai > length(alts)) 
                return(NA_character_)
            alts[ai]
        }, character(1))
    }
    read_mat <- function(f, type) {
        x <- fread(f, header = FALSE)
        if (!nrow(x)) 
            return(NULL)
        if (type == "kg") 
            setnames(x, c("chr", "pos", "ref", "alt", "aa", paste0("s", seq_len(ncol(x) - 5))))
        else setnames(x, c("chr", "pos", "ref", "alt", "gt"))
        x[, `:=`(chr, sub("^chr", "", as.character(chr), ignore.case = TRUE))]
        x[, `:=`(ref = toupper(trimws(ref)), alt = toupper(trimws(alt)))]
        if ("aa" %in% names(x)) 
            x[, `:=`(aa, toupper(trimws(aa)))]
        x
    }
    read_arch <- function(f, lab) {
        x <- read_mat(f, "arch")
        if (is.null(x)) 
            return(NULL)
        x[, `:=`(allele, arch_base(gt, ref, alt))]
        x <- x[, .(chr, pos, allele)]
        setnames(x, "allele", lab)
        x
    }
    merge_arch <- function(kg, fs, keep_arch) {
        x <- kg
        keep <- character(0)
        for (lab in keep_arch) {
            if (!lab %in% names(fs)) 
                next
            a <- read_arch(fs[[lab]], lab)
            if (is.null(a)) 
                next
            x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
            keep <- c(keep, lab)
        }
        setorder(x, pos)
        list(x = x, arch = keep)
    }
    clean_phy <- function(tr0, id0) {
        od <- file.path(phy0, tr0)
        dir.create(od, TRUE, FALSE)
        unlink(file.path(od, paste0(id0, c(".full.phy", ".full.meta.tsv", ".full.phy_phyml_tree.txt", 
            ".full.phy_phyml_stats.txt", ".main.phy", ".main.meta.tsv", ".main.phy_phyml_tree.txt", 
            ".main.phy_phyml_stats.txt"))), force = TRUE)
    }
    keys <- unique(hap[, .(trait, id, matched_archaics)])
    for (i in seq_len(nrow(keys))) {
        tr0 <- keys$trait[i]
        id0 <- keys$id[i]
        clean_phy(tr0, id0)
        keep_arch <- unlist(strsplit(keys$matched_archaics[i], ";", fixed = TRUE))
        keep_arch <- keep_arch[nzchar(keep_arch)]
        if (!length(keep_arch)) 
            next
        h <- hap[trait == tr0 & id == id0]
        if (!nrow(h)) 
            next
        d0 <- file.path(mat0, tr0, id0)
        kgf <- file.path(d0, "kg.tsv")
        fs <- arch_files(d0)
        if (!file.exists(kgf) || !length(fs)) 
            next
        kg <- read_mat(kgf, "kg")
        if (is.null(kg)) 
            next
        ma <- merge_arch(kg, fs, keep_arch)
        x <- ma$x
        keep_arch <- ma$arch
        if (!length(keep_arch)) 
            next
        samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", keep_arch))
        H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), 
            lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
        keep_modern <- apply(H, 1, function(z) {
            z <- z[!is.na(z) & z %chin% base_set]
            if (!length(z)) 
                return(FALSE)
            tab <- table(z)
            length(tab) >= 2L && min(tab) >= 2L
        })
        keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set))
        x2 <- x[keep_modern & keep_arch_called]
        if (!nrow(x2)) 
            next
        anc <- x2$aa
        anc[!anc %in% base_set] <- "N"
        od <- file.path(phy0, tr0)
        dir.create(od, TRUE, FALSE)
        make_one <- function(sub, tag) {
            if (!nrow(sub)) 
                return(invisible(NULL))
            seqs <- setNames(sub$seq, sub$hap_id)
            for (an in keep_arch) {
                z <- x2[[an]]
                z[!z %in% base_set] <- "N"
                seqs <- c(seqs, setNames(paste0(z, collapse = ""), an))
            }
            seqs <- c(seqs, Ancestral = paste0(anc, collapse = ""))
            len <- unique(nchar(seqs))
            if (length(len) != 1L) 
                return(invisible(NULL))
            lab0 <- names(seqs)
            lab12 <- substr(lab0, 1, 12)
            phyf <- file.path(od, paste0(id0, ".", tag, ".phy"))
            con <- file(phyf, "w")
            writeLines(sprintf("%d %d", length(seqs), len), con)
            for (j in seq_along(seqs)) writeLines(sprintf("%-12s%s", lab12[j], seqs[j]), con)
            close(con)
            meta <- data.table(label = lab12, original_label = lab0, type = fifelse(lab0 %in% 
                keep_arch, "archaic", fifelse(lab0 == "Ancestral", "ancestral", "modern")))
            meta <- merge(meta, sub[, .(original_label = hap_id, n, best_lineage, best_arch, 
                best_match)], by = "original_label", all.x = TRUE)
            fwrite(meta, file.path(od, paste0(id0, ".", tag, ".meta.tsv")), sep = "\t")
        }
        make_one(h, "main")
    }
}

run_make_tree <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table, ape, eoffice)
    args <- commandArgs(TRUE)
    root <- normalizePath(if (length(args)) 
        args[1]
    else Sys.getenv("BALD_RES", "D:/bald/res"), winslash = "/", mustWork = FALSE)
    phy0 <- file.path(root, "phy")
    out0 <- file.path(root, "plot")
    dir.create(out0, recursive = TRUE, showWarnings = FALSE)
    keep_existing <- Sys.getenv("MAKE_TREE_KEEP_EXISTING", "0") == "1"
    write_pptx <- Sys.getenv("MAKE_TREE_PPTX", "1") != "0"
    if (!keep_existing) 
        unlink(list.files(out0, pattern = "^s8_tree_.*\\.png$", full.names = TRUE), force = TRUE)
    ppt_full <- file.path(out0, "s8_tree_full.pptx")
    ppt_main <- file.path(out0, "s8_tree_main.pptx")
    if (write_pptx && !keep_existing && file.exists(ppt_full)) 
        file.remove(ppt_full)
    if (write_pptx && !keep_existing && file.exists(ppt_main)) 
        file.remove(ppt_main)
    min_boot <- 70
    scale_bar_len <- 0.005
    arch_col <- "#1f78b4"
    anc_col <- "#33a02c"
    hi_col <- "#d73027"
    style <- function(tag) {
        if (tag == "main") {
            list(dot_r = 0.0035, cex_tip = 0.36, cex_boot = 0.36, lwd_tree = 0.45, lwd_guide = 0.22, 
                r_lab = 1.1, r_text = 1.18, r_lim = 1.3, ppt_w = 3.6, ppt_h = 3.6)
        }
        else {
            list(dot_r = 0.0028, cex_tip = 0.27, cex_boot = 0.28, lwd_tree = 0.34, lwd_guide = 0.16, 
                r_lab = 1.08, r_text = 1.15, r_lim = 1.32, ppt_w = 4, ppt_h = 3.8)
        }
    }
    topptx_safe <- function(p, filename, width, height, append = FALSE) {
        if ("append" %in% names(formals(topptx))) {
            topptx(p, filename = filename, width = width, height = height, append = append)
        }
        else {
            if (append) 
                warning("topptx() has no append argument; output may be overwritten")
            topptx(p, filename = filename, width = width, height = height)
        }
    }
    read_phyml_tree <- function(f) {
        x <- readLines(f, warn = FALSE)
        x <- x[nzchar(trimws(x))]
        i <- grep("^\\(", x)
        if (!length(i)) 
            stop("no tree found: ", f)
        read.tree(text = x[i[1]])
    }
    desc_tips <- function(tr, node) {
        nt <- Ntip(tr)
        kid <- tr$edge[tr$edge[, 1] == node, 2]
        unlist(lapply(kid, function(k) if (k <= nt) 
            k
        else desc_tips(tr, k)), use.names = FALSE)
    }
    parent_node <- function(tr, node) {
        x <- tr$edge[tr$edge[, 2] == node, 1]
        if (length(x)) 
            x[1]
        else NA_integer_
    }
    sector_range <- function(a) {
        a <- sort((a + 2 * pi)%%(2 * pi))
        if (length(a) <= 1) 
            return(c(a, a))
        g <- c(diff(a), a[1] + 2 * pi - a[length(a)])
        j <- which.max(g)
        s <- a[(j%%length(a)) + 1]
        e <- a[j]
        if (e < s) 
            e <- e + 2 * pi
        c(s, e)
    }
    tip_type <- function(meta, labs) {
        out <- rep(NA_character_, length(labs))
        names(out) <- labs
        if (!is.null(meta) && all(c("label", "type") %in% names(meta))) {
            idx <- match(labs, meta$label)
            out[!is.na(idx)] <- as.character(meta$type[idx[!is.na(idx)]])
        }
        out[grepl("vindija|altai|chagyr|denisova|neand|archaic", labs, ignore.case = TRUE)] <- "archaic"
        out[labs == "Ancestral"] <- "ancestral"
        out[is.na(out)] <- "modern"
        out
    }
    highlight_node <- function(tr, meta) {
        nt <- Ntip(tr)
        labs <- tr$tip.label
        tp <- tip_type(meta, labs)
        arch <- labs[tp == "archaic"]
        modern <- labs[tp == "modern"]
        arch <- arch[!is.na(arch) & nzchar(arch)]
        modern <- modern[!is.na(modern) & nzchar(modern)]
        if (!length(arch) || !length(modern)) 
            return(NA_integer_)
        candidates <- function(all_arch = TRUE) {
            rbindlist(lapply((nt + 1):(nt + tr$Nnode), function(node) {
                tips <- labs[desc_tips(tr, node)]
                tips <- tips[!is.na(tips) & nzchar(tips)]
                ok_arch <- if (all_arch) 
                  isTRUE(all(arch %in% tips))
                else isTRUE(any(tips %in% arch))
                has_modern <- isTRUE(any(tips %in% modern))
                has_ancestral <- isTRUE(any(tips == "Ancestral"))
                support <- suppressWarnings(as.numeric(tr$node.label[node - nt]))
                if (length(support) != 1L) 
                  support <- NA_real_
                if (!ok_arch || !has_modern || has_ancestral || is.na(support) || support < min_boot) 
                  return(NULL)
                data.table(node = node, support = support, n_tips = length(tips), n_mod = sum(tips %in% 
                  modern))
            }), fill = TRUE)
        }
        z <- candidates(TRUE)
        if (!nrow(z)) 
            z <- candidates(FALSE)
        if (!nrow(z)) 
            return(NA_integer_)
        setorder(z, n_tips, -support, -n_mod)
        z$node[1]
    }
    draw_tree <- function(tr, meta = NULL, tag = c("full", "main")) {
        tag <- match.arg(tag)
        st <- style(tag)
        tr$tip.label <- trimws(tr$tip.label)
        if ("Ancestral" %in% tr$tip.label) {
            tr <- tryCatch(root(tr, outgroup = "Ancestral", resolve.root = TRUE), error = function(e) tr)
        }
        tr <- ladderize(tr)
        nt <- Ntip(tr)
        ia <- which(tr$tip.label == "Ancestral")
        dep <- node.depth.edgelength(tr)[1:nt]
        r_tree <- max(dep[setdiff(seq_len(nt), ia)], na.rm = TRUE)
        if (!is.finite(r_tree)) 
            r_tree <- max(dep, na.rm = TRUE)
        r_lab <- r_tree * st$r_lab
        r_text <- r_tree * st$r_text
        r_lim <- r_tree * st$r_lim
        par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, pty = "s")
        plot.phylo(tr, type = "fan", use.edge.length = TRUE, show.tip.label = FALSE, no.margin = TRUE, 
            edge.width = st$lwd_tree, x.lim = c(-r_lim, r_lim), y.lim = c(-r_lim, r_lim))
        pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx <- pp$xx
        yy <- pp$yy
        xt <- xx[1:nt]
        yt <- yy[1:nt]
        ang <- atan2(yt, xt)
        labs <- tr$tip.label
        tp <- tip_type(meta, labs)
        hi <- highlight_node(tr, meta)
        if (!is.na(hi)) {
            idx <- desc_tips(tr, hi)
            aa <- sector_range(atan2(yy[idx], xx[idx]))
            th <- seq(aa[1], aa[2], length.out = 500)
            r0 <- sqrt(xx[hi]^2 + yy[hi]^2)
            polygon(c(r0 * cos(th), rev(r_lab * cos(th))), c(r0 * sin(th), rev(r_lab * sin(th))), 
                col = adjustcolor(hi_col, alpha.f = 0.22), border = NA)
        }
        segments(xt, yt, r_lab * cos(ang), r_lab * sin(ang), lty = 3, col = "grey45", lwd = st$lwd_guide)
        symbols(xt, yt, circles = rep(r_tree * st$dot_r, nt), inches = FALSE, add = TRUE, bg = "black", 
            fg = "black", lwd = 0.15)
        lab_col <- rep("black", nt)
        lab_col[tp == "archaic"] <- arch_col
        lab_col[tp == "ancestral"] <- anc_col
        deg <- ang * 180/pi
        flip <- deg < -90 | deg > 90
        srt <- ifelse(flip, deg + 180, deg)
        for (i in seq_len(nt)) {
            text(r_text * cos(ang[i]), r_text * sin(ang[i]), labs[i], srt = srt[i], adj = if (flip[i]) 
                c(1, 0.5)
            else c(0, 0.5), cex = st$cex_tip, col = lab_col[i])
        }
        if (!is.na(hi)) {
            boot <- suppressWarnings(as.numeric(tr$node.label))
            nd <- nt + seq_len(tr$Nnode)
            show <- unique(c(parent_node(tr, hi), hi))
            idx <- match(show[show > nt], nd)
            idx <- idx[!is.na(idx) & !is.na(boot[idx])]
            if (length(idx)) 
                nodelabels(boot[idx], node = nd[idx], frame = "n", cex = st$cex_boot, col = "black")
        }
        usr <- par("usr")
        add.scale.bar(x = usr[2] - 0.14 * diff(usr[1:2]), y = usr[3] + 0.06 * diff(usr[3:4]), 
            length = scale_bar_len, lwd = 0.8, cex = 0.65)
    }
    assign("draw_tree", draw_tree, envir = .GlobalEnv)
    make_editable_plot <- function(tr, meta, tag) {
        assign(".locus_tree_tr", tr, envir = .GlobalEnv)
        assign(".locus_tree_meta", meta, envir = .GlobalEnv)
        assign(".locus_tree_tag", tag, envir = .GlobalEnv)
        assign(".locus_tree_draw_once", function() {
            draw_tree(get(".locus_tree_tr", envir = .GlobalEnv), get(".locus_tree_meta", envir = .GlobalEnv), 
                get(".locus_tree_tag", envir = .GlobalEnv))
        }, envir = .GlobalEnv)
        convertplot(.locus_tree_draw_once())
    }
    files <- list.files(phy0, pattern = "\\.main\\.phy_phyml_tree\\.txt$", recursive = TRUE, 
        full.names = TRUE)
    if (!length(files)) 
        quit(save = "no", status = 0)
    n_ppt <- c(full = 0L, main = 0L)
    for (treef in files) {
        base <- sub("\\.phy_phyml_tree\\.txt$", "", basename(treef))
        tag <- if (grepl("\\.main$", base)) 
            "main"
        else "full"
        st <- style(tag)
        metaf <- file.path(dirname(treef), paste0(base, ".meta.tsv"))
        meta <- if (file.exists(metaf)) 
            fread(metaf)
        else NULL
        tr <- tryCatch(read_phyml_tree(treef), error = function(e) NULL)
        if (is.null(tr)) 
            next
        trait <- basename(dirname(treef))
        pngf <- file.path(out0, paste0("s8_tree_", tag, "_", trait, "_", base, ".png"))
        if (!keep_existing || !file.exists(pngf) || file.info(pngf)$size == 0) {
            png(pngf, width = st$ppt_w, height = st$ppt_h, units = "in", res = 300)
            draw_tree(tr, meta, tag)
            dev.off()
        }
        p <- make_editable_plot(tr, meta, tag)
        if (write_pptx) {
            p <- make_editable_plot(tr, meta, tag)
            topptx_safe(p, if (tag == "main") 
                ppt_main
            else ppt_full, width = st$ppt_w, height = st$ppt_h, append = keep_existing || n_ppt[tag] > 
                0L)
        }
        n_ppt[tag] <- n_ppt[tag] + 1L
    }
}

args <- commandArgs(TRUE)

if (length(args) == 0) usage()

cmd <- args[1]

cmd_args <- args[-1]

switch(cmd, prep_input = run_prep_input(cmd_args), prep_archaic = run_prep_archaic(cmd_args), 
    make_hap = run_make_hap(cmd_args), filter_hap = run_filter_hap(cmd_args), make_phy = run_make_phy(cmd_args), 
    make_tree = run_make_tree(cmd_args), usage())
