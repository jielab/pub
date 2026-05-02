library(data.table)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(keys, multi = FALSE, default = NULL) {
	i <- which(args %in% keys)[1]
	if(is.na(i) || i == length(args)) return(default)
	j <- i + 1
	if(!multi) return(args[j])
	k <- j
	while(k <= length(args) && !grepl('^-', args[k])) k <- k + 1
	if(k == j) return(default)
	args[j:(k - 1)]
}

pick1 <- function(nm, cand) {
	x <- cand[cand %chin% nm]
	if(length(x)) x[1] else NA_character_
}

read_gwas <- function(f) {
	h <- names(fread(f, nrows = 0))
	cc <- c(SNP = pick1(h, c('SNP','rsid','ID')), CHR = pick1(h, c('CHR','Chr','chr')), effect_allele = pick1(h, c('EA','A1','effect_allele')), other_allele = pick1(h, c('NEA','A2','other_allele')), beta = pick1(h, c('BETA','beta','b')))
	if(anyNA(cc)) stop('Missing required GWAS columns in: ', f, '\nFound columns: ', paste(h, collapse = ', '))
	x <- fread(f, select = unname(cc), showProgress = FALSE)
	setnames(x, unname(cc), names(cc))
	x[, .(SNP = as.character(SNP), CHR = as.integer(CHR), effect_allele = toupper(as.character(effect_allele)), other_allele = toupper(as.character(other_allele)), beta = as.numeric(beta), gwas_match = TRUE)]
}

dirgwas <- get_arg(c('-dirgwas','--dirgwas'))
dirout <- get_arg(c('-dirout','--dirout'))
traits0 <- get_arg(c('-traits','--traits'), multi = TRUE)
if(is.null(dirgwas) || is.null(dirout) || is.null(traits0)) stop('Usage: Rscript 0_prepare_input.R --dirgwas DIR --dirout DIR --traits bald bald12 bald13 bald14')
dirgwas <- normalizePath(dirgwas, winslash = '/', mustWork = TRUE)
dirout <- normalizePath(dirout, winslash = '/', mustWork = FALSE)
traits <- unlist(strsplit(paste(traits0, collapse = ' '), '[, ]+'))
traits <- traits[nzchar(traits)]

dir_lead <- paste0(dirout, '/lead')
dir_risk <- paste0(dirout, '/risk')
dir.create(dir_lead, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_risk, recursive = TRUE, showWarnings = FALSE)

for(tr in traits) {
	cojo_f <- paste0(dirgwas, '/cojo/', tr, '/', tr, '.jma.cojo')
	gwas_f <- paste0(dirgwas, '/clean/', tr, '/', tr, '.gz')
	if(!file.exists(cojo_f)) stop('Missing COJO file: ', cojo_f)
	if(!file.exists(gwas_f)) stop('Missing GWAS file: ', gwas_f)

	cj <- fread(cojo_f)
	chr_col <- pick1(names(cj), c('Chr','CHR','chr'))
	bp_col <- pick1(names(cj), c('bp','BP','POS','pos'))
	if(anyNA(c(chr_col, bp_col)) || !'SNP' %chin% names(cj)) stop('COJO file must contain Chr/CHR, SNP and bp/BP columns: ', cojo_f)

	lead <- cj[, .(lead_chr = as.integer(get(chr_col)), lead_snp = as.character(SNP), lead_bp = as.integer(get(bp_col)))]
	lead <- lead[!is.na(lead_chr) & !is.na(lead_bp) & nzchar(lead_snp)]
	fwrite(lead, paste0(dir_lead, '/', tr, '.lead.tsv'), sep = '\t')

	g <- read_gwas(gwas_f)
	if(anyDuplicated(g, by = c('SNP','CHR'))) {
		g <- unique(g, by = c('SNP','CHR'))
		message('Warning: duplicated SNP+CHR rows in GWAS were de-duplicated for ', tr)
	}
	z <- merge(lead, g, by.x = c('lead_snp','lead_chr'), by.y = c('SNP','CHR'), all.x = TRUE, sort = FALSE)
	z[, fail := is.na(gwas_match) | is.na(effect_allele) | is.na(other_allele) | is.na(beta)]
	fail_f <- paste0(dir_lead, '/', tr, '.lead.fail.tsv')
	if(any(z$fail)) {
		fail_dat <- z[fail == TRUE, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, in_gwas = !is.na(gwas_match))]
		fail_dat[, fail_reason := paste(c(if(!in_gwas) 'not_in_gwas', names(.SD)[is.na(.SD) | !nzchar(as.character(.SD))]), collapse = ';'), by = .(trait, lead_chr, lead_snp, lead_bp), .SDcols = c('effect_allele','other_allele','beta')]
		fwrite(fail_dat, fail_f, sep = '\t')
		message('Warning: ', nrow(fail_dat), ' COJO lead SNPs failed for ', tr, '; written to ', fail_f)
	} else if(file.exists(fail_f)) file.remove(fail_f)
	z <- z[fail == FALSE]
	z[, `:=`(trait = tr, risk_allele = fifelse(beta >= 0, effect_allele, other_allele))]
	fwrite(z[, .(trait, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, risk_allele)], paste0(dir_risk, '/', tr, '.risk.tsv'), sep = '\t')
}

risk <- rbindlist(lapply(traits, function(tr) fread(paste0(dir_risk, '/', tr, '.risk.tsv'))), fill = TRUE)
fwrite(risk, paste0(dirout, '/risk.tsv'), sep = '\t')
message('Done: ', paste0(dirout, '/risk.tsv'), ' and ', dir_lead)
