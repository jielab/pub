library(data.table)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(keys, multi = FALSE, default = NULL){
	i <- which(args %in% keys)[1]; if(is.na(i) || i == length(args)) return(default)
	j <- i + 1; if(!multi) return(args[j]); k <- j
	while(k <= length(args) && !grepl('^-', args[k])) k <- k + 1
	if(k == j) default else args[j:(k - 1)]
}
pick1 <- function(nm, cand){ x <- cand[cand %chin% nm]; if(length(x)) x[1] else NA_character_ }
chr_int <- function(x) suppressWarnings(as.integer(sub('^chr', '', as.character(x), ignore.case = TRUE)))
fail_reason <- function(ea, oa, b, ingwas){
	x <- c(if(!isTRUE(ingwas)) 'not_in_gwas', if(is.na(ea) || !nzchar(ea)) 'effect_allele', if(is.na(oa) || !nzchar(oa)) 'other_allele', if(is.na(b)) 'beta')
	paste(x, collapse = ';')
}
read_gwas <- function(f){
	h <- names(fread(f, nrows = 0))
	cc <- c(SNP = pick1(h, c('SNP','rsid','RSID','ID')), CHR = pick1(h, c('CHR','Chr','chr','#CHROM','CHROM')), effect_allele = pick1(h, c('EA','A1','effect_allele','ALT','alt')), other_allele = pick1(h, c('NEA','A2','other_allele','REF','ref')), beta = pick1(h, c('BETA','beta','b')))
	if(anyNA(cc)) stop('Missing required GWAS columns in: ', f, '\nFound columns: ', paste(h, collapse = ', '))
	x <- fread(f, select = unname(cc), showProgress = FALSE); setnames(x, unname(cc), names(cc))
	x[, .(SNP = as.character(SNP), CHR = chr_int(CHR), effect_allele = toupper(trimws(as.character(effect_allele))), other_allele = toupper(trimws(as.character(other_allele))), beta = as.numeric(beta), gwas_match = TRUE)]
}

dirgwas <- get_arg(c('-dirgwas','--dirgwas')); dirout <- get_arg(c('-dirout','--dirout')); traits0 <- get_arg(c('-traits','--traits'), multi = TRUE)
if(is.null(dirgwas) || is.null(dirout) || is.null(traits0)) stop('Usage: Rscript prepare_input.R --dirgwas DIR --dirout DIR --traits bald bald12 bald13 bald14')
dirgwas <- normalizePath(dirgwas, winslash = '/', mustWork = TRUE); dirout <- normalizePath(dirout, winslash = '/', mustWork = FALSE)
traits <- unlist(strsplit(paste(traits0, collapse = ' '), '[, ]+')); traits <- traits[nzchar(traits)]
dir_lead <- paste0(dirout, '/lead'); dir_risk <- paste0(dirout, '/risk'); dir.create(dir_lead, TRUE, FALSE); dir.create(dir_risk, TRUE, FALSE)

for(tr in traits){
	cojo_f <- paste0(dirgwas, '/cojo/', tr, '/', tr, '.jma.cojo'); gwas_f <- paste0(dirgwas, '/clean/', tr, '/', tr, '.gz')
	if(!file.exists(cojo_f)) stop('Missing COJO file: ', cojo_f); if(!file.exists(gwas_f)) stop('Missing GWAS file: ', gwas_f)
	cj <- fread(cojo_f); chr_col <- pick1(names(cj), c('Chr','CHR','chr','#CHROM','CHROM')); bp_col <- pick1(names(cj), c('bp','BP','POS','pos'))
	if(anyNA(c(chr_col, bp_col)) || !'SNP' %chin% names(cj)) stop('COJO file must contain Chr/CHR, SNP and bp/BP columns: ', cojo_f)
	lead <- unique(cj[, .(lead_chr = chr_int(get(chr_col)), lead_snp = as.character(SNP), lead_bp = as.integer(get(bp_col)))])
	lead <- lead[lead_chr %between% c(1L, 22L) & !is.na(lead_bp) & nzchar(lead_snp)]; fwrite(lead, paste0(dir_lead, '/', tr, '.lead.tsv'), sep = '\t')
	g <- unique(read_gwas(gwas_f), by = c('SNP','CHR'))
	z <- merge(lead, g, by.x = c('lead_snp','lead_chr'), by.y = c('SNP','CHR'), all.x = TRUE, sort = FALSE)
	z[, fail := is.na(gwas_match) | is.na(effect_allele) | !nzchar(effect_allele) | is.na(other_allele) | !nzchar(other_allele) | is.na(beta)]
	fail_f <- paste0(dir_lead, '/', tr, '.lead.fail.tsv')
	if(any(z$fail)){
		fd <- z[fail == TRUE, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, in_gwas = !is.na(gwas_match))]
		fd[, fail_reason := mapply(fail_reason, effect_allele, other_allele, beta, in_gwas)]; fwrite(fd, fail_f, sep = '\t'); message('Warning: ', nrow(fd), ' lead SNPs failed for ', tr, '; written to ', fail_f)
	} else if(file.exists(fail_f)) file.remove(fail_f)
	z <- z[fail == FALSE][, `:=`(trait = tr, risk_allele = fifelse(beta >= 0, effect_allele, other_allele))]
	fwrite(z[, .(trait, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, risk_allele)], paste0(dir_risk, '/', tr, '.risk.tsv'), sep = '\t')
}
risk <- rbindlist(lapply(traits, function(tr) fread(paste0(dir_risk, '/', tr, '.risk.tsv'))), fill = TRUE)
fwrite(risk, paste0(dirout, '/risk.tsv'), sep = '\t'); 
message('Done: ', paste0(dirout, '/risk.tsv'), ' and ', dir_lead)