library(data.table)
fphe <- "/mnt/d/scripts/fphe.f.R"; if(file.exists(fphe)) source(fphe)
args <- commandArgs(TRUE)
arg <- function(k, multi = FALSE){i <- which(args %in% k)[1]; if(is.na(i) || i == length(args)) return(NULL); j <- i + 1; if(!multi) return(args[j]); z <- args[j:length(args)]; z[seq_len(which(c(grepl('^-', z[-1]), TRUE))[1])]} 
pick1 <- function(nm, z){x <- z[z %chin% nm]; if(length(x)) x[1] else NA_character_}
chr_int <- function(x){x <- toupper(sub('^CHR','',as.character(x))); x[x=='X'] <- '23'; suppressWarnings(as.integer(x))}
rd_gwas <- function(f){
	h <- names(fread(f, nrows = 0)); cc <- c(SNP=pick1(h,c('SNP','rsid','RSID','ID')), EA=pick1(h,c('EA','A1','effect_allele','ALT','alt')), NEA=pick1(h,c('NEA','A2','other_allele','REF','ref')), BETA=pick1(h,c('BETA','beta','b')))
	if(anyNA(cc)) stop('Missing GWAS columns in ', f, ': ', paste(h, collapse=','))
	x <- fread(f, select = unname(cc), showProgress = FALSE); setnames(x, unname(cc), names(cc)); unique(x[, .(SNP=as.character(SNP), EA=toupper(EA), NEA=toupper(NEA), BETA=as.numeric(BETA))], by='SNP')
}

dirgwas <- arg(c('--dirgwas','-dirgwas')); dirout <- arg(c('--dirout','-dirout')); traits <- arg(c('--traits','-traits'), TRUE)
if(is.null(dirgwas) || is.null(dirout) || is.null(traits)) stop('Usage: Rscript input_prep.R --dirgwas DIR --dirout DIR --traits bald bald12 ...')
traits <- unlist(strsplit(paste(traits, collapse=' '), '[, ]+')); traits <- traits[nzchar(traits)]
lead0 <- paste0(normalizePath(dirout, winslash='/', mustWork=FALSE), '/lead'); dir.create(lead0, TRUE, FALSE)

for(tr in traits){
	cj_f <- paste0(dirgwas, '/cojo/', tr, '/', tr, '.jma.cojo'); gw_f <- paste0(dirgwas, '/clean/', tr, '/', tr, '.gz')
	cj <- fread(cj_f); cc <- c(chr=pick1(names(cj), c('Chr','CHR','chr','#CHROM','CHROM')), bp=pick1(names(cj), c('bp','BP','POS','pos')))
	if(anyNA(cc) || !'SNP' %chin% names(cj)) stop('COJO needs Chr/SNP/bp: ', cj_f)
	lead <- unique(cj[, .(lead_chr=chr_int(get(cc['chr'])), lead_snp=as.character(SNP), lead_bp=as.integer(get(cc['bp'])))])[lead_chr %between% c(1L,23L) & !is.na(lead_bp)]
	g <- rd_gwas(gw_f); z <- merge(lead, g, by.x='lead_snp', by.y='SNP', all.x=TRUE, sort=FALSE)
	z[, fail := is.na(EA) | !nzchar(EA) | is.na(NEA) | !nzchar(NEA) | is.na(BETA)]
	fwrite(lead, paste0(lead0, '/', tr, '.lead.3col'), sep='\t')
	if(any(z$fail)) fwrite(z[fail==TRUE, .(trait=tr, lead_chr, lead_snp, lead_bp, EA, NEA, BETA, fail_reason='not_in_gwas_or_missing_allele_beta')], paste0(lead0, '/', tr, '.lead.fail.tsv'), sep='\t') else unlink(paste0(lead0, '/', tr, '.lead.fail.tsv'))
	z <- z[fail==FALSE][, .(trait=tr, lead_chr, lead_snp, lead_bp, effect_allele=EA, other_allele=NEA, beta=BETA, risk_allele=fifelse(BETA >= 0, EA, NEA))]
	fwrite(z, paste0(lead0, '/', tr, '.lead.assoc'), sep='\t')
}
fwrite(rbindlist(lapply(traits, \(tr) fread(paste0(lead0, '/', tr, '.lead.assoc'))), fill=TRUE), paste0(lead0, '/lead.assoc'), sep='\t')
