library(data.table)

args <- commandArgs(TRUE)
get_arg <- function(keys, default = NULL){
	i <- which(args %in% keys)[1]
	if(is.na(i) || i == length(args)) return(default)
	args[i + 1]
}
pick1 <- function(nm, cand){ x <- cand[cand %in% nm]; if(length(x)) x[1] else NA_character_ }
chr_int <- function(x){
	x <- as.character(x); x <- sub('^chr', '', x, ignore.case = TRUE)
	x[toupper(x) == 'X'] <- '23'; x[toupper(x) == 'Y'] <- '24'
	suppressWarnings(as.integer(x))
}
read_gwas <- function(f){
	h <- names(fread(f, nrows = 0))
	cc <- c(
		SNP  = pick1(h, c('SNP','rsid','RSID','ID')),
		CHR  = pick1(h, c('CHR','Chr','chr','#CHROM','CHROM')),
		EA   = pick1(h, c('EA','A1','effect_allele','ALT','alt')),
		OA   = pick1(h, c('NEA','A2','other_allele','REF','ref')),
		BETA = pick1(h, c('BETA','beta','b')),
		P    = pick1(h, c('P','p','pval','PVAL','p_value','P_value'))
	)
	if(anyNA(cc)) stop('Missing required GWAS columns in: ', f, '\nFound columns: ', paste(h, collapse = ', '))
	x <- fread(f, select = unname(cc), showProgress = FALSE)
	setnames(x, unname(cc), names(cc))
	x[, .(SNP = as.character(SNP), CHR = chr_int(CHR), effect_allele = toupper(trimws(as.character(EA))),
	      other_allele = toupper(trimws(as.character(OA))), beta = as.numeric(BETA), P = as.numeric(P))]
}

res1 <- get_arg(c('--res1','-res1'))
res2 <- get_arg(c('--res2','-res2'))
dirgwas <- get_arg(c('--dirgwas','-dirgwas'))
P_cut <- as.numeric(get_arg(c('--p','-p'), '1e-6'))
nmax <- as.integer(get_arg(c('--nmax','-nmax'), '50'))
traits0 <- get_arg(c('--traits','-traits'), '')
if(any(!nzchar(c(res1, res2, dirgwas)))) stop('Usage: Rscript refine_make_lead_risk.R --res1 DIR --res2 DIR --dirgwas DIR [--p 1e-6] [--nmax 50] [--traits "bald bald12"]')

mapf <- paste0(res2, '/map/map.tsv')
if(!file.exists(mapf) || file.info(mapf)$size == 0) stop('Missing or empty map file: ', mapf)
m <- unique(fread(mapf))
if(!nrow(m)) quit(save = 'no', status = 0)
m[, `:=`(trait = as.character(trait), orig_chr = chr_int(orig_chr), POS37 = as.integer(POS37), CHR = chr_int(orig_chr))]
traits <- if(nzchar(traits0)) unlist(strsplit(traits0, '[, ]+')) else unique(m$trait)
traits <- traits[nzchar(traits)]

dir.create(paste0(res2, '/lead'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(res2, '/risk'), recursive = TRUE, showWarnings = FALSE)
all_lead <- list(); all_risk <- list(); all_map <- list()

for(tr in traits){
	mt <- m[trait == tr]
	if(!nrow(mt)) next
	gwas_f <- paste0(dirgwas, '/clean/', tr, '/', tr, '.gz')
	if(!file.exists(gwas_f)) stop('Missing GWAS file: ', gwas_f)
	g <- unique(read_gwas(gwas_f), by = c('SNP','CHR'))
	z <- merge(mt, g, by = c('SNP','CHR'), allow.cartesian = TRUE)
	ldf <- paste0(res1, '/ld/', tr, '/ld.tsv')
	if(file.exists(ldf)){
		ld <- fread(ldf)
		ld[, ld_key := paste(lead_chr, lead_snp, lead_bp, snp, sep = ':')]
		z[, ld_key := paste(orig_chr, orig_snp, orig_bp, SNP, sep = ':')]
		z[, in_firstpass_ld := ld_key %in% ld$ld_key]
	} else z[, in_firstpass_ld := FALSE]

	z <- z[(P <= P_cut | in_firstpass_ld == TRUE) &
	       nchar(effect_allele) == 1 & nchar(other_allele) == 1 &
	       effect_allele %in% c('A','C','G','T') & other_allele %in% c('A','C','G','T')]
	if(!nrow(z)) next
	setorder(z, orig_chr, orig_bp, P)
	z <- z[, head(.SD, nmax), by = .(trait, orig_chr, orig_snp, orig_bp)]
	z[, `:=`(risk_allele = fifelse(beta >= 0, effect_allele, other_allele), lead_chr = orig_chr, lead_snp = SNP, lead_bp = POS37)]
	all_lead[[tr]] <- unique(z[, .(lead_chr, lead_snp, lead_bp)])
	all_risk[[tr]] <- unique(z[, .(trait, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, risk_allele)])
	all_map[[tr]]  <- unique(z[, .(trait, orig_chr, orig_snp, orig_bp, lead_chr, lead_snp, lead_bp, P, in_firstpass_ld, risk_allele)])
	fwrite(all_lead[[tr]][order(lead_chr, lead_bp)], paste0(res2, '/lead/', tr, '.lead.tsv'), sep = '\t')
	fwrite(all_risk[[tr]][order(lead_chr, lead_bp)], paste0(res2, '/risk/', tr, '.risk.tsv'), sep = '\t')
}

if(!length(all_risk)){
	fwrite(data.table(), paste0(res2, '/candidate_map.tsv'), sep = '\t')
	message('No refined candidates found.')
	quit(save = 'no', status = 0)
}

risk <- rbindlist(all_risk, fill = TRUE)
map_out <- rbindlist(all_map, fill = TRUE)
fwrite(risk[order(trait, lead_chr, lead_bp)], paste0(res2, '/risk/risk.tsv'), sep = '\t')
fwrite(risk[order(trait, lead_chr, lead_bp)], paste0(res2, '/risk.tsv'), sep = '\t')
fwrite(map_out[order(trait, orig_chr, orig_bp, P)], paste0(res2, '/candidate_map.tsv'), sep = '\t')
message('Done: ', paste0(res2, '/lead, ', res2, '/risk, ', res2, '/candidate_map.tsv'))
