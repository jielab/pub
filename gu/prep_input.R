pacman::p_load(data.table)
fphe <- "/mnt/d/scripts/f/phe.f.R"; if(file.exists(fphe)) source(fphe)

args <- commandArgs(TRUE)
arg <- function(k, multi=FALSE){i <- which(args %in% k)[1]; if(is.na(i) || i == length(args)) return(NULL); j <- i+1; if(!multi) return(args[j]); z <- args[j:length(args)]; z[seq_len(which(c(grepl('^-', z[-1]), TRUE))[1])]}
pick1 <- function(nm, z){x <- z[z %chin% nm]; if(length(x)) x[1] else NA_character_}
chr_int <- function(x){x <- toupper(sub('^CHR','',as.character(x))); x[x=='X'] <- '23'; x[x=='Y'] <- '24'; suppressWarnings(as.integer(x))}
id_mode <- function(x){x <- unique(na.omit(as.character(x))); x <- x[nzchar(x) & x!='.']; x <- head(x, 5000); if(!length(x)) return('missing'); y <- toupper(x); rs <- mean(grepl('^RS[0-9]+$', y)); cp <- mean(grepl('^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$', y)); if(rs > .8) 'rsid' else if(cp > .8) 'chrpos' else 'other'}
usable_id <- function(x){!is.na(x) & nzchar(x) & x != '.'}
allele_score <- function(ea, nea, ref, alt){
	ea <- toupper(as.character(ea)); nea <- toupper(as.character(nea))
	ref <- toupper(as.character(ref)); alt <- toupper(as.character(alt))
	pair <- ((ea == ref & nea == alt) | (ea == alt & nea == ref))
	one <- (ea == ref | ea == alt | nea == ref | nea == alt)
	fifelse(pair, 2L, fifelse(one, 1L, 0L))
}

rd_gwas <- function(f, keep_snps=NULL){
	h <- names(fread(f, nrows=0)); cc <- c(SNP=pick1(h,c('SNP','rsid','RSID','ID')), EA=pick1(h,c('EA','A1','effect_allele','ALT','alt')), NEA=pick1(h,c('NEA','A2','other_allele','REF','ref')), BETA=pick1(h,c('BETA','beta','b')))
	if(anyNA(cc)) stop('Missing GWAS columns in ', f, ': ', paste(h, collapse=','))
	if(!is.null(keep_snps)){
		keep_snps <- unique(as.character(keep_snps[!is.na(keep_snps) & nzchar(keep_snps)]))
		if(!length(keep_snps)) return(data.table(SNP=character(), EA=character(), NEA=character(), BETA=numeric()))
		tmp <- tempfile()
		fwrite(data.table(SNP=keep_snps), tmp, col.names=FALSE, sep='\t')
		idx <- match(unname(cc), h)
		cmd <- sprintf("gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{k[$1]=1; next} FNR==1{next} ($%d in k){print $%d,$%d,$%d,$%d}' %s -",
			shQuote(f), idx[1], idx[1], idx[2], idx[3], idx[4], shQuote(tmp))
		x <- fread(cmd=cmd, col.names=names(cc), showProgress=FALSE)
		unlink(tmp)
	} else {
		x <- fread(f, select=unname(cc), showProgress=FALSE)
		setnames(x, unname(cc), names(cc))
	}
	unique(x[, .(SNP=as.character(SNP), EA=toupper(EA), NEA=toupper(NEA), BETA=as.numeric(BETA))], by='SNP')
}
pvar_files <- function(dirmod, chr){
	if(chr == 23L) unique(c(paste0(dirmod, '/EUR.male.chrX.par.pvar'), paste0(dirmod, '/EUR.male.chrX.nonPar.pvar'), paste0(dirmod, '/EUR.chrX.pvar'))) else paste0(dirmod, '/EUR.chr', chr, '.pvar')
}
pvar_ids <- function(f, n=2000){
	if(!file.exists(f)) return(character())
	cmd <- sprintf("awk '/^#/{next} $3!=\".\" && $3!=\"\"{print $3; if(++n==%d) exit}' %s", as.integer(n), shQuote(f))
	system(cmd, intern=TRUE)
}
pvar_mode <- function(fs){id_mode(unlist(lapply(fs, pvar_ids), use.names=FALSE))}

rd_pvar <- function(dirmod, lead){
	out <- list(); k <- 0L
	for(chr in sort(unique(lead$lead_chr))){
		fs <- pvar_files(dirmod, chr); fs <- fs[file.exists(fs)]
		if(!length(fs)) stop('No 1000G .pvar found for chr ', chr, ' under ', dirmod)
		key <- unique(lead[lead_chr == chr, .(lead_bp, lead_snp)])
		tmp <- tempfile(); fwrite(key, tmp, col.names=FALSE, sep='\t')
		for(f in fs){
			cmd <- sprintf("awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{p[$1]=1; id[$2]=1; next} /^##/{next} /^#CHROM/{print \"CHR\",$2,$3,$4,$5; next} FNR==1 && ($1==\"CHR\" || $1==\"CHROM\"){print \"CHR\",$2,$3,$4,$5; next} (($2 in p) || ($3 in id)){gsub(/^chr/,\"\",$1); print $1,$2,$3,$4,$5}' %s %s", shQuote(tmp), shQuote(f))
			x <- fread(cmd=cmd, header=TRUE, showProgress=FALSE)
			if(nrow(x) <= 0) next
			cc <- c(CHR=pick1(names(x), c('CHR','CHROM')), POS=pick1(names(x), c('POS','BP')), ID=pick1(names(x), c('ID','SNP')), REF=pick1(names(x), c('REF')), ALT=pick1(names(x), c('ALT')))
			if(anyNA(cc)) stop('Bad pvar columns after skipping ## lines: ', f, '; names=', paste(names(x), collapse=','))
			k <- k+1L; out[[k]] <- x[, .(CHR=chr_int(get(cc['CHR'])), POS=as.integer(get(cc['POS'])), ID=as.character(get(cc['ID'])), REF=toupper(get(cc['REF'])), ALT=toupper(get(cc['ALT'])), pvar_file=f)]
		}
		unlink(tmp)
	}
	if(!length(out)) return(data.table(CHR=integer(), POS=integer(), ID=character(), REF=character(), ALT=character(), pvar_file=character()))
	unique(rbindlist(out, fill=TRUE)[!is.na(CHR) & !is.na(POS)])
}

align_1kg <- function(z, dirmod, lead0, tr){
	z[, rid := .I]
	if(is.null(dirmod)){z[, `:=`(lead_snp0=lead_snp, pvar_id=lead_snp, match_type_1kg='not_checked', pvar_file=NA_character_, pvar_bp=lead_bp)]; return(z)}
	fs <- unique(unlist(lapply(sort(unique(z$lead_chr)), \(chr) pvar_files(dirmod, chr)))); fs <- fs[file.exists(fs)]
	mode <- pvar_mode(fs); cojo_mode <- id_mode(z$lead_snp); p <- rd_pvar(dirmod, z)
	fwrite(data.table(source=c('COJO_lead','1000G_pvar'), id_mode=c(cojo_mode, mode), example=c(paste(head(z$lead_snp, 8), collapse=','), paste(head(p$ID, 8), collapse=','))), paste0(lead0, '/', tr, '.id_sanity.tsv'), sep='\t')
	if(!mode %chin% c('rsid','chrpos')) stop('ERROR: 1000G pvar ID mode is ', mode, ' for ', tr, '. Expected rsID or CHR:POS:REF:ALT. Fix the 1kg VCF/pfile IDs first; this script will not run --set-all-var-ids automatically.')
	z[, `:=`(lead_snp0=lead_snp, lead_bp0=lead_bp, pvar_id=NA_character_, pvar_bp=NA_integer_, pvar_file=NA_character_, match_type_1kg=NA_character_)]

	p[, pvar_rank := fifelse(grepl('male\\.chrX\\.(par|nonPar)\\.pvar$', pvar_file), 1L, 2L)]
	m <- data.table(rid=integer(), pvar_id=character(), pvar_bp=integer(), pvar_file=character(), match_type_1kg=character(), priority=integer())

	idm <- merge(z[usable_id(lead_snp), .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)],
		p[usable_id(ID)], by.x='lead_snp', by.y='ID', allow.cartesian=TRUE)
	idm <- idm[lead_chr == CHR]
	if(nrow(idm)){
		idm[, allele_match := allele_score(EA, NEA, REF, ALT)]
		idm <- idm[allele_match > 0]
		idm[, pos_same := as.integer(lead_bp == POS)]
		if(nrow(idm)){
			setorder(idm, rid, -allele_match, -pos_same, pvar_rank)
			idm <- idm[, .SD[1], by=rid]
			m <- rbind(m, idm[, .(
				rid, pvar_id=lead_snp, pvar_bp=POS, pvar_file,
				match_type_1kg=fifelse(pos_same == 1L, 'ID_allele_POS_same', 'ID_allele_POS_from_pvar'),
				priority=0L
			)], fill=TRUE)
		}
	}

	unmatched <- z[!(rid %in% m$rid)]
	if(nrow(unmatched)){
		pm <- merge(unmatched[, .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)],
			p, by.x=c('lead_chr','lead_bp'), by.y=c('CHR','POS'), allow.cartesian=TRUE)
		pm <- pm[(!usable_id(lead_snp)) | (!usable_id(ID))]
		if(nrow(pm)){
			pm[, allele_match := allele_score(EA, NEA, REF, ALT)]
			pm <- pm[allele_match > 0]
			if(nrow(pm)){
				setorder(pm, rid, -allele_match, pvar_rank)
				pm <- pm[, .SD[1], by=rid]
				pm[, out_id := fifelse(usable_id(ID), ID, lead_snp)]
				m <- rbind(m, pm[, .(
					rid, pvar_id=out_id, pvar_bp=lead_bp, pvar_file,
					match_type_1kg='POS_allele_missing_ID',
					priority=1L
				)], fill=TRUE)
			}
		}
	}

	h <- m
	h <- h[!is.na(pvar_id) & nzchar(pvar_id) & pvar_id!='.']
	if(nrow(h)){setorder(h, rid, priority); h <- h[, .SD[1], by=rid]; z[h, `:=`(pvar_id=i.pvar_id, pvar_bp=i.pvar_bp, pvar_file=i.pvar_file, match_type_1kg=i.match_type_1kg), on=.(rid)]}
	z[, lead_snp := fifelse(!is.na(pvar_id), pvar_id, lead_snp)]
	z[, lead_bp := fifelse(!is.na(pvar_bp), pvar_bp, lead_bp)]
	z[]
}

dirgwas <- arg(c('--dirgwas','-dirgwas')); dirout <- arg(c('--dirout','-dirout')); dirmod <- arg(c('--dirmod','-dirmod')); traits <- arg(c('--traits','-traits'), TRUE)
if(is.null(dirgwas) || is.null(dirout) || is.null(traits)) stop('Usage: Rscript input_prep.R --dirgwas DIR --dirout DIR --dirmod 1KG_DIR --traits bald bald12 ...')
traits <- unlist(strsplit(paste(traits, collapse=' '), '[, ]+')); traits <- traits[nzchar(traits)]
lead0 <- paste0(normalizePath(dirout, winslash='/', mustWork=FALSE), '/lead'); dir.create(lead0, recursive=TRUE, showWarnings=FALSE)
all_match <- list()
prep_version <- 'prep_input_v3_id_first_pos_if_missing_id'

for(tr in traits){
	cj_f <- paste0(dirgwas, '/cojo/', tr, '/', tr, '.jma.cojo'); gw_f <- paste0(dirgwas, '/clean/', tr, '/', tr, '.gz')
	cj <- fread(cj_f); cc <- c(chr=pick1(names(cj), c('Chr','CHR','chr','#CHROM','CHROM')), bp=pick1(names(cj), c('bp','BP','POS','pos')))
	if(anyNA(cc) || !'SNP' %chin% names(cj)) stop('COJO needs Chr/SNP/bp: ', cj_f)
	lead <- unique(cj[, .(lead_chr=chr_int(get(cc['chr'])), lead_snp=as.character(SNP), lead_bp=as.integer(get(cc['bp'])))])[lead_chr %between% c(1L,23L) & !is.na(lead_bp)]
	g <- rd_gwas(gw_f, lead$lead_snp); z <- merge(lead, g, by.x='lead_snp', by.y='SNP', all.x=TRUE, sort=FALSE)
	z[, fail := is.na(EA) | !nzchar(EA) | is.na(NEA) | !nzchar(NEA) | is.na(BETA)]
	fail <- z[fail==TRUE, .(trait=tr, lead_chr, lead_snp, lead_bp, effect_allele=EA, other_allele=NEA, beta=BETA, fail_reason='not_in_gwas_or_missing_allele_beta')]
	z <- z[fail==FALSE]
	if(nrow(z)) z <- align_1kg(z, dirmod, lead0, tr) else z[, `:=`(lead_snp0=lead_snp, lead_bp0=lead_bp, pvar_id=NA_character_, pvar_bp=NA_integer_, pvar_file=NA_character_, match_type_1kg=NA_character_)]
	all_match[[tr]] <- z[, .(trait=tr, lead_chr, lead_snp0, lead_bp0, lead_snp, lead_bp, EA, NEA, BETA, pvar_id, pvar_bp, pvar_file, match_type_1kg)]
	fail <- rbind(fail, z[is.na(pvar_id), .(trait=tr, lead_chr, lead_snp=lead_snp0, lead_bp=lead_bp0, effect_allele=EA, other_allele=NEA, beta=BETA, fail_reason='not_in_1000G_pvar_or_allele_mismatch')], fill=TRUE)
	z <- z[!is.na(pvar_id)][, .(trait=tr, lead_chr, lead_snp, lead_bp, effect_allele=EA, other_allele=NEA, beta=BETA, risk_allele=fifelse(BETA >= 0, EA, NEA), lead_snp0, match_type_1kg)]
	fwrite(z[, .(lead_chr, lead_snp, lead_bp)], paste0(lead0, '/', tr, '.lead.3col'), sep='\t')
	if(nrow(fail)) fwrite(fail, paste0(lead0, '/', tr, '.lead.fail.tsv'), sep='\t') else unlink(paste0(lead0, '/', tr, '.lead.fail.tsv'))
	fwrite(z, paste0(lead0, '/', tr, '.lead.assoc'), sep='\t')
}
fwrite(rbindlist(lapply(traits, \(tr) fread(paste0(lead0, '/', tr, '.lead.assoc'))), fill=TRUE), paste0(lead0, '/lead.assoc'), sep='\t')
fwrite(rbindlist(all_match, fill=TRUE), paste0(lead0, '/lead.match.tsv'), sep='\t')
writeLines(prep_version, paste0(lead0, '/prep.version'))
