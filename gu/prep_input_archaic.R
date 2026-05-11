pacman::p_load(data.table)
args <- commandArgs(TRUE)
if(length(args) < 4) {
	stop("Usage: Rscript prep_input_archaic.R kg.vcf.gz archaic.raw.vcf.gz out.vcf.gz sample_name")
}

kg_vcf  <- args[1]
arc_vcf <- args[2]
out_vcf <- args[3]
sample  <- args[4]

if(!file.exists(kg_vcf))  stop("Missing kg_vcf: ", kg_vcf)
if(!file.exists(arc_vcf)) stop("Missing arc_vcf: ", arc_vcf)

bcftools <- Sys.which("bcftools")
bgzip <- Sys.which("bgzip")
tabix <- Sys.which("tabix")
if(bcftools == "") stop("bcftools not found in PATH")
if(bgzip == "") stop("bgzip not found in PATH")
if(tabix == "") stop("tabix not found in PATH")

message("kg_vcf  = ", kg_vcf)
message("arc_vcf = ", arc_vcf)
message("out_vcf = ", out_vcf)
message("sample  = ", sample)

read_query <- function(vcf, fmt, cols, label){
	tmp <- tempfile()
	rc <- system2(bcftools, c("query", "-f", shQuote(fmt), shQuote(vcf)), stdout = tmp)
	if(rc != 0) stop("bcftools query failed for ", label, ": ", vcf)
	if(!file.exists(tmp) || file.info(tmp)$size == 0){
		unlink(tmp)
		return(as.data.table(setNames(rep(list(character()), length(cols)), cols)))
	}
	x <- fread(tmp, header = FALSE, col.names = cols, showProgress = FALSE)
	unlink(tmp)
	x
}

kg <- read_query(
	kg_vcf,
	"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n",
	c("CHR", "POS", "ID", "REF_kg", "ALT_kg"),
	"1000G template"
)

arc <- read_query(
	arc_vcf,
	"%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n",
	c("CHR", "POS", "REF_arc", "ALT_arc", "GT_arc"),
	"archaic raw"
)

if(nrow(kg) == 0) stop("kg_vcf has 0 variants: ", kg_vcf)
if(nrow(arc) == 0) warning("arc_vcf has 0 records: ", arc_vcf)

kg[, `:=`(
	idx = .I,
	CHR = as.character(CHR),
	POS = as.integer(POS),
	ID = fifelse(is.na(ID) | ID == "", ".", as.character(ID)),
	REF_kg = toupper(REF_kg),
	ALT_kg = toupper(ALT_kg)
)]

arc[, `:=`(
	CHR = as.character(CHR),
	POS = as.integer(POS),
	REF_arc = toupper(REF_arc),
	ALT_arc = toupper(ALT_arc),
	GT_arc = gsub("\\|", "/", as.character(GT_arc))
)]

x <- merge(kg, arc, by = c("CHR", "POS"), all.x = TRUE, sort = FALSE)
setorder(x, idx)

x[, GT := "./."]
x[, project_match := fifelse(is.na(REF_arc), "missing_archaic_record", "unmatched_allele")]

x[!is.na(REF_arc) & REF_arc == REF_kg & ALT_arc == ".", `:=`(
	GT = "0/0",
	project_match = "REF_arc_matches_REF_kg"
)]

x[!is.na(REF_arc) & REF_arc == ALT_kg & ALT_arc == ".", `:=`(
	GT = "1/1",
	project_match = "REF_arc_matches_ALT_kg"
)]

flip_gt <- function(gt){
	gt <- gsub("\\|", "/", as.character(gt))
	out <- rep("./.", length(gt))
	out[gt == "0/0"] <- "1/1"
	out[gt == "1/1"] <- "0/0"
	out[gt %chin% c("0/1", "1/0")] <- "0/1"
	out[gt %chin% c(".", "./.") | is.na(gt)] <- "./."
	out
}

x[!is.na(REF_arc) & REF_arc == REF_kg & ALT_arc == ALT_kg, `:=`(
	GT = GT_arc,
	project_match = "REF_ALT_exact"
)]

x[!is.na(REF_arc) & REF_arc == ALT_kg & ALT_arc == REF_kg, `:=`(
	GT = flip_gt(GT_arc),
	project_match = "REF_ALT_swapped"
)]

x[, project_priority := fcase(
	project_match == "REF_ALT_exact", 1L,
	project_match == "REF_ALT_swapped", 2L,
	project_match == "REF_arc_matches_REF_kg", 3L,
	project_match == "REF_arc_matches_ALT_kg", 4L,
	default = 9L
)]
setorder(x, idx, project_priority)
x <- x[, .SD[1], by = idx]
setorder(x, idx)

x[GT == "0", GT := "0/0"]
x[GT == "1", GT := "1/1"]
x[GT == ".", GT := "./."]
x[is.na(GT) | !(GT %chin% c("0/0", "0/1", "1/0", "1/1", "./.")), GT := "./."]

x <- x[REF_kg %chin% c("A", "C", "G", "T") & ALT_kg %chin% c("A", "C", "G", "T")]

vcf <- x[, .(
	CHR,
	POS,
	ID,
	REF = REF_kg,
	ALT = ALT_kg,
	QUAL = ".",
	FILTER = ".",
	INFO = ".",
	FORMAT = "GT",
	GT
)]

dir.create(dirname(out_vcf), recursive = TRUE, showWarnings = FALSE)

tmp_vcf <- tempfile(fileext = ".vcf")
cat("##fileformat=VCFv4.2\n", file = tmp_vcf)
cat("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", file = tmp_vcf, append = TRUE)
cat(paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample, sep = "\t"), "\n", file = tmp_vcf, append = TRUE)
fwrite(vcf, tmp_vcf, sep = "\t", append = TRUE, col.names = FALSE)

if(file.exists(out_vcf)) file.remove(out_vcf)
if(file.exists(paste0(out_vcf, ".tbi"))) file.remove(paste0(out_vcf, ".tbi"))

tmp_gz <- tempfile(fileext = ".vcf.gz")
rc1 <- system2(bgzip, c("-c", shQuote(tmp_vcf)), stdout = tmp_gz)
if(rc1 != 0) stop("bgzip failed for: ", out_vcf)
if(!file.copy(tmp_gz, out_vcf, overwrite = TRUE)) stop("failed to copy bgzip output to: ", out_vcf)
unlink(tmp_gz)
if(!file.exists(out_vcf) || file.info(out_vcf)$size == 0) stop("bgzip did not create output: ", out_vcf)

rc2 <- system2(tabix, c("-f", "-p", "vcf", shQuote(out_vcf)))
if(rc2 != 0) stop("tabix failed for: ", out_vcf)
if(!file.exists(paste0(out_vcf, ".tbi"))) stop("tabix did not create index: ", paste0(out_vcf, ".tbi"))

message("Wrote: ", out_vcf)
message("n_template=", nrow(kg))
message("n_archaic_raw=", nrow(arc))
message("n_output=", nrow(vcf))
message("GT_0_0=", sum(vcf$GT == "0/0"))
message("GT_0_1=", sum(vcf$GT %chin% c("0/1", "1/0")))
message("GT_1_1=", sum(vcf$GT == "1/1"))
message("GT_missing=", sum(vcf$GT == "./."))
message("projection_counts=", paste(names(table(x$project_match)), as.integer(table(x$project_match)), sep = ":", collapse = ","))
