#!/bin/bash

set -o pipefail

dir0=/mnt/d
export PATH="$dir0/software/bin:$PATH" 
dirmod=$dir0/refGen/1kg_phase3
dirarc=$dir0/refGen/gu
dirsoft=$dir0/software/IBDmix
dirscript=$dir0/scripts/gu
sample_info=$dir0/files/1kg.v3.sample.txt
sample_label=EUR_EAS

sample_keep=$dirmod/$sample_label.1id
dirout=$dir0/analysis/gu/ibdmix/$sample_label

chrs=({1..22} X) 
job_of_chr=2
job_in_chr=5
lod_cut=5
len_cut=1000
start_step=${START_STEP:-s1}

mkdir -p "$dirout"/{summary,report,log,tmp}
exec > >(tee "$dirout/log/ibdmix.log") 2>&1 
log(){ echo "[$(date '+%F %T')] $*"; }

for x in awk bcftools md5sum find paste gzip zcat bash Rscript; do
	command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }
done
for x in "$dirsoft/build/src/generate_gt" "$dirsoft/build/src/ibdmix" "$dirsoft/src/summary.sh" "$dirscript/ibdmix_report.R" "$sample_keep"; do
	[[ -s "$x" ]] || { log "ERROR missing required file: $x"; exit 1; }
done

check_chr=${chrs[0]} 
check_vcf=$dirmod/chr$check_chr.vcf.gz
[[ -s "$check_vcf" ]] || { log "ERROR missing check VCF: $check_vcf"; exit 1; }
check_samples=$dirout/tmp/chr$check_chr.samples
bcftools query -l "$check_vcf" > "$check_samples" || { log "ERROR: failed to query samples from $check_vcf"; exit 1; }
miss=$(awk 'NR==FNR{a[$1]=1; next} NF && !($1 in a){print $1}' "$check_samples" "$sample_keep" | head -5 | paste -sd, -)
if [[ -n "$miss" ]]; then log "ERROR: sample_keep contains IDs absent from $check_vcf, e.g. $miss"; exit 1; fi
sample_md5=$(md5sum "$sample_keep" | awk '{print $1}')
stamp=$dirout/sample_keep.md5
if [[ -s "$stamp" ]] && [[ "$(cat "$stamp")" != "$sample_md5" ]]; then
	log "ERROR: sample_keep changed for existing dirout: $dirout"
	log "Use a new sample_label or remove the old output directory."
	exit 1
fi
if [[ ! -s "$stamp" ]]; then
	old_file=$(find "$dirout" -type f \( -name '*.segments.tsv.gz' -o -name '*.burden.tsv.gz' -o -name '*.bedGraph' \) -print -quit 2>/dev/null || true)
	if [[ -n "$old_file" ]]; then log "WARNING: existing outputs have no sample_keep.md5 stamp: $old_file"; fi
fi
echo "$sample_md5" > "$stamp"
log "sample_keep: $sample_keep ($(awk 'NF{n++} END{print n+0}' "$sample_keep") samples)"
if (( job_of_chr < 1 || job_in_chr < 1 )); then
	log "ERROR: job_of_chr and job_in_chr must be positive integers"
	exit 1
fi
log "parallel jobs: job_of_chr=$job_of_chr job_in_chr=$job_in_chr max=$((job_of_chr * job_in_chr))"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: run 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_modern(){
	chr=$1; tmp_chr=$2
	g1k_vcf=$dirmod/chr$chr.vcf.gz
	outvcf=$tmp_chr/modern.chr$chr.vcf
	mlog=$dirout/log/modern.chr$chr.prep.log
	log "Prepare modern chr$chr"
	echo "[$(date '+%F %T')] bcftools view modern chr$chr" > "$mlog"

	if [[ "$chr" == "X" ]]; then
		if ! bcftools view -S "$sample_keep" -m2 -M2 -v snps "$g1k_vcf" 2>> "$mlog" |
		awk '
		BEGIN{OFS="\t"}
		/^#/{print; next}
		{
			sub(/^chr/,"",$1)
			for(i=10;i<=NF;i++){
				split($i,a,":")
				if(a[1]=="0") a[1]="0/0"
				else if(a[1]=="1") a[1]="1/1"
				else if(a[1]==".") a[1]="./."
				$i=a[1]
				for(j=2;j in a;j++) $i=$i ":" a[j]
			}
			print
		}' > "$outvcf"; then
			log "ERROR: bcftools/awk failed while preparing modern chr$chr"; return 1
		fi
	else
		if ! bcftools view -S "$sample_keep" -m2 -M2 -v snps "$g1k_vcf" 2>> "$mlog" |
		awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > "$outvcf"; then
			log "ERROR: bcftools/awk failed while preparing modern chr$chr"; return 1
		fi
	fi

	if [[ ! -s "$outvcf" ]]; then log "ERROR: failed to prepare modern chr$chr: $outvcf"; return 1; fi
}

prep_arch(){
	chr=$1; anc=$2; anc_vcf=$3; tmp_chr=$4
	outvcf=$tmp_chr/$anc.chr$chr.vcf
	alog=$dirout/log/$anc.chr$chr.prep_arch.log
	log "Prepare archaic $anc chr$chr"
	echo "[$(date '+%F %T')] bcftools view archaic $anc chr$chr" > "$alog"
	if ! bcftools view -m2 -M2 -v snps "$anc_vcf" 2>> "$alog" |
	awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > "$outvcf"; then
		log "ERROR: bcftools/awk failed while preparing archaic $anc chr$chr"; return 1
	fi
	if [[ ! -s "$outvcf" ]]; then log "ERROR: failed to prepare archaic $anc chr$chr: $outvcf"; return 1; fi
}

check_gt(){
	gt=$1
	awk 'NR==1{nf=NF; next} NF!=nf{print "bad_gt_line",NR,"NF",NF,"expected",nf,"chr",$1,"pos",$2; exit 1}' "$gt"
}

run_one(){
	chr=$1; anc=$2; anc_vcf=$3; tmp_chr=$4
	avcf=$tmp_chr/$anc.chr$chr.vcf
	mvcf=$tmp_chr/modern.chr$chr.vcf
	gt=$tmp_chr/$anc.chr$chr.gt.txt
	raw=$tmp_chr/$anc.chr$chr.raw.txt
	sum=$dirout/summary/$anc.chr$chr.segments.tsv.gz
	glog=$dirout/log/$anc.chr$chr.generate_gt.log
	ilog=$dirout/log/$anc.chr$chr.ibdmix.log
	slog=$dirout/log/$anc.chr$chr.summary.log

	log "IBDmix chr$chr $anc"
	prep_arch "$chr" "$anc" "$anc_vcf" "$tmp_chr" || return 1

	echo "[$(date '+%F %T')] generate_gt $anc chr$chr" > "$glog"
	$dirsoft/build/src/generate_gt -a "$avcf" -m "$mvcf" -o "$gt" >> "$glog" 2>&1
	if [[ ! -s "$gt" ]]; then log "ERROR: generate_gt failed: $gt"; return 1; fi
	check_gt "$gt" >> "$glog" 2>&1 || { log "ERROR: malformed gt table: $gt"; return 1; }
	rm -f "$avcf"

	echo "[$(date '+%F %T')] ibdmix $anc chr$chr" > "$ilog"
	$dirsoft/build/src/ibdmix -g "$gt" -o "$raw" -s "$sample_keep" -d 3 -t >> "$ilog" 2>&1
	if [[ ! -s "$raw" ]]; then log "ERROR: ibdmix failed: $raw"; return 1; fi

	echo "[$(date '+%F %T')] summary.sh len=$len_cut lod=$lod_cut $anc chr$chr" > "$slog"
	if ! bash "$dirsoft/src/summary.sh" "$len_cut" "$lod_cut" "$sample_label" "$raw" - 2>> "$slog" |
	awk -v anc="$anc" -v ss0="$sample_label" 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ID","chrom","start","end","length","slod","sites","positive_lods","negative_lods","sample_set","anc"; next} $1!="ID"{len=(h["length"] ? $(h["length"]) : $(h["end"])-$(h["start"])); ss=(h["pop"] ? $(h["pop"]) : ss0); pos=(h["positive_lods"] ? $(h["positive_lods"]) : ""); neg=(h["negative_lods"] ? $(h["negative_lods"]) : ""); print $(h["ID"]),$(h["chrom"]),$(h["start"]),$(h["end"]),len,$(h["slod"]),$(h["sites"]),pos,neg,ss,anc}' 2>> "$slog" |
	gzip -f > "$sum"; then
		log "ERROR: summary pipeline failed: $anc chr$chr"; return 1
	fi
	if [[ ! -s "$sum" ]]; then log "ERROR: summary failed: $sum"; return 1; fi

	rm -f "$gt" "$raw"
}

run_chr(){
	chr=$1
	log "Start chr$chr"
	tmp_chr=$dirout/tmp/chr$chr
	rm -rf "$tmp_chr"
	mkdir -p "$tmp_chr"

	prep_modern "$chr" "$tmp_chr" || { log "ERROR: failed to prepare modern chr$chr"; rm -rf "$tmp_chr"; return 1; }

	running=0
	chr_status=0
	chr_jobs=0
	for ad in "$dirarc"/*; do
		[[ -d "$ad" ]] || continue
		anc=$(basename "$ad")
		for anc_vcf in "$ad"/*chr"$chr"[_.]*.vcf.gz; do
			[[ -s "$anc_vcf" ]] || continue
			run_one "$chr" "$anc" "$anc_vcf" "$tmp_chr" &
			running=$((running+1))
			chr_jobs=$((chr_jobs+1))
			if (( running >= job_in_chr )); then
				wait -n || chr_status=1
				running=$((running-1))
			fi
			break
		done
	done
	if (( chr_jobs == 0 )); then log "ERROR: no archaic VCF found for chr$chr under $dirarc"; rm -rf "$tmp_chr"; return 1; fi
	while (( running > 0 )); do
		wait -n || chr_status=1
		running=$((running-1))
	done

	rm -f "$tmp_chr/modern.chr$chr.vcf"
	rm -rf "$tmp_chr"
	if (( chr_status != 0 )); then log "ERROR: one or more jobs failed for chr$chr"; return 1; fi
	log "Done chr$chr; temporary VCF/gt/raw removed"
}

if [[ "$start_step" == "s1" ]]; then
	running_chr=0
	s1_status=0
	for chr in "${chrs[@]}"; do
		run_chr "$chr" &
		running_chr=$((running_chr+1))
		if (( running_chr >= job_of_chr )); then
			wait -n || s1_status=1
			running_chr=$((running_chr-1))
		fi
	done
	while (( running_chr > 0 )); do
		wait -n || s1_status=1
		running_chr=$((running_chr-1))
	done
	rmdir "$dirout/tmp" 2>/dev/null || true
	if (( s1_status != 0 )); then log "ERROR: one or more chromosome jobs failed"; exit 1; fi
elif [[ "$start_step" == "s2" ]]; then
	log "Skip s1; reuse existing chromosome summary files in $dirout/summary"
else
	log "ERROR unknown START_STEP=$start_step; use s1 or s2"
	exit 1
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Summary and report
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log "Start s2: merge summary files and run report"
out=$dirout/summary/segments.tsv.gz
mapfile -t summary_files < <(find "$dirout/summary" -type f -name '*.chr*.segments.tsv.gz' ! -name '*.len*.lod*.segments.tsv.gz' | sort)
if (( ${#summary_files[@]} == 0 )); then log "ERROR: no summary files found: $dirout/summary"; exit 1; fi

if ! zcat "${summary_files[@]}" | awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} $1!="ID"' | gzip -f > "$out"; then
	log "ERROR: failed to merge summary pipeline: $out"; exit 1
fi
if [[ ! -s "$out" ]]; then log "ERROR: failed to merge summaries: $out"; exit 1; fi
log "Done: $out"

report_args=(--in "$out" --outdir "$dirout/report" --sample_keep "$sample_keep" --stats_for_locus "3:45859651-45909024;9:136130562-136150630")
if [[ -s "$sample_info" ]]; then
	if [[ "$sample_info" == *1kg* ]]; then
		by_group="super_pop"
	elif [[ "$sample_info" == *ukb* || "$sample_info" == *UKB* ]]; then
		by_group="race"
	else
		log "WARNING: unknown sample_info type, skip by-group outputs: $sample_info"
		by_group=""
	fi
	report_args+=(--sample_info "$sample_info")
	if [[ -n "$by_group" ]]; then
		report_args+=(--stats_by_group "$by_group")
	fi
else
	log "WARNING: sample_info not found, skip by-group outputs: $sample_info"
fi
if ! Rscript "$dirscript/ibdmix_report.R" "${report_args[@]}"; then
	log "ERROR: ibdmix_report.R failed"
	exit 1
fi
log "All done. Report: $dirout/report"
