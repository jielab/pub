#!/bin/bash
set -euo pipefail

dir0=/mnt/d
dirmod=$dir0/refGen/1kg_phase3 # 🏮 现代人基因数据
dirarc=$dir0/refGen/gu # 🏮 古人基因数据
sample_label=EUR_EAS # 🏮 选择分析样本
sample_keep=$dirmod/$sample_label.1id
sample_tag=$sample_label
dirout=$dir0/analysis/gu/ibdmix/$sample_label

chrs="$(seq 1 22 | paste -sd' ' -) X"
njob=8; lod_cut=5; len_cut=1000

mkdir -p $dirout/{vcf,gt,raw,summary,report,log}
log(){ echo "[$(date '+%F %T')] $*"; }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Prepare modern VCFs for sample_keep
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set -- $chrs; check_chr=$1; check_vcf=$dirmod/chr$check_chr.vcf.gz
miss=$(awk 'NR==FNR{a[$1]=1; next} NF && !($1 in a){print $1}' <(bcftools query -l $check_vcf) $sample_keep | head -5 | paste -sd, -)

sample_md5=$(md5sum $sample_keep | awk '{print $1}')
old_file=$(find $dirout -type f \( -name '*.vcf' -o -name '*.gt.txt' -o -name '*.ibdmix.noSNPLOD.txt' -o -name '*.segments.tsv.gz' \) -print -quit 2>/dev/null || true)
echo $sample_md5 > $dirout/sample_keep.md5
log "sample_keep: $sample_keep ($(awk 'NF{n++} END{print n+0}' $sample_keep) samples)"

prep_modern(){
	chr=$1; g1k_vcf=$dirmod/chr$chr.vcf.gz; outvcf=$dirout/vcf/modern.chr$chr.vcf
	[[ -s $outvcf ]] && return
	log "Prepare modern chr$chr"
	bcftools view -S $sample_keep -m2 -M2 -v snps $g1k_vcf |
	awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > $outvcf
}

for chr in $chrs; do
	prep_modern $chr &
	while [[ $(jobs -rp | wc -l) -ge $njob ]]; do sleep 2; done
done
wait


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Run IBDmix and save small segment-level summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_arch(){
	chr=$1; anc=$2; anc_vcf=$3; outvcf=$dirout/vcf/$anc.chr$chr.vcf
	[[ -s $outvcf ]] && return
	log "Prepare archaic $anc chr$chr"
	bcftools view -m2 -M2 -v snps $anc_vcf |
	awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > $outvcf
}

run_one(){
	chr=$1; anc=$2; anc_vcf=$3
	gt=$dirout/gt/$anc.chr$chr.gt.txt
	raw=$dirout/raw/$anc.chr$chr.ibdmix.noSNPLOD.txt
	sum=$dirout/summary/$anc.chr$chr.len${len_cut}.lod${lod_cut}.segments.tsv.gz
	log "IBDmix chr$chr $anc"

	prep_arch $chr $anc $anc_vcf
	[[ -s $gt ]] || generate_gt -a $dirout/vcf/$anc.chr$chr.vcf -m $dirout/vcf/modern.chr$chr.vcf -o $gt > $dirout/log/$anc.chr$chr.generate_gt.log
	[[ -s $raw ]] || ibdmix -g $gt -o $raw -s $sample_keep -d 3 -t > $dirout/log/$anc.chr$chr.ibdmix.log

	if [[ ! -s $sum ]]; then
		bash $summary_sh $len_cut $lod_cut $sample_tag $raw - |
		awk -v anc=$anc -v ss0=$sample_tag 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ID","chrom","start","end","length","slod","sites","positive_lods","negative_lods","sample_set","anc"; next} $1!="ID"{len=(h["length"] ? $(h["length"]) : $(h["end"])-$(h["start"])); ss=(h["pop"] ? $(h["pop"]) : ss0); pos=(h["positive_lods"] ? $(h["positive_lods"]) : ""); neg=(h["negative_lods"] ? $(h["negative_lods"]) : ""); print $(h["ID"]),$(h["chrom"]),$(h["start"]),$(h["end"]),len,$(h["slod"]),$(h["sites"]),pos,neg,ss,anc}' |
		gzip -f > $sum
	fi
}

for chr in $chrs; do
	for ad in $dirarc/*; do
		[[ -d $ad ]] || continue; anc=$(basename $ad)
		for anc_vcf in $ad/*chr${chr}[_.]*.vcf.gz; do
			[[ -s $anc_vcf ]] || continue
			run_one $chr $anc $anc_vcf &
			while [[ $(jobs -rp | wc -l) -ge $njob ]]; do sleep 2; done
			break
		done
	done
done
wait


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: summarize & report
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out=$dirout/report/len${len_cut}.lod${lod_cut}.segments.tsv.gz
summary_files=($dirout/summary/*.chr*.len${len_cut}.lod${lod_cut}.segments.tsv.gz)
zcat ${summary_files[@]} | awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} $1!="ID"' | gzip -f > $out
log "Done: $out"

Rscript $dir0/scripts/gu/ibdmix_report.R --in $out --outdir $dirout/report --sample_info $dir0/files/1kg.v3.sample.txt --summary_by super_pop,gender
