#!/usr/bin/bash
set -euo pipefail

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 IBDmix genome-wide scan for 1000G archaic tracts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0=/mnt/d
g1k=$dir0/data/refGen/1kg_phase3
arch0=$dir0/data/refGen/gu
ibd=$dir0/software/IBDmix
out=$dir0/analysis/gu/ibdmix
pop=EUR
sample=$g1k/$pop.sample.ids
chrs="$(seq 1 22)"
njob=6
lod_cut=5
len_cut=1000

generate_gt=$ibd/build/src/generate_gt
ibdmix=$ibd/build/src/ibdmix
summary_sh=$ibd/scripts/summary.sh
mkdir -p "$out"/{vcf,gt,raw,summary,log}

log(){ echo "[$(date '+%F %T')] $*"; }
fix_chr(){ awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/^chr/,"",$1); print}'; }

prep_modern(){
	local c=$1 kg opt=""
	kg=$g1k/ALL.chr$c.vcf.gz
	[[ -f $kg ]] || kg=$g1k/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	[[ -s $sample ]] && opt="-S $sample"
	[[ -s $out/vcf/modern.chr$c.vcf ]] || bcftools view $opt -r "$c" -m2 -M2 -v snps "$kg" | fix_chr > "$out/vcf/modern.chr$c.vcf"
}

prep_arch(){
	local c=$1 a=$2 f=$3
	[[ -s $out/vcf/$a.chr$c.vcf ]] || bcftools view -r "$c" -m2 -M2 -v snps "$f" | fix_chr > "$out/vcf/$a.chr$c.vcf"
}

run_one(){
	local c=$1 a=$2 f=$3 tag raw sum gt
	tag=$pop
	[[ -s $sample ]] || tag=ALL
	gt=$out/gt/$a.chr$c.gt.txt
	raw=$out/raw/$a.$tag.chr$c.ibdmix.txt
	sum=$out/summary/$a.$tag.chr$c.len${len_cut}.lod${lod_cut}.txt.gz

	prep_arch "$c" "$a" "$f"
	[[ -s $gt ]] || "$generate_gt" -a "$out/vcf/$a.chr$c.vcf" -m "$out/vcf/modern.chr$c.vcf" -o "$gt" > "$out/log/$a.chr$c.generate_gt.log" 2>&1

	if [[ ! -s $raw ]]; then
		if [[ -s $sample ]]; then
			"$ibdmix" -g "$gt" -o "$raw" -s "$sample" -d 3 -t --write-snps --write-lods > "$out/log/$a.chr$c.ibdmix.log" 2>&1
		else
			"$ibdmix" -g "$gt" -o "$raw" -d 3 -t --write-snps --write-lods > "$out/log/$a.chr$c.ibdmix.log" 2>&1
		fi
	fi

	[[ -s $sum ]] || "$summary_sh" "$len_cut" "$lod_cut" "$tag" "$raw" - | gzip -f > "$sum"
}

tag=$pop
[[ -s $sample ]] || tag=ALL

for c in $chrs; do
	log "Prepare 1000G chr$c"
	prep_modern "$c"
	for ad in "$arch0"/*; do
		[[ -d $ad ]] || continue
		a=$(basename "$ad")
		f=$(ls "$ad"/*chr${c}_*.vcf.gz "$ad"/*chr${c}.*.vcf.gz 2>/dev/null | head -n1)
		[[ -n ${f:-} ]] || continue
		log "IBDmix chr$c $a"
		run_one "$c" "$a" "$f" &
		while [[ $(jobs -rp | wc -l) -ge $njob ]]; do sleep 2; done
	done
done
wait

zcat "$out"/summary/*.$tag.chr*.len${len_cut}.lod${lod_cut}.txt.gz > "$out/archaic_tracts.$tag.len${len_cut}.lod${lod_cut}.tsv"
log "Done: $out/archaic_tracts.$tag.len${len_cut}.lod${lod_cut}.tsv"
