# 🏮 本流程默认使用 hg19 坐标
# 1KG: https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# Vindija: https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz
# Altai: https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz
# Chagyr: https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz	

#!/usr/bin/bash
dir0=/mnt/d
g1k=$dir0/data/refGen/1kg_phase3
vindija=$dir0/data/refGen/Vindija
altai=$dir0/data/refGen/Altai
chagyr=$dir0/data/refGen/Chagyr
denisova=$dir0/data/refGen/Denisova

dirscript=$dir0/scripts/gu
dirgwas=$dir0/data/gwas/bald
dirout=$dir0/analysis/bald
mkdir -p $dirout/ld $dirout/pick

traits="bald bald12 bald13 bald14"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩Prepare lead SNP and risk allele files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/0_prepare_input.R" --dirgwas "$dirgwas" --dirout "$dirout" --traits $traits

printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/pick/pick.tsv"

for t in $traits; do
	lead="$dirout/lead/${t}.lead.tsv"; [[ -f "$lead" ]] || { echo "missing $lead"; exit 1; }
	for c in $(seq 1 22); do
		ref="$g1k/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"; tmp="$dirout/ld/$t/chr$c"; mkdir -p "$tmp"
		[[ -f "$ref" ]] || { echo "missing $ref"; exit 1; }; [[ -f "$ref.tbi" || -f "$ref.csi" ]] || bcftools index -t "$ref"
		awk -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $1==c' "$lead" > "$tmp/lead.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"

		if [[ $(awk 'END{print NR-1}' "$tmp/lead.tsv") -gt 0 ]]; then
			awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3}' "$tmp/lead.tsv" | while IFS=$'\t' read -r snp bp0; do
				bp=$(bcftools query -f '%POS\n' -i 'ID=="'"$snp"'"' "$ref" | head -n1 || true); [[ -n ${bp:-} ]] || continue
				printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$c" "$snp" "$bp" "$c" "$bp" "$snp" >> "$tmp/ld.tsv"
				"$plink" --vcf "$ref" --double-id --allow-extra-chr --r2 --ld-snp "$snp" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0.98 --out "$tmp/$snp" >/dev/null 2>&1 || continue
				[[ -f "$tmp/$snp.ld" ]] && awk -v t="$t" -v c="$c" -v s="$snp" -v b="$bp" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="CHR_B")ic=i; else if($i=="BP_B")ip=i; else if($i=="SNP_B")is=i; else if($i=="R2")ir=i}; next} ic&&ip&&is&&ir{print t,c,s,b,$ic,$ip,$is,$ir}' "$tmp/$snp.ld" >> "$tmp/ld.tsv"
			done
			awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
		fi

		awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k])mn[k]=bp; if(!(k in mx)||bp>mx[k])mx[k]=bp; if(pos<mn[k])mn[k]=pos; if(pos>mx[k])mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1)n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
	done

	d="$dirout/ld/$t"; for x in ld block; do first=1; : > "$d/$x.tsv"; for c in $(seq 1 22); do f="$d/chr$c/$x.tsv"; [[ -f "$f" ]] || continue; [[ $first -eq 1 ]] && { cat "$f" > "$d/$x.tsv"; first=0; } || awk 'NR>1' "$f" >> "$d/$x.tsv"; done; awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"; done
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$dirout/pick/pick.tsv"
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Extract 1000G and archaic VCFs, convert to genotype matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p "$res/core_vcf"
awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$res/pick/pick.tsv" | while IFS=$'\t' read -r tr chr snp bp st en n size_bp; do
	d=$res/core_vcf/$tr/${chr}.${snp}.${bp}; lo=$st; hi=$en; [[ $bp -lt $lo ]] && lo=$bp; [[ $bp -gt $hi ]] && hi=$bp; r=${chr}:${lo}-${hi}; mkdir -p "$d"
	kg=$g1k/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; vi=$vindija/chr${chr}_mq25_mapab100.vcf.gz; al=$altai/chr${chr}_mq25_mapab100.vcf.gz; ch=$chagyr/chr${chr}.noRB.bgz.vcf.gz; de=$denisova/chr${chr}_mq25_mapab100.vcf.gz
	for f in "$kg" "$vi" "$al" "$ch" "$de"; do [[ -f $f ]] || { echo "missing $f"; exit 1; }; idx_vcf "$f"; done
	bcftools view -r "$r" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg"; tabix -f -p vcf "$d/kg.vcf.gz"
	bcftools view -r "$r" -Oz -o "$d/vindija.vcf.gz" "$vi"; tabix -f -p vcf "$d/vindija.vcf.gz"
	bcftools view -r "$r" -Oz -o "$d/altai.vcf.gz" "$al"; tabix -f -p vcf "$d/altai.vcf.gz"
	bcftools view -r "$r" -Oz -o "$d/chagyr.vcf.gz" "$ch"; tabix -f -p vcf "$d/chagyr.vcf.gz"
	bcftools view -r "$r" -Oz -o "$d/denisova.vcf.gz" "$de"; tabix -f -p vcf "$d/denisova.vcf.gz"
done

for t in $traits; do
	for d in "$res/core_vcf/$t"/*.*.*; do
		[[ -d $d ]] || continue; o=$res/mat/$t/$(basename "$d"); mkdir -p "$o"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" > "$o/kg.tsv"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/vindija.vcf.gz" > "$o/vindija.tsv"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/altai.vcf.gz" > "$o/altai.tsv"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/chagyr.vcf.gz" > "$o/chagyr.tsv"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/denisova.vcf.gz" > "$o/denisova.tsv"
	done
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Define risk haplotypes and archaic-matched loci, generate PHYLIP input files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/make_hap.R" "$res"
Rscript "$dirscript/make_phy.R" "$res"

find "$res/phy" -name '*.phy' | sort | while read -r f; do phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100; done

Rscript "$dirscript/plot_tree.R" "$res"
