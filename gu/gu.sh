#!/usr/bin/bash
dir0=/mnt/d
g1k=$dir0/data/refGen/1kg_phase3
vindija=$dir0/data/refGen/Vindija
altai=$dir0/data/refGen/Altai
chagyr=$dir0/data/refGen/Chagyr
denisova=$dir0/data/refGen/Denisova

dirscript=$dir0/scripts/gu
dirgwas=$dir0/data/gwas/bald
res=$dir0/analysis/bald
traits="bald bald12 bald13 bald14"

mkdir -p "$res" "$res/ld" "$res/pick" "$res/core_vcf" "$res/mat"
idx_vcf(){ [[ -f "$1.tbi" || -f "$1.csi" ]] || bcftools index -t "$1"; }
vcf_chr(){ bcftools index -s "$1" | awk -v c="$2" '$1==c || $1=="chr"c{print $1; exit}'; }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Prepare lead SNP, risk allele files, Define high-LD core regions in 1000G
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/prepare_input.R" --dirgwas "$dirgwas" --dirout "$res" --traits $traits

printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$res/pick/pick.tsv"
for t in $traits; do
	for c in $(seq 1 22); do
		lead=$res/lead/${t}.lead.tsv; ref=$g1k/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; tmp=$res/ld/$t/chr$c; mkdir -p "$tmp"
		[[ -f $lead && -f $ref ]] || { echo "missing $lead or $ref"; exit 1; }; idx_vcf "$ref"
		awk -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $1==c' "$lead" > "$tmp/lead.tsv"; printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"
		awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3}' "$tmp/lead.tsv" | while IFS=$'\t' read -r snp bp0; do
			bp=$(bcftools query -f '%POS\n' -i 'ID=="'"$snp"'"' "$ref" | head -n1 || true); [[ -n ${bp:-} ]] || continue
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$c" "$snp" "$bp" "$c" "$bp" "$snp" >> "$tmp/ld.tsv"
			plink --vcf "$ref" --double-id --allow-extra-chr --r2 --ld-snp "$snp" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0.98 --out "$tmp/$snp" >/dev/null 2>&1 || continue
			[[ -f "$tmp/$snp.ld" ]] && awk -v t="$t" -v c="$c" -v s="$snp" -v b="$bp" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="CHR_B")ic=i; if($i=="BP_B")ip=i; if($i=="SNP_B")is=i; if($i=="R2")ir=i}; next} ic&&ip&&is&&ir{print t,c,s,b,$ic,$ip,$is,$ir}' "$tmp/$snp.ld" >> "$tmp/ld.tsv"
		done
		awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
		awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k])mn[k]=bp; if(!(k in mx)||bp>mx[k])mx[k]=bp; if(pos<mn[k])mn[k]=pos; if(pos>mx[k])mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1)n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
	done
	d=$res/ld/$t; for x in ld block; do first=1; : > "$d/$x.tsv"; for c in $(seq 1 22); do f=$d/chr$c/$x.tsv; [[ -f $f ]] || continue; [[ $first -eq 1 ]] && { cat "$f" > "$d/$x.tsv"; first=0; } || awk 'NR>1' "$f" >> "$d/$x.tsv"; done; awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"; done
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$res/pick/pick.tsv"
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Extract 1000G and archaic VCFs, convert to genotype matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$res/pick/pick.tsv" | while IFS=$'\t' read -r tr chr snp bp st en n size_bp; do
	d=$res/core_vcf/$tr/${chr}.${snp}.${bp}; lo=$st; hi=$en; [[ $bp -lt $lo ]] && lo=$bp; [[ $bp -gt $hi ]] && hi=$bp; mkdir -p "$d"
	kg=$g1k/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; vi=$vindija/chr${chr}_mq25_mapab100.vcf.gz; al=$altai/chr${chr}_mq25_mapab100.vcf.gz; ch=$chagyr/chr${chr}.noRB.bgz.vcf.gz; de=$denisova/chr${chr}_mq25_mapab100.vcf.gz
	for f in "$kg" "$vi" "$al" "$ch" "$de"; do [[ -f $f ]] || { echo "missing $f"; exit 1; }; idx_vcf "$f"; done
	ck=$(vcf_chr "$kg" "$chr"); cv=$(vcf_chr "$vi" "$chr"); ca=$(vcf_chr "$al" "$chr"); cc=$(vcf_chr "$ch" "$chr"); cd=$(vcf_chr "$de" "$chr")
	for z in "$ck" "$cv" "$ca" "$cc" "$cd"; do [[ -n "$z" ]] || { echo "missing contig chr$chr in one VCF"; exit 1; }; done
	bcftools view -r "$ck:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg"; tabix -f -p vcf "$d/kg.vcf.gz"
	bcftools view -r "$cv:$lo-$hi" -Oz -o "$d/vindija.vcf.gz" "$vi"; tabix -f -p vcf "$d/vindija.vcf.gz"
	bcftools view -r "$ca:$lo-$hi" -Oz -o "$d/altai.vcf.gz" "$al"; tabix -f -p vcf "$d/altai.vcf.gz"
	bcftools view -r "$cc:$lo-$hi" -Oz -o "$d/chagyr.vcf.gz" "$ch"; tabix -f -p vcf "$d/chagyr.vcf.gz"
	bcftools view -r "$cd:$lo-$hi" -Oz -o "$d/denisova.vcf.gz" "$de"; tabix -f -p vcf "$d/denisova.vcf.gz"
done
for t in $traits; do for d in "$res/core_vcf/$t"/*.*.*; do [[ -d $d ]] || continue; o=$res/mat/$t/$(basename "$d"); mkdir -p "$o"; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" > "$o/kg.tsv"; for a in vindija altai chagyr denisova; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/$a.vcf.gz" > "$o/$a.tsv"; done; done; done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Define risk haplotypes and write PHYLIP files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/make_hap.R" "$res"
Rscript "$dirscript/make_phy.R" "$res"
find "$res/phy" -name '*.phy' | sort | while read -r f; do phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100; done
Rscript "$dirscript/plot_tree.R" "$res"
