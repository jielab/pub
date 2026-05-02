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
chrs="$(seq 1 22) X"
ld_r2=0.98
njob=6

mkdir -p "$res" "$res/ld" "$res/pick" "$res/coreVcf" "$res/mat"
idx_vcf(){ [[ -f "$1.tbi" || -f "$1.csi" ]] || bcftools index -t "$1"; }
vcf_chr(){ bcftools index -s "$1" | awk -v c="$2" '$1==c || $1=="chr"c{print $1; exit}'; }
log(){ echo "[$(date '+%F %T')] $*"; }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s0: Prepare lead SNP and risk allele files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/prepare_input.R" --dirgwas "$dirgwas" --dirout "$res" --traits $traits

# 加上 Positive control 🏮
if ! awk 'BEGIN{FS=OFS="\t"} NR>1 && $2=="rs35044562"{f=1} END{exit !f}' "$res/lead/bald.lead.tsv"; then
	sed -i '1a 3\trs35044562\t45859651' "$res/lead/bald.lead.tsv"
fi
if ! awk 'BEGIN{FS=OFS="\t"} NR>1 && $3=="rs35044562"{f=1} END{exit !f}' "$res/risk/bald.risk.tsv"; then
	sed -i '1a bald\t3\trs35044562\t45859651\tG\tA\t1\tG' "$res/risk/bald.risk.tsv"
fi
first=${traits%% *}
{ head -n1 "$res/risk/$first.risk.tsv"; for t in $traits; do awk 'NR>1' "$res/risk/$t.risk.tsv"; done; } > "$res/risk/risk.tsv"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Define high-LD core regions，因为并行处理，所以代码较长
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$res/pick/pick.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\treason\n" > "$res/pick/lead_1000G.fail.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tPAR\tLD_pfile\n" > "$res/pick/chrX_PAR.tsv"

run_chr(){
	local t=$1 c=$2 cn=$2 lead tmp vcf ctg nlead snp bp0 bp lo hi par pbase outvcor
	[[ $c == X ]] && cn=23
	lead=$res/lead/${t}.lead.tsv; tmp=$res/ld/$t/chr$c; mkdir -p "$tmp"
	vcf=$g1k/ALL.chr${c}.vcf.gz; [[ -f $vcf ]] || vcf=$g1k/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	idx_vcf "$vcf"; ctg=$(vcf_chr "$vcf" "$c")
	[[ -n $ctg ]] || { echo -e "$t\t$cn\tNA\tNA\tmissing_contig_chr$c" >> "$res/pick/lead_1000G.fail.tsv"; return; }

	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $1==cn || toupper($1)==c' "$lead" > "$tmp/lead.tsv"
	nlead=$(awk 'NR>1' "$tmp/lead.tsv" | wc -l); log "DO $t chr$c: $nlead lead SNPs"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"

	awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3}' "$tmp/lead.tsv" | while IFS=$'\t' read -r snp bp0; do
		lo=$((bp0-5)); [[ $lo -lt 1 ]] && lo=1; hi=$((bp0+5))
		bp=$(bcftools query -r "$ctg:$lo-$hi" -f '%POS\n' -i 'ID=="'"$snp"'"' "$vcf" | head -n1 || true)
		[[ -n ${bp:-} ]] || bp=$(bcftools query -f '%POS\n' -i 'ID=="'"$snp"'"' "$vcf" | head -n1 || true)
		[[ -n ${bp:-} ]] || { echo -e "$t\t$cn\t$snp\t$bp0\tnot_in_1000G" >> "$res/pick/lead_1000G.fail.tsv"; continue; }

		if [[ $cn -eq 23 ]]; then
			par=nonPar; ((bp>=60001 && bp<=2699520)) && par=par; ((bp>=154931044 && bp<=155260560)) && par=par
			pbase=$g1k/EUR.male.chrX.$par
			echo -e "$t\t23\t$snp\t$bp\t$par\t$pbase" >> "$res/pick/chrX_PAR.tsv"
		else
			pbase=$g1k/EUR.chr$c
		fi

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$bp" "$cn" "$bp" "$snp" >> "$tmp/ld.tsv"
		plink2 --pfile "$pbase" --allow-extra-chr --r2-unphased --ld-snp "$snp" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --out "$tmp/$snp" > "$tmp/$snp.plink2.log" 2>&1 || { echo -e "$t\t$cn\t$snp\t$bp\tplink2_ld_failed" >> "$res/pick/lead_1000G.fail.tsv"; continue; }

		outvcor=$tmp/$snp.vcor
		[[ -f $outvcor ]] && awk -v t="$t" -v c="$cn" -v s="$snp" -v b="$bp" '
			BEGIN{FS="[ \t]+";OFS="\t"}
			NR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="POS_B") ip=i; if($i=="ID_B") is=i; if($i~/R2$/) ir=i}; next}
			ip&&is&&ir{print t,c,s,b,c,$ip,$is,$ir}
		' "$outvcor" >> "$tmp/ld.tsv"
	done

	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR==1{next}
	{
		k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0
		if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp
		if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos
		u[k FS bp]=u[k FS pos]=1
	}
	END{
		print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"
		for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}
	}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
	log "DONE $t chr$c"
}

for t in $traits; do
	log "========== DO $t =========="
	for c in $chrs; do run_chr "$t" "$c" & while [[ $(jobs -rp | wc -l) -ge $njob ]]; do sleep 2; done; done
	wait

	d=$res/ld/$t
	for x in ld block; do
		first=1; : > "$d/$x.tsv"
		for c in $chrs; do
			f=$d/chr$c/$x.tsv; [[ -f $f ]] || continue
			[[ $first -eq 1 ]] && { cat "$f"; first=0; } || awk 'NR>1' "$f"
		done > "$d/$x.tsv"
		awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"
	done
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$res/pick/pick.tsv"
	log "========== DONE $t =========="
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Extract EUR 1000G and archaic VCFs, convert to genotype matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$res/pick/pick.tsv" | while IFS=$'\t' read -r tr chr snp bp st en n size_bp; do
	cl=$chr; [[ $chr -eq 23 ]] && cl=X
	lo=$st; hi=$en; [[ $bp -lt $lo ]] && lo=$bp; [[ $bp -gt $hi ]] && hi=$bp
	d=$res/coreVcf/$tr/${chr}.${snp}.${bp}; mkdir -p "$d"

	if [[ $chr -eq 23 ]]; then
		par=nonPar; ((bp>=60001 && bp<=2699520)) && par=par; ((bp>=154931044 && bp<=155260560)) && par=par
		kg=$g1k/EUR.male.chrX.$par.vcf.gz; sample_opt=""
	else
		kg=$g1k/ALL.chr${cl}.vcf.gz; [[ -f $kg ]] || kg=$g1k/ALL.chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		sample_opt="-S $g1k/EUR.sample.ids"
	fi
	vi=$vindija/chr${cl}_mq25_mapab100.vcf.gz; al=$altai/chr${cl}_mq25_mapab100.vcf.gz; ch=$chagyr/chr${cl}.noRB.bgz.vcf.gz; de=$denisova/chr${cl}_mq25_mapab100.vcf.gz
	for f in "$kg" "$vi" "$al" "$ch" "$de"; do idx_vcf "$f"; done
	ck=$(vcf_chr "$kg" "$cl"); cv=$(vcf_chr "$vi" "$cl"); ca=$(vcf_chr "$al" "$cl"); cc=$(vcf_chr "$ch" "$cl"); cd=$(vcf_chr "$de" "$cl")

	bcftools view $sample_opt -r "$ck:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg"; tabix -f -p vcf "$d/kg.vcf.gz"
	bcftools view -r "$cv:$lo-$hi" -Oz -o "$d/vindija.vcf.gz" "$vi"; tabix -f -p vcf "$d/vindija.vcf.gz"
	bcftools view -r "$ca:$lo-$hi" -Oz -o "$d/altai.vcf.gz" "$al"; tabix -f -p vcf "$d/altai.vcf.gz"
	bcftools view -r "$cc:$lo-$hi" -Oz -o "$d/chagyr.vcf.gz" "$ch"; tabix -f -p vcf "$d/chagyr.vcf.gz"
	bcftools view -r "$cd:$lo-$hi" -Oz -o "$d/denisova.vcf.gz" "$de"; tabix -f -p vcf "$d/denisova.vcf.gz"
done

for t in $traits; do
	for d in "$res/coreVcf/$t"/*.*.*; do
		[[ -d $d ]] || continue
		o=$res/mat/$t/$(basename "$d"); mkdir -p "$o"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv"
		for a in vindija altai chagyr denisova; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/$a.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{if($5=="0" || $5=="1") $5=$5"/"$5; print}' > "$o/$a.tsv"; done
	done
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Define risk haplotypes and write PHYLIP files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/make_hap.R" "$res"
Rscript "$dirscript/make_phy.R" "$res"
find "$res/phy" -name '*.phy' | sort | while read -r f; do phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100; done
Rscript "$dirscript/plot_tree.R" "$res"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: 优化、打磨一遍🔪
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
refine_p=1e-6
refine_flank=1000000
refine_nmax=50
res1=$res
res2=$res/refine
mkdir -p "$res2"/{todo,map,lead,risk}

awk -v flank="$refine_flank" 'BEGIN{FS=OFS="\t"}
NR==1{ for(i=1;i<=NF;i++) h[$i]=i
	print "trait","lead_chr","lead_snp","lead_bp","lo","hi","best_lineage"
	next }
$(h["best_lineage"])!="NA" && $(h["best_lineage"])!=""{
	lo=$(h["core_start"])-flank; if(lo<1) lo=1
	hi=$(h["core_end"])+flank
	print $(h["trait"]),$(h["lead_chr"]),$(h["lead_snp"]),$(h["lead_bp"]),lo,hi,$(h["best_lineage"])
}' "$res1/region_summary.tsv" > "$res2/todo/refine.todo.tsv"

printf "trait\torig_chr\torig_snp\torig_bp\tSNP\tPOS37\tREF37\tALT37\n" > "$res2/map/map.tsv"

awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$res2/todo/refine.todo.tsv" | while IFS=$'\t' read -r tr chr snp bp lo hi lin; do
	cl=$chr; [[ $chr -eq 23 ]] && cl=X
	if [[ $chr -eq 23 ]]; then
		par=nonPar
		((bp>=60001 && bp<=2699520)) && par=par
		((bp>=154931044 && bp<=155260560)) && par=par
		vcf=$g1k/EUR.male.chrX.$par.vcf.gz
	else
		vcf=$g1k/ALL.chr${cl}.vcf.gz
		[[ -f $vcf ]] || vcf=$g1k/ALL.chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	fi

	idx_vcf "$vcf"
	ctg=$(vcf_chr "$vcf" "$cl")

	bcftools view -r "$ctg:$lo-$hi" -m2 -M2 -v snps "$vcf" | bcftools query -f '%ID\t%POS\t%REF\t%ALT\n' | \
	awk -v tr="$tr" -v chr="$chr" -v os="$snp" -v ob="$bp" 'BEGIN{FS=OFS="\t"} $1!="."{print tr,chr,os,ob,$1,$2,$3,$4}' >> "$res2/map/map.tsv"
done

Rscript "$dirscript/rerun_lead.R" --res1 "$res1" --res2 "$res2" \
	--dirgwas "$dirgwas" --traits "$traits" --p "$refine_p" --nmax "$refine_nmax"

res=$res2
mkdir -p "$res" "$res/ld" "$res/pick" "$res/coreVcf" "$res/mat"
# 然后重新运行s1-3 🐧