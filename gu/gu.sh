#!/usr/bin/bash

dir0=/mnt/d
g1k=$dir0/data/refGen/1kg_phase3
arch0=$dir0/data/refGen/gu
dirscript=$dir0/scripts/gu
dirgwas=$dir0/data/gwas/bald
dirout=$dir0/analysis/gu
traits="bald bald12 bald13 bald14"
chrs="$(seq 1 22) X"
ld_r2=0.98
njob=8

mkdir -p "$dirout"/{lead,ld,coreVcf,mat,hap,phy}
log(){ echo "[$(date '+%F %T')] $*"; }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: COJO lead SNPs + GWAS alleles; define 1000G high-LD core blocks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/input_prep.R" --dirgwas "$dirgwas" --dirout "$dirout" --traits $traits

if ! awk 'BEGIN{FS=OFS="\t"} NR>1 && $2=="rs35044562"{f=1} END{exit !f}' "$dirout/lead/bald.lead.3col"; then
	sed -i '1a 3\trs35044562\t45859651' "$dirout/lead/bald.lead.3col"
fi
if ! awk 'BEGIN{FS=OFS="\t"} NR>1 && $3=="rs35044562"{f=1} END{exit !f}' "$dirout/lead/bald.lead.assoc"; then
	sed -i '1a bald\t3\trs35044562\t45859651\tG\tA\t1\tG' "$dirout/lead/bald.lead.assoc"
fi

first=${traits%% *}
{ head -n1 "$dirout/lead/$first.lead.assoc"; for t in $traits; do awk 'NR>1' "$dirout/lead/$t.lead.assoc"; done; } > "$dirout/lead/lead.assoc"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\treason\n" > "$dirout/lead/lead_1000G.fail.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tPAR\tLD_pfile\n" > "$dirout/lead/chrX_PAR.tsv"

run_chr(){
	local t=$1 c=$2 cn=$2 lead tmp snp bp par pbase qbase vid outvcor
	[[ $c == X ]] && cn=23
	lead=$dirout/lead/$t.lead.3col
	tmp=$dirout/ld/$t/chr$c
	mkdir -p "$tmp"

	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $1==cn || toupper($1)==c' "$lead" > "$tmp/lead.tsv"
	log "LD $t chr$c: $(awk 'NR>1' "$tmp/lead.tsv" | wc -l) lead SNPs"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"

	awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3}' "$tmp/lead.tsv" | while IFS=$'\t' read -r snp bp; do
		if [[ $c == X ]]; then
			par=nonPar
			((bp>=60001 && bp<=2699520)) && par=par
			((bp>=154931044 && bp<=155260560)) && par=par
			pbase=$g1k/EUR.male.chrX.$par
			echo -e "$t\t23\t$snp\t$bp\t$par\t$pbase" >> "$dirout/lead/chrX_PAR.tsv"
		else
			pbase=$g1k/EUR.chr$c
		fi

		qbase=$pbase.id
		if [[ ! -f $qbase.pgen ]]; then
			plink2 --pfile "$pbase" --set-all-var-ids '@:#:$r:$a' --make-pgen --out "$qbase"
		fi

		vid=$(awk -v bp="$bp" 'BEGIN{FS=OFS="\t"} NR>1 && $2==bp{print $3; exit}' "$qbase.pvar")
		[[ -n ${vid:-} ]] || { echo -e "$t\t$cn\t$snp\t$bp\tnot_in_1000G_pvar" >> "$dirout/lead/lead_1000G.fail.tsv"; continue; }

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$bp" "$cn" "$bp" "$vid" >> "$tmp/ld.tsv"
		plink2 --pfile "$qbase" --allow-extra-chr --r2-unphased --ld-snp "$vid" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --out "$tmp/$snp" > "$tmp/$snp.plink2.log" 2>&1 || { echo -e "$t\t$cn\t$snp\t$bp\tplink2_ld_failed" >> "$dirout/lead/lead_1000G.fail.tsv"; continue; }

		outvcor=$tmp/$snp.vcor
		[[ -f $outvcor ]] && awk -v t="$t" -v c="$cn" -v s="$snp" -v b="$bp" '
			BEGIN{FS="[ \t]+";OFS="\t"}
			NR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="POS_B") ip=i; if($i=="ID_B") is=i; if($i~/R2$/) ir=i}; next}
			ip&&is&&ir{print t,c,s,b,c,$ip,$is,$ir}
		' "$outvcor" >> "$tmp/ld.tsv"
	done

	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp; if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
}

for t in $traits; do
	log "========== DO $t =========="
	for c in $chrs; do
		run_chr "$t" "$c" &
		while [[ $(jobs -rp | wc -l) -ge $njob ]]; do sleep 2; done
	done
	wait

	d=$dirout/ld/$t
	for x in ld block; do
		first=1; : > "$d/$x.tsv"
		for c in $chrs; do
			f=$d/chr$c/$x.tsv
			[[ -f $f ]] || continue
			[[ $first -eq 1 ]] && { cat "$f"; first=0; } || awk 'NR>1' "$f"
		done > "$d/$x.tsv"
		awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"
	done
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$dirout/lead/pick.tsv"
	log "========== DONE $t =========="
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Extract 1000G EUR and archaic VCFs; convert to genotype matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$dirout/lead/pick.tsv" | while IFS=$'\t' read -r tr chr snp bp st en n size_bp; do
	cl=$chr
	lo=$st; hi=$en
	[[ $bp -lt $lo ]] && lo=$bp
	[[ $bp -gt $hi ]] && hi=$bp
	d=$dirout/coreVcf/$tr/${chr}.${snp}.${bp}
	mkdir -p "$d"

	if [[ $chr == 23 || $chr == X ]]; then
		cl=X
		par=nonPar
		((bp>=60001 && bp<=2699520)) && par=par
		((bp>=154931044 && bp<=155260560)) && par=par
		kg=$g1k/EUR.male.chrX.$par.vcf.gz
		bcftools view -r "X:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg"
	else
		kg=$g1k/ALL.chr$cl.vcf.gz
		[[ -f $kg ]] || kg=$g1k/ALL.chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		bcftools view -S "$g1k/EUR.sample.1id" -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg"
	fi
	tabix -f -p vcf "$d/kg.vcf.gz"

	for ad in "$arch0"/*; do
		[[ -d $ad ]] || continue
		name=$(basename "$ad")
		outname=$(echo "$name" | tr '[:upper:]' '[:lower:]')
		avcf=$(ls "$ad"/*chr${cl}_*.vcf.gz "$ad"/*chr${cl}.*.vcf.gz 2>/dev/null | head -n1)
		[[ -n ${avcf:-} ]] || continue
		bcftools view -r "$cl:$lo-$hi" -Oz -o "$d/$outname.vcf.gz" "$avcf"
		tabix -f -p vcf "$d/$outname.vcf.gz"
	done
done

for t in $traits; do
	for d in "$dirout/coreVcf/$t"/*.*.*; do
		[[ -d $d ]] || continue
		o=$dirout/mat/$t/$(basename "$d")
		mkdir -p "$o"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv"
		for f in "$d"/*.vcf.gz; do
			name=$(basename "$f" .vcf.gz)
			[[ $name == kg ]] && continue
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$f" | awk 'BEGIN{FS=OFS="\t"}{if($5=="0" || $5=="1") $5=$5"/"$5; print}' > "$o/$name.tsv"
		done
	done
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Define risk haplotypes, write PHYLIP files, run phylogeny
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rscript "$dirscript/make_hap.R" "$dirout"
Rscript "$dirscript/make_phy.R" "$dirout"
find "$dirout/phy" -name '*.phy' | sort | while read -r f; do phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100; done
Rscript "$dirscript/plot_tree.R" "$dirout"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: Refine lead SNPs around first-pass archaic loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
refine_p=1e-6
refine_flank=1000000
refine_nmax=50
res1=$dirout
dirout2=$dirout/refine
mkdir -p "$dirout2"/{todo,map,lead,ld,coreVcf,mat,hap,phy}

awk -v flank="$refine_flank" 'BEGIN{FS=OFS="\t"}
NR==1{
	for(i=1;i<=NF;i++) h[$i]=i
	print "trait","lead_chr","lead_snp","lead_bp","lo","hi","best_lineage"
	next
}
$(h["best_lineage"])!="NA" && $(h["best_lineage"])!=""{
	lo=$(h["core_start"])-flank
	if(lo<1) lo=1
	hi=$(h["core_end"])+flank
	print $(h["trait"]),$(h["lead_chr"]),$(h["lead_snp"]),$(h["lead_bp"]),lo,hi,$(h["best_lineage"])
}' "$res1/region_summary.tsv" > "$dirout2/todo/refine.todo.tsv"

printf "trait\torig_chr\torig_snp\torig_bp\tSNP\tPOS37\tREF37\tALT37\n" > "$dirout2/map/map.tsv"

awk 'BEGIN{FS=OFS="\t"} NR>1{print}' "$dirout2/todo/refine.todo.tsv" | while IFS=$'\t' read -r tr chr snp bp lo hi lin; do
	cl=$chr
	[[ $chr -eq 23 ]] && cl=X

	vcf=$g1k/ALL.chr$cl.vcf.gz
	[[ -f $vcf ]] || vcf=$g1k/ALL.chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

	bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps "$vcf" | \
		bcftools query -f '%ID\t%POS\t%REF\t%ALT\n' | \
		awk -v tr="$tr" -v chr="$chr" -v os="$snp" -v ob="$bp" 'BEGIN{FS=OFS="\t"} $1!="."{print tr,chr,os,ob,$1,$2,$3,$4}' \
		>> "$dirout2/map/map.tsv"
done

Rscript "$dirscript/rerun_lead.R" \
	--res1 "$res1" \
	--res2 "$dirout2" \
	--dirgwas "$dirgwas" \
	--traits "$traits" \
	--p "$refine_p" \
	--nmax "$refine_nmax"
	
运行完后会生成：
$dirout/refine/lead/bald.lead.3col
$dirout/refine/lead/bald.lead.assoc
$dirout/refine/lead/lead.assoc
$dirout/refine/candidate_map.tsv

一个提醒：这一步仍然依赖 map/map.tsv 里的 SNP 能和 GWAS summary 里的 SNP 匹配。如果你的 1000G VCF 的 %ID 全是 .，那么 map.tsv 会没有候选 SNP，rerun_lead.R 就无法生成 refined candidates。此时需要另外准备一个 SNP → POS37 的 map 文件。