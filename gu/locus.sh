#!/usr/bin/env bash

dir0=/mnt/d
dirmod=$dir0/refGen/1kg_phase3
arch0=$dir0/refGen/gu
dirscript=$dir0/scripts/gu
dirgwas=$dir0/data/gwas/bald
dirout=$dir0/analysis/gu/locus
traits="bald bald12 bald13 bald14"
chrs="$(seq 1 22) X"
ld_r2=0.98
njob=8
plink_threads=1

mkdir -p "$dirout"/{lead,ld,coreVcf,mat,hap,phy,plot,report,log}  # ❗这儿修改了：新增 plot/report，summary TSV 不再散落在 dirout 根目录。
log(){ echo "[$(date '+%F %T')] $*"; }
clean_msg(){ [[ -f $1 ]] && tail -n 12 "$1" | tr '\t\r\n' '   ' | sed 's/  */ /g' | cut -c1-700; }
id_mode(){ awk 'BEGIN{IGNORECASE=1} $1!="" && $1!="."{n++; if($1~/^rs[0-9]+$/) rs++; else if(toupper($1)~/^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$/) cp++; if(n==5000) exit} END{if(n==0) print "missing"; else if(rs/n>.8) print "rsid"; else if(cp/n>.8) print "chrpos"; else print "other"}'; }
check_pbase(){
	local pbase=$1 assoc=$2 label=$3 flag=$4 m1 m2 e1 e2 ns nv
	[[ -s $pbase.pgen && -s $pbase.pvar && -s $pbase.psam ]] || { echo "ERROR: missing pgen/pvar/psam for $pbase"; exit 1; }
	[[ -f $flag ]] && return
	m1=$(awk '!/^#/ && $3!=""{print $3}' "$pbase.pvar" | id_mode)
	m2=$(awk 'NR>1{print $3}' "$assoc" | id_mode)
	e1=$(awk '!/^#/ && $3!="." && $3!=""{print $3; if(++n==5) exit}' "$pbase.pvar" | paste -sd, -)
	e2=$(awk 'NR>1{print $3; if(++n==5) exit}' "$assoc" | paste -sd, -)
	ns=$(awk 'NR>1{n++} END{print n+0}' "$pbase.psam")
	nv=$(awk '!/^#/{n++} END{print n+0}' "$pbase.pvar")
	printf "label\tlead_mode\tlead_examples\tG1000_mode\tG1000_examples\tn_samples\tn_variants\tpbase\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$m2" "$e2" "$m1" "$e1" "$ns" "$nv" "$pbase" > "$flag"
	[[ $m1 == rsid || $m1 == chrpos ]] || { echo "ERROR: 1000G pvar ID mode is $m1 for $pbase.pvar. Expected rsID or CHR:POS:REF:ALT. Fix the VCF/pfile IDs first."; exit 1; }
}
write_fail(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$dirout/lead/lead_1000G.fail.tsv"; }
write_debug(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" >> "$dirout/lead/ld_debug.tsv"; }

for x in Rscript bcftools bgzip tabix plink2 phyml; do
	command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }
done
if [[ " $chrs " == *" X "* ]]; then
	for f in \
		"$dirmod/chrX.vcf.gz" "$dirmod/chrX.vcf.gz.tbi" \
		"$dirmod/EUR.chrX.pgen" "$dirmod/EUR.chrX.pvar" "$dirmod/EUR.chrX.psam" \
		"$dirmod/EUR.male.chrX.par.pgen" "$dirmod/EUR.male.chrX.par.pvar" "$dirmod/EUR.male.chrX.par.psam" \
		"$dirmod/EUR.male.chrX.par.vcf.gz" "$dirmod/EUR.male.chrX.par.vcf.gz.tbi" \
		"$dirmod/EUR.male.chrX.nonPar.pgen" "$dirmod/EUR.male.chrX.nonPar.pvar" "$dirmod/EUR.male.chrX.nonPar.psam" \
		"$dirmod/EUR.male.chrX.nonPar.vcf.gz" "$dirmod/EUR.male.chrX.nonPar.vcf.gz.tbi"; do
		[[ -s $f ]] || { log "ERROR missing required chrX file: $f"; exit 1; }
	done
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Prepare COJO lead SNPs and define high-LD core blocks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 从 *.jma.cojo 读取 COJO lead SNPs；从 GWAS summary 读取 BETA / effect allele / other allele；
# 用 1000G pvar 校正 lead SNP ID 和 b37 POS； 每个 trait每条染色体单独用 PLINK 计算 r2 >= $ld_r2 的 high-LD block；
# 主结果写入 lead/pick.tsv；只有 lead 本身的 single-SNP block 另存到 lead/pick.single.tsv。

log "s1 prep_input start"
Rscript "$dirscript/prep_input.R" --dirgwas "$dirgwas" --dirout "$dirout" --dirmod "$dirmod" --traits $traits
log "s1 prep_input done: $(wc -l < "$dirout/lead/lead.assoc") lines in lead.assoc"

# Manual extra locus. Use the 1000G/PLINK b37 POS, not the GRCh38/GWAS POS.
if ! awk 'NR>1 && $2==3 && $3=="rs35044562"{f=1} END{exit !f}' "$dirout/lead/bald.lead.assoc"; then
	awk 'BEGIN{FS=OFS="\t"} NR==1{print; print "bald",3,"rs35044562",45909024,"G","A",1,"G","rs35044562","ID_manual"; next} {print}' "$dirout/lead/bald.lead.assoc" > "$dirout/lead/.bald.lead.assoc" && mv "$dirout/lead/.bald.lead.assoc" "$dirout/lead/bald.lead.assoc"
	awk 'BEGIN{FS=OFS="\t"} NR==1{print; print 3,"rs35044562",45909024; next} {print}' "$dirout/lead/bald.lead.3col" > "$dirout/lead/.bald.lead.3col" && mv "$dirout/lead/.bald.lead.3col" "$dirout/lead/bald.lead.3col"
	log "manual add bald rs35044562 at 1000G POS 45909024"
fi

first=${traits%% *}
{ head -n1 "$dirout/lead/$first.lead.assoc"; for t in $traits; do awk 'NR>1' "$dirout/lead/$t.lead.assoc"; done; } > "$dirout/lead/lead.assoc"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\treason\tdetail\tpbase\tvid\n" > "$dirout/lead/lead_1000G.fail.tsv"
printf "trait\tlead_chr\tlead_snp\tlead_bp\tPAR\tLD_pfile\n" > "$dirout/lead/chrX_PAR.tsv"
printf "trait\tchr\tsnp\tlead_bp\tpvar_bp\tvid\tmatch\tpbase\tstatus\tn_ld_snp\tmessage\n" > "$dirout/lead/ld_debug.tsv"

run_chr(){
	local t=$1 c=$2 cn=$2 assoc tmp snp bp ea oa par pbase flag line vid ref alt refbp match n_id outpre plog outvcor rows msg cmdtxt
	[[ $c == X ]] && cn=23
	assoc=$dirout/lead/$t.lead.assoc
	tmp=$dirout/ld/$t/chr$c
	mkdir -p "$tmp"

	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $2==cn || toupper($2)==c' "$assoc" > "$tmp/lead.assoc"
	log "LD $t chr$c: $(awk 'NR>1' "$tmp/lead.assoc" | wc -l) lead SNPs"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"

	awk 'BEGIN{FS=OFS="\t"} NR>1{print $3,$4,$5,$6}' "$tmp/lead.assoc" | while IFS=$'\t' read -r snp bp ea oa; do
		[[ -n $snp && -n $bp ]] || continue
		if [[ $c == X ]]; then
			par=nonPar
			((bp>=60001 && bp<=2699520)) && par=par
			((bp>=154931044 && bp<=155260560)) && par=par
			pbase=$dirmod/EUR.male.chrX.$par
			echo -e "$t\t23\t$snp\t$bp\t$par\t$pbase" >> "$dirout/lead/chrX_PAR.tsv"
		else
			pbase=$dirmod/EUR.chr$c
		fi
		flag="$tmp/sanity.$(basename "$pbase").tsv"
		check_pbase "$pbase" "$tmp/lead.assoc" "$t chr$c" "$flag"

		line=$(awk -v snp="$snp" -v bp="$bp" 'BEGIN{FS=OFS="\t"} !/^#/ && $2==bp && $3==snp{print $2,$3,$4,$5; exit}' "$pbase.pvar")
		if [[ -n $line ]]; then read -r refbp vid ref alt <<< "$line"; match=ID_BP; else
			n_id=$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} !/^#/ && $3==snp{n++} END{print n+0}' "$pbase.pvar")
			if [[ $n_id -eq 1 ]]; then
				read -r refbp vid ref alt <<< "$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} !/^#/ && $3==snp{print $2,$3,$4,$5; exit}' "$pbase.pvar")"
				match=ID_only
				[[ $refbp != "$bp" ]] && log "WARN $t chr$c $snp: lead_bp=$bp but pvar_bp=$refbp; using pvar_bp"
			else
				line=$(awk -v bp="$bp" -v ea="$ea" -v oa="$oa" 'BEGIN{FS=OFS="\t"} !/^#/ && $2==bp{ref=toupper($4); alt=toupper($5); ea=toupper(ea); oa=toupper(oa); if((ref==ea && alt==oa) || (ref==oa && alt==ea)){print $2,$3,$4,$5; exit}}' "$pbase.pvar")
				if [[ -n $line ]]; then read -r refbp vid ref alt <<< "$line"; match=ALLELE_BP; else
					line=$(awk -v bp="$bp" 'BEGIN{FS=OFS="\t"} !/^#/ && $2==bp{id=$3; r=$4; a=$5; n++} END{if(n==1) print bp,id,r,a}' "$pbase.pvar")
					if [[ -n $line ]]; then read -r refbp vid ref alt <<< "$line"; match=UNIQUE_BP; else
						msg="not found by ID+BP, ID-only, allele+BP, or unique BP"
						write_fail "$t" "$cn" "$snp" "$bp" "not_in_1000G_pvar_or_allele_mismatch" "$msg" "$pbase" "NA"
						write_debug "$t" "$cn" "$snp" "$bp" "NA" "NA" "none" "$pbase" "pvar_miss" 0 "$msg"
						continue
					fi
				fi
			fi
		fi

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$refbp" "$cn" "$refbp" "$vid" >> "$tmp/ld.tsv"
		outpre="$tmp/$snp"; plog="$outpre.plink2.log"; outvcor="$outpre.vcor"
		cmdtxt="plink2 --pfile $pbase --allow-extra-chr --threads $plink_threads --make-founders --r2-unphased allow-ambiguous-allele --ld-snp $vid --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 $ld_r2 --out $outpre"
		printf "%s\n" "$cmdtxt" > "$outpre.cmd"
		# log "RUN $t chr$c $snp: match=$match vid=$vid lead_bp=$bp pvar_bp=$refbp"
		if plink2 --pfile "$pbase" --allow-extra-chr --threads "$plink_threads" --make-founders --r2-unphased allow-ambiguous-allele --ld-snp "$vid" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --out "$outpre" > "$plog" 2>&1; then
			if [[ -s $outvcor ]]; then
				awk -v t="$t" -v c="$cn" -v s="$snp" -v b="$refbp" '
					BEGIN{FS="[ \t]+";OFS="\t"}
					NR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="POS_B") ip=i; if($i=="ID_B") is=i; if($i~/R2$/) ir=i}; next}
					ip&&is&&ir{print t,c,s,b,c,$ip,$is,$ir}
				' "$outvcor" >> "$tmp/ld.tsv"
				rows=$(awk 'NR>1{n++} END{print n+0}' "$outvcor")
				write_debug "$t" "$cn" "$snp" "$bp" "$refbp" "$vid" "$match" "$pbase" "ok" "$rows" "plink_ok"
			else
				msg="PLINK finished but $outvcor is missing/empty"
				write_debug "$t" "$cn" "$snp" "$bp" "$refbp" "$vid" "$match" "$pbase" "ok_no_vcor" 0 "$msg"
			fi
		else
			msg=$(clean_msg "$plog")
			write_fail "$t" "$cn" "$snp" "$bp" "plink2_ld_failed" "$msg" "$pbase" "$vid"
			write_debug "$t" "$cn" "$snp" "$bp" "$refbp" "$vid" "$match" "$pbase" "plink_failed" 0 "$msg"
			log "FAIL $t chr$c $snp: $msg"
		fi
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
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7==1' "$d/block.tsv" >> "$dirout/lead/pick.single.tsv"
	log "========== DONE $t: $(awk 'NR>1{n++} END{print n+0}' "$d/block.tsv") blocks; pick_total=$(awk 'NR>1{n++} END{print n+0}' "$dirout/lead/pick.tsv"); single_total=$(awk 'NR>1{n++} END{print n+0}' "$dirout/lead/pick.single.tsv") =========="
done

log "s1 done. Check: $dirout/lead/ld_debug.tsv and $dirout/lead/lead_1000G.fail.tsv"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Extract core-region VCFs and project archaic calls to 1000G SNP templates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 根据 s1 的 lead/pick.tsv，逐个 trait-specific high-LD block 提取 1000G EUR template VCF；
# archaic VCF 先按同一区域提取 raw records，再用 prep_input_archaic.R 投影到 1000G REF/ALT 坐标；
# 这样 ALT="." 的 archaic reference-only sites 会保留为 GT=0/0，而不会原样进入 IBDmix/下游矩阵。
# 最后转成 mat/$trait/$locus/*.tsv，供 risk haplotype 和 phylogeny 使用。
# 这里只去掉完全重复的 trait-specific block，不按 SNP 全局去重，因为同一 SNP 可属于多个 trait。

while IFS=$'\t' read -r tr chr snp bp st en n size_bp; do
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
		kg=$dirmod/EUR.male.chrX.$par.vcf.gz
		[[ -s $kg ]] || { log "ERROR missing chrX 1000G VCF for $tr $snp: $kg"; exit 1; }
		bcftools view -r "X:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; exit 1; }
	else
		kg=$dirmod/chr$cl.vcf.gz
		[[ -s $kg ]] || kg=$dirmod/chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		[[ -s $kg ]] || { log "ERROR missing 1000G VCF for $tr chr$chr $snp"; exit 1; }
		bcftools view -S "$dirmod/EUR.1id" -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; exit 1; }
	fi
	tabix -f -p vcf "$d/kg.vcf.gz"

	for ad in "$arch0"/*; do
		[[ -d $ad ]] || continue
		name=$(basename "$ad")
		outname=$(echo "$name" | tr '[:upper:]' '[:lower:]')
		avcf=""
		for f in "$ad"/*chr${cl}_*.vcf.gz "$ad"/*chr${cl}.*.vcf.gz; do
			[[ -s $f ]] || continue
			avcf=$f
			break
		done
		[[ -n $avcf ]] || { log "WARN missing archaic VCF: $outname chr$cl"; continue; }

		arc_raw="$d/$outname.raw.vcf.gz"
		arc_vcf="$d/$outname.vcf.gz"
		arc_log="$d/$outname.project.log"

		# 先提取 raw archaic records，再用 1000G kg.vcf.gz 作为 template 进行 REF/ALT 投影。
		bcftools view -r "$cl:$lo-$hi" -Oz -o "$arc_raw" "$avcf" || { log "ERROR bcftools failed: archaic raw $outname $cl:$lo-$hi"; exit 1; }
		tabix -f -p vcf "$arc_raw"

		Rscript "$dirscript/prep_input_archaic.R" "$d/kg.vcf.gz" "$arc_raw" "$arc_vcf" "$outname" > "$arc_log" 2>&1 || {
			log "ERROR archaic projection failed: $outname $cl:$lo-$hi; see $arc_log"
			tail -20 "$arc_log"
			exit 1
		}

		nvar=$(bcftools view -H "$arc_vcf" | wc -l)
		dot=$(bcftools view -H "$arc_vcf" | awk 'BEGIN{FS="\t"} $5=="."{n++} END{print n+0}')
		[[ $nvar -gt 0 && $dot -eq 0 ]] || { log "ERROR bad projected archaic VCF: $arc_vcf n=$nvar ALT_dot=$dot"; exit 1; }
	done
done < <(awk 'BEGIN{FS=OFS="	"} NR>1 && !seen[$1 FS $2 FS $3 FS $4 FS $5 FS $6]++{print}' "$dirout/lead/pick.tsv")

for t in $traits; do
	for d in "$dirout/coreVcf/$t"/*.*.*; do
		[[ -d $d ]] || continue
		o=$dirout/mat/$t/$(basename "$d")
		mkdir -p "$o"
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv"
		for f in "$d"/*.vcf.gz; do
			name=$(basename "$f" .vcf.gz)
			[[ $name == kg || $name == *.raw ]] && continue
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$f" | awk 'BEGIN{FS=OFS="\t"} $4!="." && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ {if($5=="0" || $5=="1") $5=$5"/"$5; print}' > "$o/$name.tsv"
		done
	done
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Build risk haplotypes and phylogenetic trees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于 s2 的 genotype matrices 和 lead.assoc 里的 risk allele，构建 1000G risk haplotypes；
# 与 archaic reference 对比后写出 PHYLIP 文件，运行 phyml，并绘制 tree。
Rscript "$dirscript/make_hap.R" "$dirout"
find "$dirout" -maxdepth 1 -type f -name '*.tsv' -exec mv -f {} "$dirout/report/" \;
Rscript "$dirscript/make_phy.R" "$dirout"
log "s3 PHYLIP files: $(find "$dirout/phy" -name '*.phy' | wc -l)"
find "$dirout/phy" -name '*.phy' | sort | while read -r f; do phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100; done
Rscript "$dirscript/make_tree.R" "$dirout"
log "s3 tree PNG files: $(find "$dirout/plot" -maxdepth 1 -name 's8_tree_*.png' | wc -l)"
