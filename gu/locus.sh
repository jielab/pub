#!/usr/bin/env bash

# ============================================================
# 全局配置：输入/输出目录、分析 trait、染色体范围、LD 阈值和并行任务数。
# 默认只分析 bald；START_STEP=s1 从 lead SNP/LD 开始，START_STEP=s2 复用已有 LD 结果。
# ============================================================
dir0=/mnt/d
export PATH="$dir0/software/bin:$PATH"
dirmod=$dir0/refGen/1kg_phase3
arch0=$dir0/refGen/gu
dirscript=$dir0/scripts/gu
dirgwas=$dir0/data/gwas/bald
dirout=$dir0/analysis/gu/locus
sample_file=$dir0/files/1kg.v3.sample.txt
traits="bald"
read -r -a traits_arr <<< "$traits"
chrs="$(seq 1 22) X"
ld_r2=0.98
job_of_trait=1
job_in_trait=12
job_phyml=2
filter_pop=${FILTER_POP:-YRI}
filter_max_count=${FILTER_MAX_COUNT:-1}
start_step=${START_STEP:-s1}

mkdir -p "$dirout"/{lead,ld,coreVcf,mat,hap,phy,plot,report,log}
exec > >(tee "$dirout/log/locus.log") 2>&1

# ============================================================
# 通用函数：日志、文件行数、结果完整性检查、1000G ID/POS 匹配和失败记录。
# summary 也走 log，因此会同时打印到屏幕并写入 log/locus.log。
# ============================================================
log(){ echo "[$(date '+%F %T')] $*"; }
log "START_STEP=$start_step"
log "parallel jobs: job_of_trait=$job_of_trait job_in_trait=$job_in_trait job_phyml=$job_phyml"
log "haplotype population filter: $filter_pop <= $filter_max_count; sample_file=$sample_file"
clean_msg(){ [[ -f $1 ]] && tail -n 12 "$1" | tr '\t\r\n' '   ' | sed 's/  */ /g' | cut -c1-700; }
id_mode(){ awk 'BEGIN{IGNORECASE=1} $1!="" && $1!="."{n++; if($1~/^rs[0-9]+$/) rs++; else if(toupper($1)~/^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$/) cp++; if(n==5000) exit} END{if(n==0) print "missing"; else if(rs/n>.8) print "rsid"; else if(cp/n>.8) print "chrpos"; else print "other"}'; }
data_rows(){ [[ -s $1 ]] && awk 'NR>1{n++} END{print n+0}' "$1" || echo 0; }
trait_rows(){ [[ -s $1 ]] && awk -v t="$2" 'BEGIN{FS="\t"} NR>1 && $1==t{n++} END{print n+0}' "$1" || echo 0; }
count_files(){ find "$1" -name "$2" 2>/dev/null | wc -l; }
log_trait_count(){
	local label=$1 f=$2 t
	log "summary: $label (all traits; file path $f) $(data_rows "$f")"
	for t in $traits; do
		log "summary: $label (trait $t; file path $f) $(trait_rows "$f" "$t")"
	done
}
clean_report_workbooks_only(){
	find "$dirout/report" -mindepth 1 -type d -exec rm -rf {} +
	find "$dirout/report" -maxdepth 1 -type f \
		! \( -name 'all.xlsx' -o -name 'filtered.xlsx' -o -name 'selected.xlsx' \) \
		-delete
}
valid_header(){ [[ -s $1 ]] && awk 'NR==1 && NF>1{ok=1} END{exit !ok}' "$1"; }
valid_chr_ld(){
	local d=$1 nlead=$2
	[[ -s $d/ld.tsv && -s $d/block.tsv ]] || return 1
	valid_header "$d/ld.tsv" && valid_header "$d/block.tsv" || return 1
	if (( nlead == 0 )); then
		valid_header "$d/ld.tsv" && valid_header "$d/block.tsv"
	else
		[[ $(data_rows "$d/block.tsv") -gt 0 ]]
	fi
}
valid_trait_prep(){
	local t=$1
	[[ -s $dirout/lead/$t.lead.assoc && -s $dirout/lead/$t.lead.3col && -s $dirout/lead/$t.id_sanity.tsv ]] || return 1
	[[ -s $dirout/lead/prep.version ]] && grep -qx 'prep_input_v3_id_first_pos_if_missing_id' "$dirout/lead/prep.version" || return 1
	[[ $(data_rows "$dirout/lead/$t.lead.assoc") -gt 0 ]]
}
valid_mat_locus(){
	local o=$1 expected_arch=$2 n_arch
	[[ -s $o/kg.tsv && $(data_rows "$o/kg.tsv") -gt 0 ]] || return 1
	[[ -s $o/kg.samples.tsv ]] || return 1
	n_arch=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' ! -name 'kg.samples.tsv' -size +0c 2>/dev/null | wc -l)
	(( n_arch >= expected_arch ))
}
check_pbase(){
	local pbase=$1 assoc=$2 label=$3 flag=$4 m1 m2 e1 e2 ns nv
	[[ -s $pbase.pgen && -s $pbase.pvar && -s $pbase.psam ]] || { echo "ERROR: missing pgen/pvar/psam for $pbase"; exit 1; }
	[[ -f $flag ]] && return
	m1=$(awk '!/^#/ && $3!=""{print $3}' "$pbase.pvar" | id_mode)
	m2=$(awk 'NR>1{print $3}' "$assoc" | id_mode)
	e1=$(awk '!/^#/ && $3!="." && $3!=""{print $3; if(++n==5) exit}' "$pbase.pvar" | paste -sd, -)
	e2=$(awk 'NR>1{print $3; if(++n==5) exit}' "$assoc" | paste -sd, -)
	ns=$(awk 'NR>1{n++} END{print n+0}' "$pbase.psam")
	nv=NA
	printf "label\tlead_mode\tlead_examples\tG1000_mode\tG1000_examples\tn_samples\tn_variants\tpbase\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$m2" "$e2" "$m1" "$e1" "$ns" "$nv" "$pbase" > "$flag"
	[[ $m1 == rsid || $m1 == chrpos ]] || { echo "ERROR: 1000G pvar ID mode is $m1 for $pbase.pvar. Expected rsID or CHR:POS:REF:ALT. Fix the VCF/pfile IDs first."; exit 1; }
}
write_fail(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$dirout/lead/lead_1000G.fail.tsv"; }
write_debug(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" >> "$dirout/lead/ld_debug.tsv"; }
pvar_lookup(){
	local pbase=$1 lead_assoc=$2 tmp=$3 cache key
	cache="$tmp/pvar.$(basename "$pbase").lookup.tsv"
	[[ -s $cache ]] && { printf "%s\n" "$cache"; return; }
	key="$tmp/pvar.$(basename "$pbase").keys.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3,$4,$5,$6}' "$lead_assoc" > "$key"
	awk 'BEGIN{FS=OFS="\t"}
		NR==FNR{bp[$3]=1; id[$2]=1; allele[$3 FS toupper($4) FS toupper($5)]=1; allele[$3 FS toupper($5) FS toupper($4)]=1; next}
		/^#/{next}
		(($2 in bp) || ($3 in id) || (($2 FS toupper($4) FS toupper($5)) in allele)){print $2,$3,$4,$5}
	' "$key" "$pbase.pvar" > "$cache"
	rm -f "$key"
	printf "%s\n" "$cache"
}

# ============================================================
# 依赖和参考文件检查：确认 Rscript/bcftools/plink2/phyml 等命令可用。
# 如果包含 chrX，还要确认 1000G X 染色体、PAR/nonPAR 文件齐全。
# ============================================================
for x in Rscript bcftools bgzip tabix plink2 phyml; do
	command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }
done
[[ -s "$dirscript/locus.R" ]] || { log "ERROR missing R entrypoint: $dirscript/locus.R"; exit 1; }
[[ -s $sample_file ]] || { log "ERROR missing 1000G sample metadata: $sample_file"; exit 1; }
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


# ============================================================
# s1：准备 COJO/GWAS lead SNP，并按 1000G EUR LD 定义候选 loci。
# 先用 locus.R prep_input 做 ID/POS/allele 匹配，再按 trait 和染色体运行 PLINK LD。
# 多 SNP 高 LD block 写入 lead/pick.tsv；单 SNP block 写入 lead/pick.single.tsv。
# ============================================================

if [[ "$start_step" == "s1" ]]; then
need_prep=0
for t in $traits; do
	valid_trait_prep "$t" || need_prep=1
done
if (( need_prep == 0 )); then
	log "s1 prep_input skip: existing lead files are complete for traits=$traits"
else
	log "s1 prep_input start"
	Rscript "$dirscript/locus.R" prep_input --dirgwas "$dirgwas" --dirout "$dirout" --dirmod "$dirmod" --traits "${traits_arr[@]}"
	log "s1 prep_input done: $(wc -l < "$dirout/lead/lead.assoc") lines in lead.assoc"
fi

# bald 手动补入 rs35044562：这里必须使用 1000G/PLINK 的 b37 POS，不使用 GRCh38/GWAS POS。
if [[ " $traits " == *" bald "* ]] && ! awk 'NR>1 && $2==3 && $3=="rs35044562"{f=1} END{exit !f}' "$dirout/lead/bald.lead.assoc"; then
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

# 单染色体 LD：把 lead SNP 映射到 1000G pvar，按 pfile 批量运行 plink2 --r2。
# 输出每个 lead SNP 的高 LD SNP 列表 ld.tsv 和对应 block.tsv。
run_chr(){
	local t=$1 c=$2 cn=$2 assoc tmp nlead snp bp ea oa par pbase flag lookup line vid _ref _alt refbp match n_id outpre plog outvcor rows msg cmdtxt
	[[ $c == X ]] && cn=23
	assoc=$dirout/lead/$t.lead.assoc
	tmp=$dirout/ld/$t/chr$c
	mkdir -p "$tmp"

	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $2==cn || toupper($2)==c' "$assoc" > "$tmp/lead.assoc"
	nlead=$(data_rows "$tmp/lead.assoc")
	if valid_chr_ld "$tmp" "$nlead"; then
		log "LD $t chr$c: skip existing complete result ($nlead lead SNPs)"
		return 0
	fi
	log "LD $t chr$c: $nlead lead SNPs"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tvid\tpbase\tlead_bp0\tmatch\n" > "$tmp/lead.map.tsv"
	if (( nlead == 0 )); then
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$tmp/block.tsv"
		touch "$tmp/.done"
		return 0
	fi

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
		lookup=$(pvar_lookup "$pbase" "$tmp/lead.assoc" "$tmp")

		line=$(awk -v snp="$snp" -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp && $2==snp{print $1,$2,$3,$4; exit}' "$lookup")
			if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=ID_BP; else
			n_id=$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{n++} END{print n+0}' "$lookup")
			if [[ $n_id -eq 1 ]]; then
				read -r refbp vid _ref _alt <<< "$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{print $1,$2,$3,$4; exit}' "$lookup")"
				match=ID_only
				[[ $refbp != "$bp" ]] && log "WARN $t chr$c $snp: lead_bp=$bp but pvar_bp=$refbp; using pvar_bp"
			else
				line=$(awk -v bp="$bp" -v ea="$ea" -v oa="$oa" 'BEGIN{FS=OFS="\t"} $1==bp{ref=toupper($3); alt=toupper($4); ea=toupper(ea); oa=toupper(oa); if((ref==ea && alt==oa) || (ref==oa && alt==ea)){print $1,$2,$3,$4; exit}}' "$lookup")
					if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=ALLELE_BP; else
					line=$(awk -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp{id=$2; r=$3; a=$4; n++} END{if(n==1) print bp,id,r,a}' "$lookup")
					if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=UNIQUE_BP; else
						msg="not found by ID+BP, ID-only, allele+BP, or unique BP"
						write_fail "$t" "$cn" "$snp" "$bp" "not_in_1000G_pvar_or_allele_mismatch" "$msg" "$pbase" "NA"
						write_debug "$t" "$cn" "$snp" "$bp" "NA" "NA" "none" "$pbase" "pvar_miss" 0 "$msg"
						continue
					fi
				fi
			fi
		fi

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$refbp" "$cn" "$refbp" "$vid" >> "$tmp/ld.tsv"
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$t" "$cn" "$snp" "$refbp" "$vid" "$pbase" "$bp" "$match" >> "$tmp/lead.map.tsv"
	done

		awk 'BEGIN{FS=OFS="\t"} NR>1{print $6}' "$tmp/lead.map.tsv" | sort -u | while IFS= read -r pbase; do
			[[ -n $pbase ]] || continue
			outpre="$tmp/$(basename "$pbase").batch"; plog="$outpre.plink2.log"; outvcor="$outpre.vcor"
			awk -v pbase="$pbase" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $5}' "$tmp/lead.map.tsv" | sort -u > "$outpre.ld_ids"
			[[ -s $outpre.ld_ids ]] || continue
			log "LD $t chr$c: PLINK start $(basename "$pbase") with $(wc -l < "$outpre.ld_ids") lead IDs"
			cmdtxt="plink2 --pfile $pbase --allow-extra-chr --make-founders --r2-unphased allow-ambiguous-allele --ld-snp-list $outpre.ld_ids --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 $ld_r2 --out $outpre"
			printf "%s\n" "$cmdtxt" > "$outpre.cmd"
			if plink2 --pfile "$pbase" --allow-extra-chr --make-founders --r2-unphased allow-ambiguous-allele --ld-snp-list "$outpre.ld_ids" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --out "$outpre" > "$plog" 2>&1; then
				if [[ -s $outvcor ]]; then
					log "LD $t chr$c: PLINK done $(basename "$pbase"); vcor_rows=$(awk 'NR>1{n++} END{print n+0}' "$outvcor")"
					awk -v pbase="$pbase" '
					BEGIN{FS="[ \t]+";OFS="\t"}
					NR==FNR{if(NR>1 && $6==pbase) map[$5]=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4; next}
					FNR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="ID_A") ia=i; if($i=="POS_B") pb=i; if($i=="ID_B") ib=i; if($i~/R2$/) ir=i}; next}
					ia&&pb&&ib&&ir&&($ia in map){split(map[$ia],m,SUBSEP); print m[1],m[2],m[3],m[4],m[2],$pb,$ib,$ir}
				' "$tmp/lead.map.tsv" "$outvcor" >> "$tmp/ld.tsv"
				awk -v pbase="$pbase" -v dbg="$dirout/lead/ld_debug.tsv" '
					BEGIN{FS="[ \t]+";OFS="\t"}
					NR==FNR{if(NR>1 && $6==pbase) map[$5]=$0; next}
					FNR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="ID_A") ia=i}; next}
					ia&&($ia in map){n[$ia]++}
					END{for(id in map){split(map[id],m,FS); print m[1],m[2],m[3],m[7],m[4],m[5],m[8],m[6],"ok_batch",n[id]+0,"plink_ok_batch" >> dbg}}
				' "$tmp/lead.map.tsv" "$outvcor"
				else
					msg="PLINK finished but $outvcor is missing/empty"
					awk -v pbase="$pbase" -v dbg="$dirout/lead/ld_debug.tsv" -v msg="$msg" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $1,$2,$3,$7,$4,$5,$8,$6,"ok_no_vcor",0,msg >> dbg}' "$tmp/lead.map.tsv"
					log "LD $t chr$c: PLINK done $(basename "$pbase") but vcor is empty"
				fi
			else
			msg=$(clean_msg "$plog")
			awk -v pbase="$pbase" -v fail="$dirout/lead/lead_1000G.fail.tsv" -v dbg="$dirout/lead/ld_debug.tsv" -v msg="$msg" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $1,$2,$3,$7,"plink2_ld_failed",msg,$6,$5 >> fail; print $1,$2,$3,$7,$4,$5,$8,$6,"plink_failed",0,msg >> dbg}' "$tmp/lead.map.tsv"
			log "FAIL $t chr$c $(basename "$pbase") batch: $msg"
		fi
	done

	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	{
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n"
		awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp; if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n
		} > "$tmp/block.tsv"
		touch "$tmp/.done"
		log "LD $t chr$c: done; ld_rows=$(data_rows "$tmp/ld.tsv") blocks=$(data_rows "$tmp/block.tsv")"
	}

trait_ld_complete(){
	local t=$1 d=$dirout/ld/$t c cn nlead
	for c in $chrs; do
		cn=$c; [[ $c == X ]] && cn=23
		nlead=$(awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR>1 && ($2==cn || toupper($2)==c){n++} END{print n+0}' "$dirout/lead/$t.lead.assoc")
		[[ $(data_rows "$d/chr$c/lead.assoc") -eq $nlead ]] || return 1
		valid_chr_ld "$d/chr$c" "$nlead" || return 1
	done
}

# 单 trait LD 合并：合并所有染色体的 ld/block 结果。
# n>1 的 block 是 selected/high-LD loci；n==1 的 block 单独保留。
	merge_trait_ld(){
		local t=$1 d=$dirout/ld/$t x c f first pickf singlef
		pickf=$dirout/lead/$t.pick.tsv
		singlef=$dirout/lead/$t.pick.single.tsv
		for x in ld block; do
			first=1; : > "$d/$x.tsv"
			for c in $chrs; do
				f=$d/chr$c/$x.tsv
			[[ -f $f ]] || continue
			if [[ $first -eq 1 ]]; then cat "$f"; first=0; else awk 'NR>1' "$f"; fi
			done > "$d/$x.tsv"
			awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"
		done
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$pickf"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$singlef"
		awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$pickf"
		awk 'BEGIN{FS=OFS="\t"} NR>1 && $7==1' "$d/block.tsv" >> "$singlef"
		log "========== DONE $t: $(awk 'NR>1{n++} END{print n+0}' "$d/block.tsv") blocks; pick=$(data_rows "$pickf"); single=$(data_rows "$singlef") =========="
	}

# 单 trait 调度：并行运行各染色体 LD，完成后检查完整性并合并。
	run_trait_ld(){
		local t=$1 c d nlead ld_status=0
		log "========== DO $t =========="
		if trait_ld_complete "$t"; then
			log "LD $t: all chromosome results complete; merge existing results"
			merge_trait_ld "$t"
			return 0
		fi

		for c in $chrs; do
			run_chr "$t" "$c" &
			while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do
				wait -n || ld_status=1
			done
		done
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
			wait -n || ld_status=1
		done
		if (( ld_status != 0 )); then
			log "ERROR one or more LD chromosome jobs failed for $t"
			return 1
		fi

		d=$dirout/ld/$t
		for c in $chrs; do
			nlead=$(data_rows "$d/chr$c/lead.assoc")
			if ! valid_chr_ld "$d/chr$c" "$nlead"; then
				log "ERROR incomplete LD output for $t chr$c; remove $d/chr$c and rerun"
				return 1
			fi
		done
		merge_trait_ld "$t"
	}

	trait_status=0
	for t in $traits; do
		run_trait_ld "$t" &
		while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do
			wait -n || trait_status=1
		done
	done
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
		wait -n || trait_status=1
	done
	if (( trait_status != 0 )); then
		log "ERROR one or more LD trait jobs failed"
		exit 1
	fi
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
	for t in $traits; do
		[[ -s "$dirout/lead/$t.pick.tsv" ]] && awk 'NR>1' "$dirout/lead/$t.pick.tsv" >> "$dirout/lead/pick.tsv"
		[[ -s "$dirout/lead/$t.pick.single.tsv" ]] && awk 'NR>1' "$dirout/lead/$t.pick.single.tsv" >> "$dirout/lead/pick.single.tsv"
	done
	log_trait_count "selected loci" "$dirout/lead/pick.tsv"
	log_trait_count "single-SNP loci" "$dirout/lead/pick.single.tsv"

	log "s1 done. Check: $dirout/lead/ld_debug.tsv and $dirout/lead/lead_1000G.fail.tsv"
elif [[ "$start_step" == "s2" ]]; then
	for f in "$dirout/lead/pick.tsv" "$dirout/lead/lead.assoc"; do
		[[ -s "$f" ]] || { log "ERROR START_STEP=s2 requires existing file: $f"; exit 1; }
	done
	[[ $(data_rows "$dirout/lead/pick.tsv") -gt 0 ]] || { log "ERROR START_STEP=s2 requires non-empty $dirout/lead/pick.tsv; rerun START_STEP=s1"; exit 1; }
	log "Skip s1; reuse existing $dirout/lead/pick.tsv"
else
	log "ERROR unknown START_STEP=$start_step; use s1 or s2"
	exit 1
fi


# ============================================================
# s2：为 selected loci 构建现代人/古人基因型矩阵。
# 从 lead/pick.tsv 读取高 LD loci，抽取 1000G EUR core-region VCF。
# 再把 archaic VCF 投影到同一套 1000G REF/ALT 模板，写入 mat/<trait>/<locus>/*.tsv。
# 已完整生成的 kg.tsv 和 archaic TSV 会自动跳过。
# ============================================================

expected_arch=$(find "$arch0" -mindepth 1 -maxdepth 1 -type d | wc -l)
npick=$(data_rows "$dirout/lead/pick.tsv")
(( npick > 0 )) || { log "ERROR no high-LD blocks in $dirout/lead/pick.tsv; inspect $dirout/lead/lead_1000G.fail.tsv and rerun START_STEP=s1"; exit 1; }
log_trait_count "selected loci" "$dirout/lead/pick.tsv"
log "s2 start: $npick high-LD blocks from $dirout/lead/pick.tsv"

# 单 locus 矩阵：抽取 1000G VCF，逐个 archaic 样本投影 REF/ALT。
# 最终输出 kg.tsv 以及每个 archaic 样本的 TSV。
run_locus_mat(){
	local tr=$1 chr=$2 snp=$3 bp=$4 st=$5 en=$6 _n=$7 _size_bp=$8
	local cl=$chr lo=$st hi=$en d o kg par ad name outname avcf f arc_raw arc_vcf arc_log nvar dot
	[[ $bp -lt $lo ]] && lo=$bp
	[[ $bp -gt $hi ]] && hi=$bp
	d=$dirout/coreVcf/$tr/${chr}.${snp}.${bp}
	o=$dirout/mat/$tr/${chr}.${snp}.${bp}
	if valid_mat_locus "$o" "$expected_arch"; then
		log "s2 skip existing mat: $tr ${chr}.${snp}.${bp}"
		return 0
	fi
	log "s2 mat start: $tr ${chr}.${snp}.${bp} region=$chr:$lo-$hi"
	rm -rf "$d" "$o"
	mkdir -p "$d" "$o"

	if [[ $chr == 23 || $chr == X ]]; then
		cl=X
		kg=$dirmod/chrX.vcf.gz
		[[ -s $kg ]] || { log "ERROR missing chrX 1000G VCF for $tr $snp: $kg"; return 1; }
		bcftools view -r "X:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; return 1; }
	else
		kg=$dirmod/chr$cl.vcf.gz
		[[ -s $kg ]] || kg=$dirmod/chr${cl}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		[[ -s $kg ]] || { log "ERROR missing 1000G VCF for $tr chr$chr $snp"; return 1; }
		bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; return 1; }
	fi
	tabix -f -p vcf "$d/kg.vcf.gz" || { log "ERROR tabix failed: $d/kg.vcf.gz"; return 1; }
	bcftools query -l "$d/kg.vcf.gz" > "$o/kg.samples.tsv" || { log "ERROR bcftools query samples failed: $d/kg.vcf.gz"; return 1; }

	for ad in "$arch0"/*; do
		[[ -d $ad ]] || continue
		name=$(basename "$ad")
		outname=$(echo "$name" | tr '[:upper:]' '[:lower:]')
		avcf=""
		for f in "$ad"/*chr"$cl"_*.vcf.gz "$ad"/*chr"$cl".*.vcf.gz; do
			[[ -s $f ]] || continue
			avcf=$f
			break
		done
		[[ -n $avcf ]] || { log "WARN missing archaic VCF: $outname chr$cl"; continue; }

		arc_raw="$d/$outname.raw.vcf.gz"
		arc_vcf="$d/$outname.vcf.gz"
		arc_log="$d/$outname.project.log"

		# Extract raw archaic records first, then project them to the 1000G REF/ALT template.
		bcftools view -r "$cl:$lo-$hi" -Oz -o "$arc_raw" "$avcf" || { log "ERROR bcftools failed: archaic raw $outname $cl:$lo-$hi"; return 1; }
		tabix -f -p vcf "$arc_raw" || { log "ERROR tabix failed: $arc_raw"; return 1; }

		Rscript "$dirscript/locus.R" prep_archaic "$d/kg.vcf.gz" "$arc_raw" "$arc_vcf" "$outname" > "$arc_log" 2>&1 || {
			log "ERROR archaic projection failed: $outname $cl:$lo-$hi; see $arc_log"
			tail -20 "$arc_log"
			return 1
		}

		nvar=$(bcftools view -H "$arc_vcf" | wc -l)
		dot=$(bcftools view -H "$arc_vcf" | awk 'BEGIN{FS="\t"} $5=="."{n++} END{print n+0}')
		[[ $nvar -gt 0 && $dot -eq 0 ]] || { log "ERROR bad projected archaic VCF: $arc_vcf n=$nvar ALT_dot=$dot"; return 1; }
		rm -f "$arc_raw" "$arc_raw.tbi"
	done

	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv" || { log "ERROR bcftools query failed: $d/kg.vcf.gz"; return 1; }
	for f in "$d"/*.vcf.gz; do
		name=$(basename "$f" .vcf.gz)
		[[ $name == kg || $name == *.raw ]] && continue
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$f" | awk 'BEGIN{FS=OFS="\t"} $4!="." && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ {if($5=="0" || $5=="1") $5=$5"/"$5; print}' > "$o/$name.tsv" || { log "ERROR bcftools query failed: $f"; return 1; }
	done
	valid_mat_locus "$o" "$expected_arch" || { log "ERROR incomplete mat output: $o"; return 1; }
	rm -rf "$d"
	log "s2 mat done: $tr ${chr}.${snp}.${bp}; kg_rows=$(data_rows "$o/kg.tsv") archaic_files=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' | wc -l)"
}

# 单 trait 矩阵调度：对该 trait 的所有 selected loci 并行生成 mat 目录。
run_trait_mat(){
	local t=$1 tr chr snp bp st en _n _size_bp s2_status=0 n_trait
	n_trait=$(awk -v t="$t" 'BEGIN{FS="\t"} NR>1 && $1==t{n++} END{print n+0}' "$dirout/lead/pick.tsv")
	log "s2 trait start: $t ($n_trait loci)"
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		run_locus_mat "$tr" "$chr" "$snp" "$bp" "$st" "$en" "$_n" "$_size_bp" &
		while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do
			wait -n || s2_status=1
		done
	done < <(awk -v t="$t" 'BEGIN{FS=OFS="	"} NR>1 && $1==t && !seen[$1 FS $2 FS $3 FS $4 FS $5 FS $6]++{print}' "$dirout/lead/pick.tsv")
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
		wait -n || s2_status=1
	done
	if (( s2_status != 0 )); then
		log "ERROR one or more s2 locus jobs failed for $t"
		return 1
	fi
	log "s2 trait done: $t"
}

s2_trait_status=0
for t in $traits; do
	run_trait_mat "$t" &
	while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do
		wait -n || s2_trait_status=1
	done
done
while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
	wait -n || s2_trait_status=1
done
if (( s2_trait_status != 0 )); then
	log "ERROR one or more s2 trait jobs failed"
	exit 1
fi
rm -rf "$dirout/coreVcf"

# ============================================================
# s3：识别 inherited loci，构建风险单倍型，并生成系统发育树。
# locus.R make_hap 会比较风险 haplotype 与 archaic allele，输出 inherited_segments；
# locus.R make_phy 生成 PHYLIP 文件，phyml 建树，locus.R make_tree 绘制树图。
# ============================================================

# 单 trait inherited haplotype：输出 hap/site/core/region 表，并在 merge 后汇总 inherited_segments。
run_trait_hap(){
	local t=$1
	log "s3 hap start: $t"
	Rscript "$dirscript/locus.R" make_hap "$dirout" "$t" || return 1
	log "s3 hap done: $t"
}

s3_hap_status=0
for t in $traits; do
	run_trait_hap "$t" &
	while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do
		wait -n || s3_hap_status=1
	done
done
while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
	wait -n || s3_hap_status=1
done
if (( s3_hap_status != 0 )); then
	log "ERROR one or more s3 hap jobs failed"
	exit 1
fi
Rscript "$dirscript/locus.R" make_hap "$dirout" merge || { log "ERROR locus.R make_hap merge failed"; exit 1; }
inherited_tsv="$dirout/report/inherited_segments.tsv"
hap_match_tsv="$dirout/report/hap_match.tsv"
log_trait_count "inherited_segments.tsv inherited loci" "$inherited_tsv"
log "summary: loci confidence reference (file path $inherited_tsv; columns p_ils,best_lineage,matched_archaics)"
log "s3 population filter start: pop=$filter_pop max_count=$filter_max_count"
Rscript "$dirscript/locus.R" filter_hap "$dirout" "$sample_file" "$filter_pop" "$filter_max_count" || { log "ERROR locus.R filter_hap failed"; exit 1; }
filtered_hap_tsv="$dirout/report/filtered_hap_match.tsv"
filtered_inherited_tsv="$dirout/report/filtered_inherited_segments.tsv"
before_inherited_loci=$(data_rows "$inherited_tsv")
after_inherited_loci=$(data_rows "$filtered_inherited_tsv")
before_inherited_hap=$(data_rows "$hap_match_tsv")
after_inherited_hap=$(data_rows "$filtered_hap_tsv")
log_trait_count "filtered inherited loci" "$filtered_inherited_tsv"
log "FILTER SUMMARY: inherited loci before_filter=$before_inherited_loci after_filter=$after_inherited_loci criteria=${filter_pop}<=$filter_max_count"
log "FILTER SUMMARY: inherited haplotypes before_filter=$before_inherited_hap after_filter=$after_inherited_hap criteria=${filter_pop}<=$filter_max_count"
log "FILTER SUMMARY: before-filter locus counts are in $dirout/report/all.xlsx sheet inherited_loci_counts"
log "FILTER SUMMARY: before-filter haplotype counts are in $dirout/report/all.xlsx sheet inherited_haplotype_counts"
log "FILTER SUMMARY: after-filter loci/haplotypes are in $dirout/report/filtered.xlsx sheets filtered_loci and filtered_haplotypes"

# 单 trait PHYLIP：基于 inherited haplotype/archaic 匹配结果生成 phyml 输入文件。
run_trait_make_phy(){
	local t=$1
	log "s3 make_phy start: $t"
	Rscript "$dirscript/locus.R" make_phy "$dirout" "$t" || return 1
	log "s3 make_phy done: $t"
}

s3_phy_status=0
for t in $traits; do
	run_trait_make_phy "$t" &
	while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do
		wait -n || s3_phy_status=1
	done
done
while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
	wait -n || s3_phy_status=1
done
if (( s3_phy_status != 0 )); then
	log "ERROR one or more s3 make_phy jobs failed"
	exit 1
fi
log "summary: *.main.phy filtered haplotypes (file path $dirout/phy) $(count_files "$dirout/phy" '*.main.phy')"
log "summary: filtering reference (workbook $dirout/report/filtered.xlsx; criteria $filter_pop <= $filter_max_count)"
clean_report_workbooks_only

# 单 trait 建树：只对 *.main.phy 运行 phyml；已有 tree 文件时跳过。
run_trait_phyml(){
	local t=$1 f n status=0
	n=$(find "$dirout/phy/$t" -name '*.main.phy' 2>/dev/null | wc -l)
	log "s3 phyml trait start: $t ($n main.phy files; tree source=*.main.phy only)"
	while IFS= read -r f; do
		(
			if [[ -s "${f}_phyml_tree.txt" ]]; then
				log "s3 skip existing phyml tree: $f"
				exit 0
			fi
			log "s3 phyml start: $f"
			phyml -i "$f" -m HKY85 -c 4 -a e -v e -b 100 < /dev/null || exit 1
			log "s3 phyml done: $f"
		) &
		while [[ $(jobs -rp | wc -l) -ge $job_phyml ]]; do
			wait -n || status=1
		done
	done < <(find "$dirout/phy/$t" -name '*.main.phy' 2>/dev/null | sort)
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
		wait -n || status=1
	done
	if (( status != 0 )); then
		log "ERROR one or more phyml jobs failed for $t"
		return 1
	fi
	log "s3 phyml trait done: $t"
}

s3_phyml_status=0
for t in $traits; do
	run_trait_phyml "$t" &
	while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do
		wait -n || s3_phyml_status=1
	done
done
while [[ $(jobs -rp | wc -l) -gt 0 ]]; do
	wait -n || s3_phyml_status=1
done
if (( s3_phyml_status != 0 )); then
	log "ERROR one or more s3 phyml trait jobs failed"
	exit 1
fi
miss_tree=$(find "$dirout/phy" -name '*.main.phy' | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
[[ -z $miss_tree ]] || { log "ERROR missing phyml tree for: $miss_tree"; exit 1; }
Rscript "$dirscript/locus.R" make_tree "$dirout" || { log "ERROR locus.R make_tree failed"; exit 1; }
n_tree_png=$(find "$dirout/plot" -maxdepth 1 -name 's8_tree_main_*.png' | wc -l)
log "summary: loci phylogeny tree PNG (source *.main.phy; file path $dirout/plot) $n_tree_png"
[[ $n_tree_png -gt 0 ]] || { log "ERROR no tree PNG generated in $dirout/plot"; exit 1; }

# 最终 summary：运行结束后集中打印关键结果和文件路径，便于检查日志。
log_trait_count "selected loci" "$dirout/lead/pick.tsv"
log_trait_count "single-SNP loci" "$dirout/lead/pick.single.tsv"
log "summary: inherited loci before filter $before_inherited_loci"
log "summary: inherited loci after filter $after_inherited_loci"
log "summary: inherited haplotypes before filter $before_inherited_hap"
log "summary: inherited haplotypes after filter $after_inherited_hap"
log "summary: *.main.phy filtered haplotypes (file path $dirout/phy) $(count_files "$dirout/phy" '*.main.phy')"
log "summary: loci confidence reference (workbook $dirout/report/all.xlsx; sheet inherited_segments; columns p_ils,best_lineage,matched_archaics)"
log "summary: filtering reference (workbook $dirout/report/filtered.xlsx; criteria $filter_pop <= $filter_max_count)"
log "summary: loci phylogeny tree PNG (source *.main.phy; file path $dirout/plot) $n_tree_png"
clean_report_workbooks_only
log "summary: report workbooks $dirout/report/all.xlsx $dirout/report/filtered.xlsx $dirout/report/selected.xlsx"
log "ALL DONE"
