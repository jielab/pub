plink --vcf trio.vcf.gz --chr 1-22 --make-king-table --out trio.king
plink --vcf trio.vcf.gz --chr 1-22  --snps-only just-acgt  --max-alleles 2 --split-par hg19 --double-id --make-bed  --out trio
plink0 --bfile trio --chr 15 --homozyg --homozyg-snp 50 --homozyg-kb 1000 --homozyg-density 50 --homozyg-gap 1000 --out trio.roh
plink0 --bfile trio --genome full --out trio.ibd

========


module load python/anaconda3/2020.7
source activate ensembl-vep
dir=/work/sph-huangj
anodir=$dir/data/vep
grch=38
vcf=$dir/data/wgs/jiang/clean/clean.vcf.gz
samples=`bcftools query -l $vcf | tr '\n' ' '`; samples2=`echo $samples | sed 's/ /,/g'`

## merge AF and Annotations ##
bcftools annotate $vcf -a $anodir/topmed/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz -c Topmed_AF:=AF --mark-sites In_Topmed -Oz -o tmp1.vcf.gz; tabix tmp1.vcf.gz
bcftools filter tmp1.vcf.gz --exclude '(Topmed_AF >=0.01 && Topmed_AF <=0.99)' -Oz -o tmp2.vcf.gz; tabix tmp2.vcf.gz
	bcftools query -f '%CHROM\n' tmp2.vcf.gz | sort | uniq -c
# echo '##INFO=<ID=SFARI_gene,Number=1,Type=String,Description="disease region">' > SFARI.gene.h
# echo '##INFO=<ID=SFARI_cnv,Number=1,Type=String,Description="disease region">' > SFARI.cnv.h
# bcftools annotate tmp3.vcf.gz -a $adir/SFARI_genes.b$grch.bed.gz -c CHROM,FROM,TO,SFARI_gene -h SFARI.gene.h --mark-sites In_SFARI_genes -Oz -o tmp4.vcf.gz; tabix tmp4.vcf.gz
# bcftools annotate tmp4.vcf.gz -a $adir/SFARI_cnv.b$grch.bed.gz -c CHROM,FROM,TO,SFARI_cnv -h SFARI.cnv.h --mark-sites In_SFARI_cnvs -Oz -o tmp5.vcf.gz; tabix tmp5.vcf.gz
bcftools annotate tmp2.vcf.gz -a $anodir/clinvar/clinvar.b38.vcf.gz -c CLNDISDB,CLNDN,CLNREVSTAT,CLNSIG,MC,RS --mark-sites In_CLINVAR -Oz -o tmp3.vcf.gz
	bcftools query tmp3.vcf.gz -f '%CLNSIG\t%CLNREVSTAT\n' | sort | uniq -c
bcftools filter tmp3.vcf.gz --include '(CLNSIG ~ "Pathogenic")' -Ov -o final.vcf # && CLNREVSTAT ~ "expert_panel"

bcftools +split-vep tmp7.vcf.gz -c Consequence,IMPACT,SYMBOL,CANONICAL,SIFT -s worst | bcftools view -Oz -o all.final.vcf.gz; tabix all.final.vcf.gz


## add VEP prediction ##
vep -i tmp3.vcf --vcf -o jie.vcf --offline --cache --dir $dir/data/vep --cache_version 112 --species homo_sapiens --format vcf --vcf --force_overwrite --overlaps --sift b --hgvs --symbol --canonical --pick --stats_text high_conf_vep.stats --assembly GRCh$grc
vep -i tmp.vcf --vcf -o vep_out.txt --everything --protein --symbol --force_overwrite
bcftools +split-vep tmp7.vcf.gz -l
bcftools +split-vep tmp7.vcf.gz -c Consequence,IMPACT,SYMBOL,CANONICAL,SIFT -s worst | bcftools view -Oz -o all.final.vcf.gz; tabix all.final.vcf.gz
bcftools query all.final.vcf.gz -f '%IMPACT %Consequence\n' | sort | uniq -c
# !! view VEP summary !!


## For each proband, find top SNV ##
bcftools filter all.final.vcf.gz --include '(In_CLINVAR=1 && CLNSIG !~ "onflict" && CLNSIG ~ "athogenic") || (IMPACT="HIGH")' -Oz -o all.top.vcf.gz # && (In_SFARI_genes=1 || In_SFARI_cnvs=1)
bcftools view all.top.vcf.gz --samples $samples2 --no-update | \
bcftools filter --include 'GT[0] !~"\."' | \
bcftools query -H -f '%ID\t%CHROM:%POS:%REF:%ALT\t%DP\t%Topmed_AF\t%SFARI_gene\t%SFARI_cnv\t%CLNDISDB\t%CLNDN\t%CLNREVSTAT\t%CLNSIG\t%MC\t%Consequence\t%IMPACT\t%SYMBOL\t%CANONICAL\t%SIFT[\t%SAMPLE=%GT]\t%CSQ\n' > final.txt


## For each proband, find top CNV ##
#vep --plugin StructuralVariantOverlap,file=gnomad_v2_sv.sites.vcf.gz
for sample in $samples; do
#	cp $dir/$project/cnv/cnvnator/$sample.cnv.bed ./
#	awk 'BEGIN{OFS="\t"}{if ($1 !~/^@/ && $1 ~/^chr/) print $1,$2,$3,"GATK"}' $dir/$project/cnv/modelFinal/$sample.modelFinal.seg > $sample.cnv.bed
#	bcftools query $dir/$project/gcnv/$sample/genotyped-segments-*-$sample.vcf.gz -f '%CHROM %POS %INFO/END [%CN]\n' | awk '{if ($NF<=1) $NF="deletion"; else if ($NF>2) $NF="duplication"; if ($NF !=2) print $0}' | sed 's/ /\t/g' > $sample.cnv.bed
	awk -F "\t" 'NR>1 && $24 !="breakpoint" {print "chr"$1,$2,$3,$21}' $dir/$project/CNV/Annotation/$sample.final.bam_CNVs.final.hg38_multianno.xls | sed 's/ /\t/g' > $sample.cnv.bed
	bedtools intersect -a $sample.cnv.bed -b $adir/SFARI_cnv.b$grch.bed -wo > $sample.cnv.tmp1
	awk '{len1=$3-$2; len2=$7-$6; lap=$9; r1=len1/lap; r2=len2/lap; print FILENAME"\t"$0"\t|", len1, len2, r1,r2}' $sample.cnv.tmp1 | sed 's/.cnv.tmp1//' > $sample.cnv.tmp2
done
cat *.cnv.tmp2 | awk '{print $9,$1}' | sort | uniq | awk '{array[$1]=array[$1]","$2} END{for(key in array) print key,array[key]}' > vip.cnv.tmp
python $0join_file -i "vip.cnv.tmp,SPACE,0 $adir/SFARI_cnv.b$grch.bed,TAB,3" -o vip.cnv
sort -k 7nr vip.cnv > vip.cnv.txt


### check inheritance patterns (de novo vs. heritable) ##
module load genmod/3.7.2
module load gcc/7.2.0
module load denovogear/2018-05-1_6723027
dng dnm auto --ped ../vcf/sample.ped --bcf step1.final.bcf > step2.dng.out
awk 'BEGIN{print "#CHROM POS REF ALT DNG"}{if ($5 >=1 && $5<=22) print $5,$7,$9,$11,"Y"}' step2.dng.out | \
	sort -k 1,1n -k 2,2n | sed 's/ /\t/g' | bgzip > step2.dng.txt.gz
tabix -f -s 1 -b 2 -e 2 -S 1 step2.dng.txt.gz
echo '##INFO=<ID=DNG,Number=1,Type=String,Description="denovogear analysis results">' > dng.h
bcftools annotate -a step2.dng.txt.gz -c 'CHROM,POS,REF,ALT,DNG' -h dng.h -Oz -o step2a.vcf.gz step1.final.bcf
# $dir/software/from_SCC/triodenovo.0.05/bin/triodenovo --ped ../vcf/sample.ped --in_vcf step2a.vcf.gz --out_vcf step2b.tmp.vcf
genmod annotate step2a.vcf.gz -r ?? --annotate_regions | genmod models - --family_file ../vcf/sample.ped | bcftools view -Ob -o step2b.bcf
tabix step2b.bcf
bcftools filter --include 'DNG ="Y" || GeneticModels ~ "_dn"' -Ob -o step2.final.bcf step2b.bcf
bcftools query -f '%CSQ\n' step2.final.bcf | awk -F "|" '{print $2}' | sort | uniq -c 
