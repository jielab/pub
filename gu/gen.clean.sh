GRCH=37


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载数据，plink2格式
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
baseurl=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
wget -q -O index.html "$baseurl/"
grep -oP 'href="\K[^"]+' index.html | grep -v '^?' | grep -v '^/' | grep -v '../' | grep -v '/$' > files.txt
cat files.txt | while read f; do; wget -c --inet4-only --timeout=20 --tries=5 --waitretry=5 "$baseurl/$f" done

sample_file=samples_v3.20130502.ALL.panel
awk 'BEGIN{OFS="\t"} NR>1 && $3=="EUR"{print $1,$1}' $sample_file > EUR.sample.2id
awk 'BEGIN{OFS="\t"; print "#FID","IID","SEX"} NR>1 && $3=="EUR"{s=($4=="male"||$4=="M")?1:($4=="female"||$4=="F")?2:0; print $1,$1,s}' $sample_file > EUR.sex.3col
awk 'NR>1 && $3=="EUR" && tolower($4)=="male"{print $1}' $sample_file > EUR.male.sample.1id
awk 'NR>1 && $3=="EUR" && tolower($4)=="male"{print $1,$1,1}' $sample_file > EUR.male.3col

for chr in {1..22} X Y; do
	echo chr$chr
	f=$(ls *chr${chr}.*.vcf.gz | head -1)
	[ "$f" = "ALL.chr$chr.vcf.gz" ] || mv "$f" ALL.chr$chr.vcf.gz
	[ -f "$f.tbi" ] && mv "$f.tbi" ALL.chr$chr.vcf.gz.tbi; [ -f "$f.csi" ] && mv "$f.csi" ALL.chr$chr.vcf.gz.csi
	tabix -f ALL.chr$chr.vcf.gz # 🏮
	if [[ $chr == X ]]; then extra="--update-sex EUR.sex.3col --split-par b$GRCH"; elif [[ $chr == Y ]]; then extra="--update-sex EUR.sex.3col"; else extra=""; fi
	plink2 --vcf ALL.chr$chr.vcf.gz --double-id --allow-extra-chr --keep EUR.sample.2id $extra --make-pgen --out EUR.chr$chr
done


vcf=ALL.chrX.vcf.gz
if [[ $GRCH == 37 ]]; then
	par="X:60001-2699520,X:154931044-155260560"
	nonPar="X:1-60000,X:2699521-154931043,X:155260561-155270560"
	split="b37"
elif [[ $GRCH == 38 ]]; then
	par="X:10001-2781479,X:155701383-156030895"
	nonPar="X:1-10000,X:2781480-155701382,X:156030896-156040895"
	split="b38"
fi

bcftools view -S EUR.male.sample.1id -r $par -m2 -M2 -v snps -Oz -o EUR.male.chrX.par.vcf.gz "$vcf"
bcftools view -S EUR.male.sample.1id -r $nonPar -m2 -M2 -v snps -Oz -o EUR.male.chrX.nonPar.vcf.gz "$vcf"
tabix -f -p vcf EUR.male.chrX.par.vcf.gz
tabix -f -p vcf EUR.male.chrX.nonPar.vcf.gz

plink2 --vcf EUR.male.chrX.par.vcf.gz --double-id --allow-extra-chr --update-sex EUR.male.3col --split-par $split --make-pgen --out EUR.male.chrX.par
plink2 --vcf EUR.male.chrX.nonPar.vcf.gz --double-id --allow-extra-chr --update-sex EUR.male.3col --make-pgen --out EUR.male.chrX.nonPar


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载古基因数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p Vindija Altai Chagyr Denisova Denisova25 

for c in {1..22}; do
	wget -c -P Vindija	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Altai	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Chagyr	"https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz"{,.tbi}
	wget -c -P Denisova		"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Denisova25	"https://cdna.eva.mpg.de/denisova/Den25/VCF/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz"{,.tbi}
done