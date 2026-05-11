dir0=/mnt/d
source $dir0/scripts/f/0phe.f.sh

GRCH=37

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载1000G数据，plink2格式化
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
baseurl=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502 # ⚠️这个版本没有rsID🕳
baseurl=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3

wget -q -O index.html "$baseurl/"
grep -oP 'href="\K[^"]+' index.html | grep -v '^?' | grep -v '^/' | grep -v '../' | grep -v '/$' > files.txt
cat files.txt | while read f; do wget -c --inet4-only --timeout=20 --tries=5 --waitretry=5 "$baseurl/$f"; done

sample_file=$dir0/files/1kg.v3.sample.txt; # cp samples_v3.20130502.ALL.panel $sample_file
for race in EUR AFR EAS SAS AMR; do
	awk -v race=$race 'NR>1 && $3==race {print $1}' $sample_file > $race.1id
	awk -v race=$race 'NR>1 && $3==race {print $1,$1}' $sample_file > $race.2id
	awk -v race=$race 'NR>1 && $3==race {s=(tolower($4)=="male"||$4=="m")?1:(tolower($4)=="female"||$4=="f")?2:0; print $1,$1,s}' $sample_file > $race.2id_sex
	awk -v race=$race 'NR>1 && $3==race && tolower($4)=="male" {print $1}' $sample_file > $race.male.1id
	awk -v race=$race 'NR>1 && $3==race && tolower($4)=="male" {print $1,$1}' $sample_file > $race.male.2id
	awk -v race=$race 'NR>1 && $3==race && tolower($4)=="male" {print $1,$1,1}' $sample_file > $race.male.2id_sex
done


for chr in {1..22} X Y; do
	echo "DO chr$chr"
	f=$(find . -maxdepth 1 -name "*chr${chr}.*.vcf.gz" -printf "%f\n" | sort | head -1)
	[[ -n $f ]] || { echo "ERROR: missing chr$chr VCF"; exit 1; }
	if [[ $f != chr$chr.vcf.gz ]]; then
		mv "$f" chr$chr.vcf.gz
		[[ -f $f.tbi ]] && mv "$f.tbi" chr$chr.vcf.gz.tbi
		[[ -f $f.csi ]] && mv "$f.csi" chr$chr.vcf.gz.csi
	fi
	tabix -f -p vcf chr$chr.vcf.gz
	if [[ $chr == X ]]; then
		plink2 --vcf chrX.vcf.gz --double-id --allow-extra-chr --keep EUR.2id --update-sex EUR.2id_sex --split-par $split --make-pgen --out EUR.chrX
	elif [[ $chr == Y ]]; then
		plink2 --vcf chrY.vcf.gz --double-id --allow-extra-chr --keep EUR.male.2id --update-sex EUR.male.2id_sex --make-pgen --out EUR.chrY
	else
		plink2 --vcf chr$chr.vcf.gz --double-id --allow-extra-chr --keep EUR.2id --make-pgen --out EUR.chr$chr
	fi
done

Xvcf=chrX.vcf.gz
	split=b${GRCH}; par_var=X_PAR_${split}; nonpar_var=X_NONPAR_${split}
	bcftools view -S EUR.male.1id -r ${!par_var} -m2 -M2 -v snps -Oz -o EUR.male.chrX.par.vcf.gz "$Xvcf"; tabix -f EUR.male.chrX.par.vcf.gz
	bcftools view -S EUR.male.1id -r ${!nonpar_var} -m2 -M2 -v snps -Oz -o EUR.male.chrX.nonPar.vcf.gz "$Xvcf"; tabix -f -p vcf EUR.male.chrX.nonPar.vcf.gz
	plink2 --vcf EUR.male.chrX.par.vcf.gz --double-id --allow-extra-chr --update-sex EUR.male.2id_sex --split-par $split --make-pgen --out EUR.male.chrX.par
	plink2 --vcf EUR.male.chrX.nonPar.vcf.gz --double-id --allow-extra-chr --update-sex EUR.male.2id_sex --make-pgen --out EUR.male.chrX.nonPar


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载古基因数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p Vindija Altai Chagyr Denisova Denisova25 

for c in {1..22} X; do
	wget -c -P Vindija	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Altai	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Chagyr	"https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz"{,.tbi}
	wget -c -P Denisova		"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Denisova25	"https://cdna.eva.mpg.de/denisova/Den25/VCF/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz"{,.tbi}
done

# 对 Chagyr, 重新 gzip 
ls -1 *.vcf.gz | parallel -j 4 'echo {}; mv {} {}.old.gz; gunzip -c {}.old.gz | bgzip -@ 4 -c > {} && tabix -f -p vcf {}'

ls -1 */*.vcf.gz | xargs -P 8 -I {} sh -c 'echo "{}"; tabix -f -p vcf "{}"'