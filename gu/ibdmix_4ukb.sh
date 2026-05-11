#!/bin/bash

set -euo pipefail

dir0=/scratch/2026-05-10/sph-huangj
dirmod=$dir0/imp
dirarc=/data/sph-huangj/refGen/gu
dirsoft=/work/sph-huangj/software/IBDmix
dirscript=/work/sph-huangj/scripts/gu
ibdmix_sh=$dirscript/ibdmix.sh
sample_info=$dir0/files/ukb.sample_info.txt
base_out=$dir0/gu/ibdmix

batch_n=1000
queue=short
job_of_chr=7
job_in_chr=5
nthread=40
run_chrs=({1..22} X)
batch_glob=${BATCH_GLOB:-psam.*.1id}
if (( job_of_chr * job_in_chr > nthread )); then
	echo "ERROR: job_of_chr * job_in_chr exceeds nthread: $job_of_chr * $job_in_chr > $nthread" >&2
	exit 1
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Split UKB samples into 1000-sample batches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p "$base_out/sample_keep"
all_id=$base_out/sample_keep/ukb.all.1id
[[ -s "$dirmod/chr1.psam" ]] || { echo "ERROR missing psam: $dirmod/chr1.psam" >&2; exit 1; }
[[ -s "$ibdmix_sh" ]] || { echo "ERROR missing template ibdmix.sh: $ibdmix_sh" >&2; exit 1; }
[[ -s "$sample_info" ]] || echo "WARNING missing sample_info, generated reports will skip by-group outputs: $sample_info" >&2
command -v bsub >/dev/null || { echo "ERROR missing bsub command" >&2; exit 1; }
awk 'NR>1 {print $2}' "$dirmod/chr1.psam" > "$all_id"
[[ -s "$all_id" ]] || { echo "ERROR no sample IDs written: $all_id" >&2; exit 1; }

rm -f "$base_out"/sample_keep/psam.*.1id
split -d -a 3 -l "$batch_n" --numeric-suffixes=1 --additional-suffix=.1id "$all_id" "$base_out/sample_keep/psam."

nbatch=$(find "$base_out/sample_keep" -maxdepth 1 -name 'psam.*.1id' | wc -l)
echo "Total samples: $(wc -l < "$all_id")"
echo "Batch size:    $batch_n"
echo "Total batches: $nbatch"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Generate modified ibdmix.sh per batch, then bsub
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mapfile -t sample_files < <(find "$base_out/sample_keep" -maxdepth 1 -name "$batch_glob" | sort)
(( ${#sample_files[@]} > 0 )) || { echo "ERROR no batch sample files matched: $base_out/sample_keep/$batch_glob" >&2; exit 1; }

for sample_file in "${sample_files[@]:0:10}"; do
	grp=$(basename "$sample_file" | sed -E 's/^psam\.([0-9]{3})\.1id$/\1/')
	[[ "$grp" =~ ^[0-9]{3}$ ]] || { echo "ERROR bad batch filename: $sample_file" >&2; exit 1; }
	dirout=$base_out/$grp
	if [[ -d "$dirout" ]]; then
		echo "Skip existing output directory: $dirout"
		continue
	fi
	mkdir -p "$dirout"
	cp -f "$sample_file" "$dirout/psam.$grp.1id"

	awk \
		-v dir0="$dir0" \
		-v dirmod="$dirmod" \
		-v dirarc="$dirarc" \
		-v dirsoft="$dirsoft" \
		-v dirscript="$dirscript" \
		-v sample_label="$grp" \
		-v sample_keep="$dirout/psam.$grp.1id" \
		-v sample_info="$sample_info" \
		-v dirout="$dirout" \
		-v chrs="${run_chrs[*]}" \
		-v job_of_chr="$job_of_chr" \
		-v job_in_chr="$job_in_chr" '
		NR == 1 {print; print "module load python/anaconda3/2022.7; source activate; conda activate R4.5.2"; next}
		/^dir0=/         {print "dir0=" dir0; next}
		/^dirmod=/       {print "dirmod=" dirmod; next}
		/^dirarc=/       {print "dirarc=" dirarc; next}
		/^dirsoft=/      {print "dirsoft=" dirsoft; next}
		/^dirscript=/    {print "dirscript=" dirscript; next}
		/^sample_label=/ {print "sample_label=" sample_label; next}
		/^sample_keep=/  {print "sample_keep=" sample_keep; next}
		/^sample_info=/  {print "sample_info=" sample_info; next}
		/^dirout=/       {print "dirout=" dirout; next}
		/^chrs=/         {print "chrs=(" chrs ")"; next}
		/^job_of_chr=/   {print "job_of_chr=" job_of_chr; next}
		/^job_in_chr=/   {print "job_in_chr=" job_in_chr; next}
		{print}
	' "$ibdmix_sh" > "$dirout/ibdmix.sh"

	chmod +x "$dirout/ibdmix.sh"
	(
		cd "$dirout"
		bsub -q "$queue" -n "$nthread" -J "ibdmix_$grp" -o "job.$grp.LOG" -e "job.$grp.ERR" < ibdmix.sh
	)
done
