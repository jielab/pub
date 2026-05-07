#!/bin/bash

dir0=/scratch/2026-05-06/sph-huangj
dirmod=$dir0/imp
dirarc=/data/sph-huangj/refGen/gu
dirsoft=/work/sph-huangj/software/IBDmix
dirscript=/work/sph-huangj/scripts/gu
ibdmix_sh=$dirscript/ibdmix.sh
sample_info=$dir0/files/ukb.sample_info.txt
base_out=$dir0/gu/ibdmix

batch_n=1000
queue=short
njob=8
run_chrs="$(seq 1 22 | paste -sd' ' -) X"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Split UKB samples into 1000-sample batches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p "$base_out/sample_keep"
all_id=$base_out/sample_keep/ukb.all.1id
awk 'NR>1 {print $2}' "$dirmod/chr1.psam" > "$all_id"

rm -f "$base_out"/sample_keep/psam.*.1id
split -d -a 3 -l "$batch_n" --numeric-suffixes=1 --additional-suffix=.1id "$all_id" "$base_out/sample_keep/psam."

nbatch=$(find "$base_out/sample_keep" -maxdepth 1 -name 'psam.*.1id' | wc -l)
echo "Total samples: $(wc -l < "$all_id")"
echo "Batch size:    $batch_n"
echo "Total batches: $nbatch"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Generate modified ibdmix.sh per batch, then bsub
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for sample_file in "$base_out"/sample_keep/psam.00{1..3}.1id; do
	grp=$(basename "$sample_file" | sed -E 's/^psam\.([0-9]{3})\.1id$/\1/')
	dirout=$base_out/$grp
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
		-v chrs="$run_chrs" \
		-v njob="$njob" '
		/^dir0=/         {print "dir0=" dir0; next}
		/^dirmod=/       {print "dirmod=" dirmod; next}
		/^dirarc=/       {print "dirarc=" dirarc; next}
		/^dirsoft=/      {print "dirsoft=" dirsoft; next}
		/^dirscript=/    {print "dirscript=" dirscript; next}
		/^sample_label=/ {print "sample_label=" sample_label; next}
		/^sample_keep=/  {print "sample_keep=" sample_keep; next}
		/^sample_info=/  {print "sample_info=" sample_info; next}
		/^dirout=/       {print "dirout=" dirout; next}
		/^chrs=/         {print "chrs=\"" chrs "\""; next}
		/^njob=/         {print "njob=" njob; next}
		{print}
	' "$ibdmix_sh" > "$dirout/ibdmix.sh"

	chmod +x "$dirout/ibdmix.sh"
	cd "$dirout"
	bsub -q "$queue" -n "$njob" -R "span[hosts=1]" -J "ibdmix_$grp" -o "job.$grp.LOG" -e "job.$grp.ERR" < ibdmix.sh
done
