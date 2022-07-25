n1=1
n2=1
mkdir err0
for j in $(seq $n1 1 $n2);do
sort -R EUR.list | awk 'NR<70' > err0/temp_"$j".list
plink2 --pfile WGS_2020-08-05_GSA_GT_2alleles_EUR --keep err0/temp_"$j".list --export vcf --out err0/temp_"$j"

k=$RANDOM$RANDOM
for p in $(seq 1 1 4);do
        i=$(awk 'NR=='$p' {print $1}' err_list)
mkdir err"$i"
ped-sim -d 1fam.def -m refined_mf_hg38_noCHR.simmap --intf sex_av_nu_p_hg38_campbell_noCHR.tsv -i err0/temp_"$j".vcf -o err"$i"/1fam_err"$i"_"$j" --err_rate "$i" --seed "$k"
plink2 --vcf err"$i"/1fam_err"$i"_"$j".vcf --export A-transpose --out err"$i"/1fam_err"$i"_"$j"
Rscript multi_measurements_cross_single_family.R err0/1fam_err0_"$j".traw err"$i"/1fam_err"$i"_"$j".traw "$i" "$j"
done
rm err*/1fam_err*_"$j".traw
rm err*/1fam_err*_"$j".vcf
rm err0/temp_"$j".vcf
rm err0/temp_"$j".list
done
