#Downloads 1000g data from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/ 
#merge files and extract SNP only
for i in {1..22};do
nohup bcftools annotate -x INFO,^FORMAT/GT CCDG_14151_B01_GRM_WGS_2020-08-05_chr"$i".filtered.shapeit2-duohmm-phased.vcf.gz -Oz -o WGS_2020-08-05_chr"$i"_GT.vcf.gz &
done

for i in {1..22};do
nohup bcftools view --max-alleles 2 WGS_2020-08-05_chr"$i"_GT.vcf.gz -Oz -o WGS_2020-08-05_chr"$i"_GT_2alleles.vcf.gz &
done

> list
for i in {1..22};do echo WGS_2020-08-05_chr"$i"_GT_2alleles.vcf.gz >> list; done

nohup bcftools concat --naive -f list -Oz -o WGS_2020-08-05_merge_GT_2alleles.vcf.gz & 


#extract GSA SNPs (GSA_hg38.list), 503 European white peoples (EUR.list) and 22 automosomes
plink2 --vcf WGS_2020-08-05_merge_GT_2alleles.vcf.gz --extract GSA_hg38.list --keep EUR.list --autosome --make-pgen --out WGS_2020-08-05_GSA_GT_2alleles_EUR 

#simulate genotype data using 503 EUR samples as founders using ped-sim
#err_list include the error rates you want to include, one rate one row
#multi_measurements_cross_single_family.R used to calculate the measurements of each pair
#sim_table_multi_measurements_balance.txt save the pre-calculated relationship
n1=1
n2=100
mkdir err0
for j in $(seq $n1 1 $n2);do
sort -R EUR.list | awk 'NR<70' > err0/temp_"$j".list
plink2 --pfile WGS_2020-08-05_GSA_GT_2alleles_EUR --keep err0/temp_"$j".list --export vcf --out err0/temp_"$j"

k=$RANDOM$RANDOM
for p in $(seq 1 1 6);do
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

##merge all the simulated measurements files
#n is number of simulated families
n=100
miss=0
for p in $(seq 1 1 6);do
        i=$(awk 'NR=='$p' {print $1}' err_list)
        awk 'NR==1' miss0_sim_table_multi_measurements_err0_1_"$n".txt > merge_EUR_miss"$miss"_sim_table_multi_measurements_err"$i"_1_"$n".txt
        for j in $(seq 1 1 $n);do
                grep -v "individual_1" miss0_sim_table_multi_measurements_err"$i"_"$j"_"$n".txt >> merge_EUR_miss"$miss"_sim_table_multi_measurements_err"$i"_1_"$n".txt
        done
done

##Feature selection with Random Forest and SVM algorithm
##miss is missing SNP rate
##n is number of simulated families
##i is genotyping error rate
#for relationship degree
#the input data is multiple mearsurements
#the output data is the top feature of each step of forward selection 
n=100
miss=0
for p in $(seq 1 1 6);do
        i=$(awk 'NR=='$p' {print $1}' err_list)
       nohup /usr/bin/Rscript prediction_noBase_degree_RF_EUR_miss.R $i $n $miss &
done

##decide top common features selection with Random Forest across error rates
#for relationship degree
#the input data is multiple mearsurements
#the output data is the prediction accuracy after the top common features was added one by one
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
       nohup /usr/bin/Rscript best_final_prediction_sim_RF_degree.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt &
done

##Final classification with top common features selection with Random Forest
#for relationship degree
#the input data is multiple mearsurements
#the output data is prediction accuracy of each round of ten-fold cross validation
miss=0
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
nohup /usr/bin/Rscript final_prediction_sim_RF_degree.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt $miss &
done

##Hierarchical Feature selection with Random Forest and SVM algorithm for relationship type
##miss is missing SNP rate
##n is number of simulated families
##i is genotyping error rate
#the input data is multiple mearsurements
#the output data is the top feature of each step of forward selection 
n=100
miss=0
for p in $(seq 1 1 6);do
        i=$(awk 'NR=='$p' {print $1}' err_list)
        nohup /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree1.R $i $n $miss &
        nohup /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree2.R $i $n $miss &
        nohup /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree3.R $i $n $miss & 
done

##decide top common features selection with Random Forest across error rates
#for relationship type
#the input data is multiple mearsurements
#the output data is the prediction accuracy after the top common features was added one by one
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
       /usr/bin/Rscript best_final_prediction_sim_relationship_randomForest_degree1.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt
	/usr/bin/Rscript best_final_prediction_sim_relationship_randomForest_degree2.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt
	/usr/bin/Rscript best_final_prediction_sim_relationship_randomForest_degree3.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt
done

##Final classification with top common features selection with Random Forest
#for relationship type
#the input data is multiple mearsurements
#the output data (predicted relationship type) will be saved at relationship folder
miss=0
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
/usr/bin/Rscript final_prediction_sim_relationship_randomForest_degree1.R $err  merge_EUR_miss"$miss"_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss"$miss"_sim_table_multi_measurements_err"$err"_1_100.txt $miss
/usr/bin/Rscript final_prediction_sim_relationship_randomForest_degree2.R $err  merge_EUR_miss"$miss"_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss"$miss"_sim_table_multi_measurements_err"$err"_1_100.txt $miss
/usr/bin/Rscript final_prediction_sim_relationship_randomForest_degree3.R $err  merge_EUR_miss"$miss"_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss"$miss"_sim_table_multi_measurements_err"$err"_1_100.txt $miss
done


##for real family data
##to run these two script, user need download genotype data (VCF files) from GEO database
Rscript multi_measurements_realData_DNAChip.R
Rscript final_prediction_real_chip_RF_degree.R



