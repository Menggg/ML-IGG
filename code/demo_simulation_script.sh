#please install required packages
#install.packages("ranger")
#install.packages("data.table")

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
       /usr/bin/Rscript prediction_noBase_degree_RF_EUR_miss.R $i $n $miss
done

##decide top common features selection with Random Forest across error rates
#for relationship degree
#the input data is multiple mearsurements
#the output data is the prediction accuracy after the top common features was added one by one
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
       /usr/bin/Rscript best_final_prediction_sim_RF_degree.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt 
done

##Final classification with top common features selection with Random Forest
#for relationship degree
#the input data is multiple mearsurements
#the output data is prediction accuracy of each round of ten-fold cross validation
miss=0
for p in $(seq 1 1 6);do
        err=$(awk 'NR=='$p' {print $1}' err_list)
/usr/bin/Rscript final_prediction_sim_RF_degree.R $err merge_EUR_miss0_sim_table_multi_measurements_err0_1_100.txt merge_EUR_miss0_sim_table_multi_measurements_err"$err"_1_100.txt $miss 
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
        /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree1.R $i $n $miss
        /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree2.R $i $n $miss
        /usr/bin/Rscript prediction_noBase_relationship_RF_EUR_miss_degree3.R $i $n $miss 
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




