#Run command:

#python3 SVM.py data1.csv data2.csv label.csv . 

import os
from sys import argv
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
import rpy2.robjects as robjects
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

def run_SVM(Data1Path, Data2Path, LabelsPath, OutputDir):
    # read the Rdata file
    data1 = pd.read_csv(Data1Path,index_col=0,sep=',')
    data2 = pd.read_csv(Data2Path,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',')

    print(len(labels),"\n")
    
    # read the feature file
    data1 = np.log1p(data1)
    data2 = np.log1p(data2)
    
    #Classifier = LinearSVC()
    Classifier = make_pipeline( SVC(kernel = 'poly',coef0= 1))
    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    train=data1
    test=data2
    y_train=labels
    y_test=[]

    start=tm.time()
    Classifier.fit(train, y_train)
    tr_time.append(tm.time()-start)

    start=tm.time()
    predicted = Classifier.predict(test)
    ts_time.append(tm.time()-start)

    pred.extend(predicted)

    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    OutputDir = Path(OutputDir)
    pred.to_csv(str(OutputDir / Path("SVM_pred.csv")),
                index = False)
    tr_time.to_csv(str(OutputDir / Path("SVM_training_time.csv")),
                   index = False)
    ts_time.to_csv(str(OutputDir / Path("SVM_test_time.csv")),
                   index = False)

run_SVM(argv[1], argv[2], argv[3], argv[4])
