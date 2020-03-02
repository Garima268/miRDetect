import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
import joblib
import pandas as pd
import csv
#Predict precursor sequences

def predictseq(feature_file):
    RF = joblib.load("RF_mirna.pkl")
    Pred_pos = 0
    Pred_neg = 0
    sample = pd.read_csv(feature_file, header=None)
    K = sample[sample.columns[3:]]
    Pred = RF.predict(K)
    sample["Pred"] = Pred
    make_pred= sample.to_csv("features.csv", index=False, header=False)

    csv_file = csv.reader(open('features.csv', "r"), delimiter=",")
    name = []
    preseq = []
    for row in csv_file:
        if row[-1] == 'Positive':
            name.append(row[0])
            preseq.append(row[1])
            rows = zip(name, preseq)
            Pred_pos+=1
        else:
            Pred_neg+=1
    print("Novel Precursors found:",Pred_pos)   
    with open("novel-microrna.csv", "w") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
    

        
        


