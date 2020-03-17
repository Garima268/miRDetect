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
    csv_file2 = csv.reader(open('mature-seq.csv', "r"), delimiter=",")
    name = []
    preseq = []
    prelen = []
    mat_id = []
    mat_seq = []
    mat_len = []
    mfe = []
    rows = []
    Pred_pos = 0
    Pred_neg = 0
    for row in csv_file:
        if row[-1] == 'Positive':
            Pred_pos+=1
            name.append(row[0])
            preseq.append(row[1])
            prelen.append(row[3])
    for row1 in csv_file2:
        for element in list(name):
            if row1[0] == element:
                mat_id.append(row1[1])
                mat_seq.append(row1[2])
                mat_len.append(len(row1[2]))
                mfe.append(row[88])
        rows = zip(name, preseq, prelen, mat_id, mat_seq, mat_len, mfe)
    else:
            Pred_neg+=1
    print("Novel Precursors found:" ,Pred_pos)

    with open("novel-microrna.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["Pre-mirna-id", "Pre-mirna", "Pre-mirna-length", "Mature id", "Mature sequence", "Mature length", "MFE"])
        for row in rows:
            writer.writerow(row)
    
    