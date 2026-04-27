import pandas as pd
import joblib
import csv
import logging

logging.basicConfig(level=logging.INFO)


def predictseq(feature_file: str):
    """
    Predict miRNA precursor sequences using trained Random Forest model.

    Input:
        feature_file (str): Path to feature CSV

    Outputs:
        - features.csv
        - novel-microrna.csv
    """

    # ---------------------- Load Model ---------------------- #
    try:
        RF = joblib.load("RF_mirna.pkl")
    except FileNotFoundError:
        raise FileNotFoundError("RF_mirna.pkl model file not found")

    # ---------------------- Load Features ---------------------- #
    sample = pd.read_csv(feature_file, header=None)

    # Features start from column index 3
    features = sample.iloc[:, 3:]

    # ---------------------- Prediction ---------------------- #
    predictions = RF.predict(features)
    sample["Pred"] = predictions

    # Save updated features
    sample.to_csv("features.csv", index=False, header=False)

    # ---------------------- Filter Positive Predictions ---------------------- #
    positive_df = sample[sample["Pred"] == "Positive"]

    Pred_pos = len(positive_df)
    Pred_neg = len(sample) - Pred_pos

    logging.info(f"Novel Precursors found: {Pred_pos}")

    if Pred_pos == 0:
        logging.warning("No positive predictions found")
        return

    # Extract required columns
    names = positive_df.iloc[:, 0].tolist()
    preseq = positive_df.iloc[:, 1].tolist()
    prelen = positive_df.iloc[:, 3].tolist()

    # ---------------------- Load Mature Sequences ---------------------- #
    try:
        mature_df = pd.read_csv("mature-seq.csv", header=None)
    except FileNotFoundError:
        raise FileNotFoundError("mature-seq.csv not found")

    # Create lookup dictionary (O(1) instead of nested loops)
    mature_dict = {
        row[0]: (row[1], row[2], len(row[2]))
        for _, row in mature_df.iterrows()
    }

    # ---------------------- Combine Data ---------------------- #
    output_rows = []

    for idx, name in enumerate(names):
        if name in mature_dict:
            mat_id, mat_seq, mat_len = mature_dict[name]

            # Safe extraction of MFE (column 88 may not always exist)
            mfe = positive_df.iloc[idx, 88] if positive_df.shape[1] > 88 else ""

            output_rows.append([
                name,
                preseq[idx],
                prelen[idx],
                mat_id,
                mat_seq,
                mat_len,
                mfe
            ])

    # ---------------------- Write Output ---------------------- #
    with open("novel-microrna.csv", "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "Pre-mirna-id",
            "Pre-mirna",
            "Pre-mirna-length",
            "Mature id",
            "Mature sequence",
            "Mature length",
            "MFE"
        ])

        writer.writerows(output_rows)
