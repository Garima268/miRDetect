import os
import sys
import csv
import logging
from collections import Counter

from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

from config import path_vienna

# ViennaRNA setup
vienna = os.path.join(path_vienna, "interfaces/Python3/")
sys.path.append(vienna)
import RNA

logging.basicConfig(level=logging.INFO)


# ---------------------- Helper ---------------------- #
def kmer_counts(seq, k):
    """Return normalized k-mer frequencies."""
    counts = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
    total = len(seq)
    return {kmer: (count / total) * 100 for kmer, count in counts.items()}


def annotatefold(rna, foldstructure):
    """Safe annotation of RNA fold structure."""
    count_au = count_gc = tbp = count_gu = ns = nl = 0

    for i in range(len(rna)):
        if foldstructure[i] == '(':
            tbp += 1

        # safe window checks
        if i + 3 < len(foldstructure):
            if foldstructure[i:i+4] == "(((.":
                ns += 1

        if i + 4 < len(foldstructure):
            if foldstructure[i:i+5] == "(....":
                nl += 1

        if foldstructure[i] in "()":
            if rna[i] == 'C':
                count_gc += 1
            elif rna[i] == 'A':
                count_au += 1

    count_gu = tbp - (count_au + count_gc)

    return tbp, count_au, count_gc, count_gu, ns, nl


# ---------------------- Main Function ---------------------- #
def featurecount(fasta_file: str):
    """
    Extract RNA features and save to features.csv
    """

    records = list(SeqIO.parse(fasta_file, "fasta"))
    logging.info(f"Processing {len(records)} sequences")

    rows = []

    for record in records:
        name = record.name
        seq = str(record.seq)
        length = len(seq)

        if length == 0:
            continue

        # ---------------------- Basic counts ---------------------- #
        counts = Counter(seq)

        A = counts.get('A', 0)
        U = counts.get('U', 0)
        G = counts.get('G', 0)
        C = counts.get('C', 0)

        PA = (A / length) * 100
        PG = (G / length) * 100
        PC = (C / length) * 100
        PU = (U / length) * 100

        gc_ratio = (G + C) / length
        au_ratio = (A + U) / length
        ratio_gc = G / C if C != 0 else 0

        # ---------------------- k-mers ---------------------- #
        dinuc = kmer_counts(seq, 2)
        trinuc = kmer_counts(seq, 3)

        # ---------------------- Thermodynamics ---------------------- #
        try:
            Tm = mt.Tm_NN(seq)
        except:
            Tm = 0

        NTm = Tm / length if length else 0

        ss, mfe = RNA.fold(seq)
        fc = RNA.fold_compound(seq)
        _, pf = fc.pf()

        EMFE = pf
        FMFE = fc.pr_structure(ss)
        ED = fc.mean_bp_distance()

        centroid_struct, dist = fc.centroid()
        CE = fc.eval_structure(centroid_struct)
        CD = dist

        dG = mfe / length if length else 0

        # ---------------------- Structural features ---------------------- #
        tbp, aubp, gcbp, gubp, nos, nol = annotatefold(seq, ss)

        aup = aubp / length
        gup = gubp / length
        gcp = gcbp / length

        # ratios
        RBPGC = tbp / gcbp if gcbp else 0
        RBPGU = tbp / gubp if gubp else 0
        RBPAU = tbp / aubp if aubp else 0

        mfei = dG / gc_ratio if gc_ratio else 0
        mfei2 = dG / nos if nos else 0
        mfei3 = dG / nol if nol else 0
        mfei4 = dG / length if length else 0
        mfei5 = dG / au_ratio if au_ratio else 0

        # ---------------------- Row ---------------------- #
        row = [
            name, seq, ss, length,
            PA, PG, PC, PU,
            mfe, gc_ratio, au_ratio, dG,
            EMFE, mfei, mfei2, mfei3, mfei4, mfei5,
            FMFE, ED, Tm, NTm, ratio_gc,
            tbp, aup, gup, gcp,
            RBPAU, RBPGU, RBPGC,
            CE, CD
        ]

        # append k-mer features (consistent order)
        for k in sorted(dinuc):
            row.append(dinuc[k])
        for k in sorted(trinuc):
            row.append(trinuc[k])

        rows.append(row)

    # ---------------------- Write Output ---------------------- #
    with open("features.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)

    logging.info("Feature extraction completed → features.csv")
