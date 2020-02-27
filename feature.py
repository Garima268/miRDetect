import sys
import os
from config import path_vienna
vienna = os.path.join(path_vienna, "interfaces/Python3/")
sys.path.append(vienna)
import RNA
from Bio import SeqIO
import csv
import math
from Bio.SeqUtils import MeltingTemp as mt
# The RNA sequence
def featurecount(fasta_file):
    record = list(SeqIO.parse(fasta_file,"fasta"))
    names = []
    ATGC = []
    mfes= []
    sst =[]
    lgt = []
    mfei = []
    mfei2 = []
    mfei3 = []
    mfei4 = []
    amfe = []
    EMFE = []
    Tm = []
    NTm = []
    RGC = []
    FGC = []
    FAU = []
    mfei5 = []
    TBP = []
    aup = []
    gup = []
    gcp = []
    aus = []
    gcs = []
    gus = []
    FMFE = []
    ED = []
    RBPGC = []
    RBPGU = []
    RBPAU = []
    CE = []
    CD = []
    PA = []
    PG = []
    PC = []
    PU = []
    PAA = []
    PAU = []
    PAG = []
    PAC = []
    PUA = []
    PUG = []
    PUC = []
    PUU = []
    PGA = []
    PGC = []
    PGG = []
    PGU = []
    PCA = []
    PCG = []
    PCC = []
    PCU = []
    PUCA = []
    PUCC = []
    PUCG = []
    PUCU = []
    PUUC = []
    PUUU = []
    PUUA = []
    PUUG = []
    PUAC = []
    PUAU = []
    PUAA = []
    PUAG = []
    PUGC = []
    PUGU = []
    PUGA = []
    PUGG = []
    PCUA = []
    PCUC = []
    PCUG = []
    PCUU = []
    PCCA = []
    PCCC = []
    PCCG = []
    PCCU = []
    PCAC = []
    PCAU = []
    PCAA = []
    PCAG = []
    PCGA = []
    PCGC = []
    PCGG = []
    PCGU = []
    PAUA = []
    PAUC = []
    PAUU = []
    PAUG = []
    PACA = []
    PACC = []
    PACG = []
    PACU = []
    PAAC = []
    PAAU = []
    PAAA = []
    PAAG = []
    PAGC = []
    PAGU = []
    PAGA = []
    PAGG = []
    PGUA = []
    PGUC = []
    PGUG = []
    PGUU = []
    PGCA = []
    PGCC = []
    PGCG = []
    PGCU = []
    PGAC = []
    PGAU = []
    PGAA = []
    PGAG = []
    PGGA = []
    PGGC = []
    PGGG = []
    PGGU = []

    for i in range(0,len(record)):
        nm = record[i].name
        seq = record[i].seq
        seq=str(seq)
        length = len(seq)
        MeltTemp = mt.Tm_NN(seq)
        NMeltTemp = MeltTemp/length
        C_count = seq.count('C')
        G_count = seq.count('G')
        A_count = seq.count('A')
        U_count = seq.count('U')
        AA_count = seq.count('AA')
        AU_count = seq.count('AU')
        AG_count = seq.count('AG')
        AC_count = seq.count('AC')
        UA_count = seq.count('UA')
        UG_count = seq.count('UG')
        UC_count = seq.count('UC')
        UU_count = seq.count('UU')
        GA_count = seq.count('GA')
        GC_count = seq.count('GC')
        GG_count = seq.count('GG')
        GU_count = seq.count('GU')
        CA_count = seq.count('CA')
        CG_count = seq.count('CG')
        CC_count = seq.count('CC')
        CU_count = seq.count('CU')
        UCA_count = seq.count('UCA')
        UCC_count = seq.count('UCC')
        UCG_count = seq.count('UCG')
        UCU_count = seq.count('UCU')
        UUC_count = seq.count('UUC')
        UUU_count = seq.count('UUU')
        UUA_count = seq.count('UUA')
        UUG_count = seq.count('UUG')
        UAC_count = seq.count('UAC')
        UAU_count = seq.count('UAU')
        UAA_count = seq.count('UAA')
        UAG_count = seq.count('UAG')
        UGC_count = seq.count('UGC')
        UGU_count = seq.count('UGU')
        UGA_count = seq.count('UGA')
        UGG_count = seq.count('UGG')
        CUA_count = seq.count('CUA')
        CUC_count = seq.count('CUC')
        CUG_count = seq.count('CUG')
        CUU_count = seq.count('CUU')
        CCA_count = seq.count('CCA')
        CCC_count = seq.count('CCC')
        CCG_count = seq.count('CCG')
        CCU_count = seq.count('CCU')
        CAC_count = seq.count('CAC')
        CAU_count = seq.count('CAU')
        CAA_count = seq.count('CAA')
        CAG_count = seq.count('CAG')
        CGA_count = seq.count('CGA')
        CGC_count = seq.count('CGC')
        CGG_count = seq.count('CGG')
        CGU_count = seq.count('CGU')
        AUA_count = seq.count('AUA')
        AUC_count = seq.count('AUA')
        AUU_count = seq.count('AUU')
        AUG_count = seq.count('AUG')
        ACA_count = seq.count('ACA')
        ACC_count = seq.count('ACC')
        ACG_count = seq.count('ACG')
        ACU_count = seq.count('ACU')
        AAC_count = seq.count('AAC')
        AAU_count = seq.count('AAU')
        AAA_count = seq.count('AAA')
        AAG_count = seq.count('AAG')
        AGC_count = seq.count('AGC')
        AGU_count = seq.count('AGU')
        AGA_count = seq.count('AGA')
        AGG_count = seq.count('AGG')
        GUA_count = seq.count('GUA')
        GUC_count = seq.count('GUC')
        GUG_count = seq.count('GUG')
        GUU_count = seq.count('GUU')
        GCA_count = seq.count('GCA')
        GCC_count = seq.count('GCC')
        GCG_count = seq.count('GCG')
        GCU_count = seq.count('GCU')
        GAC_count = seq.count('GAC')
        GAU_count = seq.count('GAU')
        GAA_count = seq.count('GAA')
        GAG_count = seq.count('GAG')
        GGA_count = seq.count('GGA')
        GGC_count = seq.count('GGC')
        GGG_count = seq.count('GGG')
        GGU_count = seq.count('GGU')
        perc_a = float(A_count /length) * 100
        PA.append(perc_a)
        perc_g = float(G_count / length) * 100
        PG.append(perc_g)
        perc_c = float(C_count / length) * 100
        PC.append(perc_c)
        perc_u = float(U_count / length) * 100
        PU.append(perc_u)
        aa_percentage = float(AA_count / length) * 100
        PAA.append(aa_percentage)
        au_percentage = float(AU_count / length) * 100
        PAU.append(au_percentage)
        ag_percentage = float(AG_count / length) * 100
        PAG.append(ag_percentage)
        ac_percentage = float(AC_count / length) * 100
        PAC.append(ac_percentage)
        UA_percentage = float(UA_count / length) * 100
        PUA.append(UA_percentage)
        UG_percentage = float(UG_count / length) * 100
        PUG.append(UG_percentage)
        UC_percentage = float(UC_count / length) * 100
        PUC. append(UC_percentage)
        UU_percentage = float(UU_count / length) * 100
        PUU.append(UU_percentage)
        GA_percentage = float(GA_count / length) * 100
        PGA.append(GA_percentage)
        GC_percentage = float(GC_count / length) * 100
        PGC.append(GC_percentage)
        GG_percentage = float(GG_count / length) * 100
        PGG.append(GG_percentage)
        GU_percentage = float(GU_count / length) * 100
        PGU.append(GU_percentage)
        CA_percentage = float(CA_count / length) * 100
        PCA.append(CA_percentage)
        CG_percentage = float(CG_count / length) * 100
        PCG.append(CG_percentage)
        CC_percentage = float(CC_count / length) * 100
        PCC.append(CC_percentage)
        CU_percentage = float(CU_count / length) * 100
        PCU.append(CU_percentage)
        UCA_percentage = float(UCA_count / length) * 100
        PUCA.append(UCA_percentage)
        UCC_percentage = float(UCC_count / length) * 100
        PUCC.append(UCC_percentage)
        UCG_percentage = float(UCG_count / length) * 100
        PUCG.append(UCG_percentage)
        UCU_percentage = float(UCU_count / length) * 100
        PUCU.append(UCU_percentage)
        UUC_percentage = float(UUC_count / length) * 100
        PUUC.append(UUC_percentage)
        UUU_percentage = float(UUG_count / length) * 100
        PUUU.append(UUU_percentage)
        UUA_percentage = float(UUA_count / length) * 100
        PUUA.append(UUA_percentage)
        UUG_percentage = float(UUG_count / length) * 100
        PUUG.append(UUG_percentage)
        UAC_percentage = float(UAC_count / length) * 100
        PUAC.append(UAC_percentage)
        UAU_percentage = float(UAU_count / length) * 100
        PUAU.append(UAU_percentage)
        UAA_percentage = float(UAA_count / length) * 100
        PUAA.append(UAA_percentage)
        UAG_percentage = float(UAG_count / length) * 100
        PUAG.append(UAG_percentage)
        UGC_percentage = float(UGC_count / length) * 100
        PUGC.append(UGC_percentage)
        UGU_percentage = float(UGU_count / length) * 100
        PUGU.append(UGU_percentage)
        UGA_percentage = float(UGA_count / length) * 100
        PUGA.append(UGA_percentage)
        UGG_percentage = float(UGG_count / length) * 100
        PUGG.append(UGG_percentage)
        CUA_percentage = float(CUA_count / length) * 100
        PCUA.append(CUA_percentage)
        CUC_percentage = float(CUC_count / length) * 100
        PCUC.append(CUC_percentage)
        CUG_percentage = float(CUG_count / length) * 100
        PCUG.append(CUG_percentage)
        CUU_percentage = float(CUU_count / length) * 100
        PCUU.append(CUU_percentage)
        CCA_percentage = float(CCA_count / length) * 100
        PCCA.append(CCA_percentage)
        CCC_percentage = float(CCC_count / length) * 100
        PCCC.append(CCC_percentage)
        CCG_percentage = float(CCG_count / length) * 100
        PCCG.append(CCG_percentage)
        CCU_percentage = float(CCU_count / length) * 100
        PCCU.append(CCU_percentage)
        CAC_percentage = float(CAC_count / length) * 100
        PCAC.append(CAC_percentage)
        CAU_percentage = float(CAU_count / length) * 100
        PCAU.append(CAU_percentage)
        CAA_percentage = float(CAA_count / length) * 100
        PCAA.append(CAA_percentage)
        CAG_percentage = float(CAG_count / length) * 100
        PCAG.append(CAG_percentage)
        CGA_percentage = float(CGA_count / length) * 100
        PCGA.append(CGA_percentage)
        CGC_percentage = float(CGC_count / length) * 100
        PCGC.append(CGC_percentage)
        CGG_percentage = float(CGG_count / length) * 100
        PCGG.append(CGG_percentage)
        CGU_percentage = float(CGU_count / length) * 100
        PCGU.append(CGU_percentage)
        AUA_percentage = float(AUA_count / length) * 100
        PAUA.append(AUA_percentage)
        AUC_percentage = float(AUC_count / length) * 100
        PAUC.append(AUC_percentage)
        AUU_percentage = float(AUU_count / length) * 100
        PAUU.append(AUU_percentage)
        AUG_percentage = float(AUG_count / length) * 100
        PAUG.append(AUG_percentage)
        ACA_percentage = float(ACA_count / length) * 100
        PACA.append(ACA_percentage)
        ACC_percentage = float(ACC_count / length) * 100
        PACC.append(ACC_percentage)
        ACG_percentage = float(ACG_count / length) * 100
        PACG.append(ACG_percentage)
        ACU_percentage = float(ACU_count / length) * 100
        PACU.append(ACU_percentage)
        AAC_percentage = float(AAC_count / length) * 100
        PAAC.append(AAC_percentage)
        AAU_percentage = float(AAU_count / length) * 100
        PAAU.append(AAU_percentage)
        AAA_percentage = float(AAA_count / length) * 100
        PAAA.append(AAA_percentage)
        AAG_percentage = float(AAG_count / length) * 100
        PAAG.append(AAG_percentage)
        AGC_percentage = float(AGC_count / length) * 100
        PAGC.append(AGC_percentage)
        AGU_percentage = float(AGU_count / length) * 100
        PAGU.append(AGU_percentage)
        AGA_percentage = float(AGA_count / length) * 100
        PAGA.append(AGA_percentage)
        AGG_percentage = float(AGG_count / length) * 100
        PAGG.append(AGG_percentage)
        GUA_percentage = float(GUA_count / length) * 100
        PGUA.append(GUA_percentage)
        GUC_percentage = float(GUC_count / length) * 100
        PGUC.append(GUC_percentage)
        GUG_percentage = float(GUG_count / length) * 100
        PGUG.append(GUG_percentage)
        GUU_percentage = float(GUU_count / length) * 100
        PGUU.append(GUU_percentage)
        GCA_percentage = float(GCA_count / length) * 100
        PGCA.append(GCA_percentage)
        GCC_percentage = float(GCC_count / length) * 100
        PGCC.append(GCC_percentage)
        GCG_percentage = float(GCG_count / length) * 100
        PGCG.append(GCG_percentage)
        GCU_percentage = float(GCU_count / length) * 100
        PGCU.append(GCU_percentage)
        GAC_percentage = float(GAC_count / length) * 100
        PGAC.append(GAC_percentage)
        GAU_percentage = float(GAU_count / length) * 100
        PGAU.append(GAU_percentage)
        GAA_percentage = float(GAA_count / length) * 100
        PGAA.append(GAA_percentage)
        GAG_percentage = float(GAG_count / length) * 100
        PGAG.append(GAG_percentage)
        GGA_percentage = float(GGA_count / length) * 100
        PGGA.append(GGA_percentage)
        GGC_percentage = float(GGC_count / length) * 100
        PGGC.append(GGC_percentage)
        GGG_percentage = float(GGG_count / length) * 100
        PGGG.append(GGG_percentage)
        GGU_percentage = float(GGU_count / length) * 100
        PGGU.append(GGU_percentage)
        cg_percentage = float(C_count + G_count) / length
        au_percentage = float(A_count + U_count) / length
        RatioGC = float(G_count/C_count)
    # compute minimum free energy (MFE) and corresponding structure
        (ss, mfe) = RNA.fold(seq)
        fc = RNA.fold_compound(seq)
        (pp, pf) = fc.pf()
        npf = float(pf / length)
        EMFE.append(pf)
        frequency_mfe_struc = fc.pr_structure(ss)
        FMFE.append(frequency_mfe_struc)
        ensemble_diversity = fc.mean_bp_distance()
        ED.append(ensemble_diversity)
        # compute centroid structure
        (centroid_struct, dist) = fc.centroid()
        # compute free energy of centroid structure
        centroid_en = fc.eval_structure(centroid_struct)
        CE.append(centroid_en)
        CD.append(dist)
        names.append(nm)
        ATGC.append(seq)
        sst.append(ss)
        mfes.append(mfe)
        lgt.append(length)
        Tm.append(MeltTemp)
        NTm.append(NMeltTemp)
        dG = mfe/length
        amfe.append(dG)
        mfe_in5 = dG/au_percentage
        mfei5.append(mfe_in5)
        mfe_in = dG/cg_percentage
        mfe_in4 = dG/length
        mfei.append(mfe_in)
        mfei4.append(mfe_in4)
        RGC.append(RatioGC)
        FGC.append(cg_percentage)
        FAU.append(au_percentage)
        tbpc, aubp, gcbp, gubp, nos, nol = annotatefold(seq, ss)
        aul = aubp/length
        percentage_aubp = aul * 100
        gul = gubp/length
        percentage_gubp = gul * 100
        gcl = gcbp/length
        percentage_gcbp = gcl * 100
        TBP.append(tbpc)
        aup.append(aul)
        gup.append(gul)
        gcp.append(gcl)
        auso = percentage_aubp/nos
        aus.append(auso)
        guso = percentage_gubp/nos
        gus.append(guso)
        gcso = percentage_gcbp/nos
        gcs.append(gcso)
        if gcbp == 0:
            RatioBPGC = 0
        else:
            RatioBPGC = tbpc/gcbp
            RBPGC.append(RatioBPGC)
        if gubp == 0:
            RatioBPGU = 0
        else:
            RatioBPGU = tbpc / gubp
            RBPGU.append(RatioBPGU)
        if aubp == 0:
            RatioBPAU = 0
        else:
            RatioBPAU = tbpc / aubp
            RBPAU.append(RatioBPAU)
        mfe_in2 = dG/nos
        mfei2.append(mfe_in2)
        if nol == 0:
            mfe_in3 = 0
            mfei3.append(mfe_in3)
        else:
            mfe_in3 = dG / nol
            mfei3.append(mfe_in3)
        rows = zip(names, ATGC, sst, lgt, PA, PG, PC, PU, PAA, PAC, PAG, PAU, PGA, PGC, PGG, PGU, PCA, PCG, PCU, PCC, PUA,
                   PUG, PUC, PUU, PUCA, PUCC, PUCG, PUCU, PUUC, PUUU, PUUA, PUUG, PUAC, PUAU, PUAA, PUAG, PUGC, PUGU, PUGA, PUGG,
                    PCUA, PCUC, PCUG, PCUU, PCCA, PCCC, PCCG, PCCU, PCAC, PCAU, PCAA, PCAG, PCGA, PCGC, PCGG, PCGU, PAUA, PAUC,
                    PAUU, PAUG, PACA, PACC, PACG, PACU, PAAC, PAAU, PAAA, PAAG, PAGC, PAGU, PAGA, PAGG, PGUA, PGUC, PGUG, PGUU,
                    PGCA, PGCC, PGCG, PGCU, PGAC, PGAU, PGAA, PGAG, PGGA, PGGC, PGGG, PGGU, mfes, FGC, FAU, amfe, EMFE, mfei,
                   mfei2, mfei3, mfei4, mfei5, FMFE, ED, Tm, NTm, RGC, TBP,
                   aup, gup, gcp, aus, gus, gcs, RBPAU, RBPGU, RBPGC, CE, CD)
    
    with open("features.csv", "w") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)

                   
def annotatefold(rna, foldstructure):
        count_au =0
        count_gc = 0
        tbp = 0
        count_gu = 0
        ns = 0
        nl = 0
        for i in range(len(rna)):
            if foldstructure[i] == '(':
                tbp = tbp +1
            if foldstructure[i] == '(' and foldstructure[i + 1] == '(' and foldstructure[i + 2] == '(' and \
                        foldstructure[i + 3] == '.':
                ns = ns + 1
            if foldstructure[i] == '(' and foldstructure[i + 1] == '.' and foldstructure[i + 2] == '.' and \
                        foldstructure[i + 3] == '.' and foldstructure[i + 4] == '.':
                nl = nl + 1
            if rna[i] == 'C' and (foldstructure[i] == '(' or foldstructure[i] == ')'):
                count_gc = count_gc +1
            if rna[i] == 'A' and (foldstructure[i] == '(' or foldstructure[i] == ')'):
                    count_au = count_au +1
        count_gu = tbp - (count_au + count_gc)
        return (tbp, count_au, count_gc, count_gu, ns, nl)

    

