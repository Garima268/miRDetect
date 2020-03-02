import os
import sys
import Bio
from Bio import SeqIO
from config import path_db, path_blast, p_name
from shutil import copy
import subprocess
from Bio import SearchIO
from pyfasta import Fasta
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p",  action="store", dest="num", type="int", default=1, help='Input number of CPUs')
options, args = parser.parse_args()
p  = str(options.num)
path_tool = os.getcwd()
print("*********Beginning miRNA Search from Data*************")
if (int(p)>1):
    Assembly_file = sys.argv[3]
else:
    Assembly_file = sys.argv[1]
    
Assembly_file = os.path.abspath(Assembly_file)

my_file = os.path.join(path_db, "mature.ndb")
test = os.path.exists(my_file)
if test == False:
    print("*********Creating Databases*************")
    os.chdir(path_blast)
    path_mirbase = os.path.join(path_db, "mature")
    db_name = os.path.join(path_db, "mature.fa")
    make_db = subprocess.run(["./makeblastdb", "-in", db_name,"-dbtype", "nucl" , "-parse_seqids",  "-out", path_mirbase])
else:
    print("Database exists")

print("*********Running Blast for microRNAs*************")
os.chdir(path_blast)
blast_file = os.path.join(path_tool, "mirna-output.tsv")
path_mirbase = os.path.join(path_db, "mature")
run_mirna_blast = subprocess.run(["./blastn", "-query", Assembly_file, "-db", path_mirbase, "-evalue", "0.001", "-word_size", "6", "-outfmt", "6", "-max_target_seqs", "1", "-num_threads",  p , "-out", blast_file])



os.chdir(path_tool)
#Check if file is empty
if os.path.getsize(blast_file) == 0:
    print("No miRNA found")
    sys.exit()
#Code to extract IDs from TSV file 
f1 = open("mirna-output.tsv")
query_id = ""
for blast_record in SearchIO.parse(f1, "blast-tab"):
    query_id += blast_record.id + "\n"
    f = open("mirnablastid.txt", 'w')
    f.write(query_id)
    f.close()
result_file = "candidate-mirna.fasta"
noncoding_file = "mirnablastid.txt"
wanted = set()
with open(noncoding_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)
fasta_sequences = SeqIO.parse(open(Assembly_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
f.close()
result_file = os.path.join(path_tool, "candidate-mirna.fasta")
mirid = os.path.join(path_tool, "mirnablastid.txt")
xstream = 100
est_fasta = Fasta(result_file)
final_file = open("Precursor-seq.fasta", "w")
for line in open(blast_file):
    # convert to int and 0-based coords.
    qstart, qstop, sstart, sstop = [int(x) - 1 for x in line.split("\t")[6:10]]
    query, subject = line.split("\t")[:2]
    if qstart - xstream <0:
        up = 0
    else:
        up = min(0, qstart - xstream)
    down = qstop + xstream + 1
    qbegin = qstart + 1
    qcomp = qstop - 1
    part_seq = est_fasta[query][up:down]
    print(">%s" % query)
    print(part_seq)


    final_file.write(">%s" % query)
    final_file.write("\n")
    final_file.write(part_seq)
    final_file.write("\n")
final_file.close()
final_file = os.path.join(path_tool, "Precursor-seq.fasta")

my_file2 = os.path.join(path_db, "uniprot.pdb")
test2 = os.path.exists(my_file2)
if test2 == False:
    print("*********Creating Databases*************")
    os.chdir(path_blast)
    prot_db = os.path.join(path_db, p_name)
    path_uniprot = os.path.join(path_db, "uniprot")
    make_db = subprocess.run(["./makeblastdb", "-in", prot_db,"-dbtype", "prot" , "-parse_seqids",  "-out", path_uniprot])
else:
    print("Database exists")
#Code to Run Blastx 
os.chdir(path_blast)
output = os.path.join(path_tool, "output.tsv")
path_uniprot = os.path.join(path_db, "uniprot")
run_blastx = subprocess.run(["./blastx", "-query", final_file, "-db", path_uniprot, "-evalue", "0.001", "-outfmt", "6", "-max_target_seqs", "1", "-num_threads",  p ,"-out", output])

os.chdir(path_tool)

#Check if file is empty
if os.path.getsize('output.tsv') == 0:
    print("No coding sequences found")
    
print("*********Extracting IDs from TSV file*************")
#Code to extract IDs from TSV file and fasta file
if os.path.getsize('output.tsv') > 0:

    f1 = open("output.tsv")
    query_id = ""
    for blast_record in SearchIO.parse(f1, "blast-tab"):
        query_id += blast_record.id + "\n"
        f = open("blastid.txt", 'w')
        f.write(query_id)
        f.close()
    f2 = open("Precursor-seq.fasta")
    fasta_id = ""
    for seq in SeqIO.parse(f2, "fasta"):
        fasta_id += seq.id + "\n"
        d = open("fastaid.txt", 'w')
        d.write(fasta_id)
        d.close()
    #Code to get IDs of non-coding sequences
    file1 = open("fastaid.txt",'r')
    lista1 = []

    for line in file1:
            line = line.rstrip('\n')
            lista1.append(line)
    file1.close()

    file2 = open("blastid.txt",'r')
    lista2 = []

    for line in file2:
            line = line.rstrip('\n')
            lista2.append(line)
    file2.close()

    same = set(lista1) - set(lista2)

    file_out = open('noncoding-id.txt', 'w')
    for line in same:
        file_out.write(line + "\n")
    file_out.close()

    blastid = os.path.join(path_tool, "blastid.txt")
    fastaid = os.path.join(path_tool, "fastaid.txt")
    ncid = os.path.join(path_tool, "noncoding-id.txt")
    print("*********Extracting non-coding sequences*************")
    #Code to extract non-coding sequences based on IDs
    nc_file = "noncoding-seq.fasta"
    ncid_file = "noncoding-id.txt"
    wanted = set()
    with open(ncid_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)
    fasta_sequences = SeqIO.parse(open(final_file),'fasta')
    with open(nc_file, "w") as f:
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")
    f.close()
    record = list(SeqIO.parse(nc_file,"fasta"))
    my_file = open("Temporary.fasta", "w")
    for i in range(0,len(record)):
        nm = record[i].name
        seq = record[i].seq
        seq = str(seq)
        my_file.write(">")
        my_file.write(nm)
        my_file.write("\n")
        my_file.write(seq.replace("T", "U"))
        my_file.write("\n")
    my_file.close()
    Temp = os.path.join(path_tool, "Temporary.fasta")



    nc_file = os.path.join(path_tool, "noncoding-seq.fasta")




    from feature import featurecount
    print("*********Calculating features of sequences*************")
    featurecount("Temporary.fasta")
    from Predict import predictseq
    print("*********Predicting precursor sequences*************")
    predictseq("features.csv")

    feature = os.path.join(path_tool, "features.csv")
    gdx_path = os.path.join(path_tool, "noncoding-seq.fasta.gdx")
    flat_path = os.path.join(path_tool, "noncoding-seq.fasta.flat")
    mirna_gdx = os.path.join(path_tool, "candidate-mirna.fasta.gdx")
    mirna_flat = os.path.join(path_tool, "candidate-mirna.fasta.flat")
    os.remove(blastid)
    os.remove(fastaid)
    os.remove(ncid)
    os.remove(output)
    os.remove(Temp)
    os.remove(result_file)
    os.remove(feature)
    os.remove(final_file)
    os.remove(blast_file)
    os.unlink(gdx_path)
    os.unlink(flat_path)
    os.unlink(mirna_flat)
    os.unlink(mirna_gdx)
    os.remove(mirid)

else:
    record = list(SeqIO.parse(final_file,"fasta"))
    my_file = open("Temporary.fasta", "w")
    for i in range(0,len(record)):
        nm = record[i].name
        seq = record[i].seq
        seq = str(seq)
        my_file.write(">")
        my_file.write(nm)
        my_file.write("\n")
        my_file.write(seq.replace("T", "U"))
        my_file.write("\n")
    my_file.close()
    Temp = os.path.join(path_tool, "Temporary.fasta")

    from feature import featurecount
    print("*********Calculating features of sequences*************")
    featurecount("Temporary.fasta")
    from Predict import predictseq
    print("*********Predicting precursor sequences*************")
    predictseq("features.csv")
    gdx_path = os.path.join(path_tool, "candidate-mirna.fasta.gdx")
    flat_path = os.path.join(path_tool, "candidate-mirna.fasta.flat")
    feature = os.path.join(path_tool, "features.csv")
    os.remove(output)
    os.remove(Temp)
    os.remove(result_file)
    os.remove(feature)
    os.remove(final_file)
    os.remove(blast_file)
    os.unlink(gdx_path)
    os.unlink(flat_path)
    os.remove(mirid)


