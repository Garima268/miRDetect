import os
import sys
import logging
import subprocess
from optparse import OptionParser
from typing import Set, List

from Bio import SeqIO, SearchIO
from pyfasta import Fasta

from config import path_db, path_blast, p_name

# ---------------------- Logging ---------------------- #
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ---------------------- CLI ---------------------- #
parser = OptionParser()
parser.add_option("-p", action="store", dest="num", type="int", default=1,
                  help="Number of CPUs")
options, args = parser.parse_args()

p = str(options.num)
path_tool = os.getcwd()

logging.info("********* Beginning miRNA Search *************")

# ---------------------- Input Handling ---------------------- #
if int(p) > 1:
    Assembly_file = sys.argv[3]
else:
    Assembly_file = sys.argv[1]

Assembly_file = os.path.abspath(Assembly_file)

# ---------------------- Helpers ---------------------- #
def run_command(cmd: List[str], cwd: str):
    """Run subprocess safely."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def ensure_blast_db(db_file: str, db_type: str, output_name: str):
    """Ensure BLAST database exists."""
    if not os.path.exists(db_file):
        logging.info("Creating BLAST database...")
        run_command(
            ["./makeblastdb", "-in", output_name,
             "-dbtype", db_type, "-parse_seqids", "-out", db_file.replace(".ndb", "")],
            path_blast
        )
    else:
        logging.info("Database exists")


def extract_ids_from_blast(blast_file: str, output_file: str):
    """Extract query IDs from BLAST output."""
    ids = set()
    with open(blast_file) as f:
        for record in SearchIO.parse(f, "blast-tab"):
            ids.add(record.id)

    with open(output_file, "w") as f:
        for i in ids:
            f.write(i + "\n")

    return ids


def filter_fasta_by_ids(fasta_in: str, ids: Set[str], fasta_out: str):
    """Filter FASTA sequences by ID."""
    with open(fasta_out, "w") as out:
        for seq in SeqIO.parse(fasta_in, "fasta"):
            if seq.id in ids:
                SeqIO.write(seq, out, "fasta")


def convert_to_rna(input_fasta: str, output_fasta: str):
    """Convert DNA to RNA (T -> U)."""
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq).replace("T", "U")
            out.write(f">{record.id}\n{seq}\n")


# ---------------------- Step 1: miRNA BLAST ---------------------- #
ensure_blast_db(
    os.path.join(path_db, "mature.ndb"),
    "nucl",
    os.path.join(path_db, "mature.fa")
)

blast_file = os.path.join(path_tool, "mirna-output.tsv")

run_command([
    "./blastn",
    "-query", Assembly_file,
    "-db", os.path.join(path_db, "mature"),
    "-evalue", "0.001",
    "-word_size", "6",
    "-outfmt", "6",
    "-max_target_seqs", "1",
    "-num_threads", p,
    "-out", blast_file
], path_blast)

if os.path.getsize(blast_file) == 0:
    logging.info("No miRNA found")
    sys.exit()

# ---------------------- Extract Candidate Sequences ---------------------- #
mirna_ids = extract_ids_from_blast(blast_file, "mirnablastid.txt")

candidate_fasta = "candidate-mirna.fasta"
filter_fasta_by_ids(Assembly_file, mirna_ids, candidate_fasta)

# ---------------------- Extract Precursors ---------------------- #
est_fasta = Fasta(candidate_fasta)
final_file = "Precursor-seq.fasta"

with open(final_file, "w") as out, open(blast_file) as bf:
    for line in bf:
        parts = line.split("\t")
        query = parts[0]
        qstart, qstop = [int(x) - 1 for x in parts[6:8]]

        xstream = 100
        up = max(0, qstart - xstream)
        down = qstop + xstream + 1

        seq = est_fasta[query][up:down]
        out.write(f">{query}\n{seq}\n")

# ---------------------- Protein Filtering (BLASTX) ---------------------- #
ensure_blast_db(
    os.path.join(path_db, "uniprot.pdb"),
    "prot",
    os.path.join(path_db, p_name)
)

output = "output.tsv"

run_command([
    "./blastx",
    "-query", final_file,
    "-db", os.path.join(path_db, "uniprot"),
    "-evalue", "0.001",
    "-outfmt", "6",
    "-max_target_seqs", "1",
    "-num_threads", p,
    "-out", output
], path_blast)

# ---------------------- Non-coding Path ---------------------- #
def run_prediction_pipeline(input_fasta: str):
    """Shared pipeline for prediction."""
    temp = "Temporary.fasta"
    convert_to_rna(input_fasta, temp)

    from feature import featurecount
    from Predict import predictseq

    logging.info("Calculating features...")
    featurecount(temp)

    logging.info("Running prediction...")
    predictseq("features.csv")

    return temp


if os.path.getsize(output) == 0:
    logging.info("No coding sequences found")
    temp = run_prediction_pipeline(final_file)

else:
    logging.info("Filtering non-coding sequences")

    blast_ids = extract_ids_from_blast(output, "blastid.txt")

    fasta_ids = {seq.id for seq in SeqIO.parse(final_file, "fasta")}
    noncoding_ids = fasta_ids - blast_ids

    with open("noncoding-id.txt", "w") as f:
        for i in noncoding_ids:
            f.write(i + "\n")

    nc_fasta = "noncoding-seq.fasta"
    filter_fasta_by_ids(final_file, noncoding_ids, nc_fasta)

    temp = run_prediction_pipeline(nc_fasta)

# ---------------------- Cleanup ---------------------- #
def safe_remove(file):
    if os.path.exists(file):
        os.remove(file)

for f in [
    blast_file, output, "features.csv", temp,
    "mirnablastid.txt", "blastid.txt",
    "noncoding-id.txt", "noncoding-seq.fasta"
]:
    safe_remove(f)

logging.info("Pipeline completed successfully.")
