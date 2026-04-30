# valid.py
# Run with:
#   python valid.py
#
# Requirements:
#   python -m pip install biopython
#   makeblastdb -version
#   blastn -version
#
# Expected files in the same folder:
#   - H1N1_HA_single.fasta
#   - viral_HA_mixed_database.fasta

import os
import subprocess
from Bio import SeqIO

QUERY_FASTA = "H1N1_HA_single.fasta"
DATABASE_FASTA = "viral_HA_mixed_database.fasta"

# Set to True if you do NOT want the query sequence to match itself in the database.
EXCLUDE_SELF = True

# Output files
FILTERED_DB_FASTA = "viral_HA_mixed_database_no_self.fasta"
BLAST_DB_PREFIX = "viral_HA_mixed_database_no_self"
BLAST_OUT = "blast_validation_results.tsv"


def load_single_query(query_fasta):
    records = list(SeqIO.parse(query_fasta, "fasta"))

    if len(records) == 0:
        raise ValueError(f"No FASTA records found in {query_fasta}")
    if len(records) > 1:
        raise ValueError(
            f"{query_fasta} contains {len(records)} records, but this script expects exactly 1 query sequence."
        )

    query_record = records[0]
    print(f"Loaded query: {query_record.id}")
    print(f"Query length: {len(query_record.seq)}")
    return query_record


def write_database_without_self(database_fasta, query_id, output_fasta):
    db_records = list(SeqIO.parse(database_fasta, "fasta"))

    if len(db_records) == 0:
        raise ValueError(f"No FASTA records found in {database_fasta}")

    kept = [r for r in db_records if r.id != query_id]

    if len(kept) == len(db_records):
        print(f"Warning: query id {query_id} was not found in the database, so nothing was removed.")
    else:
        print(f"Removed self-hit entry {query_id} from database.")
        print(f"Database size: {len(db_records)} -> {len(kept)} records")

    SeqIO.write(kept, output_fasta, "fasta")
    return output_fasta


def make_blast_database(db_fasta, db_prefix):
    print(f"\nCreating BLAST database from {db_fasta} ...")
    subprocess.run([
        "makeblastdb",
        "-in", db_fasta,
        "-dbtype", "nucl",
        "-out", db_prefix
    ], check=True)
    print(f"BLAST database created with prefix: {db_prefix}")


def run_blastn(query_fasta, db_prefix, output_file):
    print(f"\nRunning BLASTn with Mini-BLAST parameters...")
    
    # Mapping datatypes.py BLASTN_PARAMS to NCBI command line flags
    blast_args = [
        "blastn",
        "-task", "blastn",               # Forces standard blastn instead of megablast
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", output_file,
        "-word_size", "7",               # From k_mer_size=7
        "-reward", "2",                  # From match_reward=2
        "-penalty", "-3",                # From mismatch_penalty=3 (must be negative in CLI)
        "-gapopen", "5",                 # From gap_opening=5
        "-gapextend", "2",               # From gap_extension=2
        "-xdrop_ungap", "20",            # From xdrop_ungap=20
        "-xdrop_gap", "30",              # From xdrop_gap=30
        "-xdrop_gap_final", "100",       # From xdrop_gap_final=100
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]

    subprocess.run(blast_args, check=True)
    print(f"BLAST results written to: {output_file}")


def print_top_hits(output_file, n=10):
    print("\nTop BLAST hits:")

    if not os.path.exists(output_file):
        print("No BLAST output file found.")
        return

    with open(output_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        print("No hits found.")
        return

    header = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    print("\t".join(header))

    for line in lines[:n]:
        print(line)


def main():
    if not os.path.exists(QUERY_FASTA):
        raise FileNotFoundError(f"Missing query FASTA: {QUERY_FASTA}")

    if not os.path.exists(DATABASE_FASTA):
        raise FileNotFoundError(f"Missing database FASTA: {DATABASE_FASTA}")

    query_record = load_single_query(QUERY_FASTA)

    db_fasta_to_use = DATABASE_FASTA
    db_prefix_to_use = os.path.splitext(DATABASE_FASTA)[0]

    if EXCLUDE_SELF:
        db_fasta_to_use = write_database_without_self(
            database_fasta=DATABASE_FASTA,
            query_id=query_record.id,
            output_fasta=FILTERED_DB_FASTA
        )
        db_prefix_to_use = BLAST_DB_PREFIX

    make_blast_database(
        db_fasta=db_fasta_to_use,
        db_prefix=db_prefix_to_use
    )

    run_blastn(
        query_fasta=QUERY_FASTA,
        db_prefix=db_prefix_to_use,
        output_file=BLAST_OUT
    )

    print_top_hits(BLAST_OUT, n=10)


if __name__ == "__main__":
    main()
