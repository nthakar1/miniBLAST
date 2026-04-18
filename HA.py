"""
Installation requirements: biopython, blast

Instructions for use:

To generate a database, first create a list of queries, formatted as follows:
    queries = [ 
    { "query" : "example_search_term[example_field_for_term] AND/OR other_terms[other_fields]",
      "max_results" : maximum sequence calls for this query you want to include in db,
      "label" : "whatever you want to name this query for progress tracking"
    }, ...
    ]

Then call fetch_mixed_database with the following arguments:
    queries : the list above
    dbType : "nucleotide" or "protein"
    output_filename : "example.fasta" or None if you don't want to save to a file
"""

import os
import subprocess
import time  # <-- Imported the time module for the pause
from Bio import Entrez, SeqIO

# -----------------------------
# CONFIG
# -----------------------------
Entrez.email = "nthakar@andrew.cmu.edu"   # REQUIRED by NCBI
# Entrez.api_key = "YOUR_NCBI_API_KEY"    # optional but recommended

# -----------------------------
# BLAST+ HELPERS
# -----------------------------
def make_blast_database(db_fasta, dbType="nucl"):
    """
    Convert FASTA into a local BLAST database.
    dbType must be 'nucl' or 'prot'
    """
    out_name = os.path.splitext(db_fasta)[0]

    print(f"\nCreating BLAST database from {db_fasta} ...")
    subprocess.run([
        "makeblastdb",
        "-in", db_fasta,
        "-dbtype", dbType,
        "-out", out_name
    ], check=True)

    print(f"BLAST database created with prefix: {out_name}")
    return out_name


# -----------------------------
# FETCH FUNCTIONS
# -----------------------------
def fetch_sequences_batch(query, db_type="nucleotide", max_results=20, label="query"):
    """
    Search NCBI and fetch FASTA records for one query.
    Returns a list of SeqRecord objects.
    """
    print(f"\nSearching NCBI for: {label}")

    search_handle = Entrez.esearch(db=db_type, term=query, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]
    print(f"Found {len(ids)} records for {label}")

    if not ids:
        return []

    fetch_handle = Entrez.efetch(
        db=db_type,
        id=",".join(ids),
        rettype="fasta",
        retmode="text"
    )
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    fetch_handle.close()

    # Add label to description for easier tracking
    for record in records:
        record.description = f"{label}|{record.description}"

    return records


def fetch_mixed_database(queries, dbType="nucleotide", output_filename=None):
    """
    Build a mixed FASTA database from multiple Entrez queries.
    """
    all_records = []
    seen_ids = set()

    for q in queries:
        query = q["query"]
        max_results = q.get("max_results", 20)
        label = q.get("label", "query")

        # --> Added a 1-second pause before each query to respect NCBI limits
        time.sleep(1)

        records = fetch_sequences_batch(
            query=query,
            db_type=dbType,
            max_results=max_results,
            label=label
        )

        for record in records:
            # de-duplicate by accession ID
            if record.id not in seen_ids:
                seen_ids.add(record.id)
                all_records.append(record)

    print(f"\nTotal unique records collected: {len(all_records)}")

    if output_filename:
        SeqIO.write(all_records, output_filename, "fasta")
        print(f"Saved mixed database FASTA to: {output_filename}")

    return all_records


def blast_validation(queryFile, dbFile, segment_id, dbType="nucl"):
    """
    Run BLAST+ using a query FASTA against a local FASTA database.

    queryFile: a FASTA/FNA file
    dbFile: a FASTA database file
    segment_id: sequence ID from queryFile to extract and test
    dbType: 'nucl' or 'prot'
    """
    # Build local BLAST database
    db_prefix = make_blast_database(dbFile, dbType=dbType)

    # Extract requested query sequence into temporary file
    query_records = list(SeqIO.parse(queryFile, "fasta"))
    selected = [r for r in query_records if r.id == segment_id]

    if not selected:
        available_ids = [r.id for r in query_records[:10]]
        raise ValueError(
            f"Segment ID '{segment_id}' not found in {queryFile}. "
            f"Example available IDs: {available_ids}"
        )

    temp_query = f"{segment_id}_temp_query.fasta"
    SeqIO.write(selected, temp_query, "fasta")

    output_file = f"{segment_id}_vs_{os.path.basename(db_prefix)}_blast.tsv"

    if dbType == "nucl":
        cmd = [
            "blastn",
            "-query", temp_query,
            "-db", db_prefix,
            "-out", output_file,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ]
    else:
        cmd = [
            "blastp",
            "-query", temp_query,
            "-db", db_prefix,
            "-out", output_file,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ]

    print(f"\nRunning BLAST validation for {segment_id} ...")
    subprocess.run(cmd, check=True)

    print(f"BLAST results written to: {output_file}")

    with open(output_file, "r") as f:
        print("\nTop BLAST hits:")
        for i, line in enumerate(f):
            print(line.strip())
            if i >= 9:
                break


# -----------------------------
# EXAMPLE USAGE
# -----------------------------
if __name__ == "__main__":
    
    # 1) Formatted list of queries per the usage instructions
    queries = [
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H1N1[All Fields]',
            "max_results": 25,
            "label": "H1N1_HA"
        },
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H3N2[All Fields]',
            "max_results": 25,
            "label": "H3N2_HA"
        },
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H5N1[All Fields]',
            "max_results": 25,
            "label": "H5N1_HA"
        },
        {
            "query": 'Influenza B virus[Organism] AND hemagglutinin[Title]',
            "max_results": 15,
            "label": "InfluenzaB_HA"
        }
    ]

    # 2) Call fetch_mixed_database using the defined arguments
    fetch_mixed_database(
        queries=queries,
        dbType="nucleotide",
        output_filename="viral_HA_mixed_database.fasta"
    )

    # 3) Run BLAST validation (Replace with actual query sequence and file if running)
    # blast_validation(
    #     queryFile="your_query_sequences.fasta",
    #     dbFile="viral_HA_mixed_database.fasta",
    #     segment_id="NC_002017.1",
    #     dbType="nucl"
    # )