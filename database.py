from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
import os, subprocess, time, tempfile

# Install separately: biopython, blast

Entrez.email = "mmccrear@andrew.cmu.edu"

def fetch_mixed_database(queries, output_filename):
    
    all_records = []
    
    for query_info in queries:
        print(f"\nFetching {query_info['label']}...")
        records = fetch_sequences_batch(
            query_info["query"], 
            max_results=query_info["max_results"]
        )
        
        # Add label to sequence descriptions for tracking
        for record in records:
            record.description = f"{query_info['label']}|{record.description}"
        
        all_records.extend(records)
        print(f"  Added {len(records)} sequences")
    
    # Save combined database:
    # SeqIO.write(all_records, output_filename, "fasta")
    # print(f"\nSaved {len(all_records)} total sequences to {output_filename}")
    
    # Print composition summary
    labels = {}
    for record in all_records:
        label = record.description.split('|')[0]
        labels[label] = labels.get(label, 0) + 1
    
    print("\nDatabase composition:")
    for label, count in labels.items():
        print(f"  {label}: {count} sequences")
    
    return all_records

def fetch_sequences_batch(query, db="nucleotide", max_results=500):
    """Fetch sequences in batches to respect rate limits."""
    
    # Search for matching IDs
    handle = Entrez.esearch(db=db, term=query, retmax=max_results)
    results = Entrez.read(handle)
    handle.close()
    ids = results["IdList"]
    
    if not ids:
        print(f"  No sequences found for query: {query}")
        return []
    
    # Fetch in batches
    batch_size = 50
    all_records = []
    
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        handle = Entrez.efetch(db=db, id=batch, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta-pearson"))
        handle.close()
        all_records.extend(records)
        time.sleep(0.34)  # Respect NCBI rate limits
    
    return all_records

def blastn_validation(queryFile, segment_id, dbFile):

    with tempfile.TemporaryDirectory() as tmpdir:
        # Build DB in temp dir
        db_tmp = os.path.join(tmpdir, "blast_db")
        cmd = ["makeblastdb", "-in", dbFile, "-dbtype", "nucl", "-out", db_tmp]
        subprocess.run(cmd, capture_output=True, text=True)

        # Write query to temp dir
        query_tmp = os.path.join(tmpdir, "query.fasta")
        record = next(r for r in SeqIO.parse(queryFile, "fasta") if r.id == segment_id)
        with open(query_tmp, "w") as out_handle:
            SeqIO.write(record, out_handle, "fasta")

        # Run BLASTN, output to temp dir
        blast_out = os.path.join(tmpdir, "blast_out.xml")
        blastn_cmd = ["blastn", "-query", query_tmp, "-db", db_tmp, "-outfmt", "5", "-out", blast_out]
        result = subprocess.run(blastn_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BLASTN failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")

        # Parse before tempdir is deleted
        with open(blast_out) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        print(f"Hit: {alignment.title[:60]}")
                        print(f"  E-value: {hsp.expect}")
                        print(f"  Identity: {hsp.identities}/{hsp.align_length} ({100*hsp.identities//hsp.align_length}%)")
                        print(f"  Query: {hsp.query_start}-{hsp.query_end}")
    
def main():

    ### CREATE DIFFERENT DATABASES HERE TO SAVE AS FILES, IMPORT IN MAIN.PY
    """
    queriesViral = [
        # H1N1 sequences (your target of interest)
        {
            "query": "H1N1[All Fields] AND influenza A virus[Organism]",
            "max_results": 200,
            "label": "H1N1_HA"
        },
        {
            "query": "H1N1[All Fields] AND influenza A virus[Organism]", 
            "max_results": 200,
            "label": "H1N1_NA"
        },
        
        # Other influenza subtypes (closely related background)
        {
            "query": "H3N2[All Fields] AND influenza A virus[Organism]",
            "max_results": 150,
            "label": "H3N2_HA"
        },
        {
            "query": "H5N1[All Fields] AND influenza A virus[Organism]",
            "max_results": 100,
            "label": "H5N1"
        },
        
        # Other respiratory viruses (distant background for E-value diversity)
        {
            "query": "SARS-CoV-2[Organism] AND complete genome[Title]",
            "max_results": 100,
            "label": "COVID"
        },
        {
            "query": "rhinovirus[Organism] AND complete genome[Title]",
            "max_results": 100,
            "label": "rhinovirus"
        },
        
        # Some bacterial background (adds statistical noise)
        {
            "query": "16S ribosomal RNA[Gene] AND bacteria[Organism]",
            "max_results": 150,
            "label": "bacteria_16S"
        }
    ]

    # Build the database
    mixed_h1n1_db_records = fetch_mixed_database(queriesViral, "h1n1_mixed_database.fasta")
    """

    ### BLAST VALIDATION ON CUSTOM DATABASES
    query = "human_h1n1_1934.fna"
    segment = "NC_002023.1"
    db = "h1n1_mixed_database.fasta"
    blastn_validation(query, segment, db)

if __name__ == "__main__":
    main()