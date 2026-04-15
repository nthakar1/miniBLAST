from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
import os, subprocess, time, tempfile

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

To run BLAST+ (command-based NCBI BLAST) on custome db, call blast_validation with the following arguments and view output in terminal:
    queryFile: a .fasta or .fna file
    dbFile: a .fasta file
    segment_id: for queryFile containing multiple sequences, give the string following and ">" entry up until the first space 
                (e.g., ">NC_002023.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 1, complete sequence" becomes "NC.002023.1")
    dbType: "nucl" or "prot"

"""

Entrez.email = "mmccrear@andrew.cmu.edu"

def fetch_mixed_database(queries, dbType, output_filename):
    
    if dbType not in ("nucleotide", "protein"):
        raise ValueError(f"dbType must be 'nucleotide' or 'protein', got {dbType!r}")

    all_records = []
    
    for query_info in queries:
        print(f"\nFetching {query_info['label']}...")
        records = fetch_sequences_batch(
            query_info["query"], dbType,
            max_results=query_info["max_results"]
        )
        
        # Add label to sequence descriptions for tracking
        for record in records:
            record.description = f"{query_info['label']}|{record.description}"
        
        all_records.extend(records)
        print(f"  Added {len(records)} sequences")

    # Save combined database:
    if output_filename != None:
        #SeqIO.write(all_records, output_filename, "fasta")
        #print(f"\nSaved {len(all_records)} total sequences to {output_filename}")
        with open(output_filename, "w") as f:
            for record in all_records:
                f.write(record.format("fasta"))

    # Print composition summary
    labels = {}
    for record in all_records:
        label = record.description.split('|')[0]
        labels[label] = labels.get(label, 0) + 1
    
    print("\nDatabase composition:")
    for label, count in labels.items():
        print(f"  {label}: {count} sequences")
    
    return all_records

def fetch_sequences_batch(query, dbType, max_results):
    """Fetch sequences in batches to respect rate limits."""
    
    # Search for matching IDs
    handle = Entrez.esearch(db=dbType, term=query, retmax=max_results)
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
        handle = Entrez.efetch(db=dbType, id=batch, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta-pearson"))
        handle.close()
        all_records.extend(records)
        time.sleep(0.34)  # Respect NCBI rate limits
    
    return all_records

def blast_validation(queryFile, dbFile, segment_id, dbType):

    if dbType not in ("nucl", "prot"):
        raise ValueError(f"dbType must be 'nucl' or 'prot', got {dbType!r}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Build DB in temp dir
        db_tmp = os.path.join(tmpdir, "blast_db")
        cmd = ["makeblastdb", "-in", dbFile, "-dbtype", dbType, "-out", db_tmp]
        db_result = subprocess.run(cmd, capture_output=True, text=True)
        print(f"makeblastdb stdout: {db_result.stdout}")
        print(f"makeblastdb stderr: {db_result.stderr}")

        # Write query to temp dir
        query_tmp = os.path.join(tmpdir, "query.fasta")
        record = next(r for r in SeqIO.parse(queryFile, "fasta") if r.id == segment_id)
        if record is None:
            raise ValueError(f"No record with id {segment_id!r} found in {queryFile}")
        with open(query_tmp, "w") as out_handle:
            SeqIO.write(record, out_handle, "fasta")

        # Run BLAST, output to temp dir
        blast_out = os.path.join(tmpdir, "blast_out.xml")
        if dbType == "nucl":
            blast_cmd = ["blastn", "-query", query_tmp, "-db", db_tmp, "-outfmt", "5", "-out", blast_out]
        elif dbType == "prot":
            blast_cmd = ["blastp", "-query", query_tmp, "-db", db_tmp, "-outfmt", "5", "-out", blast_out]
        result = subprocess.run(blast_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BLAST failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")

        # Parse before tempdir is deleted
        with open(blast_out) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                print(f"Alignments found: {len(record.alignments)}")
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        print(f"Hit: {alignment.title[:60]}")
                        print(f"  E-value: {hsp.expect}")
                        print(f"  Identity: {hsp.identities}/{hsp.align_length} ({100*hsp.identities//hsp.align_length}%)")
                        print(f"  Query: {hsp.query_start}-{hsp.query_end}")
    
def main():
    """
    ### CREATE DIFFERENT DATABASES HERE TO SAVE AS FILES, IMPORT IN MAIN.PY
    queriesH1N1 = [
        # H1N1 sequences (target of interest)
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
    mixed_h1n1_db_records = fetch_mixed_database(queriesH1N1, "nucleotide","h1n1_mixed_database.fasta")

    queriesRecA = [
        # E. coli RecA sequences (target of interest)
        {
            "query": "RecA[Protein Name] AND Escherichia coli[Organism]",
            "max_results": 200,
            "label": "EcoliRecA_main"
        },
        {
            "query": "RecA protein[Title] AND E. coli[Organism]", 
            "max_results": 200,
            "label": "EcoliRecA_alt"
        },
        
        # Other bacterial RecA homologs (closely related background)
        {
            "query": "RecA[Protein Name] AND Salmonella[Organism]",
            "max_results": 150,
            "label": "Salmonella_RecA"
        },
        {
            "query": "RecA[Protein Name] AND Enterobacteriaceae[Organism] NOT Escherichia coli[Organism]",
            "max_results": 100,
            "label": "Enterobacteria_RecA"
        },
        
        # Other DNA repair proteins (distant background for E-value diversity)
        {
            "query": "RadA[Protein Name] AND Archaea[Organism]",
            "max_results": 100,
            "label": "RadA_archaea"
        },
        {
            "query": "Rad51[Protein Name] AND eukaryotes[Organism]",
            "max_results": 100,
            "label": "Rad51_eukaryotes"
        },
        
        # Some unrelated proteins (adds statistical noise)
        {
            "query": "ribosomal protein L1[Protein Name] AND bacteria[Organism]",
            "max_results": 150,
            "label": "ribosomal_L1"
        }
    ]

    # Build the database
    mixed_recA_db_records = fetch_mixed_database(queriesRecA,"protein", "RecA_mixed_database.fasta")
    """

    ### BLAST VALIDATION ON CUSTOM DATABASES
    query = "Ecoli_RecA.fasta"
    segment = "sp|P0A7G6|RECA_ECOLI"
    db = "RecA_mixed_database.fasta"
    #blast_validation(query, db, segment, "nucl")
    blast_validation(query, db, segment, "prot")
    

if __name__ == "__main__":
    main()
