from Bio import Entrez, SeqIO
import time
from itertools import islice

Entrez.email = "mmccrear@nadrew.cmu.edu"


def fetch_mixed_database(queries, output_filename):
    """Fetch H1N1 + diverse viral sequences for meaningful E-values."""
    
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
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        all_records.extend(records)
        time.sleep(0.34)  # Respect NCBI rate limits
    
    return all_records


def main():

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
    mixed_db_records = fetch_mixed_database(queriesViral, "h1n1_mixed_database.fasta")

if __name__ == "__main__":
    main()