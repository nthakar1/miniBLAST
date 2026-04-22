'''
from Bio import Entrez, SeqIO

# -----------------------------
# CONFIG
# -----------------------------
# Required by NCBI so they can contact you if your script causes server issues
Entrez.email = "nthakar@andrew.cmu.edu"

def fetch_h1n1_ha(output_filename="H1N1_HA_only.fasta", max_results=25):
    """
    Searches NCBI specifically for H1N1 HA sequences and saves them to a FASTA file.
    """
    # Your exact query for H1N1 Hemagglutinin
    query = 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H1N1[All Fields]'
    
    print("Searching NCBI for H1N1 HA sequences...")
    
    # Step 1: Search for the IDs
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]
    print(f"Found {len(ids)} records. Fetching FASTA data...")

    if not ids:
        print("No records found. Exiting.")
        return

    # Step 2: Fetch the actual sequence data using those IDs
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(ids),
        rettype="fasta",
        retmode="text"
    )
    
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    fetch_handle.close()

    # Step 3: Save directly to the new file
    SeqIO.write(records, output_filename, "fasta")
    print(f"Success! Saved {len(records)} sequences straight to '{output_filename}'.")

# -----------------------------
# RUN IT
# -----------------------------
if __name__ == "__main__":
    # You can change the filename or max_results here if you want more than 25
    fetch_h1n1_ha(output_filename="H1N1_HA_only.fasta", max_results=25)
    '''

from Bio import Entrez, SeqIO

# -----------------------------
# CONFIG
# -----------------------------
Entrez.email = "nthakar@andrew.cmu.edu"

def fetch_single_h1n1_ha(output_filename="H1N1_HA_single.fasta"):
    """
    Searches NCBI specifically for H1N1 HA sequences and saves EXACTLY ONE to a FASTA file.
    """
    # Your exact query for H1N1 Hemagglutinin
    query = 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H1N1[All Fields]'
    
    print("Searching NCBI for a single H1N1 HA sequence...")
    
    # Step 1: Search, but limit the results to exactly 1 (retmax=1)
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]

    if not ids:
        print("No records found. Exiting.")
        return

    print(f"Found a record (ID: {ids[0]}). Fetching FASTA data...")

    # Step 2: Fetch the actual sequence data using that 1 ID
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=ids[0],
        rettype="fasta",
        retmode="text"
    )
    
    # Use SeqIO.read instead of parse because we expect exactly 1 record
    record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()

    # Step 3: Save directly to the new file
    SeqIO.write(record, output_filename, "fasta")
    print(f"Success! Saved 1 sequence straight to '{output_filename}'.")

# -----------------------------
# RUN IT
# -----------------------------
if __name__ == "__main__":
    fetch_single_h1n1_ha()