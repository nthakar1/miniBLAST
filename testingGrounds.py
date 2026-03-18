from bestSeeds import BestSeeds
from main import miniBLAST

def main():

    """
    # Testing seed selection
    query = "ATCG"
    ref = "GGACGGATTCCATGGATA"
    print(BestSeeds(ref, query, 3, 1, 1, 1)) # only 1 mismatch allowed
    print(BestSeeds(ref, query, 2,1,1,0)) # smaller k
    print(BestSeeds(query, query, 2,1,1,0)) # query and ref are same
    print(BestSeeds(ref, query, 3,1,1,3)) # requires exact kmer matches 
    """

    query = "ATCG"
    ref = "GGACGGATTCCATGGATA"

    print(miniBLAST(ref, query, 3, 2, 3, 5, 2, 1, 20, 1))


if __name__ == "__main__":
    main()
