from bestSeeds import BestSeeds

def main():

    # A test of seed selection. Expecting [[qStart, rStart], ...] == [[1,2], [0,6], [0,7], [1,8], [0,11], [1,12], [0,15]]
    query = "ATCG"
    ref = "GGACGGATTCCATGGATA"
    print(BestSeeds(ref, query, 3, 1, 1, 1)) # only 1 mismatch allowed


if __name__ == "__main__":
    main()
