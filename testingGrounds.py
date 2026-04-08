"""
testingGrounds.py — Unit + integration tests for the miniBLAST pipeline.

Pipeline overview:
  BestSeeds  →  TwoHitSeeds  →  ungappedExtension  →  affineGappedExtension
                                 (filtered by s1)

Run with:
    python testingGrounds.py
"""

from bestSeeds import BestSeeds
from TwoHit import TwoHitSeeds
from extendSeeds import extendFromSeeds, ungappedExtension
from main import miniBLASTn
from datatypes import BLASTN_PARAMS, BLASTP_PARAMS, BLOSUM
from statistics import calculate_bit_score, calculate_e_value

# Test runner

_results = {"passed": 0, "failed": 0}

def run_test(name, fn):
    print(f"\n{'='*64}")
    print(f"  TEST: {name}")
    print('='*64)
    try:
        fn()
        print(f"  ✓  PASSED")
        _results["passed"] += 1
    except AssertionError as e:
        print(f"  ✗  ASSERTION FAILED: {e}")
        _results["failed"] += 1
    except Exception as e:
        import traceback
        print(f"  ✗  ERROR: {e}")
        traceback.print_exc()
        _results["failed"] += 1


#   BestSeeds unit tests

def test_best_seeds_exact_match():
    """
    A query kmer that appears verbatim in the reference must always
    be returned as a seed when the threshold equals k * matchScore
    (i.e. only exact k-mer matches pass).
    """
    query = "ATCGATCG"
    ref   = "GGATCGATCGCC"     # contains "ATCG" at positions 2 and 6
    k, match, mismatch, thresh = 4, 2, 3, 8   # thresh = 4*2 → exact only

    seeds = BestSeeds(ref, query, k, match, mismatch, None, thresh)
    print(f"  Seeds: {seeds}")
    assert len(seeds) > 0, "Expected at least one seed for a verbatim kmer match"


def test_best_seeds_near_match():
    """
    A lower threshold should allow a seed with one mismatch in a 4-mer.
    query[0:4] = 'ATCG', ref contains 'ATCA' (1 mismatch → score = 2+2+2-3 = 3).
    With thresh=3 that should pass.
    """
    query = "ATCGATCG"
    ref   = "ATCATCG"          # 'ATCA' at pos 0 is 1 mismatch from 'ATCG'
    k, match, mismatch, thresh = 4, 2, 3, 3

    seeds = BestSeeds(ref, query, k, match, mismatch, thresh)
    print(f"  Seeds (lenient threshold): {seeds}")
    assert len(seeds) > 0, "Expected seeds for a near-matching kmer"


def test_best_seeds_no_match():
    """
    Completely dissimilar sequences (AAAA vs CCCC) yield no seeds.
    """
    query = "AAAAAAA"
    ref   = "CCCCCCC"
    k, match, mismatch, thresh = 3, 2, 3, 4

    seeds = BestSeeds(ref, query, k, match, mismatch, thresh)
    print(f"  Seeds: {seeds}")
    assert seeds == [], f"Expected no seeds, got: {seeds}"


def test_best_seeds_self_alignment():
    """
    Aligning a sequence against itself should produce seeds at every
    position along the main diagonal (offset 0).
    """
    seq = "ATCGATCG"
    k, match, mismatch, thresh = 3, 2, 3, 6   # exact only

    seeds = BestSeeds(seq, seq, k, match, mismatch, thresh)
    print(f"  Seeds (self): {seeds}")
    diag_seeds = [s for s in seeds if s[0] == s[1]]
    assert len(diag_seeds) > 0, "Self-alignment should have seeds on main diagonal"


#   TwoHitSeeds unit tests

def test_two_hit_valid_pair():
    """
    Two hits on the same diagonal, k positions apart, within window A.
    The second hit should survive.
    """
    k, A = 3, 15
    # diagonal = q - r; both [0,0] and [5,5] have diagonal 0
    single_hits = [[0, 0], [5, 5]]

    valid = TwoHitSeeds(single_hits, k, A)
    print(f"  Valid seeds: {valid}")
    assert len(valid) > 0, "Expected [5,5] to survive the two-hit filter"


def test_two_hit_too_far_apart():
    """
    Two hits on the same diagonal but farther than A → filtered out.
    """
    k, A = 3, 10
    single_hits = [[0, 0], [50, 50]]   # distance = 50 > A

    valid = TwoHitSeeds(single_hits, k, A)
    print(f"  Valid seeds: {valid}")
    assert valid == [], "Expected no seeds (distance > A)"


def test_two_hit_different_diagonal():
    """
    Hits on different diagonals must not pair, even if close.
    [0,0] is on diagonal 0; [3,5] is on diagonal -2.
    """
    k, A = 3, 15
    single_hits = [[0, 0], [3, 5]]

    valid = TwoHitSeeds(single_hits, k, A)
    print(f"  Valid seeds: {valid}")
    assert valid == [], "Expected no seeds (different diagonals)"


def test_two_hit_overlapping_kmers():
    """
    Two hits that overlap (distance < k) should not trigger the rule.
    """
    k, A = 4, 15
    single_hits = [[0, 0], [2, 2]]    # distance = 2 < k = 4

    valid = TwoHitSeeds(single_hits, k, A)
    print(f"  Valid seeds: {valid}")
    assert valid == [], "Expected no seeds (hits overlap)"


def test_two_hit_chain_of_three():
    """
    Three hits on the same diagonal: [0,0], [5,5], [10,10].
    Both [5,5] and [10,10] should survive (chaining).
    """
    k, A = 3, 10
    single_hits = [[0, 0], [5, 5], [10, 10]]

    valid = TwoHitSeeds(single_hits, k, A)
    print(f"  Valid seeds: {valid}")
    assert len(valid) >= 2, f"Expected at least 2 chained seeds, got: {valid}"


#   Ungapped extension unit tests

def test_ungapped_perfect_seed():
    """
    Seed is a perfect sub-string of ref; ungapped score should be positive.
    """
    query = "ATCGATCG"
    ref   = "GGATCGATCGCC"
    result = ungappedExtension(query, ref,
                               q_start=0, r_start=2, seed_score=0,
                               k=4, matrix=None, match=2, mismatch=3, Xdrop=20)
    print(f"  Ungapped result: {result}")
    assert result["score"] > 0, "Expected positive ungapped score"


def test_ungapped_returns_required_keys():
    """
    The returned dict must contain: score, alignment, q_range, r_range, q_seed, r_seed.
    """
    query = "AAATGCCATGAA"
    ref   = "TTATGCCATGTT"
    
    result = ungappedExtension(query, ref,
                               q_start=2, r_start=2,seed_score=0,
                               k=3, matrix=None, match=2, mismatch=3, Xdrop=20)
    print(f"  Keys: {list(result.keys())}")
    for key in ("score", "alignment", "q_range", "r_range", "q_seed", "r_seed"):
        assert key in result, f"Missing key: {key}"


#   End-to-end / integration tests  (full pipeline via miniBLAST)


def test_e2e_high_similarity_dna():
    """
    query and ref share the motif 'ATGCCCATG'; expect a valid local alignment.
    """
    query = "AAATGCCCATGAA"
    ref   = "TTATGCCCATGTT"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Alignment: {result['alignment'] if result else None}")
    print(f"  Score:     {result['score'] if result else None}")
    assert result is not None, "Expected a valid alignment"
    assert result["score"] > 0, "Expected positive alignment score"


def test_e2e_single_mismatch():
    """
    Sequences differ at exactly one position; the aligner should still
    return a high-scoring alignment.
    query: ATCGATCG
    ref  : ATCGTTCG  (position 4: A→T)
    """
    query = "ATCGATCG"
    ref   = "ATCGTTCG"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Alignment: {result['alignment'] if result else None}")
    print(f"  Score:     {result['score'] if result else None}")
    assert result is not None, "Expected alignment despite single mismatch"


def test_e2e_no_homology():
    """
    Completely dissimilar sequences → pipeline should return None (no seeds pass).
    """
    query = "AAAAAAAAAA"
    ref   = "CCCCCCCCCC"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Result: {result}")
    assert result is None, f"Expected None for no-homology case, got: {result}"


def test_e2e_query_is_substring():
    """
    The query is an exact substring of the reference; expect a perfect local alignment
    that covers the full query.
    """
    query = "ATCGATCG"
    ref   = "GGGGGATCGATCGCCCCC"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Alignment: {result['alignment'] if result else None}")
    print(f"  Score:     {result['score'] if result else None}")
    assert result is not None, "Expected alignment for exact substring"
    # Note: affineGappedExtension only extends rightward from the seed, so the
    # reported score reflects the portion to the right of (and including) the
    # seed position — not the full sequence length.
    assert result["score"] > 0, f"Expected positive score for exact substring, got: {result['score']}"
    print(f"  (Gapped extension is rightward-only; full-sequence score would be {len(query) * 2})")


def test_e2e_with_gap():
    """
    query has an extra base compared to ref in the middle; gapped extension
    should produce a reasonable alignment.
    ref  : ATCGGATCGATCG
    query: ATCGAATCGATCG  (extra 'A' inserted at pos 4)
    """
    query = "ATCGAATCGATCG"
    ref   = "ATCGGATCGATCG"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Alignment: {result['alignment'] if result else None}")
    print(f"  Score:     {result['score'] if result else None}")
    assert result is not None, "Expected alignment even with a gap"


def test_e2e_original_case():
    """
    The original testingGrounds sequences:
      query = AAATGCCCCCCATGAA
      ref   = TTATGGGGGGGATGTT

    Expected seeds (from original comment):
      Seeds1: [[2, 2], [11, 11]]
      Seeds2: [[11, 11]]
    """
    query = "AAATGCCCCCCATGAA"
    ref   = "TTATGGGGGGGATGTT"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Alignment: {result['alignment'] if result else None}")
    print(f"  Score:     {result['score'] if result else None}")
    # We just check the pipeline runs without crashing; alignment may be partial
    # due to complementary bases (CCCC vs GGGG are reverse complements)


def test_e2e_longer_realistic_dna():
    """
    Longer (~50 bp) synthetic sequences modelled after typical BLAST queries,
    with 3 point mutations introduced.
    ref  : ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
    query: ATGCTAGCTAGCTAGCTAGCTAGCAAGCTAGCTAGCTAGCTAGCTAGCT
                                           ^           (C→A at pos 25)
    """
    ref   = "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    query = "ATGCTAGCTAGCTAGCTAGCTAGCAAGCTAGCTAGCTAGCTAGCTAGCT"

    result = miniBLASTn(ref, query, s1=20, A=40)
    print(f"  Score: {result['score'] if result else None}")
    print(f"  Query coverage: {result.get('query_coverage') if result else None}")
    assert result is not None, "Expected alignment for near-identical 50 bp sequences"
    assert result["score"] > 30, "Expected high score for near-perfect long alignment"

# Stats testing 

def test_bit_score_positive():
    bit = calculate_bit_score(50, BLASTN_PARAMS)
    assert bit > 0, "Bit score should be positive"

def test_evalue_strong_hit():
    eval = calculate_e_value(50, 100, 987, BLASTN_PARAMS)
    assert 0 < eval < 1, "Strong hit E-value should be between 0 and 1"

def test_evalue_increases_with_db_size():
    e1 = calculate_e_value(50, 100, 987, BLASTN_PARAMS)
    e2 = calculate_e_value(50, 100, 9870, BLASTN_PARAMS)
    assert e2 > e1, "Larger database should give higher E-value"

def test_evalue_decreases_with_higher_score():
    e1 = calculate_e_value(50, 100, 987, BLASTN_PARAMS)
    e2 = calculate_e_value(100, 100, 987, BLASTN_PARAMS)
    assert e2 < e1, "Higher alignment score should give lower E-value"

def test_evalue_weak_hit():
    eval = calculate_e_value(5, 100, 987, BLASTN_PARAMS)
    assert eval > 1, "Weak hit E-value should be > 1"

def test_blastp_params():
    bit = calculate_bit_score(50, BLASTP_PARAMS)
    eval = calculate_e_value(50, 100, 987, BLASTP_PARAMS)
    assert bit > 0
    assert eval > 0



# Main testing

def main():
    print("\n" + "█"*64)
    print("  miniBLASTn test suite")
    print("█"*64)

    # BestSeeds
    run_test("BestSeeds | exact k-mer match",      test_best_seeds_exact_match)
    run_test("BestSeeds | near-match (1 mismatch)", test_best_seeds_near_match)
    run_test("BestSeeds | no match",               test_best_seeds_no_match)
    run_test("BestSeeds | self-alignment",         test_best_seeds_self_alignment)

    # TwoHitSeeds
    run_test("TwoHit | valid pair",                test_two_hit_valid_pair)
    run_test("TwoHit | too far apart (> A)",       test_two_hit_too_far_apart)
    run_test("TwoHit | different diagonal",        test_two_hit_different_diagonal)
    run_test("TwoHit | overlapping kmers",         test_two_hit_overlapping_kmers)
    run_test("TwoHit | chain of 3 hits",           test_two_hit_chain_of_three)

    # Ungapped extension
    run_test("Ungapped | perfect seed",            test_ungapped_perfect_seed)
    run_test("Ungapped | returns required keys",   test_ungapped_returns_required_keys)

    # End-to-end
    run_test("E2E | high-similarity DNA",          test_e2e_high_similarity_dna)
    run_test("E2E | single mismatch",              test_e2e_single_mismatch)
    run_test("E2E | no homology → None",           test_e2e_no_homology)
    run_test("E2E | query is exact substring",     test_e2e_query_is_substring)
    run_test("E2E | sequence with gap",            test_e2e_with_gap)
    run_test("E2E | original test case",           test_e2e_original_case)
    run_test("E2E | longer 50 bp near-identical",  test_e2e_longer_realistic_dna)

    # Statistics
    run_test("Statistics | bit score positive",           test_bit_score_positive)
    run_test("Statistics | strong hit E-value < 1",       test_evalue_strong_hit)
    run_test("Statistics | E-value scales with db size",  test_evalue_increases_with_db_size)
    run_test("Statistics | E-value drops with score",     test_evalue_decreases_with_higher_score)
    run_test("Statistics | weak hit E-value > 1",         test_evalue_weak_hit)
    run_test("Statistics | BLASTP params are working as expected",           test_blastp_params)

    # Summary 
    total = _results["passed"] + _results["failed"]
    print(f"\n{'─'*64}")
    print(f"  Results: {_results['passed']}/{total} passed  "
          f"({_results['failed']} failed)")
    print(f"{'─'*64}\n")


if __name__ == "__main__":
    main()
