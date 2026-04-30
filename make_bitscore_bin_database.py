"""
Generate a small nucleotide test database (5 sequences) where each reference
produces a hit in a specific bit-score bin when queried with a paired query.

Bit-score bins targeted (per the NCBI-style graphical summary):
    >=200, 80-200, 50-80, 40-50, <40

Scoring assumed: match=+2, mismatch=-3, blastn params lambda=0.625, K=0.41
(bit = (lambda * raw - ln K) / ln 2 ~= 0.9018 * raw + 1.29).
For an ungapped perfect match of length N, raw = 2N and bit ~= 1.80 * N + 1.29.

Design:
    Each reference is [random flank] + [exact match to a slice of the query]
    + [random flank]. The exact match length is chosen to land the bit score
    inside the target bin. Random flanks are drawn from ACGT uniformly --
    chance of long accidental secondary alignments is negligible.

Outputs (written next to this script):
    bitscore_test_query.fasta    -- single 400 bp query
    bitscore_test_database.fasta -- 5 reference sequences
"""

import os
import math
import random


def _read_first_fasta_sequence(path: str) -> str:
    """Minimal FASTA reader -- returns the first record's concatenated sequence."""
    seq_parts = []
    found_header = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if found_header:
                    break  # only want the first record
                found_header = True
                continue
            if found_header:
                seq_parts.append(line.strip())
    return "".join(seq_parts)

# Reproducible random flanks
random.seed(20260422)

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_FASTA = os.path.join(HERE, "H1N1_HA_single.fasta")

QUERY_LEN = 400  # length of the query to pull from the H1N1 HA sequence

# (bin_name, perfect-match length, offset into query where the match region
#  starts). Offsets are chosen so matches land at different positions on the
#  query for a more interesting graphical summary.
REFS = [
    # (name, perfect-match length, offset into query where the match begins)
    # match_len chosen to land comfortably inside each bin after accounting
    # for possible +2..+8 raw bonus from lucky extension into random flanks.
    ("REF01_bin1_ge200",   300,  40),   # raw 600 -> bit ~542  (well above 200)
    ("REF02_bin2_80to200",  70, 240),   # raw 140 -> bit ~127  (mid of 80-200)
    ("REF03_bin3_50to80",   36, 120),   # raw  72 -> bit  ~66  (mid of 50-80)
    ("REF04_bin4_40to50",   24, 320),   # raw  48 -> bit  ~45  (mid of 40-50)
    ("REF05_bin5_lt40",     16, 180),   # raw  32 -> bit  ~30  (comfortably <40)
]

# Flank sizes per reference (pre- and post-match random DNA). Just enough to
# make the references look like ordinary DB entries without giving chance
# alignments room to accumulate.
FLANK_LEFT  = 120
FLANK_RIGHT = 120


def random_dna(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def expected_bit_score(match_len: int, lam: float = 0.625, K: float = 0.41) -> float:
    raw = 2 * match_len
    return (lam * raw - math.log(K)) / math.log(2)


def main():
    # Pull a 400 bp chunk of real H1N1 HA as the query
    source_seq = _read_first_fasta_sequence(SOURCE_FASTA).upper().replace("N", "A")

    # Pick a window well inside the source record
    query_seq = source_seq[100:100 + QUERY_LEN]
    assert len(query_seq) == QUERY_LEN, (
        f"Source sequence too short: got {len(source_seq)} bp, need at least "
        f"{100 + QUERY_LEN} bp"
    )

    # Write the query FASTA
    query_path = os.path.join(HERE, "bitscore_test_query.fasta")
    with open(query_path, "w") as fh:
        fh.write(">TEST_QUERY_400bp_from_H1N1_HA\n")
        fh.write(query_seq + "\n")

    # Build the five reference sequences
    db_path = os.path.join(HERE, "bitscore_test_database.fasta")
    with open(db_path, "w") as fh:
        for name, match_len, start in REFS:
            assert start + match_len <= QUERY_LEN, (
                f"{name}: match window {start}..{start+match_len} exceeds "
                f"query length {QUERY_LEN}"
            )
            match_seg = query_seq[start:start + match_len]
            left  = random_dna(FLANK_LEFT)
            right = random_dna(FLANK_RIGHT)
            ref_seq = left + match_seg + right

            expected_bit = expected_bit_score(match_len)
            header = (
                f">{name} match_len={match_len} match_at_query={start}-"
                f"{start + match_len} expected_bit={expected_bit:.1f}"
            )
            fh.write(header + "\n")
            fh.write(ref_seq + "\n")

    print(f"Wrote query:    {query_path}")
    print(f"Wrote database: {db_path}")
    print()
    print(f"{'Reference':<26}{'match_len':>10}{'expected_raw':>14}{'expected_bit':>14}  bin")
    for name, match_len, _ in REFS:
        raw = 2 * match_len
        bit = expected_bit_score(match_len)
        bin_label = name.split("_", 2)[-1]
        print(f"{name:<26}{match_len:>10}{raw:>14}{bit:>14.1f}  {bin_label}")


if __name__ == "__main__":
    main()
