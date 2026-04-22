import math

def TwoHitSeeds(single_hits, k, A):
    """
    Two-hit filtration: emit a seed only when two non-overlapping k-mer hits
    share a diagonal and are within distance A of each other. Nearby valid hits
    are merged into a single representative seed.
    """
    if not single_hits:
        return []

    diagonals = {}
    for q_pos, r_pos, score in single_hits:
        diag = q_pos - r_pos
        diagonals.setdefault(diag, []).append((q_pos, r_pos, score))

    seeds = []

    for diag, hits in diagonals.items():
        hits.sort(key=lambda x: x[0])

        # Mark which hits belong to at least one valid pair
        in_valid_pair = [False] * len(hits)

        for i in range(1, len(hits)):
            q2 = hits[i][0]
            for j in range(i - 1, -1, -1):
                q1 = hits[j][0]
                gap = q2 - q1

                if gap > A:
                    break
                if gap >= k:
                    in_valid_pair[i] = True
                    in_valid_pair[j] = True
                    break

        # Merge local runs of valid hits → one seed per cluster
        i = 0
        while i < len(hits):
            if not in_valid_pair[i]:
                i += 1
                continue

            cluster = [hits[i]]
            i += 1

            while i < len(hits) and in_valid_pair[i]:
                prev_q = cluster[-1][0]
                curr_q = hits[i][0]

                if curr_q - prev_q <= A:
                    cluster.append(hits[i])
                    i += 1
                else:
                    break

            q_pos, r_pos, _ = cluster[0]
            merged_score = sum(h[2] for h in cluster)
            seeds.append((q_pos, r_pos, merged_score))

    return seeds
