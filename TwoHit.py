import math

def TwoHitSeeds(single_hits, k, A):
    """
    Merge overlapping/nearby hits on the same diagonal, while preserving
    the 2-hit rule: a merged seed is only emitted if its cluster contains
    at least two non-overlapping hits within distance A.
    
    Input:
        single_hits: list of (q_pos, r_pos, score)
        k: kmer length
        A: max allowed spacing between consecutive hits in a 2-hit group

    Output:
        list of merged seeds as (q_pos, r_pos, score), where each output
        seed represents one diagonal cluster that passed the 2-hit rule.
    """

    if not single_hits:
        return []

    # group hits by diagonal
    diagonals = {}
    for q_pos, r_pos, score in single_hits:
        diag = q_pos - r_pos
        diagonals.setdefault(diag, []).append((q_pos, r_pos, score))

    merged_valid_seeds = []

    for diag, hits in diagonals.items():
        # sort left-to-right along query
        hits.sort(key=lambda x: x[0])

        # build clusters of nearby hits on this diagonal
        cluster = [hits[0]]

        for curr_hit in hits[1:]:
            prev_hit = cluster[-1]
            distance = curr_hit[0] - prev_hit[0]

            # if close enough to still belong to same local neighborhood, keep clustering
            if distance <= A:
                cluster.append(curr_hit)
            else:
                # process finished cluster
                merged = process_two_hit_cluster(cluster, k, diag)
                if merged is not None:
                    merged_valid_seeds.append(merged)

                # start new cluster
                cluster = [curr_hit]

        # process final cluster
        merged = process_two_hit_cluster(cluster, k, diag)
        if merged is not None:
            merged_valid_seeds.append(merged)

    return merged_valid_seeds


def process_two_hit_cluster(cluster, k, diag):
    """
    Given a list of hits all on one diagonal and within a local window,
    decide whether the cluster satisfies the 2-hit rule and, if so,
    return one merged representative seed.
    """

    if len(cluster) < 2:
        return None

    # Count non-overlapping hits greedily
    nonoverlap_hits = [cluster[0]]
    last_q = cluster[0][0]

    for hit in cluster[1:]:
        q_pos = hit[0]

        # preserve the actual 2-hit condition: must be non-overlapping
        if q_pos - last_q >= k:
            nonoverlap_hits.append(hit)
            last_q = q_pos

    # Must contain at least 2 valid non-overlapping hits
    if len(nonoverlap_hits) < 2:
        return None

    # Merge the cluster into one representative seed.
    # Best choice: keep the highest-scoring hit as anchor.
    best_hit = max(cluster, key=lambda x: x[2])

    # Optional: boost score using all hits in cluster
    merged_score = sum(hit[2] for hit in cluster)

    q_pos, r_pos, _ = best_hit
    return (q_pos, r_pos, merged_score)