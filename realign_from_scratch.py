#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict
import logging

# --- Setup Logging ---
# Sets up a clear and informative log to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_sa_tag(sa_string):
    """
    Parses the SA tag from a BAM record. The SA tag lists supplementary alignments.
    Example: SA:Z:chr1,150,+,60M,60,0;...
    
    Args:
        sa_string (str): The raw string from the SA tag.
    
    Returns:
        list: A list of dictionaries, each representing a supplementary alignment.
    """
    supplementary_alignments = []
    if not sa_string:
        return supplementary_alignments
    
    parts = sa_string.split(';')
    for part in parts:
        if not part:
            continue
        fields = part.split(',')
        if len(fields) == 6:
            try:
                supp_dict = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'strand': fields[2],
                    'cigar': fields[3],
                    'mapq': int(fields[4]),
                    'nm': int(fields[5])
                }
                supplementary_alignments.append(supp_dict)
            except (ValueError, IndexError):
                # Ignore malformed SA tag fields
                continue
    return supplementary_alignments

def find_breakpoint_clusters(candidate_bam, max_dist, min_mapq):
    """
    Identifies and clusters split-read breakpoints from the candidate BAM file.
    
    Args:
        candidate_bam (str): Path to the coordinate-sorted BAM file of candidate reads.
        max_dist (int): The maximum distance between breakpoints to be considered part of the same cluster.
        min_mapq (int): The minimum mapping quality for a read to be considered.
        
    Returns:
        list: A list of clusters. Each cluster is a list of breakpoint tuples.
    """
    logging.info(f"Starting breakpoint identification from {candidate_bam}...")
    breakpoints = []
    
    with pysam.AlignmentFile(candidate_bam, "rb") as bamfile:
        for read in bamfile:
            # Skip unmapped, secondary, or low-quality reads
            if read.is_unmapped or read.is_secondary or read.mapping_quality < min_mapq:
                continue

            # A read is a split-read candidate if it has a supplementary alignment
            if read.has_tag('SA'):
                sa_tag_string = read.get_tag('SA')
                supp_alignments = parse_sa_tag(sa_tag_string)
                
                primary_chrom = read.reference_name
                primary_pos = read.reference_start + 1 # 1-based coordinate

                for supp in supp_alignments:
                    # We only care about splits on the same chromosome (intrachromosomal)
                    if supp['chrom'] == primary_chrom:
                        # Ensure we don't double-count by processing each pair once
                        # Store the breakpoint as (chrom, min_pos, max_pos)
                        breakpoint = (
                            primary_chrom,
                            min(primary_pos, supp['pos']),
                            max(primary_pos, supp['pos'])
                        )
                        breakpoints.append(breakpoint)

    logging.info(f"Found {len(breakpoints)} raw split-read events.")
    if not breakpoints:
        return []

    # Sort breakpoints to enable efficient clustering in a single pass
    breakpoints.sort()

    logging.info("Clustering sorted breakpoints...")
    clusters = []
    if not breakpoints:
        return clusters

    # Start the first cluster
    current_cluster = [breakpoints[0]]
    
    for i in range(1, len(breakpoints)):
        current_bp = breakpoints[i]
        # Get the representative breakpoint of the current cluster (the first one)
        cluster_rep = current_cluster[0]

        # Check if the current breakpoint belongs to the current cluster
        if (current_bp[0] == cluster_rep[0] and
            abs(current_bp[1] - cluster_rep[1]) <= max_dist and
            abs(current_bp[2] - cluster_rep[2]) <= max_dist):
            # It's close, add it to the current cluster
            current_cluster.append(current_bp)
        else:
            # It's too far, so the last cluster is finished. Save it.
            clusters.append(current_cluster)
            # Start a new cluster with the current breakpoint.
            current_cluster = [current_bp]
            
    # Don't forget to add the very last cluster
    clusters.append(current_cluster)

    logging.info(f"Grouped into {len(clusters)} clusters.")
    return clusters

def main():
    parser = argparse.ArgumentParser(
        description="A custom realignment algorithm to identify eccDNA breakpoints from split reads."
    )
    parser.add_argument(
        '-i', '--input_bam', required=True,
        help="Input BAM file containing coordinate-sorted candidate reads (e.g., sort_circular_read_candidates.bam)."
    )
    parser.add_argument(
        '-o', '--output_bed', required=True,
        help="Output BED file to store high-confidence eccDNA coordinates."
    )
    parser.add_argument(
        '--min_reads', type=int, default=2,
        help="Minimum number of unique split reads required to support an eccDNA call. Default: 2."
    )
    parser.add_argument(
        '--window', type=int, default=100,
        help="Maximum distance (in bp) between breakpoints to be grouped into a single cluster. Default: 100."
    )
    parser.add_argument(
        '--min_mapq', type=int, default=20,
        help="Minimum mapping quality of reads to be considered. Default: 20."
    )
    args = parser.parse_args()

    # Find clusters of breakpoints based on split reads
    clusters = find_breakpoint_clusters(args.input_bam, args.window, args.min_mapq)

    logging.info(f"Filtering clusters with at least {args.min_reads} supporting reads...")
    
    with open(args.output_bed, 'w') as bed_out:
        final_ecc_count = 0
        for i, cluster in enumerate(clusters):
            # The number of reads supporting this cluster is simply its length
            num_supporting_reads = len(cluster)

            if num_supporting_reads >= args.min_reads:
                final_ecc_count += 1
                
                # Determine a representative coordinate for the cluster (using the median)
                cluster.sort(key=lambda x: x[1])
                median_start = cluster[len(cluster) // 2][1]
                
                cluster.sort(key=lambda x: x[2])
                median_end = cluster[len(cluster) // 2][2]
                
                chrom = cluster[0][0]
                name = f"eccDNA_cluster_{final_ecc_count}"
                score = num_supporting_reads
                strand = "."
                
                # Write to BED file
                bed_out.write(f"{chrom}\t{median_start}\t{median_end}\t{name}\t{score}\t{strand}\n")

    logging.info(f"Process complete. Wrote {final_ecc_count} high-confidence eccDNAs to {args.output_bed}.")


if __name__ == '__main__':
    main()
