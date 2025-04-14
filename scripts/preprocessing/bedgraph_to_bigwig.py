import pyBigWig
from collections import defaultdict

# Path to chromosome sizes file
chrom_sizes_file = "/athena/cayuga_0083/scratch/ard4008/ref/GRCh37.chrom.sizes"

def bdg_to_bw(bedgraph_file):
    output_bigwig = bedgraph_file.replace("bedgraph", "bw")

    # Step 1: Read chromosome sizes
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)

    # Step 2: Read and sort BEDGraph data
    data = []
    with open(bedgraph_file, 'r') as f:
        for line in f:
            if line.startswith(('track', 'browser', '#')):  # Skip headers
                continue
            chrom, start, end, value = line.strip().split("\t")
            chrom = chrom.replace("chr", "") #"chr" + chrom if not chrom.startswith("chr") else chrom  # Ensure "chr" prefix
            data.append((chrom, int(start), int(end), float(value)))

    # Step 3: Check if chromosome names match
    bedgraph_chroms = {row[0] for row in data}
    chrom_sizes_chroms = set(chrom_sizes.keys())
    if not bedgraph_chroms.issubset(chrom_sizes_chroms):
        print("Warning: Chromosome names in BEDGraph do not match the sizes file!")
        # Find mismatched chromosomes
    mismatched_chroms = bedgraph_chroms - chrom_sizes_chroms
    if mismatched_chroms:
        print(f"Mismatched chromosomes in BEDGraph file: {mismatched_chroms}")

    # Step 4: Remove duplicate start positions (if any)
    grouped_data = defaultdict(list)
    for chrom, start, end, value in data:
        grouped_data[(chrom, start)].append(value)

    final_data = []
    for (chrom, start), values in grouped_data.items():
        avg_value = sum(values) / len(values)
        final_data.append((chrom, start, start + 1, avg_value))  # Assuming each entry has a span of 1

    final_data.sort(key=lambda x: (x[1], x[2]))  # Sort again

    # Step 5: Write BigWig file
    bw = pyBigWig.open(output_bigwig, "w")
    bw.addHeader(list(chrom_sizes.items()))  # Add chromosome sizes

    # Step 6: Add entries to BigWig file in batches
    #chroms, starts, ends, values = zip(*final_data)
    chroms = [f[0] for f in final_data]
    starts = [f[1] for f in final_data]
    ends = [f[2] for f in final_data]
    values = [f[3] for f in final_data]
    bw.addEntries(chroms, starts, ends=ends, values=values, span = 1)

    # Step 7: Close BigWig file
    bw.close()
    print(f"BigWig file saved as: {output_bigwig}")

# Example usage:
# bdg_to_bw("input.bedgraph")
