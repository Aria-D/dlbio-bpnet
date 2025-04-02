import pyBigWig

# Input files
#bedgraph_file = "input.bedgraph"
chrom_sizes_file = "/athena/cayuga_0083/scratch/ard4008/ref/GRCh37.chrom.sizes"  # Download from UCSC or generate
#output_bigwig = "output.bw"

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
            data.append((chrom, int(start), int(end), float(value)))

    # Sort by chromosome and start position
    data.sort(key=lambda x: (x[0], x[1]))

    # Step 3: Write BigWig file
    bw = pyBigWig.open(output_bigwig, "w")
    bw.addHeader(list(chrom_sizes.items()))  # Add chromosome sizes