from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CHROM1="NC_000001.10 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
CHROM2="NC_000002.11 Homo sapiens chromosome 2, GRCh37 primary reference assembly"
CHROM3="NC_000003.11 Homo sapiens chromosome 3, GRCh37 primary reference assembly"
CHROM4="NC_000004.11 Homo sapiens chromosome 4, GRCh37 primary reference assembly"
CHROM5="NC_000005.9 Homo sapiens chromosome 5, GRCh37 primary reference assembly"
CHROM6="NC_000006.11 Homo sapiens chromosome 6, GRCh37 primary reference assembly"
CHROM7="NC_000007.13 Homo sapiens chromosome 7, GRCh37 primary reference assembly"
CHROM8="NC_000008.10 Homo sapiens chromosome 8, GRCh37 primary reference assembly"
CHROM9="NC_000009.11 Homo sapiens chromosome 9, GRCh37 primary reference assembly"
CHROM10="NC_000010.10 Homo sapiens chromosome 10, GRCh37 primary reference assembly"
CHROM11="NC_000011.9 Homo sapiens chromosome 11, GRCh37 primary reference assembly"
CHROM12="NC_000012.11 Homo sapiens chromosome 12, GRCh37 primary reference assembly"
CHROM13="NC_000013.10 Homo sapiens chromosome 13, GRCh37 primary reference assembly"
CHROM14="NC_000014.8 Homo sapiens chromosome 14, GRCh37 primary reference assembly"
CHROM15="NC_000015.9 Homo sapiens chromosome 15, GRCh37 primary reference assembly"
CHROM16="NC_000016.9 Homo sapiens chromosome 16, GRCh37 primary reference assembly"
CHROM17="NC_000017.10 Homo sapiens chromosome 17, GRCh37 primary reference assembly"
CHROM18="NC_000018.9 Homo sapiens chromosome 18, GRCh37 primary reference assembly"
CHROM19="NC_000019.9 Homo sapiens chromosome 19, GRCh37 primary reference assembly"
CHROM20="NC_000020.10 Homo sapiens chromosome 20, GRCh37 primary reference assembly"
CHROM21="NC_000021.8 Homo sapiens chromosome 21, GRCh37 primary reference assembly"
CHROM22="NC_000022.10 Homo sapiens chromosome 22, GRCh37 primary reference assembly"

autosome_desc = [CHROM1, CHROM2, CHROM3, CHROM4, CHROM5, CHROM6, CHROM7, CHROM8, CHROM9, CHROM10,
                 CHROM11, CHROM12, CHROM13, CHROM14, CHROM15, CHROM16, CHROM17, CHROM18, CHROM19,
                 CHROM20, CHROM21, CHROM22]

new_fasta_samples = []
for record in SeqIO.parse("./GCF_000001405.13_GRCh37_genomic.fna", "fasta"):
    if record.description in autosome_desc: # if the header meets our requirement
        chrom_id = autosome_desc.index(record.description) + 1
        new_fasta_samples.append(SeqRecord(record.seq, id = "chr" + str(chrom_id), description="")) # rewrite the header of the sample

SeqIO.write(new_fasta_samples, "clean_hg19.fasta", "fasta")
    # print(record.seq)


