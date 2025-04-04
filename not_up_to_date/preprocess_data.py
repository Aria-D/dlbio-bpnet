import os
import numpy as np
import pandas as pd
from bedgraph_to_bigwig import *

#genome_dir = "/athena/cayuga_0083/scratch/ard4008/ref/GRCh37_genome.fa"
gse_to_disease_250k = {"GSE110607": "SLE",  "GSE71841": "RA"} # 
gse_to_disease_27k = {"GSE56606" : "T1D"}

outdir = "/athena/cayuga_0083/scratch/ard4008/dlbio_final_project/preprocessed_data"
illumina_metadata_path= "/athena/khuranalab/scratch/tog4002/DeepLearningProject/data/Illumina450kInfo.csv"
illumina_metadata_250k = pd.read_csv(illumina_metadata_path, sep = ",", skiprows = 7)
illumina_metadata_250k = illumina_metadata_250k.set_index("IlmnID")

illumina_metadata_path = "/athena/khuranalab/scratch/tog4002/DeepLearningProject/data/Illumina27kInfo.csv"
illumina_metadata_27k = pd.read_csv(illumina_metadata_path, sep = ",", skiprows = 150)
illumina_metadata_27k = illumina_metadata_27k.set_index("IlmnID")
illumina_metadata_27k.rename(columns={"Chr": "CHR", "MapInfo": "MAPINFO"}, inplace=True)

sample_metadata_path = "/athena/khuranalab/scratch/tog4002/DeepLearningProject/data/metadata.tsv"
sample_metadata = pd.read_csv(sample_metadata_path, sep = "\t")

data_dir = "/athena/khuranalab/scratch/tog4002/DeepLearningProject/data/"


def split_df(df, sample_name):
    df_subset = sample_metadata[sample_metadata.GSE == sample_name]
    healthy_ids, disease_ids = df_subset[df_subset.Condition == "Healthy"].Sample.tolist(), df_subset[df_subset.Condition != "Healthy"].Sample.tolist()
    return df[healthy_ids], df[disease_ids]

def get_genome_coordinates(methylation_sites, illumina_metadata):
    regions_metadata = illumina_metadata.loc[methylation_sites]
    chr, pos = ("chr" + regions_metadata.CHR.astype(str)).tolist(), regions_metadata.MAPINFO.astype(int).tolist()
    assert regions_metadata.index.tolist() == methylation_sites
    return chr, pos

def create_bed_file(df, chr, pos, scores):
    df["chr"] = chr
    df["start"] = pos
    df["end"] = np.array(pos) + 1
    df["scores"] = scores
    print(df)
    return df

for file in os.listdir(data_dir):
    if "GSE" in file and file.endswith(".tsv"):
        gse_id = file.split("_")[0].replace(".tsv", "")
        print(gse_id)
        print("Processing")
        if gse_id in gse_to_disease_250k:
            illumina_metadata = illumina_metadata_250k
            disease = gse_to_disease_250k[gse_id]
        elif gse_id in gse_to_disease_27k:
            illumina_metadata = illumina_metadata_27k
            disease = gse_to_disease_27k[gse_id]
        else: 
            continue
        
        data_tsv = pd.read_csv(f"{data_dir}/{file}", sep = "\t")
        data_tsv = data_tsv.set_index("gene")
        healthy_df, disease_df = split_df(data_tsv, gse_id)
        healthy_chr, healthy_pos = get_genome_coordinates(healthy_df.index.tolist(), illumina_metadata)
        disease_chr, disease_pos = get_genome_coordinates(disease_df.index.tolist(), illumina_metadata)

        mean_healthy_scores = healthy_df.mean(1).tolist()
        mean_disease_scores = disease_df.mean(1).tolist()
        assert len(mean_healthy_scores) == healthy_df.shape[0]
        assert len(mean_disease_scores) == disease_df.shape[0]

        healthy_bedgraph = pd.DataFrame()
        healthy_bedgraph = create_bed_file(healthy_bedgraph, healthy_chr, healthy_pos, mean_healthy_scores)
        healthy_bedgraph["chr"] = healthy_bedgraph["chr"].replace("chr22.0", "chr22") 
        healthy_bedgraph.to_csv(f"{outdir}/{disease}_healthy.bedgraph", sep = "\t", index = None, header = None)
        bdg_to_bw(f"{outdir}/{disease}_healthy.bedgraph")
        healthy_bedgraph[["chr", "start", "end"]].to_csv(f"{outdir}/{disease}_healthy_summits.bed", index = None, sep = "\t", header = None)


        disease_bedgraph = pd.DataFrame()
        disease_bedgraph = create_bed_file(disease_bedgraph, disease_chr, disease_pos, mean_disease_scores)
        disease_bedgraph["chr"] = disease_bedgraph["chr"].replace("chr22.0", "chr22") 
        disease_bedgraph.to_csv(f"{outdir}/{disease}_disease.bedgraph", sep = "\t", index = None, header = None)
        bdg_to_bw(f"{outdir}/{disease}_disease.bedgraph")
        disease_bedgraph[["chr", "start", "end"]].to_csv(f"{outdir}/{disease}_disease_summits.bed", index = None, sep = "\t", header = None)

        print(f"Processed {disease}.")
