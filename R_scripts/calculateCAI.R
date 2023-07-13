library(AnaCoDa,lib="~/R_dev")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character",default="./")
parser$add_argument("--codon_table",help="Codon Table",type="integer",default=1)

args <- parser$parse_args()
input <- args$input
codon_table <- args$codon_table

ribo <- initializeGenomeObject(paste0("Top10_phi/",input),codon_table = codon_table)
genome <- initializeGenomeObject(paste0("/data2/Labella2019/Genomes/cds_cleaned/",input),codon_table = codon_table)
cai <- getCAI(ribo,genome)
cai <- data.frame(GeneID=names(cai),Mean=unname(cai),stringsAsFactors = F)
write.csv(cai,paste0("CAI_top_10/",input),sep="\t",row.names = F,col.names=T,quote=F)
