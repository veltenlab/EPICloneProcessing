library(argparse)
library(GenomicRanges)
ap <- ArgumentParser()
ap$add_argument("-a", "--amplicons", action="store", help="The amplicon file for which genomic information is to be extracted")
ap$add_argument("-o", "--output", action="store", help="The output file name")
ap$add_argument("-s", "--assembly", action="store", help="The genome assembly version", default='mm10')
opt <- ap$parse_args()
if(opt$assembly=='mm10'){
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
}else if(opt$assembly=='hg19'){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
}else{
  stop('Unsuported genome assembly')
}
amplicon_fr <- read.table(opt$amplicons)
fasta_file <- apply(amplicon_fr, 1, function(x){
  seq_info <- genome[[x[1]]][x[2]:x[3]]
  return(paste0('>', x[4], '\n',
                seq_info))
})
writeLines(fasta_file, opt$output)