library(argparse)
library(ggplot2)
ap <- ArgumentParser()
ap$add_argument("-f", "--file", action="store", help="The barcode distribution file as input to cellfinder (txt)")
ap$add_argument("-d", "--doublet", action="store", help="The doublet detection file")
ap$add_argument("-a", "--amplicon", action="store", help="The amplicon description file")
ap$add_argument("-o", "--output", action="store", help="The output folder")
ap$add_argument("-n", "--name", action="store", help="The output file name")
opt <- ap$parse_args()
dat <- read.table(opt$file,
                  header=T)
rowinfo <- read.csv(opt$doublet)
cells <- rowinfo$Barcode
ampli_info <- read.table(opt$amplicon,
                         row.names = 1)
sel_amplis <- row.names(ampli_info)[(grepl('Non_cut',ampli_info$Type))]
reads_per_cell <- rowSums(dat[, sel_amplis])
to_plot <- data.frame(Rank=rank(-reads_per_cell),
                  Reads=reads_per_cell,
                  Type=ifelse(row.names(dat)%in%cells, 'Cells', 'Background'))
plot <- ggplot(to_plot, aes(x=log10(Rank), y=log10(Reads), color=Type))+
  xlab('log10(Cells)')+ylab('log10(Reads)')+geom_point()+theme_bw()+
  scale_color_manual(values=c(Cells='black', Background='gray80'))
ggsave(file.path(opt$output, paste0(opt$name, '_reads_per_cell.png')),
       plot)
sink(file.path(opt$output, paste0(opt$name, '_stats.txt')))
cat(paste('Aligned reads:', sum(dat), '\n'))
cat(paste('Reads in cells:', sum(dat[cells, ]), '\n'))
cat(paste('Numer of cells:', nrow(rowinfo), '\n'))
sink()
metadata_out <- file.path(opt$output, paste0(opt$name, '_metadata.csv'))
to_write <- data.frame(NA, sum(dat), NA, NA, NA, NA, nrow(rowinfo), NA, sum(dat[cells, ]), 'v2')
write.table(to_write,
    metadata_out,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE,
    sep=",")
