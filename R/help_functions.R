manhattan_options=function(){
  print("options(repr.plot.width=12, repr.plot.height=4.5,repr.plot.res = 170)")
}
regionplot_options=function(){
  print("options(repr.plot.width=11, repr.plot.height=6.5,repr.plot.res = 170)")
}


get_exons_cmd=function(chr="chr1",xmin=66716513,xmax=66716513){
  chr=gsub("chr", "", chr)
  cand_region=paste("chr",chr,":",xmin,"-",xmax,sep="")
  print(paste("gor -p ",cand_region," ref/ensgenes/ensgenes.gorz | map -c gene_symbol <(nor -h ref/ensgenes/ensgenes_exons.gorz | select gene_symbol,chromstart,chromend,exon | rename chromstart exon_chromstart | rename chromend exon_chromend) ",sep=""))
}

get_genes_cmd=function(chr="chr1",xmin=66716513,xmax=66716513){
  chr=gsub("chr", "", chr)
  cand_region=paste("chr",chr,":",xmin,"-",xmax,sep="")
  print(paste("gor -p ",cand_region," ref/ensgenes/ensgenes.gorz | select chrom,gene_start,gene_end,gene_symbol,gene_stable_id,biotype,strand", sep=""))
}