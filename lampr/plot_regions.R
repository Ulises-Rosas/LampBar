#!/usr/bin/Rscript

suppressMessages({
  
  library(data.table)
  library(dplyr)
  library(optparse)
  library(doParallel)
  library(ggplot2)
  library(reticulate)
  library(xtable)
  library(ape)

})

option_list = list(

        make_option(c("-t", "--taxon"), type="character",
                    default=NULL,
                    help="Target taxa",
                    metavar="character"),
        make_option(c("-s", "--sequences"), type="character",
                    default="curated_regions.txt",
                    help="Audited sequences from findRegion function
                    [default= %default]",
                    metavar="character"),
        make_option(c("-n", "--nproc"), type="numeric",
                    default= 2,
                    help="Number of processor used for getting diagnostic nucleotides
                    [default= %default]",
                    metavar="numeric"),
        make_option(c("-r", "--res"), type="numeric",
                    default= 75,
                    help="Resolution of plots in reports
                    [default= %default]",
                    metavar="numeric")

)

opt_parser = OptionParser(option_list=option_list, description = "

                   MinBar: Automated pipeline for finding out specie-specific primers
                   ==================================================================
                   ><{{{*> . ><(((*> . ><{{{*> . ><(((*> . ><{{{*> . ><(((*> . ><{{{*
                   ==================================================================

                              * Developed by U. Rosas and X. Velez-Zuazo
                          ")

opt = parse_args(opt_parser)

if ( is.null(opt$`taxon`) ){

        print_help(opt_parser)
        stop("--taxon argument must be supplied.\n", call.=FALSE)

}

Sys.setenv(RETICULATE_PYTHON = "path_python_executable")

use_python("path_python_executable")

urllib <- reticulate::import("urllib")

reticulate::source_python("path_ssp/minbar_py/worms.py", convert = F)
#print(opt$taxon)
# checking your taxon! there's might be a name modification
validated_name = Worms(opt$`taxon`)$get_accepted_name() %>% as.character(.)

if( !grepl("^[A-Z][a-z]+ [a-z]+$", validated_name) )
        validated_name = Worms(opt$`taxon`)$taxamatch() %>% as.character(.)

if( !grepl("^[A-Z][a-z]+ [a-z]+$", validated_name) ){
  
  writeLines("\n\n\033[0;31mCheck your taxon!\033[0m
                   ")
  stop("Name of target taxon can not be validated with WoRMS\n", call.=FALSE)
}

get_spps_names <- function(seqs,format="obitools"){

        if(format == "obitools"){

                names(seqs) %>%
                        gsub(".* organism=([A-Z][a-z]+ [a-z]+); .*", "\\1", x = .) %>%
                        return(.)
        }else if(format == "bold"){

                sapply(strsplit(names(seqs), split = "\\|"), function(x) x[2]) %>%
                        return(.)
        }
}

get_positions <- function(x, diags){
  
  pos = vector('numeric')
  
  for( node in seq_along(diags) ){
    
    item    = diags[[ node ]]
    subItem = item[ names(item) == x ][[1]]
    
    if( length(subItem) > 0 )
      pos = append( pos, values = subItem + node )
  }
  
  if( length(pos) > 0 )
    data.frame(
      
      spps  = x,
      N_pos = length(unique(pos)),
      pos   = paste(
        sort( unique(pos) ),
        collapse = " "
      )
    )
}

### spider functions
source("path_ssp/spider/R/slidingWindow.R")
source("path_ssp/spider/R/nucDiag.R")
### spider functions


# dealing with sequences right here:
#print(opt$sequences)
seq_name = strsplit(opt$`sequences`, "\\.") %>% .[[1]] %>% .[1]

# name got it from mafft software
complete_aln <- paste(seq_name, "_aln.txt", sep = "") %>% read.FASTA(.)

spps_names <- get_spps_names(seqs = complete_aln, format = "obitools")

complete_wins <- slidingWindow(complete_aln, 100, 1)

#######parallel working:
cl <- makeCluster( opt$`nproc` )
registerDoParallel(cl)

writeLines("\nGetting diagnostic nucleotides...\n")

foreach(i=1:length(complete_wins)) %dopar%
  nucDiag( complete_wins[[i]], spps_names ) -> complete_winDiag 

foreach( i =  unique(spps_names)  ) %dopar%
  get_positions( i, diags = complete_winDiag) %>%
  do.call('rbind',.) -> df_pos
stopCluster(cl)
########




if( !is.null(df_pos) )
  write.table(
    x         = df_pos,
    sep       = ",",
    file      = paste0( gsub(" ", "_",opt$taxon), "_nucDiag_from_aln.txt"),
    quote     = FALSE,
    row.names = FALSE
  )


res2 <- t( sapply( complete_winDiag, lengths ) )

all_diag_nucls <- lapply(colnames(res2), function(x){
  
        data.frame(species = x,
                   position = 1:length(res2[,x]),
                   diag_nuc = res2[,x] # bug producer
        )
  }) %>%
        do.call("rbind", .) %>%
        dplyr::mutate(species = as.factor(species))
########
#######START: Boxing area
#######
rem = paste(seq_name, "_chars_removed", sep = "") %>%
        read.table(file = ., sep = ",") %>%
        t(.) %>%
        as.data.frame(x=.) %>%
        dplyr::mutate(fact = "removed") %>%
        dplyr::mutate(chars = V1)  %>%
        dplyr::select(fact, chars)


if (  all( 1:dim(res2)[1] %in% rem$chars  ) ){
  
  bound = rem[ order( rem$chars ),  ]
  
  }else if (  all( 1:dim(res2)[1] %in% rem$chars  ) == F  ) {
    
    con = data.frame(fact= "conserved",
                     chars= c( 1:dim( res2 )[1] )[ !(1:dim(res2)[1] %in% rem$chars ) ] )
    
    bound = con[ order( con$chars ),  ]
    
  }else{
    
    con = data.frame(fact= "conserved",
                     chars= c( 1:dim( res2 )[1] )[ !(1:dim(res2)[1] %in% rem$chars ) ] )
    
    bound = rbind(rem, con) %>% .[order(.$chars),]
  
  
}

for(i in 1:nrow(bound))
        if(bound[i,1] == "conserved")
                bound[i,2] = NA

suppressWarnings({
  
  area = paste(bound$chars, collapse = " ") %>%
    strsplit(x = ., split = "NA") %>% .[[1]] %>%
    lapply(X = ., function(x){
                if(x != " "){
                        x = gsub("^[ ]+|[ ]+$", "", x)
                        data.frame(x1 = strsplit(x, " ") %>% .[[1]] %>% as.numeric(.) %>% min(.,na.rm = T) %>% {if_else( is.infinite(.), 0, .)},
                                   x2 = strsplit(x, " ") %>% .[[1]] %>% as.numeric(.) %>% max(.,na.rm = T) %>% {if_else( is.infinite(.), 0, .)},
                                   y1 = 0,
                                   y2 = Inf)
                        }
                }) %>%
        do.call("rbind", .)
})

########
######FINISH: Boxing area
#######

writeLines("\nPlotting diagnostic nucleotides...\n")

target_species = all_diag_nucls[grepl(validated_name, all_diag_nucls$species),]

gsub(" ","_",opt$taxon) %>%
        paste(., "_diag_nucls.tiff", sep = "") %>%
  tiff(filename = .,width = 8,height = 4,
       units = 'in',res = opt$`res`, compression = 'lzw+p')
ggplot(data = all_diag_nucls, aes(x = position, y = diag_nuc)) +
        geom_rect(data = area %>% dplyr::mutate(x1 = if_else(x1 > 100, x1 - 100, x1),
                                                x2 =if_else(x2 > 100, x2 - 100, x2)),
                  alpha =0.2,
                  mapping  = aes(xmin = x1,
                                 xmax = x2,
                                 ymin = y1,
                                 ymax = y2),
                  inherit.aes = F)+
        geom_line(aes(linetype=species), size =0.75)+
        scale_linetype_manual(values = rep("dotdash", length(unique(spps_names))))+
        guides(linetype = FALSE)+
        geom_line(data = target_species,
                  aes(x = position, y = diag_nuc, color = species),
                  size = 2)+
        scale_colour_brewer(palette = "Set1")+
        labs(y='Diagnostic nucleotides',
             x='Starting position of sliding windows (bp)',
             title='Diagnostic nucleotides along sequence',
             subtitle = "Complete alignment")+
        theme_bw(base_size = 15) +   #Proporción entre los caracteres y la gráfica
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 15, face = 'bold')
              #panel.background = element_rect(fill = "transparent"),
              #plot.background = element_rect(fill = "transparent", color = NA),
              #legend.background = element_rect(fill = "transparent")
              #legend.box.background = element_rect(fill = "transparent")
              )
invisible(dev.off())

writeLines("\nPlotting entropy of sequences...\n")

all_entropies <- paste(seq_name, "_aln_entropies.txt", sep = "") %>%
        data.table::fread(.) %>%
        dplyr::select(species, position, smoothed_entropy) %>%
        dplyr::mutate(position = as.numeric(position),
                      smoothed_entropy = as.numeric(smoothed_entropy),
                      species = as.factor(species))

target_species2 = all_entropies[grepl(validated_name,all_entropies$species),]

gsub(" ","_",opt$taxon) %>%
        paste(., "_entropies.tiff", sep = "") %>%
  tiff(filename = .,width = 8,height = 4,units = 'in',
       res = opt$`res`, compression = 'lzw+p')
ggplot(data = all_entropies, aes(x = position, y = smoothed_entropy)) +
        geom_point() +
        geom_rect(data = area, alpha =0.2,
                  mapping  = aes(xmin = x1,
                                 xmax = x2,
                                 ymin = y1,
                                 ymax = y2),
                  inherit.aes = F)+
        geom_point(data = target_species2,
                   aes(x = position, y = smoothed_entropy, color = species),
                   size = 2) +
        scale_colour_brewer(palette = "Set1")+
        labs(y='Smoothed entropy',
             x='Position',
             title='Entropy along sequence',
             subtitle = "Complete alignment")+
        theme_bw(base_size = 15) +   #Proporción entre los caracteres y la gráfica
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 15, face = 'bold')
              #panel.background = element_rect(fill = "transparent"),
              #plot.background = element_rect(fill = "transparent", color = NA),
              #legend.background = element_rect(fill = "transparent") 
              #legend.box.background = element_rect(fill = "transparent")
              )
invisible(dev.off())
