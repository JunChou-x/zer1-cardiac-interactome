# ============================================================
# 01_sequence_annotation.R
#
# Retrieve UniProt sequences and annotate N-terminal features
# ============================================================

library(readxl)
library(writexl)
library(dplyr)

input_file <- "data/ZER1_IPMS_significant_hits.xlsx"
output_file <- "results/01_sequence_annotation.xlsx"

df <- read_excel(input_file)

df$Nterm_2aa <- NA
df$Mature_Nt1 <- NA
df$Met_cleavage_predicted <- FALSE
df$Nt1_in_GASTC <- FALSE

cat("Retrieving UniProt sequences...\n")

for(i in 1:nrow(df)){
  
  acc <- df$Accession[i]
  
  if(!is.na(acc) && acc != ""){
    
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/",
      acc,
      ".fasta"
    )
    
    tryCatch({
      
      lines <- readLines(url,warn=FALSE)
      
      if(length(lines)>1){
        
        seq <- paste(lines[-1],collapse="")
        
        if(nchar(seq)>=2){
          
          df$Nterm_2aa[i] <- substr(seq,1,2)
          
          r2 <- substr(seq,2,2)
          
          df$Mature_Nt1[i] <- r2
          
          df$Met_cleavage_predicted[i] <-
            r2 %in% c("A","C","G","P","S","T","V")
          
          df$Nt1_in_GASTC[i] <-
            r2 %in% c("G","A","S","T","C")
        }
      }
      
    },error=function(e) NULL)
    
  }
  
  Sys.sleep(0.3)
}

# Retry failed retrievals

failed_rows <- which(
  is.na(df$Nterm_2aa) &
    !is.na(df$Accession)
)

if(length(failed_rows)>0){
  
  cat(
    sprintf(
      "Retrying %d failed UniProt queries...\n",
      length(failed_rows)
    )
  )
  
  for(i in failed_rows){
    
    acc <- df$Accession[i]
    
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/",
      acc,
      ".fasta"
    )
    
    tryCatch({
      
      lines <- readLines(url,warn=FALSE)
      
      if(length(lines)>1){
        
        seq <- paste(lines[-1],collapse="")
        
        if(nchar(seq)>=2){
          
          df$Nterm_2aa[i] <- substr(seq,1,2)
          
          r2 <- substr(seq,2,2)
          
          df$Mature_Nt1[i] <- r2
          
          df$Met_cleavage_predicted[i] <-
            r2 %in% c("A","C","G","P","S","T","V")
          
          df$Nt1_in_GASTC[i] <-
            r2 %in% c("G","A","S","T","C")
        }
      }
      
    },error=function(e) NULL)
    
    Sys.sleep(1)
  }
}

write_xlsx(df,output_file)

cat("Sequence annotation completed.\n")