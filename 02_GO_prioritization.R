# ============================================================
# 02_GO_prioritization.R
#
# GO-based prioritization of ZER1-interacting proteins
# ============================================================

library(readxl)
library(writexl)
library(dplyr)
library(org.Mm.eg.db)
library(GO.db)

input_file <- "results/01_sequence_annotation.xlsx"

df <- read_excel(input_file)

valid_genes <- df %>%
  filter(Met_cleavage_predicted) %>%
  pull(`Gene Symbol`) %>%
  na.omit()

info_df <- suppressWarnings(
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys = valid_genes,
    columns = c("GENENAME","GO"),
    keytype = "SYMBOL"
  )
)

all_go_ids <- unique(info_df$GO[!is.na(info_df$GO)])

go_terms_df <- suppressWarnings(
  AnnotationDbi::select(
    GO.db,
    keys = all_go_ids,
    columns = "TERM",
    keytype = "GOID"
  )
)

keywords <- list(
  
  ECM_Collagen =
    "extracellular matrix|extracellular structure|collagen|cell-matrix adhesion|wound healing|tissue remodeling",
  
  Muscle_Structure =
    "cardiac muscle hypertrophy|muscle tissue development|sarcomere|myofibril|actin filament|heart|cytoskeleton|Wnt",
  
  Mech_Response =
    "mechanical stimulus|response to pressure|response to stretch|shear stress",
  
  Stress_Survival =
    "apoptotic|oxidative stress|autophagy"
)

inspection_list <- list()

for(g in valid_genes){
  
  sub_info <- info_df[info_df$SYMBOL==g,]
  
  gene_go_terms <- unique(
    go_terms_df$TERM[
      go_terms_df$GOID %in% unique(sub_info$GO)
    ]
  )
  
  gene_go_terms <- gene_go_terms[!is.na(gene_go_terms)]
  
  for(cat_name in names(keywords)){
    
    hits <- grep(
      keywords[[cat_name]],
      gene_go_terms,
      ignore.case=TRUE,
      value=TRUE
    )
    
    if(length(hits)>0){
      
      inspection_list[
        [length(inspection_list)+1]
      ] <- data.frame(
        
        Gene=g,
        Category=cat_name,
        Matched_Term=hits
      )
    }
  }
}

inspection_df <- do.call(rbind,inspection_list)

final_hits_df <- inspection_df %>%
  filter(
    !grepl(
      "lymphocyte|B cell|T cell|neuron|glial",
      Matched_Term,
      ignore.case=TRUE
    )
  )

match_results <- final_hits_df %>%
  group_by(Gene) %>%
  summarize(
    Broad_Category=
      paste(unique(Category),collapse="; "),
    Pathway_Label=
      paste(unique(Matched_Term),collapse=" | "),
    .groups="drop"
  ) %>%
  mutate(Tier_Class=1)

df_all <- df %>%
  left_join(
    match_results,
    by=c("Gene Symbol"="Gene")
  ) %>%
  mutate(
    Tier_Class=
      ifelse(is.na(Tier_Class),2,Tier_Class)
  )

table_s1 <- df_all %>%
  arrange(desc(log2FC))

table_s2 <- df_all %>%
  filter(Tier_Class==1) %>%
  arrange(desc(log2FC)) %>%
  mutate(
    Rank_within_Tier1=row_number()
  )

write_xlsx(
  
  list(
    
    Table_S1_All_Candidates=table_s1,
    
    Table_S2_Tier1_Candidates=table_s2
    
  ),
  
  "results/ZER1_Supplementary_Tables.xlsx"
)

write_xlsx(
  table_s2,
  "results/Tier1_candidates.xlsx"
)

cat("GO prioritization completed.\n")