library(glue)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(splitstackshape)
library(readr)
library(purrr)
library(data.table)

#Files needed for processing
name_key <- read.delim("wormbase_name_key.txt") %>%
  select(WormBase.Gene.ID, Public.Name)
AA_SCORES <- readr::read_tsv("AA_Scores.tsv")

order_annotation <- function(df){
  df[with(df, order(CHROM, POS)), ]

} #Will reorder tables after binding rows
prn_transcript <- function(df){
  print(unique(df$TRANSCRIPT)) %>%
    toString()} #Used to list transcripts after grouping
parse_VCF <- function(df){
  parsed <- df %>%
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3","ALT" = "V4", "SAMPLE" = "V5",  "ANNOTATION" = "V6") %>%
    tidyr::separate_rows(ANNOTATION, sep=",")%>%
    tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|")
}
parse_BCSQ <- function(df){
  parsed <- df %>% rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4", "SAMPLE" = "V5", "ANNOTATION" = "V6"  ) %>% #Separate CHROM, POS, ANNOTATION
    cSplit("ANNOTATION", sep = ",", direction = "wide") #Separate Transcripts
  renamecols <- function(x){colnames(x) <- gsub("ANNOTATION", "SCRIPT", colnames(x)); x}
  pivot <- renamecols(parsed) %>%
    pivot_longer(cols = starts_with("SCRIPT"), names_to = "TRANSCRIPT", values_to = "ANNOTATION") %>% na.omit() %>%
    separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|")%>%
    as_tibble()

  return(pivot)
} #Depreciated
impact_scoring <- function(df) {  #Works for up to two multi consequences - can handle 3 but the highest may not always be resulted
  multi_con<- df %>% dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>%
    filter(.$multi_con == TRUE) %>% select(!multi_con) %>% as_tibble()  #Filter for only multiconsequnce variants

  refromated_multi <- multi_con %>%
    separate(CONSEQUENCE, sep = "&", into = c("CON1", "CON2"), remove= FALSE) %>% #Separate consequnces (Works for 2 ATM)
    mutate(CON1 = sapply(.$CON1, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    mutate(CON2 = sapply(.$CON2, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    mutate(VARIANT_IMPACT = pmax(.$CON1, .$CON2)) %>% #Select the highest impact to be variant impact
    mutate(VARIANT_IMPACT = sapply(VARIANT_IMPACT, impact_numeric_tocharacter)) %>% #Convert from numeric to character
    mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character)) %>% select(!CON1, !CON2)

  single_con <- df %>% dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>%
    filter(.$multi_con == FALSE) %>% dplyr::mutate(VARIANT_IMPACT = sapply(.$CONSEQUENCE, impact)) %>%
    select(!multi_con) %>% mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character))

  clean <- bind_rows(single_con , refromated_multi) %>%
    order_annotation()
  return(clean)
} #Uses Below functions to handle impact conversion & Multi-consequence variants
impact <- function(x) {
  dplyr::recode(x,

                "missense" ="MODERATE",
                "synonymous"="LOW",
                "stop_lost"= "HIGH",
                "stop_gained"="HIGH",
                "inframe_deletion"="MODERATE",
                "inframe_insertion"="MODERATE",
                "frameshift"="HIGH",
                "splice_acceptor"="HIGH",
                "splice_donor"="HIGH",
                "start_lost"="HIGH",
                "splice_region"="MODERATE",
                "stop_retained"="LOW",
                "5_prime_utr"="MODIFIER",
                "3_prime_utr"="MODIFIER",
                "non_coding"="MODIFIER",
                "intron"="MODIFIER",
                "intergenic"="MODIFIER",
                "inframe_altering"="MODERATE",
                "coding_sequence"="MODIFIER",
                "feature_elongation"="MODIFIER",
                "start_retained"="LOW" ,
                "*missense" ="MODERATE",
                "*synonymous"="LOW",
                "*stop_lost"= "HIGH",
                "*stop_gained"="HIGH",
                "*inframe_deletion"="MODERATE",
                "*inframe_insertion"="MODERATE",
                "*frameshift"="HIGH",
                "*splice_acceptor"="HIGH",
                "*splice_donor"="HIGH",
                "*start_lost"="HIGH",
                "*splice_region"="MODERATE",
                "*stop_retained"="LOW",
                "*5_prime_utr"="MODIFIER",
                "*3_prime_utr"="MODIFIER",
                "*non_coding"="MODIFIER",
                "*intron"="MODIFIER",
                "*intergenic"="MODIFIER",
                "*inframe_altering"="MODERATE",
                "*coding_sequence"="LOW",
                "*feature_elongation"="LOW",
                "*start_retained"="LOW")

} #Converts Consequence to Impact
impact_numeric <- function(x) {
  recode(x,
         "missense" =3,
         "synonymous"=1,
         "stop_lost"= 4,
         "stop_gained"=4,
         "inframe_deletion"=3,
         "inframe_insertion"=3,
         "frameshift"=4,
         "splice_acceptor"=4,
         "splice_donor"=4,
         "start_lost"=4,
         "splice_region"=3,
         "stop_retained"=1,
         "5_prime_utr"=2,
         "3_prime_utr"=2,
         "non_coding"=2,
         "intron"=2,
         "intergenic"=2,
         "inframe_altering"=3,
         "coding_sequence"=2,
         "feature_elongation"=2,
         "start_retained"=1 ,
         "*missense" =3,
         "*synonymous"=1,
         "*stop_lost"= 4,
         "*stop_gained"=4,
         "*inframe_deletion"=3,
         "*inframe_insertion"=3,
         "*frameshift"=4,
         "*splice_acceptor"=4,
         "*splice_donor"=4,
         "*start_lost"=4,
         "*splice_region"=3,
         "*stop_retained"=1,
         "*5_prime_utr"=2,
         "*3_prime_utr"=2,
         "*non_coding"=2,
         "*intron"=2,
         "*intergenic"=2,
         "*inframe_altering"=3,
         "*coding_sequence"=2,
         "*feature_elongation"=2,
         "*start_retained"=1)

}  #Handles multi-consequence with numneric to character
impact_numeric_tocharacter <- function(x){
  recode(x,
         "4"="HIGH",
         '3'="MODERATE",
         "2"="LOW",
         '1'="MODIFIER")
}

BCSQ_p_translate <- function(df){ #Fails for AA changes that go more than one position
  translated <- df %>%
    dplyr::mutate("AA" = stringr::str_extract(df$AMINO_ACID_CHANGE, "[A-Z]"))%>% #Grabs the first AA
    dplyr::mutate("ALT_AA" = stringr::str_sub(df$AMINO_ACID_CHANGE,-1)) %>% #Grabs the ALT AA
    dplyr::mutate("AA_POS" = stringr::str_extract(df$AMINO_ACID_CHANGE, "([0-9])+")) %>% #Grabs the AA position
    tidyr::unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|") #Unites AA, ALT
  return(translated)
}

args <- commandArgs(trailingOnly = TRUE)
myarg <- args[1]

#Read-In table
data <- read.table(glue::glue("{myarg}_BCSQ.vcf"), col.names = c("V1", "V2", "V3", "V4", "V5", "V6" ), fill = TRUE, na.strings=c("","NA"))

#Mark sites without annotation
data$V6[is.na(data$V6)] <- "No_Annotation"

#Split unannotated and annotated into two dataframes
non_annotation <- data %>% mutate(NO_ANN = str_detect(.$V6, pattern = "No_Annotation" )) %>% filter(NO_ANN == TRUE) %>% select(!NO_ANN) %>%
  rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4", "SAMPLE" = "V5", "ANNOTATION" = "V6"  )
annotated <- data %>%  mutate(NO_ANN = str_detect(.$V6, pattern = "No_Annotation" )) %>% filter(NO_ANN == FALSE) %>% select(!NO_ANN)


#Parse with new parse function
# parsed <- parse_BCSQ(data)
# #
# reformated <- parsed %>%
#   mutate(linker = str_detect(.$CONSEQUENCE, pattern = "@" )) %>%
#   filter(.$linker != TRUE | is.na(.$linker) ) %>% select(!linker) %>%
#   left_join(name_key, by = c( "GENE" = "WormBase.Gene.ID")) %>% #Convert Gene from Wormbase ID to public name
#   mutate(GENE = .$Public.Name) %>%
#   select(!Public.Name)
#

#MARK SITES WITHOUT ANNOTATIONS -

# reformated$CONSEQUENCE[is.na(reformated$V6)] <- "No_Annotation"

#SPLIT INTO TWO DF'S
# non_annotation <- reformated %>% mutate(NO_ANN = str_detect(.$CONSEQUENCE, pattern = "No_Annotation" )) %>% filter(NO_ANN == TRUE) %>% select(!NO_ANN)
# annotated <- reformated %>%  mutate(NO_ANN = str_detect(.$CONSEQUENCE, pattern = "No_Annotation" )) %>% filter(NO_ANN == FALSE) %>% select(!NO_ANN)


#ANNOTATED WORKFLOW
parsed <- parse_BCSQ(annotated) %>% mutate(linker = str_detect(.$CONSEQUENCE, pattern = "@" )) %>%
  filter(.$linker != TRUE | is.na(.$linker) ) %>% select(!linker) %>%
  left_join(name_key, by = c( "GENE" = "WormBase.Gene.ID")) %>% #Convert Gene from Wormbase ID to public name
  mutate(GENE = .$Public.Name) %>%
  select(!Public.Name)

grouped <- parsed %>% group_by(CHROM,POS,SAMPLE,REF,ALT,GENE,CONSEQUENCE,STRAND,DNA_CHANGE,AMINO_ACID_CHANGE) %>%
  nest() %>%
  mutate(TRANSCRIPTS =  as.character(map(data, prn_transcript))) %>% #Would love to remove std out from this
  select(-data) %>% as_data_frame()
#Add impact scores
impacted <- impact_scoring(grouped)
#Grantham & Blossum Scoring
G.Bscored <- BCSQ_p_translate(impacted) %>% left_join(AA_SCORES, by = "REF_ALT_AA" )



#UNANNOTATED WORKFLOW - Instead of creating consequence with no annotation make Varaint impact
step1 <- non_annotation %>% mutate(VARIANT_IMPACT = "No Annotation")



#REJOIN AND ORDER

table_ready<- bind_rows(G.Bscored , step1) %>% order_annotation() %>% select(!ANNOTATION) %>%
  select(!CON1) %>% select(!CON2) %>% select(!REF_ALT_AA) %>% select(!AA_POS)


#Work around to get Bscore and G Score Numeric
table_ready2 <- table_ready %>% mutate(GSCORE = sapply(table_ready$GSCORE, as.integer)) %>% mutate(BSCORE = sapply(table_ready$BSCORE, as.integer))

#Add in Numeric Proxy for NA - is replaced in SQL import
table_ready2$GSCORE[is.na(table_ready2$GSCORE)] <- -1
table_ready2$BSCORE[is.na(table_ready2$BSCORE)] <- 200



write.csv(table_ready2, (glue::glue("{myarg}.prepped.csv")))
