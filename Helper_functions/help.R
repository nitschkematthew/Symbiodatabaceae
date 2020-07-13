
# Convert a DNAStringSet into a data.frame to use with tidyverse functions
DNAStringSet_to_df <- function(DNAStringSet){
  seq_df <- data.frame(names = names(DNAStringSet),
                       seqs = paste(DNAStringSet))
  return(seq_df)
}

# Convert a data.frame (seq_df) into a DNAStringSet to use with Biostrings functions
df_to_DNAStringset <- function(df, seqs = "seqs", names = "names"){
  DNAstr <- DNAStringSet(df[[seqs]])
  names(DNAstr) <- df[[names]]
  return(DNAstr)
}

# FIX THIS
# Convert data.frame (seq_df) into a DNAbin object to use with ape/insect functions
# df_to_DNAbin <- function(seq_df, seqs = "seqs", names = "names"){
#   DNA <- as.character(seq_df[[seqs]])
#   names(DNA) <- seq_df[[names]]
#   return(DNA)
# }

# Convert DNAbin to DNAStringSet to use with Biostrings functions
DNAbin_to_DNAStringSet <- function(DNAbin){
  DNAstrings <- DNAbin %>% as.character() %>% lapply(.,paste0, collapse="") %>% unlist() %>% DNAStringSet()
  return(DNAstrings)
}

# Convert DNAstringset to DNAbin
DNAStringSet_to_DNAbin <- function(DNAStringSet){
  DNAbin <- as.DNAbin(DNAStringSet)
  return(DNAbin)
}

# Convert DNAbin to data.frame (seq_df) to use with tidyverse functions
DNAbin_to_df <- function(DNAbin){
  seq_df <- data.frame(names = labels(DNAbin),
                       seqs = paste(DNAbin))
  return(seq_df)
}

# Read fasta file into df
fasta_to_df <- function(filepath){
  fasta_stringset <- readDNAStringSet(filepath)
  seq_df <- DNAStringSet_to_df(fasta_stringset) 
  return(seq_df)
}

# Write seq_df to .fasta file

# Extract accession numbers (assuming each accession finishes with a .1) from DNAstringset to use in database merging
extract_accn <- function(DNAStringSet){
  names <- names(DNAStringSet)
  accn <- stringi::stri_extract(names, regex='[^.]*')
  return(accn)
}

# Dereplicate a DNAStringSet and concatenate the names together with | as the delimiter
derep_cat_names <- function(dnaSet){
  if (!is(dnaSet, "DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  if (is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), 
                           sep = "-")
  }
  
  seq_df <- DNAStringSet_to_df(dnaSet) %>%
    distinct(names, .keep_all = TRUE) %>% # Remove any duplicated IDs
    mutate(accession = str_sub(names, start = -8)) %>%
    distinct(accession, .keep_all = TRUE) %>% # Remove and duplicated accession numbers
    ungroup() %>%
    group_by(seqs) %>%
    summarise(names = paste(names, collapse = '|')) %>% # Concatenate identical sequences
    ungroup()
  
  seq_ss <- df_to_DNAStringset(seq_df, "seqs", "names")
  return(seq_ss)
}

# Function to retrieve the query from genbank nucleotide database in fasta format using accession numbers. 
# If the accession number no longer exists in the database, include the accession in the output as 404NotFound
fetch_seqs <- function(query){
  recs <- try(entrez_fetch(db = "nuccore", id = query, rettype = "fasta"))
  if(str_detect(recs[1], "Error|Fail")){
    recs <- paste0(">", query, "_404NotFound\nAAAAAAAAAA\n\n")
    return(recs)
  }
  else{
    return(recs)
  }
}

# Run virtualPCR from insect in a loop using multiple primer pairs
virtualPCR_multiplex <- function(DNAbin, fwd_list, rev_list){
  if(length(fwd_list) != length(rev_list)){message("Error: Must have equal number of primer pair combinations")}
  else{
    PCR <- list()
    for(i in 1:length(fwd_list)){
      trim <- virtualPCR(DNAbin, up = fwd_list[i], down = rev_list[i], 
                         minfsc = 50, minrsc = 50, rcdown = FALSE, trimprimers = FALSE)
      trim <- DNAbin_to_DNAStringSet(trim)
      PCR[[i]] <- trim
      message(paste("Finished round ", i))
    }
    return(PCR)
  }
}
