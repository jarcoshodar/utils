library(dplyr)
library(tidyr)
library(tibble)  # Add this for deframe()
#root with the xls files 
folder_path <- '/home/jarcos/clean/pkns/'

clean_names <- function(df) {
  names(df) <- gsub(" ", "_", names(df))  # Replace spaces with underscores
  names(df) <- gsub("[\"']", "", names(df))  # Remove quotes
  return(df)
}
filenames <- list.files(folder_path)

#the convention followed is interaction exports called int1, int2, int3 etc. 
#network expoerts called net1,net2, etc.
#presently issues with xls reading and the hyperlinks if files are unsolved
#you must open if your prefered excel equivalent and save as csv manually
allint <- data.frame()
allnet <- data.frame()
for (filename in filenames) {
  filepath <- file.path(folder_path, filename)
  
  if (startsWith(filename, 'int') && endsWith(filename, '.csv')) {
    df <- read.csv(filepath, 
                   skip = 2,
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE,
                   quote = "\"",
                   encoding = "UTF-8")
    
    df <- clean_names(df)  # Clean the column names
    allint <- bind_rows(allint, df)
    
  } else if (startsWith(filename, 'net') && endsWith(filename, '.csv')) {
    df <- read.csv(filepath,
                   skip = 2,
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE,
                   quote = "\"",
                   encoding = "UTF-8")
    
    df <- clean_names(df)  # Clean the column names
    allnet <- bind_rows(allnet, df)
  }
}
#adjust the script to use the correct column names (may change, typically can be change . or enclose in `)
temp <- data.frame(
  TF_alias = allint$`Network Object "FROM"`,
  interaction = allint$Effect,
  target = allint$`Network Object "TO"`,
  stringsAsFactors = FALSE
)
#we only keep activation and inhibition
temp <- temp[temp$interaction != 'Technical', ]

gene_dict <- allnet %>%
  group_by(Network_Object_Name) %>%      # Changed from Network.Object.Name
  summarise(Gene_Symbol = list(unique(Gene_Symbol))) %>%  # Changed from Gene.Symbol
  deframe()
#expand mappings for both columns simultaneously
expand_mappings <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      original_TF_alias = TF_alias,
      original_target = target,
      TF_alias_mappings = list(if (!is.null(gene_dict[[TF_alias]])) gene_dict[[TF_alias]] else TF_alias),
      target_mappings = list(if (!is.null(gene_dict[[target]])) gene_dict[[target]] else target)
    ) %>%
    unnest(TF_alias_mappings) %>%
    unnest(target_mappings) %>%
    ungroup() %>%
    mutate(
      TF_alias = TF_alias_mappings,
      target = target_mappings
    ) %>%
    select(-TF_alias_mappings, -target_mappings)
}

#apply the function to expand both 'TF_alias' and 'target' together
cat("\nExpanding 'TF_alias' and 'target' mappings together...\n")
before_expansion <- nrow(temp)
temp_expanded <- expand_mappings(temp)
after_expansion <- nrow(temp_expanded)
cat("Rows before expansion:", before_expansion, "\n")
cat("Rows after expansion:", after_expansion, "\n")
cat("Number of rows added due to expansion:", after_expansion - before_expansion, "\n")

#remove duplicates or they will cause trouble in grni building
temp_expanded <- temp_expanded %>% distinct()

cat("\nFinal dimensions of expanded temp dataframe:", dim(temp_expanded), "\n")

#check for multiple mappings to expand them
expanded_rows <- temp_expanded %>%
  group_by(original_TF_alias, interaction, original_target) %>%
  summarise(count = n()) %>%
  filter(count > 1)

if (nrow(expanded_rows) > 0) {
  cat("\nRows expanded due to multiple mappings:\n")
  print(expanded_rows)
}

#specific example: Check for 'AP-1'
AP1_rows <- temp_expanded %>%
  filter(
    original_TF_alias == 'AP-1' | TF_alias == 'AP-1' |
    original_target == 'AP-1' | target == 'AP-1'
  )

if (nrow(AP1_rows) > 0) {
  cat("\nRows involving 'AP-1' after mapping:\n")
  print(AP1_rows)
} else {
  cat("\nNo rows involving 'AP-1' found after mapping.\n")
}
