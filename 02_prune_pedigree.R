############################################################
# 02 - Pruning pedigree
# Strict pruning with ggroups::pruneped, then extract
# animals + parents + grandparents
############################################################

library(ggroups)
library(dplyr)
library(kinship2)

setwd("C:/Users/nisma/OneDrive/Desktop/30 herd angus")

# Read files
ped    <- read.csv("pedigree data.csv")
target <- read.csv("10herdanimals.csv")

# Make sure IDs match type
ped$animal_id <- as.character(ped$animal_id)
ped$sire_id   <- as.character(ped$sire_id)
ped$dam_id    <- as.character(ped$dam_id)

phenotyped_ids <- target$animal_id
# adjust name if needed
colnames(target)
colnames(ped)

ped_clean <- ped %>%
  mutate(
    animal_id = trimws(animal_id),
    sire_id   = trimws(sire_id),
    dam_id    = trimws(dam_id)
  ) %>%
  mutate(
    animal_id = as.character(animal_id),
    sire_id   = as.character(sire_id),
    dam_id    = as.character(dam_id)
  ) %>%
  select(animal_id, sire_id, dam_id) %>%
  rename(
    ID   = animal_id,
    SIRE = sire_id,
    DAM  = dam_id
  )

ped_pruned_strict <- pruneped(ped_clean, phenotyped_ids, mode = "strict")
head(ped_pruned_strict)
nrow(ped_pruned_strict)

########## Make sure all IDs in target are present in ped_pruned_strict$ID ##########
setdiff(target$ID, ped_pruned_strict$ID)
setdiff(target$animal_id, ped_clean$ID)

######### Sires not found in pedigree:
setdiff(na.omit(ped_clean$SIRE), ped_clean$ID)

############### anyDuplicated(ped_pruned_strict$ID) ##############
anyDuplicated(ped_pruned_strict$ID)

write.csv(ped_pruned_strict, "newpruning10herd.csv", row.names = FALSE)

PED1 <- read.csv("newpruning10herd.csv")
df   <- read.csv("10herdanimals.csv", stringsAsFactors = FALSE)

# Replace all negative numbers with NA
PED1[PED1 < 0] <- NA
head(df)
write.csv(PED1, "newpruning10herdcleanNA.csv", row.names = FALSE)

library(tidyverse)
library(kinship2)

# Read your pruned pedigree file
# Replace 'pedigree_prune_strict_NA.csv' with your actual file name
peda <- read.csv("newpruning10herdcleanNA.csv")

# Check the structure of your pedigree data
head(peda)
str(peda)

# Assuming your pedigree has columns: animal_id, sire_id, dam_id
# Extract 3-generation pedigree function

### Load focal animals (29,963 Angus)
focal <- read.csv("10herdanimals.csv", stringsAsFactors = FALSE)

### Extract the column containing the animal IDs
focal_ids <- unique(focal$animal_id)   # <-- change column name if needed

# FUNCTION: get 3 generations (animal, parents, grandparents)
# --- Build 2-generation pedigree (animals + parents + grandparents) ---
sel_ids <- unique(as.character(focal$animal_id))   # Identifying our 400 selected animals

parents <- peda %>%
  filter(ID %in% sel_ids) %>%     # find pedigree records for selected animals
  select(SIRE, DAM) %>%           # take only sire and dam columns
  unlist() %>%                    # flatten the two columns into one long vector
  unique() %>%                    # remove duplicates
  setdiff(c(NA, "", "0", "-"))    # drop blanks, missing or "0" placeholders

grandparents <- peda %>%
  filter(ID %in% parents) %>%     # find pedigree records for the parents
  select(SIRE, DAM) %>%
  unlist() %>%
  unique() %>%
  setdiff(c(NA, "", "0", "-"))

keep_ids <- unique(c(sel_ids, parents, grandparents))   # Here we combine: the selected animals, their parents
                                                        # their grandparents

ped_small <- peda %>% filter(ID %in% keep_ids)   # Filter the pedigree to just 200 animals, their parents and grandparents
cat("Trimmed pedigree size (2 generations):", nrow(ped_small), "\n")

ped_3gen <- get_3gen(peda, focal_ids)

### Check result
nrow(ped_3gen)
head(ped_3gen)

write.csv(ped_small, "10herdpedigree.csv", row.names = FALSE)
pedddd <- read.csv("10herdpedigree.csv")
