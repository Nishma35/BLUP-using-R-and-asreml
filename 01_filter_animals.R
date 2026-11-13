############################################################
# 01 - Filtering animals
# Select Angus males, split herd/year codes, sample
# 10 herds x 100 animals, and derive ADG / AMG
############################################################

library(dplyr)
library(readr)

setwd("C:/Users/nisma/OneDrive/Desktop/30 herd angus")

# Read the dataset
df <- read_csv("tyr phenotype.csv")

# Filter rows:
# 1. Remove rows with "-" in ANY column
# 2. Keep only breed_sex = 231 (Angus male)
df_angus_male <- df %>%
  filter(across(-twin, ~ !grepl("^\\-", .))) %>%   # check '-' in all columns except 'twin'
  filter(breed_sex == 231) %>%                     # select Angus male
  select(-twin)                                    # remove 'twin' and 'month'
# optional: remove the column
# select Angus male

df1 <- df_angus_male %>%
  mutate(
    hyr = substr(hyrsl, 1, nchar(hyrsl) - 2),
    sl  = substr(hyrsl, nchar(hyrsl) - 1, nchar(hyrsl))
  ) %>%
  select(-hyrsl)

df2 <- df1 %>%
  mutate(
    herd = substr(herd_year, 1, nchar(herd_year) - 2),                  # everything except last 2 digits
    year = substr(herd_year, nchar(herd_year) - 1, nchar(herd_year))    # last 2 digits
  ) %>%
  select(-herd_year)

library(dplyr)

df_summary <- df2 %>%
  group_by(herd) %>%
  summarise(count_animals = n())   # n() counts number of animals in each herd

# Count how many herds have more than 50 animals
num_large_herds <- sum(df_summary$count_animals > 100)
print(num_large_herds)

library(dplyr)

# Step 2: Identify herds with more than 100 animals
df_summary <- df2 %>%
  group_by(herd) %>%
  summarise(count_animals = n()) %>%
  filter(count_animals > 100)

set.seed(123)   # optional for reproducibility
set.seed(123)
selected_herds <- df_summary %>%
  slice_sample(n = 10)
selected_herds

df_selected <- df2 %>%
  filter(herd %in% selected_herds$herd) %>%   # keep the 10 chosen herds
  group_by(herd) %>%
  slice_sample(n = 100) %>%
  ungroup()

df_selected %>% count(herd)
nrow(df_selected)

# Step 3: Sample 100 animals from each of these herds
sampled_animals <- selected_herds %>%
  filter(herd %in% df_summary$herd) %>%   # keep only the herds with >100 animals
  group_by(herd) %>%
  slice_sample(n = 100)                   # randomly select 100 animals per herd

# Step 4: View the result
head(df_selected)

# Optional: save the sampled dataset to a new CSV
write.csv(df_selected, "10herdanimals.csv", row.names = FALSE)

# Read your data
df <- read.csv("10herdanimals.csv", header = TRUE)

# Convert age in months to days (approx: 1 month = 30 days)
df$age_days <- df$sl_age * 30

# Calculate Average Daily Gain (ADG)
df$ADG <- df$sl_w / df$age_days   # kg per day

# Calculate Average Monthly Gain (AMG)
df$AMG <- df$sl_w / df$sl_age     # kg per month

write.csv(df, "10herdanimalsadg$amg.csv", row.names = FALSE)
