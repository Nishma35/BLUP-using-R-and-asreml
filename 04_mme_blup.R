############################################################
# 04 - Mixed model equations (MME) by hand
# Build X, Z, G, R, solve the MME for EBVs, compute
# PEV / SE / accuracy, and check model fit
############################################################

# Set working directory and load libraries
setwd("C:/Users/nisma/OneDrive/Desktop/30 herd angus")
library(pedigreemm)
library(Matrix)

# Read and prepare pedigree
ped <- read.csv("10herdpedigree.csv", header = TRUE, na.strings = "")
ped$SIRE <- as.character(ped$SIRE)
ped$DAM  <- as.character(ped$DAM)
ped$ID   <- as.character(ped$ID)

pd_edited <- editPed(sire = ped$SIRE, dam = ped$DAM, label = ped$ID)
pd_short  <- pedigreemm::pedigree(sire  = pd_edited$sire,
                                  dam   = pd_edited$dam,
                                  label = pd_edited$label)

# Get relationship matrix
A <- getA(pd_short)

# Read and prepare phenotype data
phen <- read.csv("10herdanimals.csv")
phen$herd_y    <- interaction(phen$herd, phen$year)
phen$animal_id <- as.character(phen$animal_id)
phen$month     <- as.factor(phen$month)
phen$sl        <- as.factor(phen$sl)

# Remove NA values from response variable
phen <- phen[!is.na(phen$sl_cl), ]

# ---- CORRECT Variance Components ----
sig2A <- 0.2731423   # vm(animal_id, ainv) from ASReml
sig2H <- 0.4470498   # herd_y!herd_y from ASReml
sig2E <- 0.4838579   # units!units from ASReml <- ONLY THIS!

# ---- CORRECT Design Matrices ----
# Response
y <- as.matrix(phen$sl_cl)
X <- model.matrix(~ month + sl + sl_age, data = phen)

# Animal matrix (Z1) - CORRECTED
animal_levels <- colnames(A)
Z1 <- matrix(0, nrow(phen), ncol(A))
for (i in 1:nrow(phen)) {
  idx <- which(animal_levels == phen$animal_id[i])
  if (length(idx) > 0) Z1[i, idx] <- 1
}

# Herd-year matrix (Z2) - CORRECTED
herd_years <- unique(phen$herd_y)
Z2 <- matrix(0, nrow(phen), length(herd_years))
for (i in 1:nrow(phen)) {
  idx <- which(herd_years == phen$herd_y[i])
  Z2[i, idx] <- 1
}

# Combine Z
Z <- cbind(Z1, Z2)

######### checking if Z1 is correct #################
setdiff(phen$animal_id, animal_levels)
which(rowSums(Z1) == 0)
table(rowSums(Z1))
head(animal_levels)
head(pd_edited$label)
i <- 20
phen$animal_id[i]
which(Z1[i, ] == 1)
animal_levels[which(Z1[i, ] == 1)]

# ---- CORRECT G matrix ----
# G = [A*σ²_A,  0      ]
#     [ 0,      I*σ²_H ]
G_animal <- sig2A * A
G_herd   <- sig2H * diag(length(herd_years))   # NOT m2!
G <- bdiag(G_animal, G_herd)

# ---- CORRECT R matrix ----
R <- sig2E * diag(nrow(phen))

# ---- CORRECT Inverses ----
R_inv <- (1 / sig2E) * diag(nrow(phen))   # Simple inverse for diagonal R
G_inv <- solve(G)

# ---- Build MME ----
# Left side: [X'R⁻¹X   X'R⁻¹Z      ]
#            [Z'R⁻¹X   Z'R⁻¹Z + G⁻¹]
C11 <- t(X) %*% R_inv %*% X
C12 <- t(X) %*% R_inv %*% Z
C21 <- t(Z) %*% R_inv %*% X
C22 <- t(Z) %*% R_inv %*% Z + G_inv
C <- rbind(cbind(C11, C12),
           cbind(C21, C22))

# Right side: [X'R⁻¹y]
#             [Z'R⁻¹y]
rhs1 <- t(X) %*% R_inv %*% y
rhs2 <- t(Z) %*% R_inv %*% y
rhs  <- rbind(rhs1, rhs2)

# ---- Solve ----
solution <- solve(C, rhs)

# ---- Extract results ----
n_fixed      <- ncol(X)
n_animals    <- ncol(A)
n_herd_years <- length(herd_years)
hist
fixed_effects  <- solution[1:n_fixed]
animal_effects <- solution[(n_fixed + 1):(n_fixed + n_animals)]
herd_effects   <- solution[(n_fixed + n_animals + 1):(n_fixed + n_animals + n_herd_years)]
hist(animal_effects)

# Convert vector to data frame
animal_effects_df <- data.frame(AnimalEBV = animal_effects)

# Save to CSV
write.csv(animal_effects_df, "animal_effects.csv", row.names = FALSE)

# 1. Calculate Accuracies from your MME results
# PEV = diagonal of inverse MME matrix for animal effects
C_inv <- solve(C)

# Positions of animal effects in the inverse
animal_start <- n_fixed + 1
animal_end   <- n_fixed + n_animals
animal_pev   <- diag(C_inv)[animal_start:animal_end]

############ SE #############
animal_se <- sqrt(animal_pev)

# Calculate accuracy: sqrt(1 - PEV/σ²_A)
accuracy <- sqrt(1 - (animal_pev / sig2A))

# 2. Create simple EBV table
ebv_table <- data.frame(
  Animal_ID = animal_levels,
  EBV       = animal_effects,
  SE        = animal_se,
  Accuracy  = round(accuracy, 3)
)

# Save full table
write.csv(ebv_table, file = "EBV_table_full.csv", row.names = FALSE)

# Create a smaller table with just ID and EBV
ebv_simple <- ebv_table[, c("Animal_ID", "EBV")]

# Save it
write.csv(ebv_simple, file = "EBV_table_simple.csv", row.names = FALSE)

# 1. Calculate corrected phenotypes
corrected_pheno <- phen$sl_cl - X %*% fixed_effects
y <- phen$sl_cl

# 2. Match EBVs to phenotypes (same order as phen data)
ebv_for_pheno <- numeric(nrow(phen))
for (i in 1:nrow(phen)) {
  animal_idx <- which(animal_levels == phen$animal_id[i])
  ebv_for_pheno[i] <- animal_effects[animal_idx]
}
df <- ebv_for_pheno

df_ebv <- data.frame(
  Animal_ID       = phen$animal_id,   # from your phen data
  Phenotype       = phen$sl_cl,       # original phenotype
  Corrected_Pheno = corrected_pheno,  # phenotype adjusted for fixed effects
  EBV             = ebv_for_pheno     # EBVs aligned to phenotypes
)

# 3. Create the plot
plot(ebv_for_pheno, corrected_pheno,
     xlab = "Estimated Breeding Value (EBV)",
     ylab = "y",
     main = "EBVs vs Corrected Phenotypes",
     pch = 19, col = "blue", cex = 0.8)

# Add regression line
abline(lm(corrected_pheno ~ ebv_for_pheno), col = "red", lwd = 2)
cor_y <- cor(y, ebv_for_pheno)

# Add correlation
cor_value <- cor(ebv_for_pheno, corrected_pheno)
legend("topleft", legend = paste("r =", round(cor_y, 3)),
       bty = "n", cex = 1.2)

# 4. Save the plot
dev.copy(png, "EBV_vs_CorrectedPheno.png", width = 800, height = 600, res = 100)
dev.off()
hist(ebv_table$EBV)

#################################################################
table(rowSums(Z2))
which(rowSums(Z2) == 0)
set.seed(1)
rows <- sample(1:nrow(phen), 5)
cbind(
  phen$herd_y[rows],
  herd_years[apply(Z2[rows, ], 1, which.max)]
)
dim(Z)   # "nrow(phen) × (n_animals + n_herd_years)"
diag(G)[1:10]
diag(G)[n_animals + 1:n_animals + 10]
diag(G_inv)[1:10]
diag(G_inv)[(n_animals + 1):(n_animals + 10)]
max(abs(C - t(C)))
dim(C11)
dim(C12)
dim(C21)
dim(C22)
length(rhs)
nrow(C)
cor(ebv_for_pheno, corrected_pheno)

n_fixed   <- ncol(X)
n_random  <- ncol(Z)
n_animals <- length(animal_levels)   # 165? 300? whatever your dataset uses
n_hy      <- length(herd_years)

n_fixed
n_random
n_animals
n_hy

ebv_animal    <- solution[(n_fixed + 1):(n_fixed + n_animals)]
ebv_for_pheno <- Z1 %*% ebv_animal
cor(ebv_for_pheno, corrected_pheno)

# --- 1. Predicted phenotypes ---
# Extract herd-year effects
herd_year_effects <- solution[(n_fixed + n_animals + 1):(n_fixed + n_animals + n_hy)]

# Predicted components
Xb <- X %*% fixed_effects          # fixed effects
Zu <- Z1 %*% ebv_animal            # animal genetic effect
Zh <- Z2 %*% herd_year_effects     # herd-year effect

# Full predicted phenotype
y_hat <- Xb + Zu + Zh

# --- 2. Correlation between observed y and predicted y_hat ---
y <- phen$sl_cl
cor_y_yhat <- cor(y, y_hat)
cat("Correlation between y and y_hat:", round(cor_y_yhat, 3), "\n")

# --- 3. BLUP residuals and correlation with y_hat ---
residuals_blup <- y - y_hat
cor_yhat_resid <- cor(y_hat, residuals_blup)
cat("Correlation between y_hat and BLUP residuals:", round(cor_yhat_resid, 3), "\n")

# --- Predicted phenotypes (y_hat) ---
herd_year_effects <- solution[(n_fixed + n_animals + 1):(n_fixed + n_animals + n_hy)]
Xb <- X %*% fixed_effects
Zu <- Z1 %*% ebv_animal
Zh <- Z2 %*% herd_year_effects
y_hat <- Xb + Zu + Zh

# Observed phenotypes
y <- phen$sl_cl

# BLUP residuals
residuals_blup <- y - y_hat

# --- 1. Scatter plot: y vs y_hat ---
plot(y_hat, y,
     xlab = "Predicted Phenotype (y_hat)",
     ylab = "Observed Phenotype (y)",
     main = "Observed vs Predicted Phenotypes",
     pch = 19, col = "blue", cex = 0.8)
abline(a = 0, b = 1, col = "red", lwd = 2)   # 45-degree line
cor_y_yhat <- cor(y, y_hat)
legend("topleft", legend = paste("r =", round(cor_y_yhat, 3)), bty = "n", cex = 1.2)

# --- 2. Histogram of BLUP residuals ---
hist(residuals_blup,
     breaks = 30,
     main = "Histogram of BLUP Residuals",
     xlab = "Residuals",
     col = "lightblue", border = "white")

# --- 3. QQ plot of residuals ---
qqnorm(residuals_blup,
       main = "QQ Plot of BLUP Residuals")
qqline(residuals_blup, col = "red", lwd = 2)

# --- Optional: Correlation y_hat vs residuals (should be near 0) ---
cor_yhat_resid <- cor(y_hat, residuals_blup)
cat("Correlation between y_hat and BLUP residuals:", round(cor_yhat_resid, 3), "\n")

# Scatter plot
plot(y_hat, residuals_blup,
     xlab = "Predicted values (y_hat)",
     ylab = "BLUP residuals",
     main = "BLUP residuals vs Predicted values",
     pch = 16)

# Horizontal reference line at 0
abline(h = 0, lty = 2)

# Add correlation text
cor_val <- cor(y_hat, residuals_blup)
legend("topright",
       legend = paste("Correlation =", round(cor_val, 3)),
       bty = "n")
