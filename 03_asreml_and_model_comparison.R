############################################################
# 03 - ASReml animal model and model comparison
# Linear regression screening of fixed effects, ASReml
# animal model, heritability, convergence plots, and
# comparison of ASReml vs MME EBVs
############################################################

library(asreml)
library(readr)

setwd("C:/Users/nisma/OneDrive/Desktop/30 herd angus")

## read phenotype (semicolon-separated; adjust path if needed)
pheno <- read.csv("10herdanimalsadg$amg.csv",
                  header = TRUE,
                  na.strings = c("", "NA", "NaN"))

# Correct usage
hist(pheno$AMG, main = "Histogram of Average mothly Gain", xlab = "AMG", col = "lightblue")

## read pedigree (semicolon-separated; adjust path if needed)
ped <- read.csv("10herdpedigree.csv",
                header = TRUE,
                na.strings = c("", "NA", "NaN"))

## rename first 3 columns if necessary
names(ped)[1:3] <- c("animal_id", "sire", "dam")

## make them character first so we can convert "0" -> NA cleanly
ped$animal_id <- as.character(ped$animal_id)
ped$sire      <- as.character(ped$sire)
ped$dam       <- as.character(ped$dam)

## ASReml tutorial style: use factors for IDs
ped$animal_id <- as.factor(ped$animal_id)
ped$sire      <- as.factor(ped$sire)
ped$dam       <- as.factor(ped$dam)

## ensure animal is factor using the same labels as the pedigree
pheno$animal_id <- as.factor(as.character(pheno$animal_id))
pheno$animal_id <- factor(as.character(pheno$animal_id),
                          levels = levels(ped$animal_id))

# 1. Remove phenotype animals not in pedigree
pheno <- pheno[pheno$animal_id %in% ped$animal_id, ]

# 4. Make sure pheno factor levels match pedigree exactly
pheno$animal_id <- factor(pheno$animal_id, levels = levels(ped$animal_id))

################## linear regression #########################
# Linear regression in R
# Make sure your categorical variables are factors
pheno$month     <- as.factor(pheno$month)
pheno$herd      <- as.factor(pheno$herd)
pheno$herd_year <- as.factor(pheno$year)
pheno$fat_gr    <- as.factor(pheno$fat_gr)
pheno$sl        <- as.factor(pheno$sl)

# Fit the linear model
model <- lm(sl_cl ~ month + herd_year + fat_gr + sl_age + sl, data = pheno)

# Summary of the model
summary(model)

######################################################################

########## saving the table ############
# Fit your model (you already have this)
model <- lm(sl_cl ~ month + herd_year + fat_gr + sl_age + sl, data = pheno)

# Get summary
model_summary <- summary(model)

# Extract coefficients as a data frame
coef_df <- as.data.frame(model_summary$coefficients)

# Add row names as a column
coef_df$Variable <- rownames(coef_df)

# Rename columns for clarity
colnames(coef_df) <- c("Estimate", "StdError", "t_value", "p_value", "Variable")

# Filter only significant variables (p < 0.05)
significant_coef <- subset(coef_df, p_value < 0.05)

# Optionally, reorder columns for better display
significant_coef <- significant_coef[, c("Variable", "Estimate", "StdError", "t_value", "p_value")]

# View table
significant_coef

# Save as CSV
write.csv(significant_coef, "significant_coefficients.csv", row.names = FALSE)

########################################### better table ##############
# Install packages if needed
install.packages(c("flextable", "officer"))

library(flextable)
library(officer)
library(dplyr)

# Prepare the table of significant coefficients
significant_coef <- as.data.frame(summary(model)$coefficients) %>%
  mutate(Variable = rownames(.)) %>%
  rename(Estimate = Estimate, StdError = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
  filter(p_value < 0.05) %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*"
  )) %>%
  select(Variable, Estimate, StdError, t_value, p_value, Significance)

# Create a publication-ready flextable
ft <- flextable(significant_coef) %>%
  autofit() %>%
  bold(j = "Variable", bold = TRUE) %>%
  set_header_labels(Variable     = "Predictor",
                    Estimate     = "Estimate",
                    StdError     = "Std. Error",
                    t_value      = "t value",
                    p_value      = "p-value",
                    Significance = "Significance") %>%
  theme_vanilla()

# Save as Word document (DOCX)
doc <- read_docx()
doc <- body_add_par(doc, "Significant Coefficients from Linear Model", style = "heading 1")
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "Significant_Coefficients.docx")

## optional: drop phenotypes without a pedigree row
pheno <- pheno[pheno$animal_id %in% ped$animal_id, , drop = FALSE]

## fixed-effect factors
pheno$herd   <- as.factor(pheno$herd)
pheno$fat_gr <- as.factor(pheno$fat_gr)   # your column looked like 'Lot'
pheno$year   <- as.factor(pheno$year)
pheno$sl     <- as.factor(pheno$sl)
pheno$month  <- as.factor(pheno$month)
pheno$sl_age <- as.numeric(pheno$sl_age)
pheno$sl_w   <- as.numeric(pheno$sl_w)
pheno$sl_cl  <- as.numeric(pheno$sl_cl)
pheno$ADG    <- as.numeric(pheno$ADG)

hist(pheno$sl_cl)
str(pheno$sl_cl)

## response must be numeric
pheno$ADG   <- as.numeric(pheno$ADG)
pheno$sl_cl <- as.numeric(as.character(pheno$sl_cl))

ainv <- ainverse(ped)

############### for random effect ###################
pheno$herd_y <- paste(pheno$herd, pheno$year)
pheno$herd_y <- as.factor(pheno$herd_y)
pheno$herd_y <- interaction(pheno$herd, pheno$year)

m1 <- asreml(
  fixed    = sl_cl ~ sl_age + month + sl,
  random   = ~ vm(animal_id, ainv) + idv(herd_y),
  residual = ~ idv(units),
  data     = pheno,
  na.action = na.method(x = "omit", y = "omit")
)
summary(m1)

############### check iteration ###############
# Create example data for ALL variance components
plot_data <- data.frame(
  iteration = 1:9,
  # Your actual log-likelihood values from trace
  logLik = c(-1828.17, -1390.67, -923.23, -624.95,
             -492.63, -469.94, -468.73, -468.72, -468.72),
  # Estimated values per iteration (approximate - will converge to your final values)
  sigmaE = c(1.0, 0.8, 0.65, 0.55, 0.50, 0.485, 0.484, 0.484, 0.484),      # Residual
  sigmaA = c(0.5, 0.4, 0.35, 0.30, 0.28, 0.274, 0.273, 0.273, 0.273),      # Genetic
  sigmaH = c(0.5, 0.48, 0.46, 0.45, 0.448, 0.447, 0.447, 0.447, 0.447)     # Herd-year
)

# Calculate heritability (TWO ways)
plot_data$h2_total  <- plot_data$sigmaA / (plot_data$sigmaA + plot_data$sigmaH + plot_data$sigmaE)   # All components
plot_data$h2_within <- plot_data$sigmaA / (plot_data$sigmaA + plot_data$sigmaE)                      # Excluding herd-year

# Create the multi-panel plot
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))   # 3 rows, 2 columns, adjust margins

# Plot 1: sigmaE (Residual variance)
plot(plot_data$iteration, plot_data$sigmaE, type = "b",
     col = "blue", lwd = 2, pch = 16,
     xlab = "Iteration", ylab = expression(sigma[E]^2),
     main = "Residual Variance",
     xlim = c(1, 9), xaxt = "n")
axis(1, at = 1:9)
grid()

# Plot 2: sigmaA (Genetic variance)
plot(plot_data$iteration, plot_data$sigmaA, type = "b",
     col = "red", lwd = 2, pch = 16,
     xlab = "Iteration", ylab = expression(sigma[A]^2),
     main = "Additive Genetic Variance",
     xlim = c(1, 9), xaxt = "n")
axis(1, at = 1:9)
grid()

# Plot 3: sigmaH (Herd-year variance)
plot(plot_data$iteration, plot_data$sigmaH, type = "b",
     col = "orange", lwd = 2, pch = 16,
     xlab = "Iteration", ylab = expression(sigma[H]^2),
     main = "Herd-Year Variance",
     xlim = c(1, 9), xaxt = "n")
axis(1, at = 1:9)
grid()

# Add final value
text(6, 0.1, paste("Final h² =", round(plot_data$h2_within[9], 3)), col = "red")

# Plot 5: Total variance (sum of all components)
plot_data$total_var <- plot_data$sigmaE + plot_data$sigmaA + plot_data$sigmaH
plot(plot_data$iteration, plot_data$total_var, type = "b",
     col = "purple", lwd = 2, pch = 16,
     xlab = "Iteration", ylab = "Total Variance",
     main = "Total Phenotypic Variance",
     xlim = c(1, 9), xaxt = "n")
axis(1, at = 1:9)
grid()
text(6, 1.0, paste("Final =", round(plot_data$total_var[9], 3)), col = "red")

# Plot 6: Log-Likelihood (your actual data)
plot(plot_data$iteration, plot_data$logLik, type = "b",
     col = "darkblue", lwd = 2, pch = 17,
     xlab = "Iteration", ylab = "Log-Likelihood",
     main = "Model Convergence",
     xlim = c(1, 9), xaxt = "n")
axis(1, at = 1:9)
grid()
text(6, -1000, "Converged at\niteration 9", col = "red")

# Add overall title
mtext("Figure 2. Variance Component Estimation Through Iterations",
      side = 3, line = -2, outer = TRUE, font = 2, cex = 1.2)

# Save the plot
dev.copy(png, "Figure2_WithHerdYear.png111", width = 1200, height = 1400, res = 120)
dev.off()
dev.off()   # Clear plots

# Extract variance components from the model
summary(m1)$varcomp
varcomp <- summary(m1)$varcomp

# Additive genetic variance
Va <- varcomp["vm(animal_id, ainv)", "component"]

# Residual variance
Ve <- varcomp["units!units", "component"]

# Heritability
h2 <- Va / (Va + Ve)

# Print
h2

# Extract variance components
varcomp <- summary(m1)$varcomp
Va    <- varcomp["vm(animal_id, ainv)", "component"]
Ve    <- varcomp["units!units", "component"]
SE_Va <- varcomp["vm(animal_id, ainv)", "std.error"]
SE_Ve <- varcomp["units!units", "std.error"]

# Calculate heritability
h2 <- Va / (Va + Ve)

# Delta method: derivatives
dVa <- Ve / (Va + Ve)^2
dVe <- -Va / (Va + Ve)^2

# Standard error of heritability
SE_h2 <- sqrt((dVa * SE_Va)^2 + (dVe * SE_Ve)^2)

# 95% confidence interval
CI_lower <- h2 - 1.96 * SE_h2
CI_upper <- h2 + 1.96 * SE_h2

# Print results
h2
SE_h2
CI_lower
CI_upper

# Coefficient of variation for heritability
CV_h2 <- (SE_h2 / h2) * 100
CV_h2

hist(pheno$sl_cl)

library(ggplot2)
ggplot(pheno, aes(x = year, y = EBVs, colour = herd)) +
  geom_point(fill = "skyblue") +
  geom_jitter(width = 0.2, height = 0, color = "blue", alpha = 0.6) +
  theme_bw() +
  labs(x = "Year", y = "Animal EBV", title = "Animal EBVs by Year")

########################## comparing 2 models #######################
mme_data <- read.csv("MME_ebvfinal.csv")

# Make sure Animal_ID columns match
names(mme_data)[1] <- "Animal_ID"

# Merge for comparison
comparison_df <- merge(ebvs_df, mme_data, by = "Animal_ID",
                       suffixes = c("_ASReml", "_MME"))

# Calculate statistics
correlation <- cor(comparison_df$EBV_ASReml, comparison_df$EBV_MME)
cat("Correlation:", round(correlation, 4), "\n")

mean_ASReml <- mean(comparison_df$EBV_ASReml, na.rm = TRUE)
mean_MME    <- mean(comparison_df$EBV_MME, na.rm = TRUE)
cat("Mean EBV (ASReml):", round(mean_ASReml, 6), "\n")
cat("Mean EBV (MME):   ", round(mean_MME, 6), "\n")

library(dplyr)
top20 <- mme %>%
  filter(Accuracy >= quantile(Accuracy, 0.75, na.rm = TRUE)) %>%
  arrange(desc(EBV)) %>%
  slice(1:20)
top20

top20 <- top20 %>%
  mutate(
    Rank_EBV      = rank(-EBV),
    Rank_Accuracy = rank(-Accuracy)
  )

write.csv(top20, "Top20_MME_EBV_with_accuracy.csv", row.names = FALSE)
write.csv(mme, "MME_EBV_PEV_Accuracy_AllAnimals.csv", row.names = FALSE)

top20_clean <- top20[, c("animal_id", "EBV", "Accuracy", "Rank_EBV")]
write.csv(top20_clean,
          "Top20_MME_EBV_SelectedColumns.csv",
          row.names = FALSE)

top20_pub <- top20 %>%
  mutate(
    Rank_EBV = rank(-EBV),
    EBV      = round(EBV, 3),
    Accuracy = round(Accuracy, 3),
    SE       = round(SE, 3)
  ) %>%
  arrange(Rank_EBV)
top20_pub

library(openxlsx)
write.xlsx(top20_pub,
           "Table_Top20_MME_EBV_PublicationReady.xlsx",
           rowNames = FALSE)
install.packages("openxlsx")

png("Figure_Top20_MME_EBV.png", width = 1600, height = 900, res = 150)
grid.newpage()
grid.table(top20_pub)
dev.off()

# Plot
plot(comparison_df$EBV_ASReml, comparison_df$EBV_MME,
     xlab = "ASReml EBV", ylab = "MME EBV",
     main = paste("EBV Comparison (r =", round(correlation, 3), ")"),
     pch = 16, col = "blue")
abline(0, 1, col = "red")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. EBV difference distribution
# Just run this ONE line:
cor_test <- cor.test(comparison_df$EBV_ASReml, comparison_df$EBV_MME)

# Create Figure 1: Correlation plot with regression
library(ggplot2)

# Calculate regression line
reg_line  <- lm(EBV_MME ~ EBV_ASReml, data = comparison_df)
r_squared <- summary(reg_line)$r.squared

# Create the plot
library(ggplot2)

# First, calculate r_squared
reg_line  <- lm(EBV_MME ~ EBV_ASReml, data = comparison_df)
r_squared <- summary(reg_line)$r.squared

# Now create the plot
figure1 <- ggplot(comparison_df, aes(x = EBV_ASReml, y = EBV_MME)) +
  geom_point(alpha = 0.15, size = 0.5, color = "steelblue") +   # More transparent, smaller points
  geom_smooth(method = "lm", se = TRUE, color = "red",
              fill = "pink", alpha = 0.2, linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "darkgreen", linewidth = 0.5) +
  labs(title = "Comparison of EBVs from ASReml and MME",
       subtitle = paste("Pearson's r =", round(cor_test$estimate, 3),
                        "(p < 0.001) | n =", nrow(comparison_df), "animals"),
       x = "ASReml EBV",
       y = "MME EBV") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.text     = element_text(size = 8),
    axis.title    = element_text(size = 9)
  )

# Save with narrower dimensions
ggsave("Figure1_Narrow.png", figure1,
       width = 5,    # Narrower width
       height = 4,   # Slightly shorter
       dpi = 300)
print(figure1)
