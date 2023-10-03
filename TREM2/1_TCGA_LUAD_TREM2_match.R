### 1. PSM
rm(list=ls()); library(MatchIt); outdir="./Results/"

### 2. Call dataset
df <- read.csv('./Results/Merged_TCGA_LUAD_TREM2_L67_H33.csv', header = TRUE, row.names = 'case_submitter_id')
df <- df[,-c(1)]
str(df)

### 3. Before Match confirm the variables as numeric
sum(df$age); range(df$age)
sum(df$gender); range(df$gender)
sum(df$Stage); range(df$Stage)
sum(df$pack_years_smoked); range(df$pack_years_smoked)
dim(df)

### 4. Matching

## 4.1 Nearest neighbor matching - default method - 1:1 matching
m.1 <- matchit(Group ~ age + gender + Stage + pack_years_smoked, data=df, method="nearest", ratio=1)
summary(m.1, standardize=TRUE)
head(m.1$match.matrix)

## 4.2 Visualize
# Graph 1; Jitter plot
png(filename = "./Results/m1_jitter.png", width = 400, height = 300, units = "px", pointsize = 12)
plot(m.1, type="jitter")
dev.off()

# Graph 2; Histogram
plot(m.1, type="hist")

# Graph 3; eQQ
png(filename = "./Results/m1_eQQ.png", width = 300, height = 600, units = "px", pointsize = 12)
plot(m.1)
dev.off()

# Graph 4; Summary
plot(summary(m.1, standardize=TRUE))

## 4.3 Identify pairs
matched_dat <- match.data(m.1)
table(matched_dat$Group)
head(matched_dat)
attach(matched_dat)

### 5. Survival analysis
library(survival)

# 5.1 Log-rank test
survdiff(Surv(days_to_death, Status==1) ~ Group, data=matched_dat)

# 5.2 KM Plot 1
mfit = survfit(Surv(days_to_death, Status==1) ~ Group, data=matched_dat)
plot(mfit, xlab = 'Time(Days)', ylab = "Cumulative Survival Rate", lty = c(1, 2), conf.int = F, mark = 3)

# 5.3 KM Plot 2
library(packHV)
png(filename = "./Results/KM_TCGA_LUAD_L67_H33_TREM2.png", width = 500, height = 400, units = "px", pointsize = 12)
plot_km(Surv(days_to_death, Status) ~ Group, data=matched_dat, 
        conf.int = FALSE, xlab = "Time (Days)", ylab = "Cumulative Survival Rate", lwd = 1,
        lty = c(1, 2), mark = 3, main='Kaplan-Meier Curves')
dev.off()

# 5.4 Save matched data
write.csv(matched_dat,"./Results/Matched_TCGA_LUAD_TREM2_L67_H33.csv", row.names = TRUE)
