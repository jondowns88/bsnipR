###############################################################################
# PROJECT SETUP
###############################################################################
#Load libraries
library(tidyverse)
library(MASS)
library(ggplot2)
library(ggthemes)
library(GGally)
library(sva)
library(caret)

#Set output folders
data_dir <- 'C:\\Users\\jondo\\Documents\\PyR\\Upwork\\BSNIP\\data\\'
output_dir<- 'C:\\Users\\jondo\\Documents\\PyR\\Upwork\\BSNIP\\outputs\\'

###############################################################################
# READ IN DATA AND PREPARE
###############################################################################
#Read in dataframe (also, first column was a weird character...let's fix it)
df <- read.csv("fa_49bundles.csv") %>%
  rename(subject = 1)
colnames(df) <- iconv(colnames(df), 'latin1', 'ASCII', sub = "")

#Inspect the data
head(df)
df %>% group_by(site) %>% summarize(n = n())

#Pull only numeric fields
fa <- df[,5:ncol(df)]

#Mean imputation of the FA df
fa_imputed <- fa
for(i in 1:ncol(fa)){
  fa_imputed[,i][is.na(fa_imputed[,i])] <- mean(fa_imputed[,i], na.rm = TRUE)
}

#Write imputed DF to csv file (column names are preserved)
write.csv(fa_imputed, paste0(data_dir, "df_NAN_R.csv"))

#Specify number of components for LDA and PCA, respectively
n_comp_LDA <- 1
n_comp_PCA <- 5

#Center and scale the entire dataframe (slight differences in rounding r/t Python)
#This seems to be the biggest source of difference relative to Python
scale_mod <- caret::preProcess(fa_imputed)
x <- predict(scale_mod, fa_imputed) %>%
  as.matrix()

#Pull the site vector
site_vector <- df$site

###############################################################################
# PERFORM LDA (pre)
###############################################################################
#Make chart labels
pre_plots_name <- "BSNIP data before ComBat - effect of site"
post_plots_name <- "BSNIP data after ComBat - effect of site"

#Fit the LDA model and make predictions
lda <- lda(site_vector ~ x) #Fit the model
linearDiscriminant <- lda %>% predict(data.frame(x)) #Make predictions

#Convert to df, add true values
linearDiscriminant <- data.frame(linearDiscriminant) %>%
  mutate(site = site_vector)

#Get predictive power (%)
lda_pred_power <- sum(linearDiscriminant$class == linearDiscriminant$site)/
  nrow(linearDiscriminant)

#Plot LDA results (with fill for the 'true' groups. 
pre_lda_plot_lab <- paste0("LDA - ", pre_plots_name, "\n", "Predict: ", round(lda_pred_power, 2))
pre_lda_plot <- ggplot(linearDiscriminant, aes(x = LD1, group = site, fill = site)) +
  geom_density(alpha = 0.7) +
  ggtitle(pre_lda_plot_lab) +
  theme_economist_white()

#Save the plot
ggsave(filename = paste0(output_dir, "LDA_", 
                         gsub("(\\s|-){1,3}", "_", pre_plots_name), ".png"),
       plot = pre_lda_plot)

###############################################################################
# PERFORM PCA (pre)
###############################################################################
#Run PCA, convert to DF, add site as a vector
pca <- princomp(x)$scores[,1:5] %>%
  data.frame() %>%
  mutate(site = site_vector)

#Generate the pairs plot
pre_pca_plot_lab <- paste0("PCA - ", pre_plots_name)
pre_pca_plot <- ggpairs(pca[,1:5],
        title = pre_pca_plot_lab,
        diag = list(continuous = 'barDiag'),
        upper = list(continuous = 'points'),
        mapping = ggplot2::aes(color = pca$site),
        legend = c(1,1))

#Save the plot
ggsave(filename = paste0(output_dir, "PCA_", 
                         gsub("(\\s|-){1,3}", "_", pre_plots_name), ".png"),
       plot = pre_pca_plot)

###############################################################################
# PERFORM ComBat
###############################################################################
#Run the combat model
new <- t(ComBat(dat = t(fa_imputed),
              batch = site_vector)) %>%
  as.matrix()

#We want to scale the data with the same model we used originally.
x_combat <- predict(scale_mod, new)

###############################################################################
# PERFORM LDA (post)
###############################################################################
#Fit the LDA model and make predictions
lda_combat <- lda(site_vector ~ x_combat) #Fit the model
linearDiscriminant_combat <- lda_combat %>% predict(data.frame(x_combat)) #Make predictions

#Convert to df, add true values
linearDiscriminant_combat <- data.frame(linearDiscriminant_combat) %>%
  mutate(site = site_vector)

#Get predictive power (%)
lda_pred_power_combat <- sum(linearDiscriminant_combat$class == linearDiscriminant_combat$site)/
  nrow(linearDiscriminant_combat)

#Plot LDA results (with fill for the 'true' groups. 
post_lda_plot_lab <- paste0("LDA - ", post_plots_name, "\n", "Predict: ", round(lda_pred_power_combat, 2))
post_lda_plot <- ggplot(linearDiscriminant_combat, aes(x = LD1, group = site, fill = site)) +
  geom_density(alpha = 0.7) +
  ggtitle(post_lda_plot_lab) +
  theme_economist_white()

#Save the plot
ggsave(filename = paste0(output_dir, "LDA_", 
                         gsub("(\\s|-){1,3}", "_", post_plots_name), ".png"),
       plot = post_lda_plot)

###############################################################################
# PERFORM PCA (post)
###############################################################################
#Run PCA, convert to DF, add site as a vector
pca_combat <- princomp(x_combat)$scores[,1:5] %>%
  data.frame() %>%
  mutate(site = site_vector)

#Generate the pairs plot
post_pca_plot_lab <- paste0("PCA - ", post_plots_name)
post_pca_plot <- ggpairs(pca_combat[,1:5],
                        title = post_pca_plot_lab,
                        diag = list(continuous = 'barDiag'),
                        upper = list(continuous = 'points'),
                        mapping = ggplot2::aes(color = pca_combat$site),
                        legend = c(1,1))

#Save the plot
ggsave(filename = paste0(output_dir, "PCA_", 
                         gsub("(\\s|-){1,3}", "_", post_plots_name), ".png"),
       plot = post_pca_plot)
