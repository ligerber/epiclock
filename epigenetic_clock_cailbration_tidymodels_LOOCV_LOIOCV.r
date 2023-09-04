#load required packages

library(tidymodels)
library(GEOquery)
library(tidyverse)
library(skimr)

###
#obtain DNA methylation data

setwd("C:/Users/GER094/OneDrive - CSIRO/Legacy Publication")
#epi <- read.csv("./Data/normalized_betasN93_N139_365confidenceNoDups_AgeInc_MappedTAdu_BowtieQ10.csv", check.names = FALSE)
epi <- read.csv("./Data/normalized_betasN93_N139_greyExcluded_DupIncluded_AgeInc_IDs_mappedTAdu_BowtieQ10.csv", check.names = FALSE)
geo = setNames(data.frame(t(epi[,-1])), epi[,1])

###
#Data summary & preprocessing
###
geo$Age <- log(geo$Age +1)
Age <- geo$Age
#reduce dataset for code testing purposes
set.seed(123)
keep <- sample(1:ncol(geo), 900)
geo <- geo[,keep]

geo <-cbind(geo,Age)


###
#separate data into training and testing sets
###

methyl_dolphin_split <- initial_split(geo,prop = 0.8,strata=Age) #80% of data is used as training set, strata = age indicated that the distribution of 'age' is similar in the training and test data

methyl_dolphin_train <- training(methyl_dolphin_split) #extract the training set from the split dataset
methyl_dolphin_test <- testing(methyl_dolphin_split) # extract the test set 


###
#create preprocessing recipe
###

methyl_recipe <- #creates new recipe using the training data 
  recipe(methyl_dolphin_train) %>%
  update_role(everything()) %>% #all columns aredesignated as predictors
  update_role(Age,new_role = "outcome")  %>% #the column age is is the target variable to predict
  step_center(all_predictors()) %>% #centering of the predictor variable
  step_scale(all_predictors()) #scaling of predictor variable

#setp of a linear regression model with Lasso regularization
glmn_fit <- 
  linear_reg( mixture = 0.5, penalty = tune()) %>% #linear regression using Lasso, set to 0.5 to combne Lasso and Ridge regression. Tune indicates that the penalty term will be tuned during model selection
  set_engine("glmnet") 

#Cross-validatioon
folds <- vfold_cv(methyl_dolphin_train, v = 10, strata = Age, breaks= 2) #v indicates the number of folds per cross-validation, breaks specifies the number of times the cross-validation process will be repeated with different random splits

#Workflow setup
methyl_wf <- workflow() %>%
  add_model(glmn_fit) %>%
  add_recipe(methyl_recipe)

lasso_grid <- tibble(penalty = 10^seq(-3, 0, length.out = 50))

library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

#Hyperparamete tuning
methyl_res <- methyl_wf %>% 
  tune_grid(resamples = folds,
            grid = lasso_grid,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(rmse))


#Model evaluation
autoplot(methyl_res)

#Selection of best model and workflow finalisation
best_mod <- methyl_res %>% select_best("rmse")
best_mod

#Model fitting and prediction (fit the final model to the training data and 
#generate predictions for the testing data using augment())
final_fitted <- finalize_workflow(methyl_wf, best_mod) %>%
  fit(data = methyl_dolphin_train)

#Evaluation metrics RMSE and Pearson correlation coefficient
methyl_dolphin_aug <- augment(final_fitted, methyl_dolphin_test)
rmse(methyl_dolphin_aug,truth = Age, estimate = .pred)

methyl_dolphin_aug$tAge <- exp(methyl_dolphin_aug$Age)-1
methyl_dolphin_aug$tAgePred <- exp(methyl_dolphin_aug$.pred)-1

plot(methyl_dolphin_aug$tAgePred,methyl_dolphin_aug$tAge)

r <- cor(methyl_dolphin_aug$tAge, methyl_dolphin_aug$tAgePred, method = "pearson")

rmse(methyl_dolphin_aug,truth = tAge, estimate = tAgePred)

error <- methyl_dolphin_aug$tAgePred - methyl_dolphin_aug$tAge
MAE <- sqrt(median(error^2))
MAE

### Inference of importance of individual CpGs and their visualisation

library(vip)
library(cowplot)

#get the importance from glmnet using the select penalty
importance <- final_fitted %>%
  extract_fit_parsnip() %>%
  vi(lambda = best_mod$penalty) %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  )

#how many CpGs are retained
table(importance$Importance>0)

#plot the top 10 CpGs
importance %>% slice_max(Importance,n=10) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) + theme_cowplot()

#helper function to plot a CpG beta values against age
plotCpG <- function(cpg,dat){
  
  ggplot(dat,aes(x=!!sym(cpg),y=exp(Age)-1))+ #exp -1 to change from log age back to age in years
    geom_point() +
    theme_cowplot()
  
}

#plot the most important CpGs
importance %>% 
  slice_max(Importance,n=4) %>%
  pull(Variable) %>% 
  as.character() %>% 
  map(plotCpG,geo) %>%
  plot_grid(plotlist = .)

### LOOCV
geo$Individual_ID <- row.names(geo)
geo$Individual_ID <- as.character(geo$Individual_ID)

#predictedAgeID <- data.frame()

# Define the LOOCV setup
loocv_folds <- loo_cv(geo, strata = Age)

# Initialize a vector to store correlation coefficients
pred_values <- numeric(nrow(geo))
age_values <- numeric(nrow(geo))
IDs <- numeric(nrow(geo))

# Perform LOOCV loop
for (i in 1:nrow(loocv_folds)) {
  # Get the current fold
  current_fold <- loocv_folds$splits[[i]]
  
  # Train the model on all samples except the left-out one
  current_train <- analysis(current_fold)
  current_model <- final_fitted %>%
    fit(data = current_train)
  
  # Predict the left-out sample's age
  current_test <- assessment(current_fold)
  current_predictions <- predict(current_model, new_data = current_test) %>%
    bind_cols(current_test)
  
  # Calculate the correlation between predicted and actual ages
  pred_values[i] <- current_predictions$.pred
  age_values[i] <- current_predictions$Age
  IDs[i] <- current_predictions$Individual_ID
  print(i)
}


# Evaluation metrics RMSE and Pearson correlation coefficient
# Calculate the correlation between predicted and actual ages

age_values <- exp(age_values)-1
pred_values <- exp(pred_values)-1

cor_LOOCV <- cor(age_values, pred_values, method = "pearson")
cor_LOOCV


plot(age_values,pred_values)

error <- pred_values - age_values
MAE <- sqrt(median(error^2))
MAE

LOOCV_results <- data.frame(IDs,pred_values,age_values)

#write.csv(LOOCV_results, "FittedValues365_LOOCV_tidymodels_LOOCV_MappedTAdu_BowtieQ10.csv")


# Define the LOIOCV setup - not yet tested!

geo$Individual_ID <- row.names(geo)
geo$Individual_ID <- as.character(geo$Individual_ID)
geo$Individual_ID <- str_replace_all(geo$Individual_ID, '.1|.2|.3', '')

#predictedAgeID <- data.frame()

unique_individuals <- unique(geo$Individual_ID)

# Initialize vectors to store predicted and actual ages
all_pred_values <- numeric()
all_age_values <- numeric()
all_IDs <- numeric()

# Loop over unique individual IDs
for (i in 1:length(levels(as.factor(geo$Individual_ID)))){
  
  # Filter data to exclude the current individual
  current_data <-geo[geo$Individual_ID!=levels(as.factor(geo$Individual_ID))[i],]
  left_out_individual <- geo[geo$Individual_ID==levels(as.factor(geo$Individual_ID))[i],]
  
  # Fit model on filtered data
  current_fitted <- final_fitted %>% fit(data = current_data)
  
  # Predict the age for left-out individual
  current_prediction <- augment(current_fitted, new_data = left_out_individual)
  
  all_pred_values <- c(all_pred_values, current_prediction$.pred)
  all_age_values <- c(all_age_values, left_out_individual$Age)
  all_IDs <- c(all_IDs, left_out_individual$Individual_ID)
  
  print(paste("Iteration", i, "/", length(unique_individuals), "completed"))
}

# Evaluation metrics RMSE and Pearson correlation coefficient
# Calculate the correlation between predicted and actual ages

all_age_values <- exp(all_age_values)-1
all_pred_values <- exp(all_pred_values)-1

cor_LOIOCV <- cor(all_age_values, all_pred_values, method = "pearson")
print("Mean Correlation:")
cor_LOIOCV


plot(all_age_values,all_pred_values)

error <- all_pred_values - all_age_values
MAE <- sqrt(median(error^2))
MAE

LOIOCV_results <- data.frame(all_IDs,all_pred_values,all_age_values)

#write.csv(LOOCV_results, "FittedValues365_LOOCV_tidymodels_LOOCV_MappedTAdu_BowtieQ10.csv")
