#load required packages

library(tidymodels)
library(GEOquery)
library(tidyverse)
library(skimr)
set.seed(123)
###
#obtain DNA methylation data

#epi <- read.csv("./Data/normalized_betasN93_N139_365confidenceNoDups_AgeInc_MappedTAdu_BowtieQ10.csv", check.names = FALSE)
epi <- read.csv("./normalized_betasN93_N139_greyExcluded_DupIncluded_AgeInc_IDs_mappedTAdu_BowtieQ10_allOver35y.csv", check.names = FALSE)
geo = setNames(data.frame(t(epi[,-1])), epi[,1])

geo <- geo[!(row.names(geo) %in% c("PUC.2")),]

###
#Data summary & preprocessing
###
geo$Age <- log(geo$Age +1)
Age <- geo$Age

#reduce dataset for code testing purposes - run these three lines

#keep <- sample(1:ncol(geo), 900)
#geo <- geo[,keep]
#geo <-cbind(geo,Age)


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
            metrics = metric_set(rmse),
            )


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
r

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
  
  ggplot(dat,aes(x=!!sym(cpg),y=exp(Age)-1)) +
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


## apply the model above to all samples

# read in file containing all epigenetic data (of known and unkown individuals)
epiAll <- read.csv("./normalized_betasN93_N139_allSamples_MappedTAdu_BowtieQ10.csv", check.names = FALSE)
methyl_dolphin_all = setNames(data.frame(t(epiAll[,-1])), epiAll[,1])

## no individual should be excluded, but in case there is a good reason, modify the line below
#geo <- geo[!(row.names(geo) %in% c("PUC.2")),]

###
#Data summary & preprocessing
###
#generate predictions for the testing data using augment())

#Evaluation metrics RMSE and Pearson correlation coefficient
methyl_dolphin_aug_all <- augment(final_fitted, methyl_dolphin_all)

methyl_dolphin_aug_all$tAgePred <- exp(methyl_dolphin_aug_all$.pred)-1

hist(methyl_dolphin_aug_all$tAgePred)

PredAge <- cbind(row.names(methyl_dolphin_all),methyl_dolphin_aug_all$tAgePred)
  
write.csv(PredAge, 'AgePredictedTidyModels_MappedTAdu_BowtieQ10.csv')
#If needed, LOOCV of the best fitted model can be gained from LOOCV script