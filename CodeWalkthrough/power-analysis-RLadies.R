### Power simulation script for the language learning experiment
### By Marianne de Heer Kloots
### with help from Laurel Brehm
### edited for the R-Ladies workshop 06-06-2019

library(readr)    # for loading data
library(dplyr)    # for data manipulation
library(lme4)     # for glmer computation
library(simr)     # for changing glmer effect sizes
library(ggplot2)  # for plots

sims <- 1    # number of simulations per setting
pp_per_lang <- seq(5, 15)    # numbers of participants to simulate
condition_effect_size <- c(0.1)    # group size effect sizes to simulate
structure_effect_scale <- c(0.15)    # factors to scale the effect of structure with

power_simulation <- function(number_of_subjects_per_language, group_size_effect, structure_effect_scale){
  #' Power simulation
  #' This function runs one simulation computing results for a given number of participants,
  #' condition (group size) effect and linguistic structure effect.
  
  # predict learnability (average accuracy) scores per condition for the given effect sizes, based on the pilot model
  predicted_accs <- predict_accuracies(pilot_test_acc_regression_linear, group_size_effect, structure_effect_scale)

  # create the given number of participants to simulate results for
  total_number_of_subjects <- number_of_subjects_per_language * number_of_languages
  
  simulated_participants <- data.frame(Language.ID = rep(experiment_language_info$Language.ID, each = number_of_subjects_per_language), 
                                       Participant.ID = 1:total_number_of_subjects)
  
  simulated_data <- experiment_languages %>%
    left_join(., simulated_participants, by = 'Language.ID') %>%
    arrange(Participant.ID)
  
  # simulate all accuracy scores per condition (= per language)
  simulated_data$Accuracy <- NA
  for (lang in unique(simulated_data$Language.ID)){
    n <- nrow(simulated_data[simulated_data$Language.ID == lang,])
    acc <- predicted_accs[predicted_accs$Language.ID == lang,]$Accuracy
    simulated_data[simulated_data$Language.ID == lang,]$Accuracy <- rbinom(n=n, size=1, prob=acc)
  }

  # glmer models on simulated results
  # comparing model1 & model2 will give us the effect of structure
  # comparing model1 & model3 will give us the effect of group size
  model1 <- glmer(Accuracy ~ Condition + LinguisticStructure + (1|Participant.ID) + (1|Item.ID), 
                  family="binomial", 
                  data=simulated_data, 
                  control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  model2 <- glmer(Accuracy ~ Condition + (1|Participant.ID) + (1|Item.ID), 
                  family="binomial", 
                  data=simulated_data,
                  control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  model3 <- glmer(Accuracy ~ LinguisticStructure + (1|Participant.ID) + (1|Item.ID), 
                  family="binomial", 
                  data=simulated_data,
                  control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  
  return(c(round(anova(model2, model1)[2,8], 5), round(anova(model3, model1)[2,8], 5)))
}

predict_accuracies <- function(model, group_size_effect, structure_effect_scale){
  #' Average accuracy score prediction
  #' This function computes the learnability (average accuracy) predictions per experimental condition (group size + structure score):
  #' the probability of success parameter that we use to sample from binomial distributions in the simulation.
  #' This function returns a dataframe of predictions per condition.
  
  # rescale the effect of LinguisticStructure (possible through simr)
  fixef(model)[2] <- structure_effect_scale*fixef(model)[2]
  
  group_size_scaling = 1 + group_size_effect

  # create dataframe to put the learnability predictions in
  pp_ids <- data.frame(Language.ID = factor(experiment_language_info$Language.ID), 
                       Participant.ID = 1:number_of_languages)
  acc_predictions <- experiment_languages %>%
    left_join(., pp_ids, by = 'Language.ID') %>%
    arrange(Participant.ID)
  
  # predict the learnability (average accuracy) scores per condition
  acc_predictions$Accuracy <- predict(model, newdata=acc_predictions, type="response", re.form=NA)
  
  # scale the learnability predictions for big groups by the given group effect size (maxed at .99 learnability)
  predicted_accs <- acc_predictions %>%
    mutate(Accuracy = ifelse(Condition == 'big', 
                             ifelse(group_size_scaling*Accuracy > .99, .99, 
                                    group_size_scaling*Accuracy), 
                             Accuracy)) %>%
    group_by(Language.ID, LinguisticStructure, Condition) %>%
    summarise(Accuracy = mean(Accuracy)) %>%
    arrange(LinguisticStructure)
  
  return(predicted_accs)
}

# load the languages that we will use in the experiment (2x5 languages; for 2 group sizes and 5 structure levels)
experiment_language_info <- read_delim('selected_languages.csv', delim = ';') %>%
  select('Language.ID', 'LinguisticStructure') %>%
  mutate(LinguisticStructure = as.numeric(scale(LinguisticStructure, center=TRUE, scale=FALSE)))
experiment_languages <- read_csv('all_languages.csv') %>%
  filter(Language.ID %in% experiment_language_info$Language.ID) %>%
  select('Language.ID', 'Condition', 'Item.ID') %>%
  left_join(experiment_language_info, ., by='Language.ID')
number_of_languages <- nrow(experiment_language_info)

# glmer for the pilot data (only small languages on 3 structure levels), that we base our simulations on
test_trials <- read_csv('pilot_data.csv') %>%
  filter(Task == 'MemorizationTest') %>%
  mutate(LinguisticStructure = as.numeric(scale(LinguisticStructure, center=TRUE, scale=FALSE))) %>%
  select('LinguisticStructure', 'Condition', 'Participant.ID', 'Item.ID', 'Accuracy') %>%
  arrange(LinguisticStructure)

pilot_test_acc_regression_linear <- glmer(Accuracy ~ LinguisticStructure + (1 | Item.ID) + (1 | Participant.ID), 
                                          family = "binomial", 
                                          data = test_trials)

StructureEffectSizes <- data.frame(ScaleFactor=structure_effect_scale,
                                   Linear=structure_effect_scale*fixef(pilot_test_acc_regression_linear)[2])

write.csv(StructureEffectSizes, 'structure_effect_sizes.csv', row.names=FALSE)

# create dataframe for the simulation results
nrows <- length(pp_per_lang)*length(condition_effect_size)*length(structure_effect_scale)*sims
powerdf <- data.frame(Participants=numeric(nrows),
                      LinguisticStructure=numeric(nrows), 
                      Condition=numeric(nrows),
                      ConditionEffectSize=numeric(nrows),
                      StructureEffectSize=numeric(nrows))

# loop over all settings to run the power simulation
row = 1
for (n in pp_per_lang){
  for (s in structure_effect_scale){
    for (c in condition_effect_size){
      
      # compute results for this setting (for the given number of simulations)
      cat(paste("\n\n", n, " ", s, " ", c, " \n"))
      pb <- txtProgressBar(max=sims)
      out <- replicate(sims, { setTxtProgressBar(pb, getTxtProgressBar(pb) + 1); power_simulation(n, c, s) })
      
      # put results in dataframe
      powerdf[row:(row+sims-1),]$Participants <- rep(n*number_of_languages, sims)
      powerdf[row:(row+sims-1),]$LinguisticStructure <- out[1,]
      powerdf[row:(row+sims-1),]$Condition <- out[2,]
      powerdf[row:(row+sims-1),]$ConditionEffectSize <- rep(c, sims)
      powerdf[row:(row+sims-1),]$StructureEffectScale <- rep(s, sims)
      
      # write results to csv
      write.csv(powerdf, 'power_analysis_results.csv', row.names=FALSE)
      
      # update row number
      row <- row + sims
    }
  }
}

## POWER PLOTS
# (from the demo simulation)
demo_structure_power <- powerdf %>%
  group_by(Participants) %>%
  summarize(DemoStructurePower=mean(LinguisticStructure < 0.05))

ggplot(demo_structure_power, aes(x=Participants, y=DemoStructurePower)) + 
  geom_line(colour="blue") +
  ggtitle("Power for effect of structure") +
  ylim(0, 1)

demo_condition_power <- powerdf %>%
  group_by(Participants) %>%
  summarize(DemoConditionPower=mean(Condition < 0.05))

ggplot(demo_condition_power, aes(x=Participants, y=DemoConditionPower)) + 
  geom_line(colour="red") +
  ggtitle("Power for effect of condition") +
  ylim(0, 1)

# (from my real power analysis)
real_power_results <- read_csv('real_power_analysis_results.csv')

structure_power <- real_power_results %>%
  group_by(Participants) %>%
  summarize(StructurePower=mean(LinguisticStructure < 0.05))

ggplot(structure_power, aes(x=Participants, y=StructurePower)) +
  geom_line(colour="blue") +
  ggtitle("Power for effect of structure") +
  ylim(0, 1)

condition_power <- real_power_results %>%
  group_by(Participants) %>%
  summarize(ConditionPower=mean(Condition < 0.05))

ggplot(condition_power, aes(x=Participants, y=ConditionPower)) +
  geom_line(colour="red") +
  ggtitle("Power for effect of condition") +
  ylim(0, 1)