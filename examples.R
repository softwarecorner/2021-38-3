################################################### -
## Title: example for article
## Author: Ryan Peterson
## Date Created: Tue Jun 29 11:48:04 2021
################################################### -

library(tidyverse)
library(sparseR)
library(glinternet)
library(kableExtra)

data("iris")

head(iris)

me <- sparseR(Sepal.Width ~ ., data = iris, k = 0, seed = 333)
me

srl <- sparseR(Sepal.Width ~ ., data = iris, k = 2, seed = 333)
srl

summary(srl, at = "cv1se")


par(mfcol = c(2,2), mar = c(4, 4, 3, 1))
plot(srl)
effect_plot(srl, "Petal.Width", by = "Species", at = "cvmin")
effect_plot(srl, "Petal.Width", by = "Species", at = "cv1se")

par(mfcol = c(2,2), mar = c(4, 4, 3, 1))
plot(srl)
effect_plot(srl, "Petal.Width", by = "Species", at = "cvmin", 
            plot.args = list(ylim = c(1.5, 5), legloc = "topright"), 
            main = "Min CV model")
effect_plot(srl, "Petal.Width", by = "Species", at = "cv1se", 
            plot.args = list(ylim = c(1.5, 5), legloc = "topright"), 
            main = "1SE CV model")


### glinternet example

set.seed(321)
X <- iris %>% 
  select(-Sepal.Width) %>% 
  mutate(Species = as.numeric(Species) - 1)

cv_fit <- glinternet.cv(X, Y = iris$Sepal.Width, numLevels = c(1,1,1,3))
plot(cv_fit)

p_me <- predict(me)
p_srl <- predict(srl)
p_gln <- as.vector(predict(cv_fit, X))

gln_res <- tibble(p_gln, y = iris$Sepal.Width) %>% 
  yardstick::metrics(y, p_gln) %>% 
  rename("glinternet"= .estimate) 
srl_res <- tibble(p_srl, y = iris$Sepal.Width) %>% 
  yardstick::metrics(y, p_srl) %>% 
  rename("SRL"= .estimate) 
me_res <- tibble(p_me, y = iris$Sepal.Width) %>% 
  yardstick::metrics(y, p_me) %>% 
  rename("Main effects only"= .estimate) 

gln_res %>% 
  bind_cols(srl_res[,3]) %>% 
  bind_cols(me_res[,3]) %>% 
  rename("Metric" = .metric) %>% 
  mutate(Metric = toupper(Metric)) %>% 
  select(-.estimator) %>% 
  kable(digits = 2) %>% 
  kable_styling(c("striped"), full_width = FALSE)
  
