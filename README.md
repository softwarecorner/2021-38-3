# R packages for selecting important interactions via regularization

Ryan Peterson

Department of Biostatistics and Informatics, 
Colorado School of Public Health, 
University of Colorado Anschutz Medical Campus

_Editor note: This article was originally published in the [Biometric Bulletin (2021) Volume 38 Issue 3](https://www.biometricsociety.org/publications/biometric-bulletin). The example code is included as `examples.R` file._

Have you ever presented null results to disappointed researchers, and then been asked the question “but what about interactions; are any of those significant?” I have heard this question from clinicians and researchers from many fields of science. While usually asked in earnest, **this question is a dangerous one**; the sheer number of interactions can greatly inflate the number of false discoveries in the interactions, resulting in difficult-to-interpret models with many unnecessary interactions. Still, there are times when these expeditions are necessary and fruitful. Thankfully, useful tools are now available to help with the process. This article discusses two regularization-based approaches: Group-Lasso INTERaction-NET (glinternet) and the Sparsity-Ranked Lasso (SRL). The glinternet method implements a hierarchy-preserving selection and estimation procedure, while the SRL is a hierarchy-preferring regularization method which operates under ranked sparsity principles (in short, ranked sparsity methods ensure interactions are treated more skeptically than main effects *a priori*).

## Useful package #1: ranked sparsity methods via **sparseR**

While currently in a beta-phase, the **sparseR** package has been designed to make dealing with interactions and polynomials much more analyst-friendly. Building on the recipes package, **sparseR** has many built-in tools to facilitate the prepping of a model matrix with interactions and polynomials; these features are presented in the package website located at https://petersonr.github.io/sparseR/. The simplest way to implement the SRL in **sparseR** is via a single call to the `sparseR()` function, here demonstrated with Fisher’s `iris` data set: 

```r
(srl <- sparseR(Sepal.Width ~ ., data = iris, k = 1, seed = 1))
```

```
Model summary @ min CV:
-----------------------------------------------------
  lasso-penalized linear regression with n=150, p=18
  (At lambda=0.0015):
    Nonzero coefficients: 10
    Cross-validation error (deviance): 0.07
    R-squared: 0.62
    Signal-to-noise ratio: 1.64
    Scale estimate (sigma): 0.267

  SR information:
             Vartype Total Selected Saturation Penalty
         Main effect     6        4      0.667    2.45
 Order 1 interaction    12        6      0.500    3.46


Model summary @ CV1se:
-----------------------------------------------------
  lasso-penalized linear regression with n=150, p=18
  (At lambda=0.0070):
    Nonzero coefficients: 7
    Cross-validation error (deviance): 0.08
    R-squared: 0.57
    Signal-to-noise ratio: 1.33
    Scale estimate (sigma): 0.285

  SR information:
             Vartype Total Selected Saturation Penalty
         Main effect     6        3      0.500    2.45
 Order 1 interaction    12        4      0.333    3.46
```

```r
summary(srl, at = "cv1se")
```

```
lasso-penalized linear regression with n=150, p=18
At lambda=0.0070:
-------------------------------------------------
  Nonzero coefficients         :   7
  Expected nonzero coefficients:   1.38
  Average mfdr (7 features)    :   0.198

                                  Estimate       z     mfdr Selected
Species_setosa                    0.810513 17.9513  < 1e-04        *
Sepal.Length                      0.191210  9.3371  < 1e-04        *
`Petal.Length:Petal.Width`        0.119640  5.0379  < 1e-04        *
`Petal.Width:Species_versicolor`  0.275341  3.1640 0.055680        *
`Sepal.Length:Petal.Length`      -0.052711 -3.2466 0.078121        *
`Sepal.Length:Species_setosa`     0.062782  2.5978 0.251076        *
Species_versicolor               -0.001653 -0.8052 1.000000        *
```

We see (via print and summary functions) that two models are displayed by default corresponding to two “smart” choices for the penalization parameter λ. The first model printed refers to the model where λ is set to minimize the cross-validated error, while the second one refers to a model where λ is set to a value such that the model is as sparse as possible while still being within 1 SD of the minimum cross-validated error. Visualizations are also available via sparseR that can help visualize both the solution path and the resulting model (interactions can be very challenging to interpret without a good figure!) 

```r
plot(srl)
effect_plot(srl, "Petal.Width", by = "Species", at = "cvmin")
effect_plot(srl, "Petal.Width", by = "Species", at = "cv1se")
```

<img width="381" alt="Picture 1" src="https://user-images.githubusercontent.com/2189134/131218514-d7ce360a-2d7f-4b66-9461-5940486a9ccc.png">

Note that while ranked sparsity principles were motivated by the estimation of the lasso (in a paper currently under review), they can also be implemented with MCP, SCAD, or elastic net and for binary, normal, and survival data. Finally, sparseR includes some functionality to perform forward-stepwise selection using a sparsity-ranked modification of BIC, as well as post-selection inferential techniques using sample splitting and bootstrapping.

## Useful package #2: hierarchy-preserving regularization via **glinternet**

Some argue that when it comes to interactions, hierarchy is very important (i.e., an interaction shouldn’t be included in a model without its constituent main effects). While ranked sparsity methods do *prefer* hierarchical models, they can often still produce non-hierarchical ones. The **glinternet** package and the function of the same name uses regularization for model selection under hierarchy constraint, such that all candidate models are hierarchical. **Glinternet** can handle both continuous and categorical predictors, but requires pre-specification of a numeric model matrix. It can be performed as follows:  

```r
X <- iris %>% 
  select(-Sepal.Width) %>% 
  mutate(Species = as.numeric(Species) - 1)

set.seed(321)
cv_fit <- glinternet.cv(X, Y = iris$Sepal.Width, numLevels = c(1,1,1,3))
```

The `cv_fit` object contains necessary information from the cross-validation procedure and the fits themselves stored in a series of lists. A more in-depth tutorial to extract coefficients (and facilitate a model interpretation) using the **glinternet** package can be found at https://strakaps.github.io/post/glinternet/. Importantly, both the **glinternet** and **sparseR** methods have associated predict methods which can yield predictions on new (or the training) data, shown below. For comparison, we also fit a “main effects only” model with **sparseR** by setting `k = 0`. 

```r
me <- sparseR(Sepal.Width ~ ., data = iris, k = 0, seed = 333)
p_me <- predict(me)
p_srl <- predict(srl)
p_gln <- as.vector(predict(cv_fit, X))
```

With a little help from the **yardstick** package’s `metrics()` function, we can compare the accuracy of each model’s predictions using root-mean-squared error (RMSE), R-squared (RSQ), and mean absolute error (MAE); see table below. Evidently, **glinternet** and SRL are similar in terms of their predictive performance. However, both outperform the main effects model considerably, suggesting interactions among other variables do have signal worth capturing when predicting `Sepal.Width`. 

Metric | glinternet | SRL | Main effects only
------|-----------|---------|------------
RMSE  | 0.24 |	0.24 |	0.26
RSQ   | 0.69 |	0.70 |	0.63
MAE   | 0.19 |	0.18 |	0.20

## Other packages worth mentioning: ncvreg, hierNet, visreg, sjPlot

The SRL and other sparsity-ranked regularization methods implemented in **sparseR** would not be possible without the **ncvreg** package, which performs the heavy-lifting in terms of model fitting, optimization, and cross-validation. The **hierNet** package is another hierarchy-enforcing procedure that may yield better models than **glinternet**, however the latter is more computationally efficient especially for situations with a medium-to-large number of covariates. Finally, when interactions or polynomials are included in models, figures are truly worth a thousand words, and packages such as **visreg** and **sjPlot** have great functionality for plotting the effects of interactions. 

## References

- Bien J and Tibshirani R (2020). hierNet: A Lasso for Hierarchical Interactions. R package version 1.9. https://CRAN.R-project.org/package=hierNet 
- Breheny P and Burchett W (2017). Visualization of Regression Models Using visreg. The R Journal, 9: 56-71. 
- Breheny P and Huang J (2011). Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection. Ann. Appl. Statist., 5: 232-253.
- Kuhn M and Vaughan D (2021). yardstick: Tidy Characterizations of Model Performance. R package version 0.0.8. https://CRAN.R-project.org/package=yardstick
- Lim M and Hastie T (2020). glinternet: Learning Interactions via Hierarchical Group-Lasso Regularization. R package version 1.0.11. https://CRAN.R-project.org/package=glinternet 
- Lüdecke D (2021). sjPlot: Data Visualization for Statistics in Social Science. R package version 2.8.8. https://CRAN.R-project.org/package=sjPlot
- Peterson R (2021). sparseR: Variable selection under ranked sparsity principles for interactions and polynomials. https://github.com/petersonR/sparseR/. 
- Peterson R and Cavanaugh J (2021+). Ranked Sparsity: A Cogent Regularization Framework for Selecting and Estimating Feature Interactions and Polynomials. [arXiv:2107.07594](https://arxiv.org/abs/2107.07594) 

