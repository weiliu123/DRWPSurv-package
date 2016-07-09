# DRWPSurv-package
 This package implements the DRWPSurv method which predicts survival outcome using topologically inferred pathway activities.
# Details

Package: DRWPSurv

Type: Package

Title: Predicting survival outcome using pathway activities

Version: 1.0

Date: 2016-05-02

Author: Wei Liu

Maintainer: Wei Liu <30330590@qq.com>

Depends: igraph, Matrix, survival, glmnet

Description: This package implements the DRWPSurv method which predicts survival outcome using topologically infe

License: GPL(>=2)
# Index of help topics:

DRW:       Directed Random Walk

DRWPSurv-package:       Predicting survival outcome using pathway activities

dGMMirGraph:      The global pathway graph

fit.DRWPSurv:      Fit a Lasso-Cox model using DRWPSurv

getPathActivity:      Inferring pathway activity

getW:      Calculating the weights of genes

mRNA_matrix:      The expression data

pathSet:      Pathway set

predict.DRWPSurv:      Make predictions from a "DRWPSurv" object

survData:      Survival data
# Examples

data(dGMMirGraph)

data(pathSet)

data(mRNA_matrix)

data(survData)

trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[2], floor(4/5*dim(mRNA_matrix)[2]))

testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[2], trainSmpl.Idx)

trainSmpl <- mRNA_matrix[ ,trainSmpl.Idx]

testSmpl <- mRNA_matrix[ ,testSmpl.Idx]
fit <- fit.DRWPSurv(x.mRNA = t(trainSmpl), y = survData[trainSmpl.Idx,], DEBUG=TRUE,
standardize=TRUE, globalGraph = dGMMirGraph, pathSet = pathSet,
Gamma=0.7, alpha= 1, nfolds = 5)

predict.DRWPSurv(object = fit, newx.mRNA = t(testSmpl), type="link",s="lambda.min")
