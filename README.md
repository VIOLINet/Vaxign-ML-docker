# Vaxign-ML-docker

## Vaxign-ML

Vaxign-ML is a supervised machine learning classification to predict protective antigens. To identify the best machine learning method with optimized conditions, 5 machine learning algorithms (logistic regression, support vector machine, k-nearest neighbors, random forest, and extreme gradient boosting) were tested with biological and physiochemical features extracted from the Protegen database. Nested five-fold cross-validation and leave-one-pathogen-out validation were used to ensure unbiased performance assessment and the capability to predict vaccine candidates for a new emerging pathogen. The best performing model, Vaxign-ML (extreme gradient boosting trained on all Protegen data), was compared to three publicly available reverse vaccinology programs with a high-quality benchmark dataset, and showed superior performance in predicting protective antigens.

## Installation 
(Docker version >=1.13.1, API version >=1.26)

$ docker pull e4ong1031/vaxign-ml:latest

$ wget https://raw.githubusercontent.com/VIOLINet/Vaxign-ML-docker/master/VaxignML.sh

$ chmod a+x VaxignML.sh

$ ./VaxignML.sh [INPUT_FASTA] [OUTPUT_DIRECTORY] [ORGANISM_TYPE]

(You may need root privilege to run docker commands)

## Docker
https://hub.docker.com/r/e4ong1031/vaxign-ml
