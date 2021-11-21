# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Summary       : Script to load thromboSeq processed HTSeq, perform data QCs, differential expression 
#                 of splice junction analysis, and PSO-enhanced development of thromboSeq classification
#                 algorithms
# Date          : 14th of November 2018
# Revision      : 19th of November 2018
# Revision v1.1 : 28th of November 2018
# Revision v1.2 : 29th of December 2018
# Revision v1.3 : 20th of March 2019
# Revision v1.4 : 29th of August 2019
# Revision v1.5 : 17th of April 2021
# req. packages : methods       [general functions]
#                 reshape       [visualize data in a heatmap]
#                 affy          [visualize data in a heatmap]
#                 ggplot2       [visualize data in a heatmap]
#                 dendroextras  [visualize data in a heatmap]
#                 e1071         [SVM-algorithm and classification]
#                 ROCR          [generate ROC curves and calculate AUC-values]
#                 caret         [calculate highly correlated transcripts during spliced RNA biomarker panel selection]
#                 pROC          [calculate ROC 95% confidence intervals]
#                 ppso          [run particle-swarm optimization algorithms, 
#                                included as ppsoThromboSeq.tar in the Supplemental Code (/bin/ppsoThromboSeq.tar); 
#                                install via: install.packages("ppsoThromboSeq.tar", type="source", repos=NULL)]
#                 edgeR         [store data in a data object called a DGEList, calculate normalization factor, 
#                                perform ANOVA Likelihood-ratio analysis]
#                 doSNOW/doMC   [parallel computing in order to enable efficient data analysis]
#                 foreach       [enable for repeated executions]
#                 RUVSeq        [perform RUVg data correction]
#
# note: If you do not have one of these packages installed, install by source("https://bioconductor.org/biocLite.R") and then biocLite("packageName")

# open R by running in terminal the following command: R

# continue in R
# load the required functions by running the following commands: 
source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO.R')

# Collect the intron-spanning read counts outputted by the mapping pipeline 
# by running the following command:
counts <- collect.read.counts(inputDir = "symlinksHTSEQ/")

# Prepare a DGE object for subsequent data analysis.
# For this, run the following command:
dge <- prepare.dge.object(sample.info = "sampleInfo.csv",
                          gene.info = "bin/dgeGenesEnsembl75.RData",
                          read.counts = counts)

# Filter the dataset for low abundant RNAs.
# For this, run the following command:
dge <- filter.for.platelet.transcriptome.TrainEval.group(dge)

# Perform thromboSeqQC pipeline.
# For this, run the following command:
predictedGroup = 'Malignant'
controlGroup = 'nonMalignant'
predictedGroupColor = '#a8ddb5'
controlGroupColor = '#43a2ca'
trainingEvaluationColor = '#95A5A6'
validationColor = '#212F3C'
analysis.settings <- data.frame(predictedGroup,controlGroup,predictedGroupColor,controlGroupColor,trainingEvaluationColor,validationColor,stringsAsFactors=FALSE)
training.evaluation.samples = NA ## here, provide the colnames of samples assigned to training and evaluation

dgeIncludedSamples <- thromboSeqQC.TrainEval(dge,
                                             training.series.only.samples = training.evaluation.samples,
                                             TrainEvalBased = TRUE)

# Perform thromboSeq ANOVA differential expression analysis of splice junctions
# For this, run the following command:
thromboSeq.anova <- thromboSeqANOVA(dgeIncludedSamples)

# Perform PSO-enhanced thromboSeq classifier development.
# For this, run the following command:
thromboPSO <- thromboSeqPSO(dge = dgeIncludedSamples, 
                            n.particles = 100,
                            n.iterations = 10,
                            number.cores = 8
                            )

# Validate the developed thromboSeq algorithm.
# For this, run the following command:
thromboPSOreadout <- thromboSeqPSO.readout(dge = dgeIncludedSamples,
                                           replace.counts.validation = 0, 
                                           filter.clinical.characteristics.validation = NA,
                                           filter.clinical.characteristics.group = NA,
                                           filter.clinical.characteristics.specified = NA,
                                           readout.training = T, 
                                           readout.evaluation = T, 
                                           readout.validation = T, 
                                           apply.rule.in.readout = F, 
                                           rule.in.setting = NA, 
                                           apply.rule.out.readout = F, 
                                           rule.out.setting = NA,
                                           additional.dge.to.predict = NA,
                                           number.cores = 8, 
                                           clinical.data.in.output = NA
                                           )

# Run LOOCV analysis on all samples
thromboSeqLOOCV <- thromboSeq.LOOCV(dge = dge)

# Perform control experiments for thromboSeq classifier development.
# For this, run the following command:
control <- thromboSeqPSO.controls(dge = dgeIncludedSamples,
                                  filter.clinical.characteristics.validation = NA,
                                  filter.clinical.characteristics.group = NA,
                                  filter.clinical.characteristics.specified = NA,
                                  thromboSeqPSO.shuffled = T,
                                  thromboSeqPSO.iterations = T,
                                  number.cores = 8
                                  )

# Perform digitalSWARM analysis.
# For this, run the following command:
digitalSWARM.output <- digitalSWARM( dge = dgeIncludedSamples,
                                     variable.to.assess = c("lib.size"),
                                     variable.to.assess.thresholds = c(0.0, 1.0),
                                     select.biomarker.FDR = FALSE,
                                     FDR.value.threshold = 0.1, 
                                     percentage.for.training = 30,
                                     percentage.for.evaluation = 30,
                                     percentage.for.verification = 20,
                                     percentage.for.validation = 20,
                                     training.samples.provided = NULL,
                                     evaluation.samples.provided = NULL,
                                     verification.samples.provided = NULL,
                                     validation.samples.provided = NULL,
                                     n.training.iterations = 200,
                                     minimum.n.transcripts.biomarkerpanel = 2,
                                     n.particles = 100,
                                     n.iterations = 50,
                                     number.cores = 8, 
                                     verbose = TRUE
)

# Perform control experiments for digitalSWARM analysis.
# For this, run the following command:
digitalSWARMshuffled <- digitalSWARM.shuffled(dge = dgeIncludedSamples,
                                              FDR.value.threshold = 0.1,
                                              variable.to.assess = c("lib.size"),
                                              variable.to.assess.thresholds = c(0.0, 1.0),
                                              select.biomarker.FDR = FALSE,
                                              percentage.for.training = 30,
                                              percentage.for.evaluation = 30,
                                              percentage.for.verification = 20,
                                              percentage.for.validation = 20,
                                              training.samples.provided = NULL,
                                              evaluation.samples.provided = NULL,
                                              verification.samples.provided = NULL,
                                              validation.samples.provided = NULL,
                                              n.training.iterations = 200,
                                              minimum.n.transcripts.biomarkerpanel = 2,
                                              n.particles = 100,
                                              n.iterations = 50,
                                              n.shuffled = 100,
                                              number.cores = 8, 
                                              verbose = TRUE,
                                              matrix = digitalSWARM.output)

# Optionally save and store the R-session
sessionInfo <- sessionInfo()
save.image(paste(timestamp(),".RData",sep=""))
# End