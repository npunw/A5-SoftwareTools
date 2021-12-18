library(rentrez)
library(Biostrings)
library(tidyverse)
library(e1071)
library(caret)
library(Boruta)
library(pROC)
########DATA ACQUISITION###########

#Search for COI sequences for our two species using entrez
mellifera.COI.search <- entrez_search(db="nuccore", term = "Apis mellifera[ORGN] AND 10:1000[SLEN] AND COI[GENE]", retmax =1000, use_history = T)
cerana.COI.search <- entrez_search(db="nuccore", term = "Apis cerana[ORGN] AND 10:1000[SLEN] AND COI[GENE]", retmax =1000, use_history = T)


m_summ <- entrez_summary(db = "nuccore", id = mellifera.COI.search$ids)
c_summ <- entrez_summary(db = "nuccore", id = cerana.COI.search$ids)

m_fetch <- entrez_fetch(db = "nuccore", rettype = "fasta", web_history = mellifera.COI.search$web_history)
c_fetch <- entrez_fetch(db = "nuccore", rettype = "fasta", web_history = cerana.COI.search$web_history)

#We will write our sequences to a fasta file for each species
write(m_fetch, "mellifera.fasta", sep = "\n")
write(c_fetch, "cerana.fasta", sep = "\n")

#Read in our fasta files
m_seq <- readDNAStringSet("mellifera.fasta")
c_seq <- readDNAStringSet("cerana.fasta")

#Build a dataframe for each species' sequence set
df.mellifera <- cbind(as.data.frame(m_seq)$x, as.data.frame(rep("mellifera", length(m_seq))))

df.cerana <- cbind(as.data.frame(c_seq)$x, as.data.frame(rep("cerana", length(c_seq))))

colnames(df.mellifera) <- c("Sequence", "Species")
colnames(df.cerana) <- c("Sequence", "Species")

########DATA EXPLORATION###########

#lets make a function that will generate a histogram given an Apis species' dataset. This returns a statement of the sequence length median and plots our histogram in our "Plots" tab. You can specify colour of the bars as well.
generate.hist <- function (df.taxa, col ="grey"){
  taxa_str <- deparse(substitute(df.taxa))
  taxa_name <- tail(str_split(taxa_str, "\\.")[[1]], n=1)
  histo <- hist(nchar(df.taxa$Sequence),
                    xlab = "Sequence Length (bp)",
                    main = paste("COI Sequence Length Distribution \nfor Apis", taxa_name),
                    col = col)
  taxa_med <- median(nchar(df.taxa$Sequence))
  abline(v=taxa_med, col="red", lwd=3)
  print (paste("Median sequence length: ",taxa_med))
}

#Here are our histograms
generate.hist(df.mellifera, "yellow")
generate.hist(df.cerana, "black")

#Lets join our species data into one dataframe
df.all <- rbind(df.mellifera, df.cerana)

#Lets only accept sequences between 400 and 900bp
df.all <- subset(df.all, nchar(Sequence) >400 & nchar(Sequence) <900)

#If we want to take at look at overall sequence length distribution
hist(nchar(df.final$Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "COI Sequence Length Distribution Apis")

max(str_count(df.all$Sequence,"N" ))
# at most, our sequences have 4 "N"s, no need to remove any sequences 

########DATA ANALYSIS###########

#Before analysis, add some features for the SVM
df.all$label <- df.all$Species

df.all[df.all$label == "mellifera",]$label <- "-1"
df.all[df.all$label == "cerana",]$label <- "1"
df.all$label <- as.factor(df.all$label)

df.all <- cbind(df.all, as.data.frame(letterFrequency(DNAStringSet(df.all$Sequence), letters = c("A", "C","G", "T"))))


#Monomer frequencies
df.all$Aprop <- (df.all$A) / (df.all$A + df.all$T + df.all$C + df.all$G)
df.all$Tprop <- (df.all$T) / (df.all$A + df.all$T + df.all$C + df.all$G)
df.all$Gprop <- (df.all$G) / (df.all$A + df.all$T + df.all$C + df.all$G)
df.all$Cprop <- (df.all$C) / (df.all$A + df.all$T + df.all$C + df.all$G)


#Di-mer features added
df.all <- cbind(df.all, as.data.frame(dinucleotideFrequency(DNAStringSet(df.all$Sequence), as.prob = TRUE)))

#Tri-mer features added as well
df.all <- cbind(df.all, as.data.frame(trinucleotideFrequency(DNAStringSet(df.all$Sequence), as.prob = TRUE)))

#Tetra-mer features added
df.all <- cbind(df.all, as.data.frame(oligonucleotideFrequency(DNAStringSet(df.all$Sequence), width = 4, as.prob = TRUE)))



#Before anything else, lets split our dataset into train-test sets 75/25
set.seed(1234)
train_ind <- sample(seq_len(nrow(df.all)), size = round(0.75 * (nrow(df.all))))

df.Train <- df.all[train_ind, ]
df.Test <- df.all[-train_ind, ]

#We have 344 feature columns. This is quite a lot; we want to have a conservative amount of features. We should not choose features arbitrarily; we need a method of choosing features meaningfully.We will use the package "Boruta" for feature selection.


#Perhaps we should consider tetra-mers.
df.features <- df.all[, 3:347]

#Due to the many features considered (>300), this may take a while.
boruta_output <- Boruta(label ~., data=na.omit(df.features), doTrace=0)
#boruta_signif <- getSelectedAttributes(boruta_output, withTentative = F)
#print(boruta_signif)


imps <- attStats(boruta_output)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
head(imps2[order(-imps2$meanImp), ],5)  
# Pick the 5 features with the highest importance score; i.e., the highest ranked features according to boruta

#ACCT, AC, CCTG, CTAT, and ACTA, seem to be important tetra-mers.

#Now lets fit our SVM with our selected features
svmfit = svm(label ~ ACCT + AC + CCTG + CTAT + ACTA,
             data = df.Train,
             kernel = "linear",
             cost = 10,
             scale = FALSE)


#Lets look at how our SVM predicts against the test set
table(df.Test$label, predict(svmfit, newdata=df.Test))
#Not good. The -1,-1 and 1,1, cells in the table represent correct classifications. While all mellifera sequences were correctly classified, none of the cerana were correctly classified.

#perhaps lets try different kernels
svmfit = svm(label ~ ACCT + AC + CCTG + CTAT + ACTA,
             data = df.Train,
             kernel = "polynomial",
             cost = 10,
             scale = FALSE)
table(df.Test$label, predict(svmfit, newdata=df.Test))
#No change

svmfit = svm(label ~ ACCT + AC + CCTG + CTAT + ACTA,
             data = df.Train,
             kernel = "radial",
             cost = 10,
             scale = FALSE)
table(df.Test$label, predict(svmfit, newdata=df.Test))
#No change


#Ive noticed many tetra-mers were selected by feature selection. Lets look more closely at monomers, dimers, and trimers and exclude tetramers for features selection.
df.features <- df.all[, 3:91]


boruta_output <- Boruta(label ~., data=na.omit(df.features), doTrace=0)
imps <- attStats(boruta_output)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
head(imps2[order(-imps2$meanImp), ],5)

#looks like ACC, AC, CAC, TAC, and CTA are important features.
#Now lets fit our SVM with our selected features

svmfit = svm(label ~ ACC + AC + CAC + TAC + CTA,
             data = df.Train,
             kernel = "linear",
             cost = 10,
             scale = FALSE)

#Lets look at how our SVM predicts against the test set
table(df.Test$label, predict(svmfit, newdata=df.Test))
#Better, but still not good. Most of the cerana sequences are misclassified.


#There are no monomer proportions in our selected features. Perhaps we should consider adding those into the model.


svmfit = svm(label ~ ACC + AC + CAC + CTA + ACT +Aprop,
             data = df.Train,
             kernel = "linear",
             cost = 10,
             scale = FALSE)

#Lets look at how our SVM predicts against the test set
table(df.Test$label, predict(svmfit, newdata=df.Test))
#Better!

#Lets add another monomer
svmfit = svm(label ~ ACC + AC + CAC + CTA + ACT +Aprop +Gprop,
             data = df.Train,
             kernel = "linear",
             cost = 10,
             scale = FALSE)
table(df.Test$label, predict(svmfit, newdata=df.Test))
#No improvement from last, better only include one monomer. Lets try with just T

svmfit = svm(label ~ ACC + AC + CAC + TAC + CTA + Tprop,
             data = df.Train,
             kernel = "linear",
             cost = 10,
             scale = FALSE)
table(df.Test$label, predict(svmfit, newdata=df.Test))
#Not bad! Lets stop here.

#We can graph all the features and their importance from the Boruta output
df.features <- select(df.all, label, ACC, AC, CAC, TAC, CTA, Tprop)
boruta_output <- Boruta(label ~., data=na.omit(df.features), doTrace=0)
plot(boruta_output, cex.axis=.7, las=2, xlab="k-mer", main="Relative Variable Importance from Final Model")  


#Lets cross-validate our SVM

#Here we separate our dataset into 5 subsets or "folds"
#this variable will separate the data by index into 5 groups
folds = createFolds(df.all$label, k = 5)

#We will  now cross validate our SVM. 
cv = lapply(folds, function(x) {
#our training fold will be the four folds that is not our test fold
training_fold = df.all[-x, ] 
test_fold = df.all[x, ]
classifier = svmfit
y_pred = predict(classifier, newdata = test_fold[-3])

cm = table(test_fold[, 3], y_pred)
#This calculation for accuracy is (TP + TN) / (TP + TN + FP + FN)
accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
return(accuracy)
})

overall.accuracy = mean(as.numeric(cv))
overall.accuracy
#Nice! Overall accurary of 0.910

#Lets make a ROC for our SVM
ROC1 <- roc(df.Test$label, as.numeric(predict(svmfit, newdata=df.Test)), plot=T)

#USe ggplot for an aesthetic ROC
ggroc(ROC1, linetype=1, colour="blue", size=1, legacy.axes = T) +
  ggtitle ("ROC for SVM \nClassifer of Apis species") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="red", linetype="dashed") +
  xlab("False Positive Rate (1-Specificity)") +
  ylab("True Positive Rate (Sensitivity)") +
  coord_fixed(ratio=1)+
  annotate("text", x=0.75, y=0.25, label= paste("AUC: ", round(ROC1$auc, 3)))+



#Lets now compare with a logistic regression model. We will use the same features from the boruta output.
GeneMod1 <- glm(label ~ ACC + AC + CAC + TAC + CTA + Tprop, data = df.Train, family = binomial)


log.cv = lapply(folds, function(x) {
  #our training fold will be the foour folds that is not our test fold
  training_fold = df.all[-x, ] 
  test_fold = df.all[x, ]
  classifier = GeneMod1
  y_pred = predict(classifier, newdata = test_fold[-3])
  
  cm = table(test_fold[, 3], y_pred)
  #This calculation for accuracy is (TP + TN) / (TP + TN + FP + FN)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

log.overall.accuracy <- mean(as.numeric(log.cv))
log.overall.accuracy
#Not bad. The logistic regression model performs well according to the cross validation. We get an overall accuracy of 0.738, slightly lower than our SVMs 0.910



