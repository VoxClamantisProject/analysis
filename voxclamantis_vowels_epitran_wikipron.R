# Analyze VoxClamantis v1.0 Vowel Formants
# Created by Eleanor Chodroff
# 1 June 2020

# This script takes as input the  Epitran and Wikipron midpoint_f1_f2 formant text files, the per_utt_mcd text files (scores), and reading_info.csv

# The script returns the numbers reported in the paper (# of languages, # of language families, counts of vowels, correlation tables for f1 and f2) and the figures

# If SAVE_SURPRISAL_INPUT is set to TRUE, it also saves the a text file for each reading that serves as input to the Dispersion analysis in the wilderness-analysis-master Python scripts

require(tidyverse)
require(reshape2)
require(xtable)
require(data.table)
require(jsonlite)
require(ggthemes)
require(ggExtra)
require(ggforce)

#################
### CHANGE ME ###
#################

INPUT_SUFFIX <- "_formants_mid.csv"
DATASET <- "epiwiki"
SAVE_SURPRISAL_INPUT <- TRUE
epiDir <- "./epitran_formants/"
wikiDir <- "./wikipron_formants/"
scoresDir <- "./per_utt_mcd/"
surprisalDir <- "./epiwiki_processed/"
outputDir <- "~/Desktop/"
inventory <- read_csv("./reading_info.csv")

#################
### FUNCTIONS ###
#################

# convert hertz to ERB
hz2erb <- function(x) {
    output <- (.0043 * x) + 1
    output <- 21.4 * log10(output)
    return(output)
}

# read in file
readFiles <- function(fileDir, file_list) {
	myFile <- read_csv(paste0(fileDir, file_list[i]),col_names=T, col_types=cols(
  		file = col_character(),
  		vowel = col_character(),
  		prec = col_character(),
  		foll = col_character(),
  		dur = col_double(),
  		f1_mid = col_double(),
  		f2_mid = col_double()))
  	return(myFile)
}

# get number of vowels by row count
getNumVowels <- function(x, newcolname) {
	count <- sapply(x, nrow)
	count <- enframe(count)
	colnames(count) <- c("lang", newcolname)
	return(count)
}

# exclude utterances with MCD score below some value (6 used below)
removeUttMcd <- function(x, cutoff) {
	x <- filter(x, mcd < cutoff)
}

# get rid of non-vocalic labels: only necessary for unitran
removeNonVowels <- function(x) {
  	x <- na.omit(x)
  	x <- subset(x, !vowel %in% c("QUESTIONMARK", "ONE", "ASTERISK", "QM", "AVAGRAHA"))
  	x <- subset(x, !grepl("p", vowel))
}

# remove outliers
removeOutliers <- function(x) {
 	x <- filter(x, count > 50)
  	x <- filter(x, erbF1 < meanf1plusSD)
  	x <- filter(x, erbF1 > meanf1minSD)
  	x <- filter(x, erbF2 < meanf2plusSD)
  	x <- filter(x, erbF2 > meanf2minSD)
  	x <- filter(x, dur < 0.3)
}

# populate dataset with mean values for each vowel (col) and reading (row)
vowelValues <- function(dataset, measure) {
	# reshape the data
	x0 <- dataset %>% select(c(lang, vowel, {{measure}})) %>% spread(vowel, {{measure}})
}

# get all correlations of means between vowel pairs
getCors <- function(vowel_values, count_cutoff) {
	x <- vowel_values
	x <- cor(x[,-1], use="pairwise.complete.obs")
	x <- reshape2::melt(x)
	names(x) <- c("vowel.x","vowel.y","value")

	# create dataset with vowel pair cors and counts
	x <- merge(x, vowelCooccur, by=c('vowel.x','vowel.y'))
	
	# get rid of NAs (where a vowel pair occurs only once)
	x <- na.omit(x)
	
	# keep only pairs with some minimum number of readings
	x <- subset(x, count > count_cutoff & vowel.x!=vowel.y)
	
	# rearrange by correlation magnitude and get rid of duplicates
	x <- x[order(x$value, decreasing=T),]
	x$keep <- rep(c("TRUE","FALSE"), nrow(x)/2)
	x <- subset(x, keep=="TRUE")
	x$keep <- NULL
	
	# rearrange by correlation count
	x <- x[order(x$count, decreasing=T),]
	return(x)
}

# get p-values for each correlation
getPvals <- function(corTable, vowel_values) {
	x0 <- vowel_values
	for (i in 1:nrow(corTable)) {
    	vowel1 <- corTable$vowel.x[i]
    	vowel2 <- corTable$vowel.y[i]
    	input1 <- as.numeric(unlist(x0[,which(colnames(x0) == vowel1)]))
    	input2 <- as.numeric(unlist(x0[,which(colnames(x0) == vowel2)]))
    	mycor <- cor.test(~input1 + input2)
    	pval <- mycor$p.value
    	corTable$pval[i] <- pval
    }
    return(corTable)
}

# make correlation table nice for publication / latex
latexCors <- function(x, cutoff) {
	x <- subset(x, vowel.x!="4")
	x$vowel.x <- paste0("phon{", x$vowel.x, "}")
	x$vowel.y <- paste0("phon{", x$vowel.y, "}")
	x$value <- round(x$value, 2)
	x$pval <- round(x$pval, 4)
	x <- subset(x, count > cutoff)
	x <- x %>% arrange(desc(value)) %>% select(c(vowel.x, vowel.y, count, value, pval))
	colnames(x) <- c("V1", "V2", "# Readings", "$r$", "$p$")
	return(x)
}

# prep the dataset for the correlation scatterplot
prepFigure <- function(dataset_with_means, vowel_values, sd_measure, vowel1, vowel2) {
	x0 <- vowel_values
	x0.sd <- dataset_with_means %>% select(c(lang, vowel, {{sd_measure}})) %>% spread(vowel, {{sd_measure}})
	x <-cbind(x0[,colnames(x0)%in%c('lang', vowel1,vowel2)], x0.sd[,colnames(x0.sd) %in% c('lang',vowel1, vowel2)])
	names(x) <- c('lang',vowel1,vowel2, 'lang2', paste0("sd", "_", vowel1,), paste0("sd", "_", vowel2))
	x$lang2 <- NULL
	x <- data.frame(x)
	x <- merge(x, inventory, by ="lang")
	x$MCD <- x$MCD1
	return(x)
}

# make correlation scatterplots
makeFigure <- function(prep_figure_data, label1, label2, corInfo, uniformityFigure) {
	x <- prep_figure_data[,2]
	y <- prep_figure_data[,3]
	sd_x <- prep_figure_data[,4]
	sd_y <- prep_figure_data[,5]
	minX <- round(min(x, na.rm=T))
	maxX <- ceiling(max(x, na.rm=T))
	minY <- round(min(y, na.rm=T))
	maxY <- ceiling(max(y, na.rm=T))

	if (uniformityFigure==TRUE) {
		minX <- minY <- min(minX, minY)
		maxX <- maxY <- max(maxX, maxY)
	}
	
	diffX <- maxX-minX
	diffY <- maxY-minY
	
	p <- ggplot(prep_figure_data, aes(x, y)) + geom_point(color = "white", fill = "white", size = 0) + 
		geom_ellipse(aes(x0 = x, y0 = y, a = sd_x*0.1, b = sd_y*0.1, angle = 0), 
		color = "slategray", fill = "slategray", alpha = 0.04) + 
		geom_smooth(method = lm, size = 0.75, se = TRUE, color = "black") + 
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) + 
		theme_few(20) + xlab("ERB") + ylab("ERB") + 
		scale_x_continuous(limit = c(minX,maxX)) + scale_y_continuous(limit = c(minY,maxY)) + 
		annotate("text", family="serif", x = minX + 0.8*diffX, y = minY, size = 7, label = corInfo, fontface = "italic") + 
		annotate("text", family="serif", x = minX + (diffX/2), y = maxY, size = 9, label = label1) + 
		annotate("text", family="serif", x = maxX, y = minY + (diffY/2), size = 9, label = label2, angle=-90) + 
		theme(text=element_text(family="serif", size=22))

	ggMarginal(p, type="histogram", xparams = list(colour="slategray", fill="ghostwhite"), 
	yparams = list(colour="ghostwhite", fill="slategray"))
}


##############################
######## READ IN DATA ########
##############################

epi_list <- list.files(path=epiDir, pattern="*_mid.csv")
wiki_list <- list.files(path=wikiDir, pattern="*_mid.csv")
scores_list <- list.files(path=scoresDir, pattern="*.scores")

# READ IN PER-UTT-MCD SCORES

scores <- list()
for (i in 1:length(scores_list)) {
	scoresi <- read_table2(paste0(scoresDir, scores_list[i]), col_names = FALSE)
	colnames(scoresi) <- c("file", "mcd")
	rec_lang <- substring(scores_list[i], 1, 6)
	scores[[rec_lang]] <- scoresi
}
	
scores[["FINV38"]] <- scores[["FIN38V"]]
scores[["FINV38"]]$file <- gsub("FIN38VN", "FINV38N", scores[["FINV38"]]$file)


### READ IN EPITRAN AND WIKIPRON DATA

langDFs <- list()

# READ IN EPITRAN FILES
for (i in 1:length(epi_list)) {
	langID <- substr(epi_list[i], 1, 6)
  	langDFs[[langID]] <- readFiles(epiDir, epi_list)
 }

# READ IN WIKIPRON FILES
for (i in 1:length(wiki_list)) {
	langID <- substr(wiki_list[i], 1, 6)
  	langDFs[[langID]] <- readFiles(wikiDir, wiki_list)
}

# PROCESS EPITRAN AND WIKIPRON FILES

for (i in 1:length(langDFs)) {
  	# convert hz to erb
  	langID <- names(langDFs[i])
  	langDFs[[langID]]$lang <- langID
  	langDFs[[langID]]$erbF1 <- hz2erb(langDFs[[langID]]$f1_mid)
  	langDFs[[langID]]$erbF2 <- hz2erb(langDFs[[langID]]$f2_mid)

  	# merge formants with utterance-level scores
  	langDFs[[langID]] <- left_join(langDFs[[langID]], scores[[langID]], by=c("file"="file"))
  	
  	# merge formants with wilderness inventory data (reading_info.csv)
    langDFs[[langID]] <- left_join(langDFs[[langID]], inventory, by="lang")
}

# get rid of vowel length which is non-contrastive in URDWTC labeling
langDFs[["URDWTC"]]$vowel <- gsub(":", "", langDFs[["URDWTC"]]$vowel)

# copy dataset
langDFs_orig <- langDFs 

###################################
######## OUTLIER EXCLUSION ########	     
###################################

# count number of vowels
numVowels <- getNumVowels(langDFs, "orig")

# prep readings with MCD > 8 for removal
removeme <- inventory$lang[which(inventory$MCD1 > 8)]
removeme <- removeme[removeme %in% names(langDFs)]

# remove files above file-mcd cutoff
for (i in 1:length(removeme)) {
	langDFs[[removeme[i]]] <- NULL
}

# count remaining vowels
count <- getNumVowels(langDFs, "file_mcd")
numVowels <- left_join(numVowels, count, by="lang") 

# remove utterances below mcd cutoff
langDFs <- lapply(langDFs, removeUttMcd, cutoff=6)

# count remaining vowels
count <- getNumVowels(langDFs, "utt_mcd")
numVowels <- left_join(numVowels, count, by="lang") 

# get rid of non vowels (Unitran only)
#langDFs <- lapply(langDFs, removeNonVowels)

# get mean, sd, count, and mean+/-2sd per vowel
vowelInfo <- list()

for (i in 1:length(langDFs)) {
  	vowelInfo[[i]] <- langDFs[[i]] %>% group_by(vowel) %>% summarise(
  	meanf1 = mean(erbF1, na.rm=T), sdf1=sd(erbF1, na.rm=T), 
  	meanf2 = mean(erbF2, na.rm=T), sdf2=sd(erbF2, na.rm=T), 
  	count=length(vowel))
  	
  	vowelInfo[[i]]$meanf1plusSD <- vowelInfo[[i]]$meanf1+2*vowelInfo[[i]]$sdf1
  	vowelInfo[[i]]$meanf1minSD <- vowelInfo[[i]]$meanf1-2*vowelInfo[[i]]$sdf1
  	vowelInfo[[i]]$meanf2plusSD <- vowelInfo[[i]]$meanf2+2*vowelInfo[[i]]$sdf2
  	vowelInfo[[i]]$meanf2minSD <- vowelInfo[[i]]$meanf2-2*vowelInfo[[i]]$sdf2
}

# merge these numbers with the main list
for (i in 1:length(vowelInfo)) {
 	langDFs[[i]] <- left_join(langDFs[[i]], vowelInfo[[i]], by="vowel")
}

# remove outliers	
langDFs <- lapply(langDFs, removeOutliers)

### REMOVE READINGS WITH NO DATA
nodata <- c()
for (i in 1:length(langDFs)) {
	if (nrow(langDFs[[i]]) < 50) {
		nodata <- append(nodata, names(langDFs[i]))
	}
}
	
for (i in 1:length(nodata)) {
	langDFs[[nodata[i]]] <- NULL
}

# count remaining vowels
count <- getNumVowels(langDFs, "outlier_after")
numVowels <- left_join(numVowels, count, by="lang")

numVowels$retained_mcd <- numVowels$utt_mcd/numVowels$orig
numVowels$retained_total <- numVowels$outlier_after/numVowels$orig
numVowels$retained_mcd <- ifelse(numVowels$retained_mcd < 0.01, "NA", numVowels$retained_mcd)
numVowels$retained_mcd <- as.numeric(as.character(numVowels$retained_mcd))

numVowelsTable <- summary(numVowels)
summary(numVowels)

############################
### SAVE RETENTION STATS ###
############################

print(xtable(numVowelsTable), file=paste0(outputDir, DATASET, "-retention-vowels.tex"))
write_csv(numVowels, paste0(outputDir, DATASET, "-retention-vowels-all.csv"))

##########################################
######## GET VOWEL FREQUENCY DATA ########
##########################################

# update number of tokens per vowel / remove old vowels
vowelFreqs <- list()

for (i in 1:length(langDFs)) {
  	langID <- names(langDFs[i])
  	vowelFreqs[[langID]] <- langDFs[[langID]] %>% group_by(vowel) %>% summarise(count=length(vowel))
  	vowelFreqs[[langID]]$total <- sum(vowelFreqs[[langID]]$count)
  	vowelFreqs[[langID]]$prop <- vowelFreqs[[langID]]$count/vowelFreqs[[langID]]$total
  	vowelFreqs[[langID]]$logprop <- log(vowelFreqs[[langID]]$prop)
  	vowelFreqs[[langID]]$lang <- langID
}

df.vowelFreqs <- do.call(rbind.data.frame, vowelFreqs)
    
nVowels <- df.vowelFreqs %>% group_by(lang) %>% summarise(count=length(vowel))
	
for (i in 1:length(langDFs)) {
  	removeCols <- which(names(langDFs[[i]]) %in% c('count', 'meanf1plusSD', 'meanf1minSD', 'meanf2plusSD','meanf2minSD'))
  	langDFs[[i]] <- langDFs[[i]][,-removeCols]
  	langDFs[[i]] <- left_join(langDFs[[i]], vowelFreqs[[i]], by="vowel")
}

###########################################################
### SAVE PROCESSED FILES FOR PYTHON DISPERSION ANALYSIS ###
###########################################################

if (SAVE_SURPRISAL_INPUT == TRUE) {
    for (i in 1:length(langDFs)) {
        lang <- langDFs[[i]]
        langid <- names(langDFs[i])
        langtmp <- lang[,-which(colnames(lang)%in%c("start", "end", "language", "code", "MCD0", "MCD1", "genus", "family", "macroarea"))]
        write.csv(langtmp, paste0(surprisalDir, langid, "_formants_mid_fin.csv"), quote=F, row.names=F)
    	}
}

#############################
### GET DESCRIPTIVE STATS ###
#############################

meansPerVowelPerLang <- list()

for (i in 1:length(langDFs)) {
  langID <- names(langDFs[i])  
  meansPerVowelPerLang[[langID]] <- langDFs[[langID]] %>% group_by(vowel) %>% summarise(
  mean_f1 = mean(erbF1), mean_f2 = mean(erbF2), mean_dur = mean(dur), 
  sd_f1 = sd(erbF1), sd_f2 = sd(erbF2), sd_dur = sd(dur), lang=langID, family=family[1])
}

myMeans <- do.call(rbind.data.frame, meansPerVowelPerLang)
myMeans <- left_join(myMeans, df.vowelFreqs, by=c("lang", "vowel"))

##############################
### SAVE IMPORTANT NUMBERS ###
##############################

formantNumbers <- matrix(ncol=1, nrow=11)
formantNumbers[1,1] <- length(langDFs) # number readings
formantNumbers[2,1] <- sum(vowelFreqs$count) # total number of vowels

# of vowels per language + average
vowelsPerLang <- vowelFreqs %>% group_by(lang) %>% summarise(count=sum(count))
summary(vowelsPerLang)

formantNumbers[3,1] <- mean(vowelsPerLang$count)
formantNumbers[4,1] <- median(vowelsPerLang$count)
formantNumbers[5,1] <- min(vowelsPerLang$count)
formantNumbers[6,1] <- max(vowelsPerLang$count)
formantNumbers[7,1] <- quantile(vowelsPerLang$count, 0.25)
formantNumbers[8,1] <- quantile(vowelsPerLang$count, 0.75)

# range of vowel inventory sizes per language
summary(nVowels)
formantNumbers[9,1] <- median(nVowels$count)
formantNumbers[10,1] <- min(nVowels$count)
formantNumbers[11,1] <- max(nVowels$count)

formantNumbers <- data.frame(formantNumbers)
formantNumbers$description <- c("# langs", "total # vowels", "mean # per lang", "median # per lang", "min # per lang", "max # per lang", "1q # per lang", "3q # per lang", "median # categories per lang", "min # categories per lang", "max # categories per lang")
formantNumbers <- formantNumbers %>% select(description, formantNumbers)

# counts of languages per vowel category
counts <- myMeans %>% group_by(vowel) %>% summarise(count=length(lang)) %>% arrange(desc(count))
head(counts)

# language family count
langFamilies <- inventory %>% filter(lang %in% myMeans$lang) %>% group_by(family) %>% count(family) %>% arrange(desc(n))

# language count
langs <- names(langDFs)
inv <- filter(inventory, lang %in% langs)

write.csv(nVowels, paste0(outDir, "nVowelCategories_", DATASET, ".csv"), quote=F, row.names=F)
write.table(formantNumbers, paste0(outDir, "formantNumbers_", DATASET, ".txt"), quote=F, sep="\t")
write.csv(vowelsPerLang, paste0(outDir, "nVowelTokens_", DATASET, ".csv"), quote=F, row.names=F)
write.csv(counts, paste0(outDir, "vowelFrequencies_", DATASET, ".csv"), quote=F, row.names=F)
write.csv(langFamilies, paste0(outDir, "langFamilies_", DATASET, ".csv"), quote=F, row.names=F)
write.csv(inv, paste0(outDir,"vowel_langs.csv"), row.names=F)

##################################
### GET CORRELATIONS BTW MEANS ###
##################################

# populate matrix with mean measure for each vowel (col) and language (row)
f1 <- vowelValues(myMeans, mean_f1)
f2 <- vowelValues(myMeans, mean_f2)

# count up the number of times each vowel pair occurs across readings
vowelCount1 <- myMeans %>% group_by(lang, vowel) %>% summarise(count=length(vowel))
vowelCount2 <- vowelCount1
vowelCooccur <- merge(vowelCount1, vowelCount2, by="lang")
vowelCooccur <- vowelCooccur %>% group_by(vowel.x, vowel.y) %>% summarise(count=length(vowel.x))

# get correlations
f1.cors <- getCors(f1, 3)
f2.cors <- getCors(f2, 3) 

# get p-values
f1.cors <- getPvals(f1.cors, f1)
f2.cors <- getPvals(f2.cors, f2) 

# make table nice for publication
f1.cors.latex <- latexCors(f1.cors, 9)
f2.cors.latex <- latexCors(f2.cors, 9) 

###############################
### SAVE CORRELATION TABLES ###
###############################

# still saving as csv because of mismatching number of significant digits; post-process elsewhere
write.csv(f1.cors.latex, file=paste0(outDir, "f1cors_", DATASET,".csv"), row.names=F, quote=F)
write.csv(f2.cors.latex, file=paste0(outDir, "f2cors_", DATASET,".csv"), row.names=F, quote=F)

####################
### MAKE FIGURES ###
####################

f1.iu <- prepFigure(myMeans, f1, sd_f1, "i", "u")

cor.test(~i + u, f1)

pdf(paste0(outputDir, "voxclamantis_iu_F1_ellipse_ERB_", DATASET, ".pdf"))
makeFigure(f1.iu, "F1 /i/", "F1 /u/", "r = 0.79, p < 0.001", TRUE) 
dev.off()








