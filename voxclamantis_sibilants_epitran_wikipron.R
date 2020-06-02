# Analyze VoxClamantis v1.0 Sibilants
# Created by Eleanor Chodroff
# 1 June 2020

# This script takes as input the Epitran and Wikipron sibilant *info.csv and *sibilants.csv text files, the per_utt_mcd text files (scores), and reading_info.csv

# The script returns the numbers reported in the paper (# of languages, # of language families, counts of sibilants, correlation of mid-frequency peak for /s/ and /z/) and the figure

library(tidyverse)
library(ggExtra)
library(ggthemes)
library(ggforce)
library(jsonlite)
library(data.table)
library(xtable)

#################
### CHANGE ME ###
#################

DATASET <- "epiwiki"
epiDir <- "./epitran_sibilants/"
wikiDir <- "./wikipron_sibilants/"
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
readSibFiles <- function(fileDir, file_list) {
	d <- read_csv(paste0(fileDir, file_list[i]))
	d$erbfreq <- hz2erb(d$peak3k7k)
	d$stim <- NULL
    
    sibinfo <- gsub("sibilants.csv", "sibInfo.csv", file_list[i])
    d2 <- read_csv(paste0(fileDir, sibinfo))
    d <- left_join(d, d2, by = c("file", "trial", "sib"))
    d <- left_join(d, inventory, by = "lang")
    d <- d %>% select(c(file, sib, lang), everything())  
    return(d)  
}

# remove non-sibilants: only necessary for unitran
removeNonSibs <- function(x) {
	x <- filter(x, !sib %in% c('SIL', 'SHADDA', 'SOFTSIGN', 'SEMICOLON', 'SIX'))
	return(x)
}

# only keep 's' and 'z' sibilants
keepSZ <- function(x) {
	x <- filter(x, sib %in% c('s','z'))
	return(x)
}

# get number of sibilant categories per reading
getNumSibCats <- function(x) {
	numSibCats <- length(unique(x$sib))
	return(numSibCats)
}

# get number of vowels by row count
getNumSibs <- function(x, newcolname) {
	count <- sapply(x, nrow)
	count <- enframe(count)
	colnames(count) <- c("lang", newcolname)
	return(count)
}

# exclude utterances with MCD score below some value (6 used below)
removeUttMcd <- function(x, cutoff) {
	x <- filter(x, mcd < cutoff)
}

# remove outliers
removeOutliers <- function(x) {
	x <- filter(x, count > 50)
	x <- subset(x, erbfreq < meanfreqmplusSD)
	x <- subset(x, erbfreq > meanfreqmminSD)
	x <- subset(x, dur < meandurplusSD)
	x <- subset(x, dur > meandurminSD)
}

####################
### READ IN DATA ###
####################

epi_list <- list.files(path = epiDir, pattern = "*sibilants.csv")
wiki_list <- list.files(path = wikiDir, pattern = "*sibilants.csv")
scores_list <- list.files(path = scoresDir, pattern = "*.scores")

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
	langDFs[[langID]] <- readSibFiles(epiDir, epi_list)
	langDFs[[langID]] <- left_join(langDFs[[langID]], scores[[langID]], by = c("file" = "file"))
}

# READ IN WIKIPRON FILES
for (i in 1:length(wiki_list)) {
	langID <- substr(wiki_list[i], 1, 6)
	langDFs[[langID]] <- readSibFiles(wikiDir, wiki_list)
	langDFs[[langID]] <- left_join(langDFs[[langID]], scores[[langID]], by = c("file" = "file"))
}

langDFs_orig <- langDFs

###################################
######## OUTLIER EXCLUSION ########	     
###################################

# only keep 's' and 'z' sibilants
langDFs <- lapply(langDFs, keepSZ)

# remove languages w/o /s/ *and* /z/
numSibCats <- sapply(langDFs, getNumSibCats)
lessThanTwo <- names(numSibCats[which(numSibCats < 2)])
	
for (i in 1:length(lessThanTwo)) {
	langDFs[[lessThanTwo[i]]] <- NULL
	}

# count number of sibilants
numSibs <- getNumSibs(langDFs, "orig")

# prep readings with MCD > 8 for removal
removeme <- inventory$lang[which(inventory$MCD1 > 8)]
removeme <- removeme[removeme %in% names(langDFs)]

# remove files above file-mcd cutoff
for (i in 1:length(removeme)) {
	langDFs[[removeme[i]]] <- NULL
}

# count remaining sibilants
count <- getNumSibs(langDFs, "file_mcd")
numSibs <- left_join(numSibs, count, by = "lang") 

# remove utterances below mcd cutoff
langDFs <- lapply(langDFs, removeUttMcd, cutoff = 6)

# count remaining sibilants
count <- getNumSibs(langDFs, "utt_mcd")
numSibs <- left_join(numSibs, count, by = "lang") 
 
# get mean, sd, count, and mean+/-2sd per sibilant for mid-frequency peak and duration
sibInfo <- list()

for (i in 1:length(langDFs)) {
	langID <- names(langDFs[i])
    sibInfo[[langID]] <- langDFs[[langID]] %>% group_by(lang, sib) %>% 
    summarise(meanfreqm = mean(erbfreq), sdfreqm = sd(erbfreq), meandur = mean(dur), sddur = sd(dur), count = length(sib))
    
    	sibInfo[[i]]$meanfreqmplusSD <- sibInfo[[i]]$meanfreqm + 2 * sibInfo[[i]]$sdfreqm
	sibInfo[[i]]$meanfreqmminSD <- sibInfo[[i]]$meanfreqm - 2 * sibInfo[[i]]$sdfreqm
	sibInfo[[i]]$meandurplusSD <- sibInfo[[i]]$meandur + 2 * sibInfo[[i]]$sddur
	sibInfo[[i]]$meandurminSD <- sibInfo[[i]]$meandur - 2 * sibInfo[[i]]$sddur
	
    langDFs[[langID]] <- left_join(langDFs[[langID]], sibInfo[[langID]], by = c('lang', 'sib'))
}

# remove outliers
langDFs <- lapply(langDFs, removeOutliers)

### REMOVE READINGS WITH NO DATA
numSibCats <- sapply(langDFs, getNumSibCats)
lessThanTwo <- names(numSibCats[which(numSibCats < 2)])
	
for (i in 1:length(lessThanTwo)) {
	langDFs[[lessThanTwo[i]]] <- NULL
}

# count remaining sibilants
count <- getNumSibs(langDFs, "outlier_after")
numSibs <- left_join(numSibs, count, by="lang") 

numSibs$retained_mcd <- numSibs$utt_mcd / numSibs$orig
numSibs$retained_total <- numSibs$outlier_after / numSibs$orig
numSibs$retained_mcd <- ifelse(numSibs$retained_mcd < 0.01, "NA", numSibs$retained_mcd)
numSibs$retained_mcd <- as.numeric(as.character(numSibs$retained_mcd))

numSibsTable <- summary(numSibs)
summary(numSibs)

############################
### SAVE RETENTION STATS ###
############################

print(xtable(numSibsTable), file = paste0(outputDir, DATASET, "-retention-sibilants.tex"))
write_csv(numSibs, paste0(outputDir, DATASET, "-retention-sibilants-all.csv"))

#############################
### GET DESCRIPTIVE STATS ###
#############################

meansPerSibPerLang <- list()

for (i in 1:length(langDFs)) {
  langID <- names(langDFs[i])  
  meansPerSibPerLang[[langID]] <- langDFs[[langID]] %>% group_by(sib) %>% summarise(
  meanpeak = mean(erbfreq), sdpeak = sd(erbfreq), count = length(sib), lang = langID, family = family[1])
}

myMeans <- do.call(rbind.data.frame, meansPerSibPerLang)
mymeanpeak <- myMeans %>% select(-c(sdpeak, count)) %>% pivot_wider(names_from = sib, values_from = meanpeak)
mysdpeak <- myMeans %>% select(-c(meanpeak, count)) %>% pivot_wider(names_from = sib, values_from = sdpeak)
mycount <- myMeans %>% select(-c(meanpeak, sdpeak)) %>% pivot_wider(names_from = sib, values_from = count)

peak.sz <- merge(mymeanpeak, mysdpeak, by = c("lang", "family"))
peak.sz <- merge(peak.sz, mycount, by = c("lang", "family"))
colnames(peak.sz) <- c("lang", "family", "s", "z", "sd_s", "sd_z", "count_s", "count_z")

cor.test(~s + z, peak.sz)

##############################
### SAVE IMPORTANT NUMBERS ###
##############################

sibNumbers <- matrix(ncol=1, nrow=12)
sibNumbers[1,1] <- length(langSZ)
sibNumbers[2,1] <- sum(subset(nSibs, sib == "s")$count)
sibNumbers[3,1] <- sum(subset(nSibs, sib == "z")$count)
sibNumbers[4,1] <- mean(nSibsLang$count)
sibNumbers[5,1] <- median(nSibsLang$count)
sibNumbers[6,1] <- min(nSibsLang$count)
sibNumbers[7,1] <- max(nSibsLang$count)
sibNumbers[8,1] <- quantile(nSibsLang$count, 0.25)
sibNumbers[9,1] <- quantile(nSibsLang$count, 0.75)
sibNumbers[10, 1] <- nrow(subset(peak.sz, !is.na(z)))
sibNumbers[11, 1] <- cor.test(~s + z, mymeanpeak)$estimate
sibNumbers[12, 1] <- cor.test(~s + z, mymeanpeak)$p.value

sibNumbers <- data.frame(sibNumbers)
colnames(sibNumbers) <- c("value")
rownames(sibNumbers) <- c("# languages", "# phon{s}", "# phon {z}", "mean tokens", "median tokens", "min tokens", "max tokens", "1q tokens", "3q tokens", "# languages with both phon{s} and phon{z}", "sz cor", "sz pval")

langs <- unique(filter(peak.sz, !is.na(z))$lang)
inv <- filter(inventory, lang %in% langs)

lang_families <- inv %>% group_by(family) %>% summarise(count = n()) %>% arrange(desc(count))

write.table(sibNumbers, paste0(outputDir, "sibNumbers_", DATASET, ".txt"), sep = "\t", quote = F)
write.csv(inv, paste0(outputDir, "sibilant_langs.csv"), row.names = F)
write.csv(lang_families, paste0(outputDir, "sibilant_lang_families.csv"), row.names = F)

###################
### MAKE FIGURE ###
###################

corLabel <- "r = 0.87, p < 0.001"
scale_min <- 25
scale_max <- 31
diff <- scale_max - scale_min

p <- ggplot(peak.sz, aes(x = s, y = z)) + geom_point(size = 0, color = "slategray", alpha = 0) + 
	geom_ellipse(aes(x0 = s, y0 = z, a = sd_s * 0.1, b = sd_z * 0.1, angle = 0), color = "slategray",fill = "slategray",alpha = 0.04) + 
	geom_smooth(method = lm, size = 0.75, se = TRUE, color = "black") + 
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) + 
	theme_few(20) + xlab("ERB") + ylab("ERB") + 
	scale_x_continuous(limit = c(scale_min, scale_max)) + scale_y_continuous(limit = c(scale_min, scale_max)) + 
	annotate("text", family = "serif", x = scale_min + 0.8 * diff, y = scale_min, size = 7, label = corLabel, fontface = "italic") + 
	annotate("text", family = "serif", x = scale_min + (diff / 2), y = scale_max, size = 9, label = "mid-freq peak /s/") + 
	annotate("text", family = "serif", x = scale_max, y = scale_min + (diff / 2), size = 9, label = "mid-freq peak /z/", angle = -90) + 
	theme(legend.position = c(0.15,0.87), text = element_text(family = "serif", size = 22))

ggMarginal(p, type = "histogram", xparams = list(color = "slategray", fill = "ghostwhite"), yparams = list(color = "ghostwhite", fill = "slategray"))

quartz.save(type = "pdf", file = paste0(outputDir, "voxclamantis_sz_ellipse_", DATASET, "_erb.pdf"))

