# analysis

voxclamantis_vowels_epitran_wikipron.R: analysis script for vowel midpoint F1 and F2

This script takes as input the Epitran and WikiPron midpoint F1 and F2 files, the per-utt-mcd scores files, reading_info.csv (inventory), and returns counts of tokens, families, and languages presented in the paper, correlation tables, and the correlation scatterplot. It optionally saves text files of the midpoint F1 and F2 following outlier exclusion that serves as input to the Python dispersion analysis.

voxclamantis_sibilants_epitran_wikipron.R: analysis script for sibilant mid-frequency peak

This script takes as input the Epitran and WikiPron sibilant info.csv and sibilant.csv files, the per-utt-mcd scores files, reading_info.csv (inventory), and returns counts of tokens, families, and languages presented in the paper, the correlation of mean mid-frequency peak /s/ and /z/, and the correlation scatterplot. 
