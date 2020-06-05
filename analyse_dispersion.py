from os import listdir, path
import pathlib
import argparse
import math
import pandas as pd
import numpy as np
import scipy.stats


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--info-file', type=str,
        help='The file where to save results.')
    parser.add_argument(
        '--formants-type', type=str,
        default='erb', choices=['erb', 'hz'],
        help='The unit of measurement for formants.')
    args = parser.parse_args()
    args.entropy_col = 'entropy' if args.formants_type == 'erb' else 'entropy_hertz'
    return args


def read_csv(fname, entropy_col):
    df = pd.read_csv(fname, sep='\t', index_col=0)
    df = df[df.n_samples > 50]
    df['entropy_prob'] = df[entropy_col] * df['prob']
    df.dropna(inplace=True)
    return df


def get_correlation(df, x, y):
    spearman, spearman_p = scipy.stats.spearmanr(df[x], df[y])
    pearson, pearson_p = scipy.stats.pearsonr(df[x], df[y])

    return spearman, spearman_p, pearson, pearson_p


def main():
    args = get_args()
    df = read_csv(args.info_file, args.entropy_col)

    df_lang = df.groupby('lang').agg('mean')
    df_lang['n_vowels'] = df.groupby('lang').agg('count')['vowel']
    df_lang['lang_entropy'] = df.groupby('lang').agg('sum')['entropy_prob']

    corrs_entropies = get_correlation(df_lang, 'lang_entropy', 'n_vowels')

    print('Correlations between dispersion and vowel inventory size:')
    print('\tspearman:', corrs_entropies[0])
    print('\tspearman_p:', corrs_entropies[1])
    print('\tpearson:', corrs_entropies[2])
    print('\tpearson_p:', corrs_entropies[3])


if __name__ == '__main__':
    main()
