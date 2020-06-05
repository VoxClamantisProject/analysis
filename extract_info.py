from os import listdir, path
import pathlib
import argparse
import math
import pandas as pd
import numpy as np
from tqdm import tqdm


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--src-path', type=str, required=True,
        help='The path to the source files.')
    parser.add_argument(
        '--tgt-file', type=str, required=True,
        help='The file where to save results.')
    return parser.parse_args()


def get_filenames(filepath):
    filenames = [path.join(filepath, f)
                 for f in listdir(filepath)
                 if path.isfile(path.join(filepath, f))]
    fnames = sorted(filenames)
    fnames = [x for x in fnames if '_formants_mid_fin.csv' in x]
    return fnames


def read_csv(fname):
    df = pd.read_csv(fname, sep=',')

    df.erbF1 = df.erbF1.apply(float)
    df.erbF2 = df.erbF2.apply(float)
    df.f1_mid = df.f1_mid.apply(float)
    df.f2_mid = df.f2_mid.apply(float)

    return df


def get_lang_name(fname):
    return fname.split('/')[-1][:6]


def get_cov(matrix):
    return np.cov(matrix, rowvar=False)


def get_gaussian_entropy(cov, base=2):
    # entropy = multivariate_normal.entropy(cov=cov)
    _, logdet = np.linalg.slogdet(2 * np.pi * np.e * cov)
    return 0.5 * logdet / math.log(base)


def get_vowel_info(df, df_full, lang):
    df = df.copy()
    n_samples = df.shape[0]

    if n_samples > 2:
        formants_erb = df[['erbF1', 'erbF2']].values
        cov = get_cov(formants_erb)
        entropy = get_gaussian_entropy(cov)
        cov = cov.tolist()

        formants_hertz = df[['f1_mid', 'f2_mid']].values
        cov_hertz = get_cov(formants_hertz)
        entropy_hertz = get_gaussian_entropy(cov_hertz)
        cov_hertz = cov_hertz.tolist()
    else:
        cov = None
        entropy = None
        cov_hertz = None
        entropy_hertz = None

    vowel = df.vowel.iloc[0]
    prob = n_samples / df_full.shape[0]

    return {
        'vowel': vowel,
        'lang': lang,
        'prob': prob,
        'cov': cov,
        'entropy': entropy,
        'cov_hertz': cov_hertz,
        'entropy_hertz': entropy_hertz,
        'n_samples': n_samples,
    }


def main():
    args = get_args()
    fnames = get_filenames(args.src_path)
    vowel_info = []
    for fname in tqdm(fnames, desc='Extracting files'):
        tqdm.write('Extracting file %s' % (fname))

        lang = get_lang_name(fname)
        df = read_csv(fname)

        for vowel in df.vowel.unique():
            df_vowel = df[df.vowel == vowel]
            vowel_info += [get_vowel_info(df_vowel, df, lang)]

        df_results = pd.DataFrame(vowel_info)
        df_results = df_results[df_results.n_samples > 1]
        df_results.to_csv(args.tgt_file, sep='\t')


if __name__ == '__main__':
    main()
