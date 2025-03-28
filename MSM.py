import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from scipy import stats


def load_DEG_profile(path: str, gene_index_col: str, DEG_format: str):
    # imports .csv data and indexes based on gene index (user can make it ensembl id, gene symbol, etc)
    df = pd.read_csv(path, low_memory=False, index_col=gene_index_col)[[DEG_format]]
    return df


def map_to_binary(profile, metric: str, threshold: float, path: str):
    # maps a metric to binary depending on the input threshold
    profile['Binary_Notation'] = np.where(profile[metric] > threshold, 1,
                                          np.where(profile[metric] < -threshold, -1, 0))

    profile = profile['Binary_Notation']
    profile.name = path.split(".")[0]

    # returns the proptions of each map
    map_proportions = (profile.value_counts() / len(profile)) * 100
    map_proportions.name = path.split(".")[0]

    return profile, map_proportions


def prepare_DEG_Profiles(DEG_profiles_path: list, DEG_formats: list, gene_index_col: str):
    DEG_Profiles = []
    map_data = []
    for path, DEG_format in zip(DEG_profiles_path, DEG_formats):
        df = load_DEG_profile(path=path, gene_index_col=gene_index_col, DEG_format=DEG_format[0])
        profile, map_proportions = map_to_binary(profile=df, metric=DEG_format[0], threshold=DEG_format[1], path=path)

        DEG_Profiles.append(profile)
        map_data.append(map_proportions)

    DEG_Profiles = pd.concat(DEG_Profiles, axis=1)
    map_data = pd.concat(map_data, axis=1)

    # computes basic statistics
    map_data['Mean'] = map_data.mean(axis=1)
    map_data['Std_Dev'] = map_data.drop(columns=['Mean']).std(axis=1)
    # map_data.to_csv(f'{datetime.now().strftime("%Y_%m_%d_%H_%M_%S")}_MAPPING_PROPORTIONS.csv')

    return DEG_Profiles, map_data


def do_MSM(DEG_Profiles):
    # Overlap is when both profiles have DEG != 0
    DEG_Profiles['Overlap'] = np.where((DEG_Profiles.iloc[:, 0] != 0) & (DEG_Profiles.iloc[:, 1] != 0), True, False)

    DEG_Overlap = (DEG_Profiles['Overlap'].value_counts()[True] / len(DEG_Profiles))

    # Concordant is when overlap DEG is in the same direction
    DEG_Profiles['Concordant'] = np.where(
        (DEG_Profiles.iloc[:, 0] != 0) & (DEG_Profiles.iloc[:, 1] != 0) &
        (DEG_Profiles.iloc[:, 0] == DEG_Profiles.iloc[:, 1]),
        True,
        False
    )

    # Disconcordant is when overlap DEG is in the opposite direction
    DEG_Profiles['Disconcordant'] = np.where(
        (DEG_Profiles.iloc[:, 0] != 0) & (DEG_Profiles.iloc[:, 1] != 0) &
        (DEG_Profiles.iloc[:, 0] != DEG_Profiles.iloc[:, 1]),
        True,
        False
    )

    Concordance = (
            (DEG_Profiles['Concordant'].value_counts()[True] - DEG_Profiles['Disconcordant'].value_counts()[True]) /
            DEG_Profiles['Overlap'].value_counts()[True])

    return Concordance, DEG_Overlap


def do_MSM_permutation(DEG_Profiles, n_permutations: int, test_column_index: int):
    # saves the observed concordance and overlap
    observed_concordance, DEG_Overlap = do_MSM(DEG_Profiles)

    # creates a distribution of concordance values by shuffling the test vector
    permutation_matrix = DEG_Profiles.copy()
    permuted_concordance_values = np.array([observed_concordance])
    for iteration in range(n_permutations):
        permutation_matrix.iloc[:, test_column_index] = np.random.permutation(DEG_Profiles.iloc[:, test_column_index])

        permuted_concordance, _ = do_MSM(permutation_matrix)
        permuted_concordance_values = np.append(permuted_concordance_values, permuted_concordance)

    # converts into a standard normal distribution and returns observed p-value (always first bc this method only appends)
    permuted_concordance_values = (
                                              permuted_concordance_values - permuted_concordance_values.mean()) / permuted_concordance_values.std()
    p_value = 2 * (1 - stats.norm.cdf(abs(permuted_concordance_values[0])))

    return DEG_Overlap, observed_concordance, p_value


def main(DEG_profiles_path: list, gene_index_col: str,n_permutations: int,
         test_column_index: int, DEG_formats: list):

    df, map_data = prepare_DEG_Profiles(DEG_profiles_path=DEG_profiles_path,
                              DEG_formats=DEG_formats,
                              gene_index_col=gene_index_col,
                              )

    DEG_Overlap, observed_concordance, p_value = do_MSM_permutation(DEG_Profiles=df,
                                                                    n_permutations=n_permutations,
                                                                    test_column_index=test_column_index)

    return DEG_Overlap, observed_concordance, p_value, map_data
