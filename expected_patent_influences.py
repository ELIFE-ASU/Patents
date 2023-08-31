import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
from sklearn.linear_model import LinearRegression
import scipy.stats as stats
import ast
import math
import datetime
from dateutil.relativedelta import relativedelta
import random


def classification_to_list(classification):
    """ Turns a string of classifications into a list

    Args:
        classification (str): USPTO classifications, guarenteed to be in a list

    Returns:
        list: list of different USPTO classifications which are applied to a particular patent
    """
    try:
        return ast.literal_eval(classification)
    except ValueError:
        return []
    except SyntaxError:
        return classification.strip('][').split(',')


def filter_classification(classification):
    """ Removes subheadings of patents

    Args:
        classification (str): a particular USPTO classification, in the form "classification/subheading"

    Returns:
        str: USPTO classification, without the subheading or backslash
    """
    try:
        return classification.split("/")[0]
    except AttributeError:
        return ""


def test_results_df(results_df):
    """ Asserts results dataframe expanded as expected

    Args:
        results_df (pandas dataframe): holds patent data regarding authors, assignees, and classifications
    """

    assert len(results_df) == 1247658, "Classification Filter failed"

    print("Classification results expanded", end="\n\n")


def get_results_df(fp):
    """ Load & expand results dataframe for analysis

    Args:
        fp (str): filepath to authors/assignee/classification results

    Returns:
        pandas dataframe: expanded pandas dataframe containing A/A/C results
    """
    tqdm.pandas()

    results_df = pd.read_csv(fp)

    print(f"Expanding classification results...")
    results_df["classification"] = results_df["classification"].progress_apply(
        classification_to_list)

    results_df = results_df.explode("classification")

    results_df["filtered_classification"] = results_df["classification"].apply(
        filter_classification)

    test_results_df(results_df)

    return results_df


def patent_stats_per_assignee(results_df):
    """ Get average and standard deviation of patents per assignee

    Args:
        results_df (pandas dataframe): full expanded results dataframe  

    Returns:
        float: average patents per assignee
        float: standard deviation of patents per assignee distribution
    """
    patents_per_assignee = []
    assignee_lengths = []

    time_format = "%Y-%m-%d"

    print("Calculating patents per assignee statistics & time per patents...")
    for assignee in tqdm(list(results_df["assignees"].unique())):
        sub_df = results_df[results_df["assignees"] == assignee]
        if len(sub_df) > 10:
            patents_per_assignee.append(len(sub_df))

            start_date = datetime.datetime.strptime(min(sub_df["date"]),
                                                    time_format)
            end_date = datetime.datetime.strptime(max(sub_df["date"]),
                                                  time_format)

            assignee_lengths.append((end_date.year - start_date.year) * 12 +
                                    end_date.month - start_date.month)

    avg_patents_per_assignee = np.mean(patents_per_assignee)
    stdev_patents_per_assignee = np.std(patents_per_assignee)

    print(
        f"Average patents per assignee: {round(avg_patents_per_assignee,3)}"
        f", stdev patents per assignee: {round(stdev_patents_per_assignee,3)}\n"
    )

    avg_months_per_assignee = np.mean(assignee_lengths)
    stdev_months_per_assignee = np.std(assignee_lengths)

    print(f"Average time per assignee: {round(avg_months_per_assignee,3)}"
          f", stdev time per assignee: {round(stdev_months_per_assignee,3)}\n")

    return avg_patents_per_assignee, stdev_patents_per_assignee, avg_months_per_assignee, stdev_months_per_assignee


def time_per_assignee(results_df):
    """ Calculating the average time an assignee exists for

    Args:
        results_df (pandas dataframe): full expanded results dataframe  

    Returns:
        float: average time (in months) an assignee sticks arount
        float: standard deviation of time distribution
    """
    assignee_lengths = []

    print("Calculating average time per assignee...")
    for assignee in tqdm(list(results_df["assignees"].unique())):
        sub_df = results_df[results_df["assignees"] == assignee]
        if len(sub_df) > 10:
            format = "%Y-%m-%d"
            start_date = datetime.datetime.strptime(min(sub_df["date"]), format)
            end_date = datetime.datetime.strptime(max(sub_df["date"]), format)

            assignee_lengths.append((end_date.year - start_date.year) * 12 +
                                    end_date.month - start_date.month)

    avg_months_per_assignee = np.mean(assignee_lengths)
    std_months_per_assignee = np.std(assignee_lengths)

    print(f"Average time per assignee: {round(avg_months_per_assignee,3)}"
          f", stdev time per assignee: {round(std_months_per_assignee,3)}\n")

    return avg_months_per_assignee, std_months_per_assignee


def build_profiles(avg_patents_per_assignee, stdev_patents_per_assignee,
                   num_profiles):
    """ Build random assignee profiles, each profile having x patents in it

    Args:
        avg_patents_per_assignee (float): average patents per assignee
        stdev_patents_per_assignee (float): standard deviation of patents per assignee distribution
        num_profiles (int): number of profiles to build

    Returns:
        list: list of n profiles, each with a size based on previously calculated statistics
    """
    print("Generating random profiles...")
    profile_patent_size = []

    lower = 1
    upper = 10000  #upper value set to be arbitarily high

    profile_patent_size = stats.truncnorm.rvs(
        (lower - avg_patents_per_assignee) / stdev_patents_per_assignee,
        (upper - avg_patents_per_assignee) / stdev_patents_per_assignee,
        loc=avg_patents_per_assignee,
        scale=stdev_patents_per_assignee,
        size=num_profiles)

    profile_patent_size = [int(math.floor(x)) for x in profile_patent_size]

    print("Random profile length:", len(profile_patent_size), ";",
          profile_patent_size[0:10], "\n")

    return profile_patent_size


def random_dates2(start, end, n, unit='D', seed=None):
    """
    Calculates random dates in YYYY-MM-DD format 
    (from https://stackoverflow.com/questions/50559078/generating-random-dates-within-a-given-range-in-pandas)
    """
    ndays = (end - start).days + 1
    return start + pd.to_timedelta(np.random.randint(0, ndays, n), unit=unit)


def get_starting_time(months):
    """ Generate a random starting time for a given "profile" - this is the date of the first patent

    Samples are guaranteed to be between Jan 1976 & (Dec 2022 - number of months)

    Args:
        months (int): number of months that a particular profile exists for

    Returns:
        (list of pyDateTime objects): randomly generated samples, all in pyDateTime format with year, month, and date
    """
    low_range = datetime.datetime(year=1976, month=1, day=1)
    high_range = datetime.datetime(year=2022, month=12, day=31) - relativedelta(
        months=months, day=0)

    #Generate a random starting date, between the low & high range
    start_date = random_dates2(low_range, high_range, 1,
                               seed=random.randint).to_pydatetime()[0]

    ## Randomly generate num_patents numbers of timestamps
    return start_date


def get_timestamps(avg_months_per_assignee, stdev_months_per_assignee,
                   num_profiles):
    ### For each profile, calculate starting and ending time

    # List of list data structure: [[t_start, t_end], [t_start, t_end], ..., [t_start, t_end]],
    # t is a timestamp, start/end is the starting and ending timestamp of that particular "company" profile

    print("Generating a set of timestamps for each profile...")
    profile_timestamps = []

    #Sample number of months
    lower = 1
    upper = 556  #upper value set to be arbitarily high

    num_months = stats.truncnorm.rvs(
        (lower - avg_months_per_assignee) / stdev_months_per_assignee,
        (upper - avg_months_per_assignee) / stdev_months_per_assignee,
        loc=avg_months_per_assignee,
        scale=stdev_months_per_assignee,
        size=num_profiles)

    num_months = [math.floor(x) for x in num_months]

    for months in tqdm(num_months):
        #Sample number of compounds per patent
        starting_timestamp = get_starting_time(months)
        ending_timestamp = starting_timestamp + relativedelta(months=months,
                                                              day=0)

        profile_timestamps.append((starting_timestamp, ending_timestamp))

    print("Random profile timestamps:", len(profile_timestamps), ";",
          profile_timestamps[0:10], "\n")

    return profile_timestamps


def get_all_profile_patent_timestamps(random_profile_sizes,
                                      random_profile_timestamps):
    ### For each profile, generate timestamps of all patents

    # Data structure: list of lists, each sublist with n timestamps (n from profile_patent_size)
    print("Building timestamps for all patents...")
    profile_patent_timestamps = []

    # Loop through patent sizes (number of patents per company), and starting & ending point
    for i in tqdm(range(len(random_profile_sizes))):
        #Generate a list of timestamps for patents, uniformly sampled from within the starting & ending point of the "company"
        timestamps = random_dates2(random_profile_timestamps[i][0],
                                   random_profile_timestamps[i][1],
                                   random_profile_sizes[i],
                                   seed=random.randint).to_pydatetime()

        profile_patent_timestamps.append(timestamps)

    print(f"All {len(profile_patent_timestamps)} patent timestamps built.")
    print(f"0th profile length: {len(profile_patent_timestamps[0])}, {profile_patent_timestamps[0][0:10]}\n")

    return profile_patent_timestamps


def get_avg_MA_per_timestamped_patent(timestamp, MA_month_avg_dict,
                                      MA_month_std_dict, n_cpds):
    """ Calculates a list of expected MAs from a given timestamp

    Args:
        timestamp (datetime.datetime object): Specific date of a random "patent"
        MA_df (pandas dataframe): holds MA avg & std dev values at each month
        n_cpds (int): number of compounds in a given patent
    """
    #Get year & month of timestamp
    month = timestamp.month
    if month < 10:
        month = "0" + str(timestamp.month)
    else:
        month = str(timestamp.month)

    MA_avg = MA_month_avg_dict[str(timestamp.year) + "-" + month]
    MA_std = MA_month_std_dict[str(timestamp.year) + "-" + month]

    lower = 1
    upper = 1000  #upper value set to be arbitarily high

    MAs = stats.truncnorm.rvs((lower - MA_avg) / MA_std,
                              (upper - MA_avg) / MA_std,
                              loc=MA_avg,
                              scale=MA_std,
                              size=n_cpds)
    #MAs = np.random.normal(avg_cpds_per_patent, stdev_cpds_per_patent, n_cpds)

    return np.mean(MAs)


def get_MAs_of_patents(profile_patent_timestamps, avg_cpds_per_patent,
                       stdev_cpds_per_patent, MA_month_avg_dict,
                       MA_month_std_dict):
    ### Get timestamps of compounds - multiply each timestamp by the number of compounds within a patent

    # Data structure - list of lists, each sublist with n * [c0, c1, c2, ...] timestamps
    # (n = number of patents, [c0, c1, ..] = number of compounds per patent)

    profile_avg_MAs = []
    profile_cpds_per_patent = []

    print("Finding MA values per profile...")
    for profile in tqdm(profile_patent_timestamps):
        avg_MAs = []
        n_cpds_per_profile = []
        for timestamp in profile:

            #Sample number of cpds
            lower = 1
            upper = 10000  #upper value set to be arbitarily high

            n_cpds = stats.truncnorm.rvs(
                (lower - avg_cpds_per_patent) / stdev_cpds_per_patent,
                (upper - avg_cpds_per_patent) / stdev_cpds_per_patent,
                loc=avg_cpds_per_patent,
                scale=stdev_cpds_per_patent,
                size=1)[0]
            n_cpds_per_profile.append(math.floor(n_cpds))

            avg_MAs.append(
                get_avg_MA_per_timestamped_patent(timestamp, MA_month_avg_dict,
                                                  MA_month_std_dict,
                                                  math.floor(n_cpds)))
            
        # print(f"\tTesting: profile size = {len(profile)}, avg MAs[0:10] = {avg_MAs[0:10]}")
        # print(f"\tn_cpds_per_profile = {len(n_cpds_per_profile)}, avg_cpds_per_patent: {avg_cpds_per_patent}")

        profile_avg_MAs.append(avg_MAs)
        profile_cpds_per_patent.append(n_cpds_per_profile)

    print(f"Found {len(profile_avg_MAs)} MA values: {profile_avg_MAs[0][0:10]}\n")

    return profile_avg_MAs, profile_cpds_per_patent


def get_deltaMA(profile_avg_MAs, profile_patent_timestamps):
    ### Get deltaMA profile of expected data

    print("Fnding deltaMA slopes...")
    linear_regressor = LinearRegression()

    profile_MA_slopes = []

    for i in tqdm(range(len(profile_avg_MAs))):
        X = np.array([x.toordinal() for x in profile_patent_timestamps[i]
                     ]).reshape(-1, 1)
        Y = np.array(profile_avg_MAs[i]).reshape(-1, 1)
        #try:
        reg = linear_regressor.fit(X, Y)
        # Y_pred = linear_regressor.predict(X)

        # r2 = reg.score(X, Y)

        #Multiply the occurances of the slope by the number of compounds per "assignee"
        profile_MA_slopes.extend([reg.coef_[0][0] * (max(X) - min(X))[0]
                                 ])  #*sum(profile_cpds_per_patent[i]))

    print(
        f"Found {len(profile_MA_slopes)} deltaMA values; {profile_MA_slopes[0:10]}\n"
    )

    return profile_MA_slopes


def get_positive_slopes(profile_MA_slopes):
    """ Find the number of positive slopes from a particular sample

    Args:
        profile_MA_slopes (list): list of deltaMA slopes
    """
    positive_count = 0
    for MA_slope in profile_MA_slopes:
        if MA_slope > 0:
            positive_count += 1

    positive_slope_percentage = float(positive_count) / len(profile_MA_slopes)

    print(f"Percent positive slopes: {positive_slope_percentage}")

    return positive_slope_percentage


def main():
    #Load results dataframe and MA data
    results_df = get_results_df("Data/Patents/patent_MA_results.csv")
    MA_df = pd.read_csv("Data/AssemblyValues/MA_df_months_FULL.csv")
    MA_month_avg_dict = dict(zip(MA_df.month, MA_df.avg))
    MA_month_std_dict = dict(zip(MA_df.month, MA_df.stdev))

    #Precalculated (see deltaMA_factors.ipynb)
    avg_cpds_per_patent = 108.41867561794324
    stdev_cpds_per_patent = 272.51858323850996

    avg_patents_per_assignee, stdev_patents_per_assignee, avg_months_per_assignee, stdev_months_per_assignee = patent_stats_per_assignee(
        results_df)

    print(f"Avg Patents Per Assignee: {avg_patents_per_assignee}\n"
            f"Stdev Patents Per Assignee: {stdev_patents_per_assignee}\n"
            f"Avg Months Per Assignee: {avg_months_per_assignee}\n"
            f"Stdev Months Per Assignee: {stdev_months_per_assignee}")
    print("\n" + " - " * 20 + "\n")

    num_profiles = 769

    random_profile_sizes = build_profiles(avg_patents_per_assignee,
                                          stdev_patents_per_assignee,
                                          num_profiles)

    print(f"Random profile sizes: {random_profile_sizes[0:10]}")
    print("\n" + " - " * 20 + "\n")

    random_profile_timestamps = get_timestamps(avg_months_per_assignee,
                                               stdev_months_per_assignee,
                                               num_profiles)

    print(f"Random profile timestamps: {random_profile_timestamps[0:10]}")
    print("\n" + " - " * 20 + "\n")

    profile_patent_timestamps = get_all_profile_patent_timestamps(
        random_profile_sizes, random_profile_timestamps)

    profile_avg_MAs, profile_cpds_per_patent = get_MAs_of_patents(
        profile_patent_timestamps, avg_cpds_per_patent, stdev_cpds_per_patent,
        MA_month_avg_dict, MA_month_std_dict)

    profile_MA_slopes = get_deltaMA(profile_avg_MAs, profile_patent_timestamps)

    percent_positive = get_positive_slopes(profile_MA_slopes)


if __name__ == "__main__":
    main()
