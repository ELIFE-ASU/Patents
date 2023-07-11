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


def main():
    #Load & expand results dataframe
    results_df = get_results_df("Data/Patents/patent_MA_results.csv")


if __name__ == "__main__":
    main()
