import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from sklearn.linear_model import LinearRegression
import ast
from random import sample
import scipy.stats as stats


#Matplotlib defaults (eventually)
plt.style.use("default")
plt.rcParams.update({'font.size': 22})


#Read in results dataframe
print("Reading in results dataframe")
results_df = pd.read_csv("Data/Patents/patent_MA_results.csv")

def classification_to_list(classification):
    try:
        return ast.literal_eval(classification)
    except ValueError:
        return []
    except SyntaxError:
        return classification.strip('][').split(',')
    
tqdm.pandas()

print("Classification column building")
results_df["classification"] = results_df["classification"].progress_apply(classification_to_list)
#Explode classifications
print("Exploding dataframe")
results_df = results_df.explode("classification")

## Remove backslash and two numbers after it (and put them in a different column, in case they are useful later...)
def filter_classification(classification):
    try:
        return classification.split("/")[0]
    except AttributeError:
        return ""
    
print("Filtering classification")
results_df["filtered_classification"] = results_df["classification"].progress_apply(filter_classification)

plt.figure(figsize=(16,12))

linear_regressor = LinearRegression()

class_sizes = []
class_rsqs = []
class_deltaMAs = []
class_slopes = []

#Alternative graphing (finding slope change in MA)
linear_regressor = LinearRegression()

MA_slopes = []

print("Finding slopes of unique classifications")
for c in tqdm(list(results_df["filtered_classification"].unique())):
    sub_df = results_df[results_df["filtered_classification"] == c]
    if len(sub_df) > 10:

        X = sub_df["date_ordinal"].values.reshape(-1,1)
        Y = sub_df["MA_avg"].values.reshape(-1,1)
        try:
            reg = linear_regressor.fit(X, Y)
            Y_pred = linear_regressor.predict(X)

            #Calculate deltaMA (change in MA of linear regression, taking into consideration negative slopes)
            MA_slopes.append([reg.coef_[0][0] * (max(X) - min(X))[0]]*len(sub_df))
                        
        except ValueError as e:
            pass

del(results_df)

#Sample 80% of assignees, flatten 2d lists

dropout_MA_slopes = []

print("Finding dropout slopes")
for i in tqdm(range(1000)):
    dropout_MA_slopes.append(list(np.concatenate(sample(MA_slopes, round(len(MA_slopes) * 0.8))).flat))

del(MA_slopes)

#Various statistics (mean/median, skewness via skewtest)
def calculate_skewness_stats(dropout_MA_slopes):
    avgs = []
    medians = []
    skews = []
    skewtests = []

    for s in tqdm(dropout_MA_slopes):
        avgs.append(np.mean(s))
        medians.append(np.median(s))
        skews.append(stats.skew(s))
        skewtests.append(stats.skewtest(s))

    return avgs, medians, skews, skewtests

print("Calculating skewness stats")
avgs, medians, skews, skewtests = calculate_skewness_stats(dropout_MA_slopes)

del(dropout_MA_slopes)

for i in range(1,10):
    print("Mean:", round(avgs[i],2), " Median:", round(medians[i],2), 
          " Skew:", skews[i], " Skewtest:", skewtests[i])
    
plt.figure(figsize=(8,8))

plt.plot(np.arange(0, len(avgs), 1), avgs, label="Mean")
plt.plot(np.arange(0, len(avgs), 1), medians, label="Median")
plt.plot(np.arange(0, len(avgs), 1), skews, label="Skew")

plt.legend() #bbox_to_anchor=(1.05,1.025))

plt.xlabel("Sample")
plt.ylabel("MA")

plt.savefig("classification_dropout.png", format="png")
