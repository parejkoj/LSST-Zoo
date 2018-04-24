# ZooMeR -- Zooniverse to Machine LeaRning

#Code Contained in this file takes the output of a zooniverse classification run and
#changes it to a form usable by ML algorithms

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import Imputer
from sklearn.cross_validation import train_test_split
from sklearn.metrics import mean_squared_error
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy
import json

def ZooMeR(stack_features,zoo_classifications):

    #Retrieve Parameters
    params = pd.read_csv(stack_features)
    del params['coord_dec'] #Irrelevant (Hopefully) and poorly formatted for ML
    del params['coord_ra']
    del params['id']
    del params['parent']

    #Retrieve Training Set Classifications
    csv_path = "../classifications/"+zoo_classifications
    clsfn = pd.read_csv(csv_path)
    subset = clsfn.loc[clsfn['workflow_name'] == "Difference Imaging Classifier"]
    np.shape(subset)
    im_class = []
    for item, row in subset.iterrows():
        s_data = json.loads(row.subject_data) #Subject Data
        s_data = s_data.get(s_data.keys()[0])
        im=s_data.get(s_data.keys()[0])[47:-4]
        a_data = json.loads(row.annotations)[0] #Annotations
        classification = a_data['value']
        im_class.append([int(im),classification])
    im_class = sorted(im_class)

    #Generate Data Frame that holds classifications and respective measured quantities
    #This is done as a dataframe is a more cohesive object for data analysis
    d = {}
    for key in params.keys():
        d[key] = []
    d["Classification"]=[]
    im_data = params.set_index('image').T.to_dict()
    for im in im_class:
        d["Classification"].append(im[1])
        d["image"].append(im[0])
        dat_dict = im_data.get(im[0])
        for key in dat_dict:
            d[key].append(dat_dict.get(key))
    df = pd.DataFrame(d)
    #Read in the features in the data frame, filter for columns relevant for ML
    features = df.columns.tolist()

    features = [c for c in features if c not in ["image", "Classification"]]
    target = "Classification"  #Predict on Classifications

    #Cut 'useless' features and problematic features
    for feature in features:
        if df[feature].isnull().all():
            del df[feature]
        elif np.mean(df[feature]) == np.inf:
            del df[feature]
        elif "flag" in feature: #Flags don't contribute to ML based on initial testing
            del df[feature]

    features = df.columns.tolist()

    features = [c for c in features if c not in ["image", "Classification"]]

    #Imputer
    imp = Imputer(missing_values='NaN', strategy='median', axis=0, verbose = 1)
    imp.fit(df[features])
    features_imp = imp.transform(df[features])
    return features_imp, df
