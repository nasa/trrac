# imports <- data handling
import pandas as pd
import numpy as np

# imports <- path navigation and data accession
import RNASeqUtility
import os
import itertools
import json
from pathlib import PurePath, Path

# imports <- visualization
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# imports <- data set partioning
from random import sample
from random import seed
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.covariance import GraphicalLassoCV

# imports <- performance metrics
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.inspection import permutation_importance

# imports <- data preprocessing
from sklearn.preprocessing import StandardScaler

# imports <- models
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from utility.py_mrmr import mrmr_classif
from utility.py_mrmr import random_forest_classif
from utility import figure_generator

# imports <- runtime
import argparse

from platform import python_version


# set the default plotting dimensions 7" by 5" with 200 dots per inch 
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':200})

"""RNASeqModel
Object representation of actions necessary to load, concatenate, and train on
transcriptomics data sets extracted from NASA GeneLab repository for the 
purposes of the Transalational Radiation Research and Countermeasures (TRRaC)
program for the NASA Space Radiation Element (SRE).
"""
class RNASeqModel:
    def __init__(self):
        # empty dictionary to hold RNA Seq counts data and respective metadata
        self.count_dict = {}
        self.roc_dict = {}
        

    ### DATA SET CURATION ###
    """add_study
    adds GeneLab data set (GLDS) id to the counts dictionary
    
    study_id <- string containing the GLDS ID
    """
    def add_study(self, study_id):
        self.count_dict[study_id] = {}
    
    """add_study_counts
    adds GLDS counts associated with GLDS ID
    
    study_id <- string containing the GLDS ID
    study_counts <- DataFrame containing RNA-Seq counts
    """
    def add_study_counts(self, study_id, study_counts):
        self.count_dict[study_id]['counts'] = study_counts
    
    """add_study_metadata
    adds GLDS study metadata associated with GLDS ID
    
    study_id <- string containing the GLDS ID
    study_metadata <- DataFrame containing RNA-Seq metadata
    """    
    def add_study_metadata(self, study_id, study_metadata):
        self.count_dict[study_id]['metadata'] = study_metadata
    
    """build_combined_df
    concatenates all GLDS studies into a single DataFrame
    
    custom_concat <- string filepath for concatenated DataFrame
    """
    def build_combined_df(self, custom_concat=None, verbose=True):
        # create concatenated DataFrame if one is not provided
        if custom_concat is None:
            # inform the user the order that counts are concatenated
            # this is included to ensure consistent matching with metadata
            if verbose:
                print("Order that counts are concatenated: {}".format(self.count_dict.keys()))
            # iterate over counts dictionary to produce concatenated counts
            self.concat_df = pd.concat([x['counts'] for x in self.count_dict.values()])
        # if externally defined concat dataframe is provided, do not use loaded counts
        else:
            # inform user that concatenated dataframe will not be based off loaded counts
            # in this iteration
            print("Overwriting concatenated dataframe with file {}".format(custom_concat))
            # counts need to be transposed to match orientation of samples along 0-axis (rows)
            self.concat_df = pd.read_csv(custom_concat, index_col=0).transpose()
        
        # filter the DataFrame to remove ERCC spike ins which are only present in a subset of studies
        self.concat_df = self.concat_df.loc[:,self.concat_df.isna().sum(axis=0) == 0]
        # this is a temporary method of encoding that works given the metadata
        self.concat_df['target'] = -1
        # any sample name that contains "F" has corresponded to Space Flight
        # TODO: implement this better based upon a metadata file match
        self.concat_df.loc[self.concat_df.index.str.contains("F"), 'target'] = 1
    
    """load_filenames
    list directory contents from provided directory
    
    data_dir <- string filepath with data contents as flat files
    ftype <- string indicating the type of file being loaded
    verbose <- boolean for printing and troubleshooting
    """    
    def load_filenames(self, data_dir, ftype='', verbose=False):
        p = Path.cwd()
        # generate full file paths
        for d in data_dir:
            p = p / d
        fpaths = sorted([x for x in p.iterdir()])
        fpaths = list(filter(lambda x: (x.name[-4] == '.'), fpaths))
        # long form versus short hand file load outputting
        if verbose:
            for f in fpaths:
                print("Loading the {} file: {}".format(ftype, f.name))
        else:
            print("Loaded files:", ",".join([f.name for f in fpaths]))
        return fpaths
    """feature_subset
    filter counts dataframe by provided list of features
    
    count_dfs <- list of counts data
    feature_list <- string filename of features
    """
    def feature_subset(self, count_dfs, feature_list, verbose=False):
        features = pd.read_csv(feature_list).iloc[:,-1] 
        count_dfs = [x.loc[:,features] for x in count_dfs]
        if verbose:
            print("The subset dimensions of the count dataframes are: ", [x.shape for x in count_dfs])
        return count_dfs
    
    """scale_data
    apply the standard scale to a log transformed input
    
    count_dfs <- list of count data
    """
    def scale_data(self, count_dfs, exclude_test):
        for i in range(len(count_dfs)):
            scaler = StandardScaler()
            df = count_dfs[i]
            if exclude_test:
                test_ids = []
                for k in exclude_test.keys():
                    test_ids.append(exclude_test[k]['id']['SF'])
                    test_ids.append(exclude_test[k]['id']['GC'])
                test_ids = list(itertools.chain(*test_ids))
                valid_ids = count_dfs[i].index.isin(test_ids)
                train_ids = np.arange(df.index.values.shape[0])[~valid_ids]
                df = df.iloc[train_ids]
                
            scaler.fit(np.log2(df+1))
            
            count_dfs[i] = pd.DataFrame(scaler.transform(np.log2(count_dfs[i]+1)), columns=count_dfs[i].columns, index=count_dfs[i].index)
        return count_dfs    
    
    """load_counts
    load in the unnormalizard counts data
    
    data_dir <- array of strings with filepath to data contents as flar files
    metadata_dir <- array of strings with filepath to metadata contents as flat file
    verbose <- boolean for printing and troubleshooting
    scale <- boolean for applying transformation to loaded counts data
    feature_list <- string filepath to feature variables for subsetting
    """
    def load_counts(self, data_dir=['data', 'raw_counts'], metadata_dir=['data', 'single_study_metadata'], verbose=False, scale=False, exclude_test=None, feature_list=None):
        # call load_filenames for import of flat text files
        count_files = self.load_filenames(data_dir, ftype='counts', verbose=verbose)
        metadata_files = self.load_filenames(metadata_dir, ftype='metadata', verbose=verbose)
        
        # remove count files that do not have corresponding metadata file
        keep_files = []

        # perform check of counts against metadata files
        # do not include counts if metadata is not available
        rem_files = []
        for x in count_files:
            if x.name[:-4] in [y.name[:-4] for y in metadata_files]:
                keep_files.append(x)
            else:
                if verbose:
                    print("Removing {} counts due to lack of corresponding metadata".format(x.name))
                else:
                    rem_files.append(x.name)
        if not verbose:
            print("Removed following counts:", ",".join(rem_files))
        count_files = keep_files
        
        # read in the files 
        count_dfs = [pd.read_csv(x, index_col=0).transpose() for x in count_files]
        metadata_dfs = [pd.read_csv(x, delimiter="\t") for x in metadata_files]
        
        ### APPLY OPTIONAL BLOCK PARAMETERS ###
        
        # apply subsetting to feature list if provided
        if feature_list:
            count_dfs = self.feature_subset(count_dfs, feature_list, verbose=verbose)
        
        # fit scale and transform data
        if scale:
            count_dfs = self.scale_data(count_dfs, exclude_test=exclude_test)

        fnames = [x.name.split('.')[0] for x in count_files]
       
       # save counts into object model
        for id, cdf, mdf in zip(fnames, count_dfs, metadata_dfs):
            self.add_study(id)
            self.add_study_counts(id, cdf)
            self.add_study_metadata(id, mdf)
    
    """treatment_filter
    filter the counts data by the treatments per the metadata file
    
    id_name <- string indicating the study id
    factor <- string indicating the column used for filtering
    valid_levels <- array of strings indicating valid values for factor column
    """
    def treatment_filter(self, id_name, factor, valid_levels, verbose=False):
        for k, v in self.count_dict.items():
            if verbose:
                print("Handling filtering of study {}".format(k))
            counts_df = v['counts']
            metadata_df = v['metadata']
            # get valid sample names 
            re = "|".join(valid_levels)
            valid_ids = metadata_df.loc[metadata_df[factor].str.contains(re), id_name]
            valid_ids = np.intersect1d(counts_df.index.values, valid_ids.tolist())
            self.count_dict[k]['counts'] = counts_df.loc[valid_ids, :]

    ### DATA SET PROFILING ###

    """profile_data
    generate summary statistics for counts data
    """
    def profile_data(self):
        profile_dict = {}
        kwargs = dict(axis=0)
        for id in self.count_dict.keys():
            study_counts = self.count_dict[id]['counts']
            means = study_counts.mean(**kwargs)
            stdevs = study_counts.std(**kwargs)
            vars = study_counts.var(**kwargs)
            kurts = study_counts.kurt(**kwargs)
            profile_dict[id] = {'mean':means, 
                                'std':stdevs, 
                                'var':vars,
                                'kurt':kurts}
        self.profile_dict = profile_dict
        
    """plot_profile
    metric <- string key for profile dict
    generate_csv <- boolean that saves output
    
    generates plots from data profiling metrics
    """
    def plot_profile(self, metric, generate_csv=False):
        # statistics generated from profile_data function call
        stats = [[k, v[metric]] for k, v in self.profile_dict.items()]
        # optional save statistics to csv
        if generate_csv:
            pd.concat([v[metric] for v in self.profile_dict.values()], axis=1).dropna(axis=0).to_csv('./data/stats/{}.csv'.format(metric))

        # define histogram parameters
        kwargs = dict(alpha=0.33, bins=250)
        
        # generate histogram of output        
        for i in stats:
            x = i[1].to_numpy()
            plt.hist(x[np.isfinite(x)], **kwargs, label=i[0]) # kurts
        
        # set plotting params and save histogram output            
        plt.title('Frequency Histogram of {} at Gene Level'.format(metric)); plt.ylabel('Frequency'); plt.legend()
        plt.savefig('{}.png'.format(metric))
    
    ### FEATURE SELECTION ### 
    """compute_mrmr
    X <- dataframe counts matrix 
    y <- series label vector
    K <- integer hyperparameter with number of iterations
    
    executes the minimum redundancy maximum relevance framework
    """
    def compute_mrmr(self, X, y, K):
        results = mrmr_classif(X, y, K)
        return results
    
    ### MODEL FITTING ###
    
    """update_roc
    trues <- array of labels 
    scores <- array of scores
    model_name <- name for fitted model
    
    saves the receiver operator characteristic curve fit results into a dictionary
    """
    def update_roc(self, trues, scores, model_name):
        fpr, tpr, thresholds = roc_curve(trues, scores)
        auc = roc_auc_score(trues, scores)
        self.roc_dict[model_name] = {}
        self.roc_dict[model_name]['fpr'] = fpr
        self.roc_dict[model_name]['tpr'] = tpr
        self.roc_dict[model_name]['thresholds'] = thresholds
        self.roc_dict[model_name]['auc'] = auc
        return self.roc_dict[model_name]
    
    """fit_model
    models <- list of strings with models selected
    mrmr <- boolean flag that indicates mrmr feature subset
    roc <- boolean flag that indicates roc curve fitting
    pfi <- boolean flag that indicates permutation feature importance calculation
    test_set <- dictionary with test set samples
    accuracy_block <- boolean flag that indicates accuracy testing calculation
    verbose <- boolean flag for function monitoring
    
    fits statistical models to data
    """
    def fit_model(self, models, mrmr=False, roc=False, pfi=False, test_set=None, accuracy_block=False, verbose=False):
        model_params = {
            "rf": {
                    "n_estimators": 100, 
                    "max_depth": 5, 
                    "random_state" : 0, 
                    "min_samples_leaf": 1,
                    "min_samples_split":2,
                    "oob_score":True,
                    "class_weight":"balanced"
                },
            "svm": {
                    "penalty":'l2',
                    "loss":'squared_hinge',
                    "C":1.0,
                    "class_weight":'balanced',
                    "random_state":12345
            },
            "lda": {
                "solver":"eigen",
                "shrinkage":"auto"
            },
            "lda-svd" : {
                "solver":"svd"
                },
            "LDA_lsqr_params" :  {"solver":"lsqr", "covariance_estimator":GraphicalLassoCV(assume_centered=True, n_jobs=-1)},
            "glm" : {
                "random_state":0
            }
        }
        
        model_list = []
        model_names = []
        if "rf" in models:
            model_list.append(RandomForestClassifier(**model_params["rf"]))
            model_names.append('rf')
        if "svm" in models: 
            model_list.append(LinearSVC(**model_params["svm"]))
            model_names.append('svm')
        if "lda" in models:
            model_list.append(LDA(**model_params["lda"]))
            model_names.append('lda')
        if "glm" in models:
            model_list.append(LogisticRegression(**model_params["glm"]))
            model_names.append('glm')

        
        # perform stratified kfold splitting of the data sets
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=12345)
        
        X = self.concat_df.iloc[:,:-1]; y = self.concat_df.iloc[:,-1]
        test_ids = []
        if test_set:
            for k in test_set.keys():
                test_ids.append(test_set[k]['id']['SF'])
                test_ids.append(test_set[k]['id']['GC'])
            test_ids = list(itertools.chain(*test_ids))
            valid_ids = X.index.isin(test_ids)
            train_ids = np.arange(X.index.values.shape[0])[~valid_ids]
            test_ids = np.arange(X.index.values.shape[0])[valid_ids]            
        
        if mrmr:
            X = X.loc[:,self.Genes['Gene']]
        
        dict_frame = {"accuracy":[], "model_coefs":[], "importance":[]}
        results = {m : dict_frame.copy() for m in models}
        acc_list = []; model_coef = []; roc_y_scores = []; roc_y_trues = []; 
        imp_arr = []
        for (model, model_name) in zip(model_list, model_names):
            train_test = skf.split(X, y)
            if test_set:
                train_test = zip([train_ids], [test_ids])
            if accuracy_block:
                train_test=skf.split(X.iloc[train_ids,:], y.iloc[train_ids])
                
            for train, test in train_test:
                model = model.fit(X.iloc[train,:], y.iloc[train])
                accuracy = accuracy_score(model.predict(X.iloc[test,:]), y.iloc[test])
                acc_list.append(accuracy)
                ### MODEL SPECIFIC ACTIONS ###
                if "rf" == model_name:
                    model_coef.append(model.feature_importances_)
                else:
                    model_coef.append(model.coef_.T)
                
                if roc:
                    if "rf" == model_name:
                        y_scores = model.predict_proba(X.iloc[test])[:, 1]
                    elif "lda" == model_name:
                        y_scores = model.predict_proba(X.iloc[test])[:, 1]
                    else:
                        y_scores = model.decision_function(X.iloc[test])
                    roc_y_scores.append(y_scores); roc_y_trues.append(y.iloc[test])
                            
                if pfi:
                    importances = self.permutation_scoring(model, X.iloc[test], y.iloc[test], sort=False)
                    imp_arr.append(importances)
                
                ### MODEL SPECIFIC VERBOSITY ###
                if verbose:
                    print("The average accuracy was {:.3f}".format(np.mean(acc_list)))
                    if "rf" == model_name:
                        print('train -  {}   |   test -  {}   |   accuracy - {:.3f}'.format(np.bincount(y.iloc[train]+1), np.bincount(y.iloc[test]+1), accuracy))
                        print("OOB Score is {}".format(model.oob_score_))
                        print("Number of features seen is {}".format(model.n_features_in_))
                    if "lda" == model_name:
                        print("The number of features seen during model fitting: {}".format(model.n_features_in_))
                        print("The shape of the coefficient matrix from the model: {}".format(model.coef_.shape))

            if roc: 
                roc_dict = self.update_roc(np.concatenate(roc_y_trues, axis=0), np.concatenate(roc_y_scores, axis=0), model_name)
                results[model_name]['roc'] = pd.DataFrame(roc_dict)
                results[model_name]['predictions'] = pd.concat([pd.DataFrame(np.concatenate(roc_y_trues, axis=0), columns=['truth']), pd.DataFrame(np.concatenate(roc_y_scores, axis=0), columns=['probability'])], axis=1)
                roc_y_scores.clear(); roc_y_trues.clear()
                
            if pfi:
                results[model_name]['importance'] = pd.concat(imp_arr, axis=0, ignore_index=True); imp_arr.clear()
                
            results[model_name]["accuracy"] = pd.DataFrame(acc_list, columns=["Accuracy"]); acc_list.clear()
            results[model_name]["model_coefs"] = pd.DataFrame(np.mean(model_coef, axis=0), columns=['Score'], index=X.columns.values); model_coef.clear()
        
        return results
    
    """permutation_scoring
    estimator <- fitted model object
    X <- dataframe matrix of counts data
    Y <- Series vector of predicted labels
    sort <- boolean indicator to sort values by importance
    seed <- integer to constrain randomness
    
    helper function to run permutation feature importance scoring 
    """
    def permutation_scoring(self, estimator, X, Y, sort=True, seed=12345):
        # run the permutation feature importance
        if python_version()[:3] == "3.9":
            result = permutation_importance(estimator, X.to_numpy(), Y.to_numpy(), n_repeats=100, random_state=seed, scoring='roc_auc', n_jobs=-1)#, max_samples=10)
        else:
            result = permutation_importance(estimator, X, Y, n_repeats=100, random_state=seed, n_jobs=-1)
        # run a sort on the features
        sorted_importances_idx = result.importances_mean.argsort()
        # save to table
        # print(result.importances_mean[sorted_importances_idx])
        if sort:
            importances = pd.DataFrame(
                result.importances[sorted_importances_idx].T,
                columns=X.columns[sorted_importances_idx]
            )
        else:
            importances = pd.DataFrame(
                result.importances.T,
                columns=X.columns
            )
        return importances
                
    """generate_directionality
    uses the direction of average counts to create a directionality vector for features
    """
    def generate_directionality(self):
        data = self.concat_df.iloc[:,:-1]
        labels = self.concat_df.iloc[:,-1]
        
        sf_mask = labels == 1
        gc_mask = labels == -1
        
        sf = data.loc[sf_mask,:]
        gc = data.loc[gc_mask,:]
        
        print(sf.mean(axis=0))
        print(gc.mean(axis=0))
        print(np.greater(sf.mean(axis=0), gc.mean(axis=0)))
        
        dir_mask = np.greater(gc.mean(axis=0), sf.mean(axis=0))
        
        dir_vector = np.ones(dir_mask.shape[0])
        dir_vector[dir_mask] = -1
        
        print(dir_vector)
        
        dir_df = pd.DataFrame(dir_vector, columns=["Direction"], index=data.columns)
        dir_df.to_csv("directional.csv", header=True)
    
    """plot_pca
    mrmr <- boolean flag that indicates if the mrmr feature subset should be used
    
    plots a PCA of the data based on classification target
    """
    def plot_pca(self, mrmr=False):
        pca = PCA(n_components=2)
        
        X = self.concat_df.iloc[:,:-1]
        y = self.concat_df.iloc[:,-1]
        target_names = ['Ground Control', 'Space Flight']
        target_names = ['GC', 'SF']

        labels = np.empty(X.shape[0], dtype=object)
        for k in self.count_dict.keys():
            study_mask = np.isin(X.index.values, self.count_dict[k]['metadata']['Sample Name'])
            labels[study_mask] = k
        
        if mrmr:
            X = X.loc[:,self.Genes['Gene']]
        
        X_r = pca.fit(X).transform(X)
        colors = ["navy", "darkorange"]

        figure_generator.plot_pca(pca, X_r, y, target_names, colors, labels, mrmr)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments for MRMR run')

    parser.add_argument("-K", "--K", help="Specifies the number of MRMR features to select.", default=30, type=int)
    parser.add_argument("-O", "--outfile_name", help="Name for output MRMR csv.", default="mrmr")
    parser.add_argument("-F", "--feature_list", help="csv with gene set filtering for concatenated matrix.", default=None)
    parser.add_argument("-C", "--custom_concat", help="Provide custom concatenated data file.", default=None)
    parser.add_argument("-G", "--gene_list", help="csv with gene lists for model fitting.", default=None)
    parser.add_argument("-B", "--block_flags", help="json with flag specifications", default=None)
    parser.add_argument("-V", "--verbose", help="toggle output from functions", action="store_true")
    
    args = parser.parse_args()
        
    kwargs = { 
              "data_dir":['data', 'norm_counts'],
              "feature_list": args.feature_list 
              }
    
    run_flags = {}
    if args.block_flags is not None:
        with open(args.block_flags, 'r') as f:
            json_str = f.read()
            run_flags = json.loads(json_str)    
    
    data_process_block = run_flags.get("data_process_block", True)
    save_concat_df_block = run_flags.get("save_concat_df_block", False)
    load_concat_df_to_clipboard_block = run_flags.get("load_concat_df_to_clipboard_block", False)
    profile_data_block = run_flags.get("profile_data_block", False)
    generate_pickle_block = run_flags.get("generate_pickle_block", False)
    mrmr_fitting_block = run_flags.get("mrmr_fitting_block", False)
    rf_fitting_block = run_flags.get("rf_fitting_block", False)
    svm_fitting_block = run_flags.get("svm_fitting_block", False)
    lda_fitting_block = run_flags.get("lda_fitting_block", False)
    accuracy_block = run_flags.get("accuracy_block", False)
    pca_block = run_flags.get("pca_block", False)
    roc_block = run_flags.get("roc_block", True)
    pfi_block = run_flags.get("pfi_block", False)
    random_block = run_flags.get("random_block", False)
    test_block = run_flags.get("test_block", False)

    gene_list_flag = args.gene_list is not None
    
    if test_block:
        test_labels = {
        '47': {
            'id': {
                'SF':["Mmus_C57-6T_LVR_FLT_Rep1_F1"], 
                'GC':["Mmus_C57-6T_LVR_GC_Rep3_G5"]
                }
            },
        '168_rr1': {
            'id':{
                'SF':["Mmus_C57-6J_LVR_RR1_FLT_wERCC_Rep1_M25"],
                'GC':["Mmus_C57-6J_LVR_RR1_GC_wERCC_Rep5_M40"]
                }
            },
        '168_rr3': {
            'id': {
                    'SF':["Mmus_BAL-TAL_LVR_RR3_FLT_wERCC_Rep1_F1"], 
                    'GC':["Mmus_BAL-TAL_LVR_RR3_GC_wERCC_Rep4_G5"]
                    }
                },
        '242': {
            'id': {
                'SF':["Mmus_C57-6J_LVR_FLT_C1_Rep1_F1"], 
                'GC':["Mmus_C57-6J_LVR_GC_C2_Rep1_G1","Mmus_C57-6J_LVR_CC_C1_Rep1_C1-1","Mmus_C57-6J_LVR_CC_C2_Rep1_C2-1"]
                }
            },
        '245': {
            'id': {
                'SF':["Mmus_C57-6T_LVR_FLT_LAR_Rep1_F1","Mmus_C57-6T_LVR_FLT_ISS-T_Rep4_F10","Mmus_C57-6T_LVR_FLT_ISS-T_Rep2_F8"],
                'GC':["Mmus_C57-6T_LVR_GC_ISS-T_Rep1_G4","Mmus_C57-6T_LVR_GC_ISS-T_Rep2_G9","Mmus_C57-6T_LVR_GC_LAR_Rep1_G6","Mmus_C57-6T_LVR_GC_LAR_Rep2_G4"]
                }
            },
        '379': {
            'id': {
                'SF': ["RR8_LVR_FLT_ISS-T_YNG_FI7","RR8_LVR_FLT_ISS-T_YNG_FI8","RR8_LVR_FLT_ISS-T_YNG_FI9"
                        ,"RR8_LVR_FLT_ISS-T_OLD_FI11","RR8_LVR_FLT_ISS-T_OLD_FI12","RR8_LVR_FLT_ISS-T_OLD_FI13"],
                'GC':["RR8_LVR_GC_ISS-T_YNG_GI8","RR8_LVR_GC_ISS-T_YNG_GI9","RR8_LVR_GC_ISS-T_YNG_GI10"
                      ,"RR8_LVR_GC_ISS-T_OLD_GI11","RR8_LVR_GC_ISS-T_OLD_GI12","RR8_LVR_GC_ISS-T_OLD_GI13"]
                }
            },
    }
    else:
        test_labels = None
    
    model = RNASeqModel()

    # this block creates the combined dataframe 
    if data_process_block:
        model.load_counts(verbose=args.verbose, scale=True, exclude_test=test_labels,**kwargs)    
        model.treatment_filter(id_name = "Sample Name", factor='Factor Value[Spaceflight]', valid_levels=['Space Flight', 'Ground Control', 'Cohort Control #1', 'Cohort Control #2'], verbose=args.verbose)
        model.build_combined_df(args.custom_concat)
    
    print(model.concat_df.shape)
    
    if random_block:
        seed(12345)
        ids = sample(list(np.arange(model.concat_df.shape[1]-1)),60)
        ids.append(-1)
        print(f"Random gene IDs {ids[:5]}")
        model.concat_df = model.concat_df.iloc[:,ids]

    if args.gene_list is not None:
        gene_list = pd.read_csv(args.gene_list, index_col=0)#.iloc[:100,]
        model.Genes = gene_list
    
    if pca_block:
        model.plot_pca(mrmr=args.gene_list)
        

    # this block is used to export the concatenated data set
    if save_concat_df_block:
        pd.DataFrame(model.concat_df).transpose().to_csv("concat_df.csv", header=True)
        model.generate_directionality()
    # this block saves the concatenated dataframe to the clipboard for troubleshooting    
    if load_concat_df_to_clipboard_block:
        model.concat_df.transpose().to_clipboard()
        
    # this block profiles the mean and std of the data 
    if profile_data_block:
        model.profile_data()
        model.plot_profile('mean')
    
    # # this block computes MRMR and exports 
    if mrmr_fitting_block:
        mrmr = model.compute_mrmr(model.concat_df.iloc[:,:-1], model.concat_df.iloc[:,-1], K=args.K)
        mrmr_df = pd.DataFrame(mrmr, columns=["Gene"])
        mrmr_df.to_csv("{}.csv".format(args.outfile_name), header=True)
    
    model_names = ['rf', 'svm', 'lda']
    chosen_models = [x for x, y in zip(model_names, [rf_fitting_block, svm_fitting_block, lda_fitting_block]) if y]
    print("Models chosen are: ", chosen_models)    
    # this block loads a gene list and computes performance
    if args.gene_list is not None:
        gene_list = pd.read_csv(args.gene_list, index_col=0)#.iloc[:100,]
        model.Genes = gene_list
    

    results = model.fit_model(chosen_models, mrmr=gene_list_flag, roc=roc_block, pfi=pfi_block, test_set=test_labels)

    if accuracy_block and gene_list_flag:
        dim = gene_list.shape[0]
        print(f"dim is {dim}")
        acc_array = np.zeros([dim, 3])
        for i in range(2,dim+1):
            model.Genes = gene_list.iloc[:i,]
            print("Gene {} added: ID {}".format(i, model.Genes.iloc[i-1,0]))
            results = model.fit_model(model_names, mrmr=gene_list_flag, roc=False, pfi=False, test_set=test_labels, accuracy_block=accuracy_block)
            for k, m in enumerate(model_names):
                acc_array[i-1,k] = results[m]['accuracy'].mean()
        col_names = model_names
        acc_cols = [col_names[x] for x,i in enumerate(model_names) if i]
        pd.DataFrame(acc_array, columns=col_names).to_csv("{}.csv".format(args.outfile_name))
        
    print("Done.")
    
