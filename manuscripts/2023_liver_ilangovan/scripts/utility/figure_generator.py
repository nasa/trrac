from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mtick
import numpy as np
import seaborn as sns
from upsetplot import plot, generate_counts

mpl.rcParams['font.size'] = 18
mpl.rcParams['legend.title_fontsize'] = 16
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16


def roc_plot(roc_dict, title="ROC"):
    fig, ax = plt.subplots(figsize=(6,4), dpi=300)
    lw=1

    full_name = {'rf':"RF", 'svm':"SVM", 'lda':'LDA', 'glm':'GLM'}
    marker_map = {'rf':"o", 'svm':"X", 'lda':'s', 'glm':'v'}
    for k,v in roc_dict.items():
        ax.plot(v['fpr'], v['tpr'], lw=lw, label=f"{full_name[k]} (AUC = {v['auc'].mean():.2f})", marker=marker_map[k])
    # ax.plot(df['fpr'], df['tpr'], color='yellow', lw=1, label=f"ROC Curve (area = {df['auc'].mean():.2f})")
    ax.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend(loc="lower right")

def accuracy_plot(df, add_label=True,elbow=None, right=None):
    if elbow is None:
        elbow=df['num_mrmr_features'].max()/2
    if right is None:
        right=df['num_mrmr_features'].max()
    fig, ax = plt.subplots(figsize=(6,4), dpi=300)
    full_name = {'rf':"RF", 'svm':"SVM", 'lda':'LDA'}
    df = df.replace({"model":full_name})
#     df['accuracy'] = df['accuracy']*100
    sns.lineplot(data=df, x="num_mrmr_features", y="accuracy", hue="model", ax=ax, markers=True, style="model",markevery=5)
    min_acc = df[['num_mrmr_features', 'accuracy']].groupby(by='num_mrmr_features').min()
    
    ax.set_title("Reduced Gene Set Classifier Accuracy")
    ax.set_xlabel("MRMR feature count")
    ax.set_ylabel("Mean Accuracy (%)")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels)
    
    ax.set_xticks(np.arange(0,df['num_mrmr_features'].max(), step=20))
    top=0.95; bottom=0.7
    ax.set_xlim(left=1, right=right); ax.set_ylim(bottom=0.7, top=0.92); 
    print(min_acc.loc[elbow,'accuracy'])
    top=0.92; bottom=0.7
    ax.axvline(x=elbow, ymax=(min_acc.loc[elbow,'accuracy']-bottom)/(top-bottom), lw=1, color=sns.color_palette("deep")[0], linestyle='dashed')    
    
    
    annotations = (df['accuracy'].round(2)*100).astype(int).astype(str)

    plt.annotate("Elbow Point",  
        (elbow, 
        min_acc.loc[elbow,'accuracy']*(0.94)), 
        c=sns.color_palette("deep")[0])

    for i, label in enumerate(annotations):
        if i % 10 == 0 and i >df['num_mrmr_features'].max() and add_label:
            plt.annotate(label,  
                         (df.iloc[i,0], 
                          df.iloc[i,2]*(0.97)), 
                         c=sns.color_palette("deep")[3])
    def format_fn(tick_val, tick_pos):
        return str(int(tick_val*100))
    ax.yaxis.set_major_formatter(format_fn)

def generate_ranking_from_permutation_importance(df):
    df[df < 0] = 0
    df_ranking = df.mean(axis=0).rank(method="max", ascending=False).astype(int).values
    df_ordered_by_ranking = df.columns[df_ranking-1]
    return (df_ranking, df_ordered_by_ranking)

def upset_plot():
    print("Generating Upset")
    
"""generate a pca plot using supplied labels and features
pca <- principal component analysis object generated from PCA function
X_r <- transformed data matrix using fitted PCA as data frame
y <- response variable as a series
target_names <- labels associated with the target response as a list of strings 
colors <- the colors for the target_names given as a list of strings
labels <- the labels used to encode each sample as a series 
mrmr <- boolean indicator if the provided features were from mrmr
"""
def plot_pca(pca, X_r, y, target_names, colors, labels, mrmr):
    # generate a pyplot object that will hold two plots, one for the target encoding and the other encoding based on the provided labels
    fig, ax = plt.subplots(figsize=(6,8), dpi=300, nrows=2, ncols=1, sharex=True, sharey=True)
    lw = 2
    
    # extract the variance explained from the pca object for axis labeling
    var_explained = pca.explained_variance_ratio_        
#         print(X.index.values[X_r[:, 0] > 29040000])
    # generate a mapping for the markers to use in the plot based on the target names provided
    marker_map = {key:mpl.markers.MarkerStyle.filled_markers[i] for i,key in enumerate(target_names)}
    # use the provided colors to  generate the color encoding for the scatter plot
    for color, i, target_name in zip(colors, [-1, 1], target_names):
        ax[0].scatter(
            X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=0.8, lw=lw, label=target_name, marker=marker_map[target_name]
        )

    ax[0].legend(bbox_to_anchor=(1.05, 1),
                         loc='upper left', borderaxespad=0.)
    title = "PCA of Concatenated Study Matrix"
    if mrmr:
        title += " (MRMR)"
    else:
        title += " (All Genes)"
    ax[0].set_title(title, pad=20)
    # set PC1 axis label
#     ax[0].set_xlabel(f"PC1 ({np.round(var_explained[0], 2)*100}%)")
    # set PC2 axis label
    
#         plt.savefig('PCA.png')

    # here we set up a generalizable section that can be used to provide alternative labels for the pca plot
    # this is the mission/lab based mapping
    mission_map = {'47':'RR1(C)','168_rr1':'RR1(N)','168_rr3':'RR3','242':'RR9','245':'RR6','379':'RR8'}
    order_map = {'RR1(C)':0,'RR1(N)':1,'RR3':2,'RR9':5,'RR6':3,'RR8':4}
    
    # this is the strain based mapping
    mission_map = {'47':'C57BL/6NTac','168_rr1':'C57BL/6J','168_rr3':'BALB/c','242':'C57BL/6J','245':'C57BL/6NTac','379':'BALB/cAnNTac'}
    order_map = {'C57BL/6J':0,'C57BL/6NTac':1,'BALB/c':2,'BALB/cAnNTac':3}

    # this is the age based mapping
    mission_map = {'47':'32','168_rr1':'16','168_rr3':'12','242':'10','245':'32','379':'10-32'}
    order_map = {'10':0,'12':1,'16':2,'10-32':3, '32': 4}
    
    # this is the sex based mapping
    mission_map = {'47':'F','168_rr1':'F','168_rr3':'F','242':'M','245':'F','379':'F'}
    order_map = {'M':0,'F':1}

    mission_labels = np.array([mission_map[i] for i in labels])    
    rgb_values = sns.color_palette("Set2", len(set(mission_labels)))
    color_map = dict(zip(set(mission_labels), rgb_values))
    marker_map = dict(zip(set(mission_labels), ['o', 'v', '^', '<', '>', 's'][:len(set(mission_labels))]))
    

#     marker_map = dict(zip(set(labels), ['o', 'v', '^', 'o', '^', 'v'][:len(set(labels))]))
#     order_map = {'C57BL/6NTac':0,'C57BL/6J':1,'BALB/c':2,'C57BL/6J':5,'C57BL/6NTac':3,'BALB/cAnNTac':4}
#     rgb_values = sns.color_palette("Set2", len(set(order_map.keys())))
#     color_map = dict(zip(set(order_map.keys()), rgb_values))
    for study_id in set(mission_labels):
        ax[1].scatter(
            X_r[mission_labels == study_id, 0], X_r[mission_labels == study_id, 1], alpha=0.8, lw=lw, color=color_map[study_id], label=study_id, marker=marker_map[study_id]
        )
#     plt.legend(loc="best", shadow=False, scatterpoints=1)
#     ax[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=12, mode="expand", borderaxespad=0.)
    # reordering labels
    legend_handles, legend_labels = ax[1].get_legend_handles_labels()
    order = [order_map[o] for o in legend_labels]
    
    
    print([legend_labels[i] for i in order])
    ax[1].legend([legend_handles[i] for i in np.argsort(order)], [legend_labels[i] for i in np.argsort(order)], bbox_to_anchor=(1.05, 1),
                         loc='upper left', borderaxespad=0.)
    title = "PCA of Concatenated Study Matrix"
    if mrmr:
        title += " (MRMR)"
    else:
        title += " (All Genes)"
#     ax[1].set_title(title)
    # set PC1 axis label
    ax[1].set_xlabel(f"PC1 ({np.round(var_explained[0], 2)*100}%)")
    # set PC2 axis label
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.ylabel(f"PC2 ({np.round(var_explained[1], 2)*100}%)       ")
