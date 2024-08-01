import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import numpy as np
import matplotlib
from scipy import stats
import sys
from pathlib import Path
from scipy.stats import ranksums
from sklearn.metrics.pairwise import cosine_similarity

def get_stats(df):
    if len(df.columns) == 2:
        #spearman rank correlation
        spcor = stats.spearmanr(df, nan_policy='raise')
        #cosine similarity
        cosi = cosine_similarity(df.T)
        cosidf = pd.DataFrame(cosi)
        cosidf.columns = df.columns
        cosidf.index = df.columns
        cosidf1 = cosidf.where(np.triu(np.ones(cosidf.shape)).astype(np.bool_))
        np.fill_diagonal(cosidf1.values, np.nan) 
        cosidf2 = cosidf1.stack().reset_index()
        cosidf2.columns = ["cell_pair1", "cell_pair2", "cosine_similarity"]
        cosidf2['spearman_correlation'] = spcor.statistic
    else:
        #spearman rank correlation
        spcor = stats.spearmanr(df, nan_policy='raise')
        spcor1 = pd.DataFrame(spcor.statistic)
        spcor1.columns = df.columns
        spcor1.index = df.columns
        spcor2 = spcor1.where(np.triu(np.ones(spcor1.shape)).astype(np.bool_))
        np.fill_diagonal(spcor2.values, np.nan)
        spcor3 = spcor2.stack().reset_index()
        spcor3.columns = ["cell_pair1", "cell_pair2", "spearman_correlation"]
        #cosine similarity
        cosi = cosine_similarity(df.T)
        cosidf = pd.DataFrame(cosi)
        cosidf.columns = df.columns
        cosidf.index = df.columns
        cosidf1 = cosidf.where(np.triu(np.ones(cosidf.shape)).astype(np.bool_))
        np.fill_diagonal(cosidf1.values, np.nan) 
        cosidf2 = cosidf1.stack().reset_index()
        cosidf2.columns = ["cell_pair1", "cell_pair2", "cosine_similarity"]
        cosidf2['spearman_correlation'] = spcor3["spearman_correlation"]
    return(cosidf2)


PBMC_PATH=Path("PBMC_ifnb_data")
PROCESS_PATH = PBMC_PATH / "TRIAGE_Cluster"
PLOT_PATH = PBMC_PATH / "corr_cos"



bwall = []
for i in range(11):
    bwall.append(round(0.1+i/100, 3))
j=0.3
while j <=5:
    bwall.append(round(j,3))
    j+=0.1

#integrated data
intdf = pd.read_csv(PBMC_PATH / 'pbmc_ifnb_Seurat_CCAint_data_HVG.csv', index_col=0)
bw_list = list()
for bw in bwall:
    pme = pd.read_csv(PROCESS_PATH / str("TRIAGE_Cluster_bw{:.2f}_metadata.csv".format(bw)), index_col=0)
    peakls = pme["Peak"].unique()
    peak_list = list()
    for peak in peakls:
        intdf1 = intdf[pme[pme["Peak"] == peak]["cell_name"]]
        if len(intdf1.columns) > 1:
            dfstat = get_stats(intdf1)
            dfstat1=[dfstat.loc[:, 'cosine_similarity'].mean(), dfstat.loc[:, 'spearman_correlation'].mean(),'peak_' + str(peak)]
            peak_list.append(dfstat1)
    peakstat = pd.DataFrame(peak_list, columns=["cosine_similarity", "spearman_correlation","peak"])
    peakstat['bandwidth'] = 'bw' + str(bw)
    peakstat['clust_num'] = peakstat["peak"].nunique()
    peakstat["method"] = "TRIAGE-Cluster"
    bw_list.append(peakstat)
bwstat = pd.concat(bw_list, ignore_index=True) 
bwstat.to_csv(PBMC_PATH / 'pbmc_ifnb_integrated_TRIAGECluster_bw01_5_peak_mean_correlation_similarity.csv')

seupca = pd.read_csv(PBMC_PATH /'pbmc_ifnb_integrated_Seurat_within_cluster_mean_correlation_similarity.csv', index_col=0)
peakpca = pd.read_csv(PBMC_PATH / 'pbmc_ifnb_integrated_TRIAGECluster_bw01_5_peak_mean_correlation_similarity.csv', index_col=0)

#the range of total clusters from low to high is different for Seurat and TRIAGE-Cluster
#end the plot at same number of clusters on x-axis
max_cls = seupca["clust_num"].max()
peakpca_3 = peakpca.query("clust_num <= @max_cls")
pcacorall = pd.concat([peakpca_3, seupca], ignore_index=True) 

#plot cosine similarity
t1=pcacorall.loc[pcacorall['method'] == "TRIAGE-Cluster"]["cosine_similarity"].values
t2=pcacorall.loc[pcacorall['method'] == "Seurat"]["cosine_similarity"].values
t3=ranksums(t1, t2)
fig, ax = plt.subplots(figsize=(12,12))
sns.lineplot(data= pcacorall, x="clust_num", y="cosine_similarity", hue="method", ax=ax)
ax.set(ylim=(0,1))
ax.text(35, 0.05, f'pvalue {str(t3[1])}')
plt.savefig(PLOT_PATH / "pbmc_ifnb_integrated_data_all_cosine_similarity.pdf", dpi=300)

#plot spearman correlation
t1=pcacorall.loc[pcacorall['method'] == "TRIAGE-Cluster"]["spearman_correlation"].values
t2=pcacorall.loc[pcacorall['method'] == "Seurat"]["spearman_correlation"].values
t3=ranksums(t1, t2)
fig, ax = plt.subplots(figsize=(12,12))
sns.lineplot(data= pcacorall, x="clust_num", y="spearman_correlation", hue="method", ax=ax)
ax.set(ylim=(0,1))
ax.text(35, 0.05, f'pvalue {str(t3[1])}')
plt.savefig(PLOT_PATH / "pbmc_ifnb_integrated_data_all_spearman_correlation.pdf", dpi=300)


#PCA data, same as above integrated data, only change the input for "intdf"
intdf = pd.read_csv(PBMC_PATH / 'pbmc_ifnb_Seurat_CCAint_pca.csv', index_col=0)
intdf = intdf.T
bw_list = list()
for bw in bwall:
    pme = pd.read_csv(PROCESS_PATH / str("TRIAGE_Cluster_bw{:.2f}_metadata.csv".format(bw)), index_col=0)
    peakls = pme["Peak"].unique()
    peak_list = list()
    for peak in peakls:
        intdf1 = intdf[pme[pme["Peak"] == peak]["cell_name"]]
        if len(intdf1.columns) > 1:
            dfstat = get_stats(intdf1)
            dfstat1=[dfstat.loc[:, 'cosine_similarity'].mean(), dfstat.loc[:, 'spearman_correlation'].mean(),'peak_' + str(peak)]
            peak_list.append(dfstat1)
    peakstat = pd.DataFrame(peak_list, columns=["cosine_similarity", "spearman_correlation","peak"])
    peakstat['bandwidth'] = 'bw' + str(bw)
    peakstat['clust_num'] = peakstat["peak"].nunique()
    peakstat["method"] = "TRIAGE-Cluster"
    bw_list.append(peakstat)
bwstat = pd.concat(bw_list, ignore_index=True) 
bwstat.to_csv(PBMC_PATH / 'pbmc_ifnb_PCA_TRIAGECluster_bw01_5_peak_mean_correlation_similarity.csv')

seupca = pd.read_csv(PBMC_PATH /'pbmc_ifnb_PCA_Seurat_within_cluster_mean_correlation_similarity.csv', index_col=0)
peakpca = pd.read_csv(PBMC_PATH / 'pbmc_ifnb_PCA_TRIAGECluster_bw01_5_peak_mean_correlation_similarity.csv', index_col=0)
max_cls = seupca["clust_num"].max()
peakpca_3 = peakpca.query("clust_num <= @max_cls")
pcacorall = pd.concat([peakpca_3, seupca], ignore_index=True) 

#plot cosine similarity
t1=pcacorall.loc[pcacorall['method'] == "TRIAGE-Cluster"]["cosine_similarity"].values
t2=pcacorall.loc[pcacorall['method'] == "Seurat"]["cosine_similarity"].values
t3=ranksums(t1, t2)
fig, ax = plt.subplots(figsize=(12,12))
sns.lineplot(data= pcacorall, x="clust_num", y="cosine_similarity", hue="method", ax=ax)
ax.set(ylim=(0,1))
ax.text(35, 0.05, f'pvalue {str(t3[1])}')
plt.savefig(PLOT_PATH / "pbmc_ifnb_PCA_data_all_cosine_similarity.pdf", dpi=300)

#plot spearman correlation
t1=pcacorall.loc[pcacorall['method'] == "TRIAGE-Cluster"]["spearman_correlation"].values
t2=pcacorall.loc[pcacorall['method'] == "Seurat"]["spearman_correlation"].values
t3=ranksums(t1, t2)
fig, ax = plt.subplots(figsize=(12,12))
sns.lineplot(data= pcacorall, x="clust_num", y="spearman_correlation", hue="method", ax=ax)
ax.set(ylim=(0,1))
ax.text(35, 0.05, f'pvalue {str(t3[1])}')
plt.savefig(PLOT_PATH / "pbmc_ifnb_PCA_data_all_spearman_correlation.pdf", dpi=300)

