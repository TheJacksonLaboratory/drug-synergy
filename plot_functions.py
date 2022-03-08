import numpy as np
import pandas as pd
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.patches import Wedge
from sklearn.cluster import KMeans

from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import pdist


def plotDendrogramWithKnownPairs(Z, dfC, dfKS):

    fig, ax = plt.subplots(2, 1, figsize=(10, 5), gridspec_kw={'height_ratios':[1,3]})

    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.5
    n_clusters = 10
    cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1))
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    D = dendrogram(Z, ax=ax[0], color_threshold = (Z[-n_clusters,2] + Z[-n_clusters + 1,2]) / 2, above_threshold_color='k', orientation='top', no_labels=True)
    leaves = D['leaves']
    names = dfC.index[leaves]
    ax[0].axis("off")

    maxco = dfC.stack().max()

    d = []
    for drug1, drug2 in dfKS.values[:]:
        if (drug1 in dfC.index) and (drug2 in dfC.index):
            p1, p2 = np.where(names == drug1)[0][0] / len(dfC), np.where(names == drug2)[0][0] / len(dfC)
            spr = np.abs(p1 - p2) / 2
            co = dfC.loc[drug1, drug2]
            w = Wedge((min(p1, p2) + spr, 1), spr, 180, 0, width=0.001, color=cm.hot_r(0.25 + (co / (2. * maxco)))) # coolwarm
            ax[1].add_artist(w)

            d.append(spr)
    d = np.array(d)

    ax[1].set_ylim([0.5, 1.0])
    ax[1].axis("off")

    plt.subplots_adjust(hspace=0)

    hierarchy.set_link_color_palette(None)
    matplotlib.rcParams['lines.linewidth'] = origLineWidth

    return


def plotH2(dfC, Z, pvalues):

    fig, ax = plt.subplots(2, 2, figsize=(7, 7), gridspec_kw={'height_ratios':[1,3], 'width_ratios':[3,1]})

    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.5
    n_clusters = 10
    cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1))
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    dendrogram(Z, ax=ax[1, 1], color_threshold = (Z[-n_clusters,2] + Z[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='right', no_labels=True)
    D = dendrogram(Z, ax=ax[0, 0], color_threshold = (Z[-n_clusters,2] + Z[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='top', no_labels=True)

    leaves = D['leaves']
    names = dfC.index[leaves]
    ax[1, 1].axis("off")
    ax[0, 0].axis("off")

    ax[1, 0].axis("off")

    se = pd.Series(index=pd.MultiIndex.from_tuples(pvalues), data=1)
    se = pd.concat([se, se.reorder_levels([0, 2, 1])])
    se = se.loc[~se.index.duplicated()]
    se = se.droplevel(0)
    se = se.groupby(level=[0, 1]).sum()
    se = se.loc[se.index.to_frame()[0] != se.index.to_frame()[1]]
    dfDr = se.unstack(0).fillna(0.).reindex(dfC.index, axis=0).reindex(dfC.index, axis=1).fillna(0.).astype(int)

    dfDr = dfDr.replace(0, np.nan)
    dfDr = dfDr.iloc[leaves[::-1], leaves]
    masked_M = np.ma.array(dfDr.values, mask=np.isnan(dfDr.values))

    #dfCc = dfC.iloc[leaves[::-1], leaves]
    #masked_M = np.ma.array(dfCc.values, mask=np.isnan(dfCc.values))
    cmap = copy.copy(plt.cm.bwr)
    cmap.set_bad('lightgrey')
    vmin, vmax = None, None

    im = ax[1, 0].imshow(masked_M, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='None', 
                   extent=(-0.5, masked_M.shape[0] - 0.5, masked_M.shape[1] - 0.5, -0.5))

    ax[0, 1].axis("off")
    clb = fig.colorbar(im, ax=ax[0, 1], fraction=0.5, shrink=0.85, orientation='horizontal', label='Count')
    clb.ax.tick_params(labelsize=10)

    plt.subplots_adjust(hspace=0.001)
    fig.tight_layout()

    hierarchy.set_link_color_palette(None)
    matplotlib.rcParams['lines.linewidth'] = origLineWidth

    return


def plotH3(dfC, Z, pvalues):
    
    se = pd.Series(index=pd.MultiIndex.from_tuples(pvalues), data=1)
    se = pd.concat([se, se.reorder_levels([0, 2, 1])])
    se = se.loc[~se.index.duplicated()]
    df = se.groupby(level=[0,1]).sum().unstack(0).fillna(0.)
    
    #return df.T

    fig, ax = plt.subplots(2, 2, figsize=(7, 7), gridspec_kw={'height_ratios':[1,3], 'width_ratios':[3,1]})

    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.5
    n_clusters = 8
    cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1))
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    Z1 = linkage(pdist(df.values, metric='euclidean'), method='ward', optimal_ordering=True)
    D1 = dendrogram(Z1, ax=ax[0, 0], color_threshold = (Z1[-n_clusters,2] + Z1[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='top', no_labels=True)
    
    Z2 = linkage(pdist(df.values.T, metric='euclidean'), method='ward', optimal_ordering=True)
    D2 = dendrogram(Z2, ax=ax[1, 1], color_threshold = (Z2[-n_clusters,2] + Z2[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='right', no_labels=True)

    leaves1 = D1['leaves']
    leaves2 = D2['leaves']
    ax[1, 1].axis("off")
    ax[0, 0].axis("off")

    ax[1, 0].axis("off")

    print(df.shape, df.shape[0] ** 2)
    
    df = df.replace(0., np.nan)
    
    df = df.iloc[leaves1[::-1], leaves2]

    masked_M = np.ma.array(df.values.T, mask=np.isnan(df.values.T))
    cmap = copy.copy(plt.cm.bwr)
    cmap.set_bad('lightgrey')
    vmin, vmax = None, None

    im = ax[1, 0].imshow(masked_M, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='None', 
                   extent=(-0.5, masked_M.shape[0] - 0.5, masked_M.shape[1] - 0.5, -0.5))

    ax[0, 1].axis("off")
    clb = fig.colorbar(im, ax=ax[0, 1], fraction=0.5, shrink=0.85, orientation='horizontal', label='Count')
    clb.ax.tick_params(labelsize=10)

    plt.subplots_adjust(hspace=0.001)
    fig.tight_layout()

    hierarchy.set_link_color_palette(None)
    matplotlib.rcParams['lines.linewidth'] = origLineWidth

    return


def plotH4(dfC, Z, pvalues, seed=0):
    
    se = pd.Series(index=pd.MultiIndex.from_tuples(pvalues), data=1)
    se = pd.concat([se, se.reorder_levels([0, 2, 1])])
    se = se.loc[~se.index.duplicated()]
    #df = se.groupby(level=[0,1]).sum().unstack(0).fillna(0.)
    df = se.unstack(0).fillna(0.)
    print(df.shape)
    
    cl = None
    if df.shape[0] > 10 ** 3:
        np.random.seed(seed)
        cellClusterIndex = KMeans(n_clusters=200).fit(df.values).labels_.astype(int)
        cl = pd.Series(data=df.index.values, index=cellClusterIndex)
        df.index = cellClusterIndex
        df = df.groupby(level=0).mean()
        print(df.shape)
    
    fig, ax = plt.subplots(2, 2, figsize=(7, 7), gridspec_kw={'height_ratios':[1,3], 'width_ratios':[3,1]})

    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.5
    n_clusters = 8
    cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1))
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    Z1 = linkage(pdist(df.values, metric='euclidean'), method='ward', optimal_ordering=True)
    D1 = dendrogram(Z1, ax=ax[0, 0], color_threshold = (Z1[-n_clusters,2] + Z1[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='top', no_labels=True)
    
    Z2 = linkage(pdist(df.values.T, metric='euclidean'), method='ward', optimal_ordering=True)
    D2 = dendrogram(Z2, ax=ax[1, 1], color_threshold = (Z2[-n_clusters,2] + Z2[-n_clusters + 1,2]) / 2, 
                   above_threshold_color='k', orientation='right', no_labels=True)

    leaves1 = D1['leaves']
    leaves2 = D2['leaves']
    ax[1, 1].axis("off")
    ax[0, 0].axis("off")

    ax[1, 0].axis("off")

    df = df.replace(0., np.nan)
    
    df = df.iloc[leaves1[::-1], leaves2]

    masked_M = np.ma.array(df.values.T, mask=np.isnan(df.values.T))
    cmap = copy.copy(plt.cm.bwr)
    cmap.set_bad('lightgrey')
    vmin, vmax = None, None

    im = ax[1, 0].imshow(masked_M, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='None', 
                   extent=(-0.5, masked_M.shape[0] - 0.5, masked_M.shape[1] - 0.5, -0.5))

    ax[0, 1].axis("off")
    clb = fig.colorbar(im, ax=ax[0, 1], fraction=0.5, shrink=0.85, orientation='horizontal', label='Count')
    clb.ax.tick_params(labelsize=10)

    plt.subplots_adjust(hspace=0.001)
    fig.tight_layout()

    hierarchy.set_link_color_palette(None)
    matplotlib.rcParams['lines.linewidth'] = origLineWidth

    return cl

