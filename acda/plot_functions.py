import os
import copy
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

    fig, ax = plt.subplots(2, 1, figsize=(15, 7), gridspec_kw={'height_ratios':[1,3]})

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

    return fig

def plotHeatmapPredictedSynergy(dfC, Z, pvalues):

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

    return fig

def makeBarplotSingleDatasets(df_res_single, figsize=(10, 8), c=['green', 'gold', 'navy', 'grey', 'crimson'], width=0.15, labelsAbove=False, saveName=None, dpi=300):

    oneax = None
    fig, ax = plt.subplots(figsize=figsize)
    for pos in np.arange(len(df_res_single)):
        df_temp = df_res_single.iloc[pos].unstack()
        bpos = np.array([-width*8/5, -width*4/5, width*0/5, width*4/5, width*8/5])
        #bpos = np.array([-width*6/4, -width*2/4, width*2/4, width*6/4])
        bars = ax.bar(pos + bpos, df_temp['avg'].values, width, label=df_res_single.index[pos], yerr=df_temp['sem'].values, color=c, 
                      align='center', alpha=1.0, ecolor='black', capsize=2, edgecolor='w', linewidth=0.25)
        if oneax is None:
            oneax = bars

        if labelsAbove:
            t0 = ax.text(pos, df_temp['avg'].max() + 0.08, df_res_single.index.values[pos][0].replace('_', ' '), ha='center')
            t1 = ax.text(pos, df_temp['avg'].max() + 0.05, df_res_single.index.values[pos][1], ha='center')

            t0.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
            t1.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])

    if not labelsAbove:
        ax.set_xticks(range(len(df_res_single)))
        ax.set_xticklabels([v[0].replace('_', ' ').capitalize() + '\n' + v[1].replace('_', ' ') for v in df_res_single.index.values], rotation=90, va='top')
    else:
        ax.set_xticks([])

    ax.set_ylabel('Pearson corr. coef.', fontsize=16)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.legend(oneax, df_temp.index, frameon=False, loc='upper right', fontsize=12)
    ax.set_ylim([-0.1, 1.095])
    ax.tick_params(axis='y', labelsize=16)

    fig.tight_layout()
    
    if not saveName is None:
        fig.savefig(saveName, facecolor='w', dpi=dpi)
    
    return fig

def makeBarplotCrossDatasets(df_res_cross, figsize=(12, 7), c=['green', 'gold', 'navy', 'grey', 'crimson'], width=0.15, labelsAbove=False, saveName=None, dpi=300):

    oneax = None
    fig, ax = plt.subplots(figsize=figsize)
    for pos in np.arange(len(df_res_cross)):
        df_temp = df_res_cross.iloc[pos].unstack()
        bpos = np.array([-width*8/5, -width*4/5, width*0/5, width*4/5, width*8/5])
        #bpos = np.array([-width*6/4, -width*2/4, width*2/4, width*6/4])
        bars = ax.bar(pos + bpos, df_temp['avg'].values, width, label=df_res_cross.index[pos], yerr=df_temp['sem'].values, color=c, 
                      align='center', alpha=1.0, ecolor='black', capsize=2, edgecolor='w', linewidth=0.25)
        if oneax is None:
            oneax = bars

        if labelsAbove:
            t0 = ax.text(pos, df_temp['avg'].max() + 0.08, df_res_cross.index.values[pos][0].replace('_', ' '), ha='center')
            t1 = ax.text(pos, df_temp['avg'].max() + 0.05, df_res_cross.index.values[pos][1], ha='center')
            t2 = ax.text(pos, df_temp['avg'].max() + 0.02, df_res_cross.index.values[pos][2], ha='center')

            t0.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
            t1.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
            t2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])

    if not labelsAbove:
        ax.set_xticks(range(len(df_res_cross)))
        ax.set_xticklabels([v[0].replace('_', ' ').capitalize() + '\n' + v[1].replace('_', ' ') + '-' + v[2].replace('_', ' ') for v in df_res_cross.index.values], rotation=90, va='top')
    else:
        ax.set_xticks([])

    ax.set_ylabel('Pearson corr. coef.', fontsize=16)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.legend(oneax, df_temp.index, frameon=False, loc='upper right', fontsize=12)
    ax.set_ylim([-0.3, 0.9])
    ax.tick_params(axis='y', labelsize=16)

    fig.tight_layout()
    
    if not saveName is None:
        fig.savefig(saveName, facecolor='w', dpi=dpi)
    
    return fig

def drawOneDownsampled(usedf, ax, panel, a, b, c, v, loc='lower right', xlabel='Training set size', ylabel='Pearson corr. coef.', col='val'):

    gb = usedf.groupby(level=[0, 1], axis=0)
    dftemp = gb.mean()[col].unstack(0)
    dftemp.plot(ax=ax)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlim([a - c, b + c])
    ax.set_ylim([0, v])

    dftemp_err = gb.sem()[col].unstack(0)

    for m in dftemp.columns:
        ax.errorbar(dftemp.index, dftemp[m], yerr=dftemp_err[m], fmt='.k', capsize=3);

    l = ax.legend(loc=loc)   
    xdif = ax.get_xlim()[1] - ax.get_xlim()[0]
    ydif = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax.text(ax.get_xlim()[0] - 0.15*xdif, ax.get_ylim()[1] - 0.05*ydif, panel, fontsize=16)

    return
