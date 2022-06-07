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

import plotly.express as px
import plotly.graph_objects as go
from plotly.offline import plot as plot_offline
from plotly.offline import plot_mpl


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


def makeViolinPlot(df_sel, genes, dimPanels, dimCategories, delimiterIn = '|', delimiterOut = ' ', panelWidth = 5, panelHeight = 5, title = '{name} {gene}', exportData = True, xlabel = '$log(count+1)$', ylabel = '', addPoints = True, linesColor = 'black', linesWidth = 1.0, cmap = cm.jet, fontsize = 10, showMedians = True, showExtrema = True, showFractions = True, showMean = True, meanMarkerSize = 7., meanMarkerColor = 'white', meanMarkerShape = 'o', excludeZeroValues = False, violinWidths = 0.85, violinAlpha = 0.9, pointsColor = 'black', pointsSize = 1.0, pointsAlpha = 0.5, pointsPushBack = True, sharex = True, sharey = True, dpi = 300, extension = 'png', saveDir='', colorFactor=None, **kwargs):
    
    ''' Exloratory analysis of the numeric values distributions using matplotlib violinplot.
        Parameters:
            df_sel: pandas.DataFrame
                Table where rows are unique object identifiers, columns are [dimPanels, dimCategories, gene1, gene2, ...]
                    Numeric columns should be without any missing values
            genes: list
                List of genes names to plot, these should be a (sub)set of the df_sel columns
            dimPanels: str
                Name of the categorical variable is for saparation into panels. Option 'All' can be used too
            dimCategories: str
                Name of the categorical variable is for saparation into categories within a panel. Option 'All' can be used too
            panelWidth: float, Default 5
                Width of a panel, including the tick labels
            panelHeight: float, Default 5
                Height of a panel, including the tick labels
            title: str, Default '{name} {gene}'
                Template for panel names
            exportData: float, Default True
                Whether to export data summary into an excel file
            xlabel: str, Default '$log(count+1)$'
                x-axis label
            ylabel: str, Default ''
                y-axis label
            addPoints: boolean, Default True
                Whehter to include scattered points on violins
            linesColor: str, Default 'black' 
                Line color
            linesWidth: float, Default 1.0
                Line width
            cmap: matplotlib.colormap or callable, Default cm.jet
                Colormap or its string name
            fontsize: float, Default 10
                Size of labels font
            showMedians: boolean, Default True
                Whehter to display median
            showExtrema: boolean, Default True
                Whehter to display max and min
            excludeZeroValues: boolean, Default False
                If True then zeros and missing values are not used in calculation of the probability densities
            violinWidths: float, Default 0.85
                Relative violin widths
            violinAlpha: float, Default 0.7
                Transparency of the violins
            pointsColor: str, Default 'black'
                Color of the points
            pointsSize: float, Default 1.0
                Size of the points
            pointsAlpha: float, Default 0.7
                Transparency of the points
            pointsPushBack: boolean, Default True
                If False then points will be drawn in front of all other objects
            sharex: boolean, Default True
                    Whehter to share x-axis
            sharey: boolean, Default True
                    Whehter to share y-axis
            dpi: float, Default 300
                Resolution of the figure
            extension: str, Default 'png'
                Format extension of the figure
        Returns:
            None
        
        Usage:
            DCS.makeViolinPlot(data, ['Numeric 1', 'Numeric 2'], dimPanels='Property A', dimCategories='Property B')
    '''
    
    if dimPanels == 'All' or dimCategories == 'All':
        df_sel['All'] = ['All'] * df_sel.shape[0]

    for dim in [dimPanels, dimCategories]:
        if not dim in df_sel.columns:
            if delimiterIn in dim:
                cols = dim.split(delimiterIn)

                for col in cols:
                    if not col in df_sel.columns:
                        print('Column %s not found' % col)

                        return

                df_sel = df_sel.astype({col: str for col in cols})

                df_sel[dim] = df_sel[cols[0]].copy()
                for col in cols[1:]:
                    df_sel[dim] += delimiterOut + df_sel[col]
            else:
                print('Column %s not found' % dim)

                return

    df_sel = df_sel.astype({dimPanels: str, dimCategories: str})

    df_sel = df_sel.fillna(0.)
        
    panels = np.unique(df_sel[dimPanels].values)
    allCategories = np.sort(df_sel[dimCategories].value_counts().index.values)
    allCategories = np.array(allCategories, dtype=str)[::-1]
    
    n_rows, n_cols = len(genes), len(panels)

    vmin, vmax = df_sel[genes].values.ravel().min(), 1.05 * df_sel[genes].values.ravel().max()
    vmin -= 0.05 * (vmax - vmin)
    
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(panelWidth * n_cols, panelHeight * n_rows), sharex=sharex, sharey=sharey)

    for igene, gene in enumerate(genes):
        for ind, panel in enumerate(panels):

            if n_rows == 1 and n_cols == 1:
                axt = ax
            elif n_rows == 1:
                axt = ax[ind % n_cols]
            elif n_cols == 1:
                axt = ax[igene]
            else:
                axt = ax[igene, ind % n_cols]

            data = df_sel[df_sel[dimPanels] == panel].set_index(dimCategories)[gene].groupby(level=0).agg(list).reindex(allCategories).fillna(0.)
            vdata = [v if type(v) is list else [v] for v in data.values.tolist()]

            if excludeZeroValues:
                pvdata = [np.array(v)[np.array(v)!=0].tolist() for v in vdata]
                pvdata = [v if len(v)>0 else [0] for v in pvdata]
            else:
                pvdata = vdata

            parts = axt.violinplot(pvdata, vert=False, showmedians=showMedians, showextrema=showExtrema, widths=violinWidths)

            if addPoints:
                for i, v in enumerate(vdata):
                    try:
                        axt.scatter(v, 0.75 * (np.random.rand(len(v)) - 0.5) + 1 + i, 
                                    marker='o', color=pointsColor, 
                                    s=pointsSize, zorder=-np.inf if pointsPushBack else np.inf, alpha=pointsAlpha)
                    except:
                        pass

            if showFractions:
                for i, v in enumerate(vdata):
                    try:
                        v = np.array(v)
                        nz = len(v[v!=0])

                        if nz > 0:
                            axt.text(0.975*vmax, 1 + i, '%s%%\n%s' % (np.round(100.*nz/len(v), 1), nz), va='center', ha='right', fontsize=fontsize).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
                    except:
                        pass

            if showMean:
                for i, v in enumerate(vdata):
                    try:
                        m = np.array(v).mean()

                        if m > 0:
                            axt.plot([m], [1 + i], meanMarkerShape, ms=meanMarkerSize, markerfacecolor=meanMarkerColor, markeredgecolor='black')
                    except:
                        pass

            for obj in list(parts):
                try:
                    parts[obj].set(color=linesColor, linewidth=linesWidth)
                except:
                    pass

            try:
                for ipc, pc in enumerate(parts['bodies']):
                    if not colorFactor is None:
                        color = cmap(divmod(ipc, colorFactor)[0] / (len(parts['bodies']) / colorFactor))
                    else:
                        color = cmap(ipc / len(parts['bodies']))
                    pc.set_facecolor(color)
                    pc.set_edgecolor(linesColor)
                    pc.set_alpha(violinAlpha)
            except Exception as exception:
                print(exception)
                pass
                
            axt.set_xlim([vmin, vmax])

            if ind % n_cols == 0 or not sharey:
                axt.tick_params(axis='y', labelsize=fontsize, rotation=0)
                axt.get_yaxis().set_tick_params(direction='out')
                axt.yaxis.set_ticks_position('left')
                axt.set_yticks(np.arange(1, len(allCategories) + 1))
                axt.set_yticklabels(allCategories)
                axt.set_ylim(0.25, len(allCategories) + 0.75)
                
                if ylabel != '':
                    axt.set_ylabel(ylabel, fontsize=fontsize)
                
            if (igene == (len(genes) - 1)) or not sharex:
                axt.tick_params(axis='x', labelsize=fontsize)          

                if xlabel != '':
                    axt.set_xlabel(xlabel, fontsize=fontsize)
                
            axt.set_title(title)
            #axt.set_title(title.format(name=panel, gene=gene) if panel!='All' else gene)

    fig.tight_layout()
    
    saveName = dimPanels + ' ' + dimCategories
    saveName = saveName.replace(delimiterIn, '_')
    saveFigure(fig, saveDir, saveName, extension=extension, dpi=dpi, **kwargs)

    return fig

def makeSankeyDiagram(df, colormapForIndex = None, colormapForColumns = None, linksColor = 'rgba(100,100,100,0.6)', title = '', attemptSavingHTML = False, quality = 4, width = 400, height = 400, border = 20, nodeLabelsFontSize = 15, nameAppend = '_Sankey_diagram', saveDir=''):

    '''Make a Sankey diagram, also known as 'river plot' with two groups of nodes
    Parameters:
        df: pandas.DataFrame 
            With counts (overlaps)
        colormapForIndex: dictionary, Default None
            Colors to use for nodes specified in the DataFrame index
        colormapForColumns: dictionary, Default None
            Colors to use for nodes specified in the DataFrame columns
        linksColor: str, Default 'rgba(100,100,100,0.6)'
            Color of the non-overlapping links
        title: str, Default ''
            Title to print on the diagram
        interactive: boolean , Default False
            Whether to launch interactive JavaScript-based graph
        quality: int, Default 4
            Proportional to the resolution of the figure to save
        nodeLabelsFontSize: int, Default 15
            Font size for node labels
        nameAppend: str, Default '_Sankey_diagram'
            Name to append to the figure file
    Returns:
        None
        
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()
        DCS.makeSankeyDiagram(df)
    '''

    try:
        temp_index = pd.MultiIndex.from_arrays([df.index, [colormapForIndex[item] for item in df.index]], names=['label', 'color'])
        temp_columns = pd.MultiIndex.from_arrays([df.columns, [colormapForColumns[item] for item in df.columns]], names=['label', 'color'])
        df.index = temp_index
        df.columns = temp_columns
    except Exception as exception:
        print('Using default node colors')
        colormapForIndex = None
        colormapForColumns = None

    if (colormapForIndex is None) or (colormapForColumns is None):
        nodeColors = ['rgba(150,0,10,0.8)'] * len(df.index) + ['rgba(10,0,150,0.8)'] * len(df.columns)
        nodeLabels = df.index.to_list() + df.columns.to_list()
    else:
        nodeLabels = df.index.get_level_values('label').to_list() + df.columns.get_level_values('label').to_list()
        nodeColors = df.index.get_level_values('color').to_list() + df.columns.get_level_values('color').to_list()

    sources, targets, values, labels = [], [], [], []
    for i, item in enumerate(df.index):
        sources.extend([i] * len(df.loc[item]))
        targets.extend(list(range(len(df.index), len(df.index) + len(df.loc[item]))))
        values.extend([j for j in df.loc[item].values])
        if type(item) is tuple:
            labels.extend([str(item[0]) + ' -> ' + str(jtem[0]) for jtem in df.loc[item].index])
        else:
            labels.extend([str(item) + ' -> ' + str(jtem) for jtem in df.loc[item].index])

    colorscales = [dict(label=label, colorscale=[[0, linksColor], [1, linksColor]]) for label in labels]

    if not nodeColors is None:
        for i in range(len(sources)):
            if nodeColors[sources[i]] == nodeColors[targets[i]]:
                newColor = ','.join(nodeColors[sources[i]].split(',')[:3] + ['0.6)'])
                colorscales[i] = dict(label=labels[i], colorscale=[[0, newColor], [1, newColor]])

    fig = go.Figure(data=[go.Sankey(valueformat = '', valuesuffix = '', textfont = dict(color = 'rgb(0,0,0)', size = nodeLabelsFontSize, family = 'Arial'),
        node = dict(pad = 20, thickness = 40, line = dict(color = 'white', width = 0.0), label = nodeLabels, color = nodeColors,), # hoverlabel=dict(bordercolor = 'yellow')
        link = dict(source = sources, target = targets, value = values, label = labels, colorscales = colorscales, hoverinfo='all'),)],) #line ={'color':'rgba(255,0,0,0.8)', 'width':0.1}

    if not title is None:
        fig.update_layout(title_text=title, font_size=10)

    fig.update_layout(margin=dict(l=border, r=border, t=border, b=border))

    try:
        fig.write_image(os.path.join(saveDir, nameAppend + '.png'), width=width, height=height, scale=quality)

    except Exception as exception:
        print('Cannot save static image (likely due to missing orca). Saving to interactive html')
        attemptSavingHTML = True

    if attemptSavingHTML:
        fig.update_layout(margin=dict(l=200, r=200, t=100, b=100))
        plot_offline(fig, filename=os.path.join(saveDir, nameAppend + '.html'), auto_open=False)

    return fig

def alignSeries(se1, se2, tagForMissing):

    '''Align two pandas.Series
    Parameters:
        se1: pandas.Series
            Series with the first set of items
        se2: pandas.Series 
            Series with the second set of items
        tagForMissing: str, Default 'Missing'
            Label to assign to non-overlapping items
    Returns:
        pandas.DataFrame
            Contains two aligned pandas.Series
        
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()
        df = DCS.alignSeries(pd.Index(['A', 'B', 'C', 'D']).to_series(), pd.Index(['B', 'C', 'D', 'E', 'F']).to_series())
    '''
        
    se1.index.name = 'index'
    se2.index.name = 'index'

    append = lambda se1, se2: pd.concat([se1, pd.Series(index=se2.index.difference(se1.index), data=[tagForMissing] * len(se2.index.difference(se1.index)))], axis=0, sort=False)

    se1 = append(se1, se2)
    se2 = append(se2, se1)

    se1.name = 'se1'
    se2.name = 'se2'

    return pd.concat((se1, se2.loc[se1.index]), axis=1, sort=True)

def getCountsDataframe(se1, se2, tagForMissing = 'N/A'):

    '''Get a pandas.DataFrame with cross-counts (overlaps) between two pandas.Series
    Parameters:
        se1: pandas.Series
            Series with the first set of items
        se2: pandas.Series 
            Series with the second set of items
        tagForMissing: str, Default 'N/A'
            Label to assign to non-overlapping items
    Returns:
        pandas.DataFrame
            Contains counts
        
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()
        df = DCS.getCountsDataframe(se1, se2)
    '''

    df = alignSeries(se1, se2, tagForMissing)

    counts = {group[0]:{k:len(v) for k, v in group[1].groupby(by='se1').groups.items()} for group in df.reset_index(drop=True).set_index('se2').groupby('se2')}

    df = pd.DataFrame.from_dict(counts).fillna(0.0).astype(int)

    moveTag = lambda df: pd.concat([df.iloc[np.where(df.index != tagForMissing)[0]], df.iloc[np.where(df.index == tagForMissing)[0]]], axis=0, sort=False) if tagForMissing in df.index else df

    return moveTag(moveTag(df.T).T)

def saveFigure(fig, saveDir, label = 'Figure', extension = 'png', dpi = 300, close = True, attemptSavingHTML = False, verbose=0):

    '''Function used internally to save and close figures
    Parameters:
        saveDir: str
            Path of directories to save the object to
        label: str, Default 'Figure'
            Name of the figure to save
        extension: str, Default '.png'
            Path of directories to save the object to
            
        dpi: int, Default 300
            Figure resolution if rasterized
        close: boolean: Default True
            Whether to close the figure after saving
    Returns:
        None
    Usage:
        saveFigure(fig, saveDir, label, extension, dpi)
    '''

    if saveDir != os.path.join('') and not os.path.exists(saveDir):
        os.makedirs(saveDir)

    try:
        if not extension[0] == '.':
            extension = ''.join(['.', extension])
    except Exception as exception:
        if verbose >= 1:
            print(exception)
            print('Figure extension/format error')
            print('Example of acceptable extension: \".png\"')

        return

    if extension in ['.png', '.jpeg', '.tiff']:
        try:
            fig.savefig(os.path.join(saveDir, label + extension), dpi=dpi)
        except Exception as exception:
            if verbose >= 1:
                print(exception)

    elif extension in ['.svg', '.eps', '.pdf']:
        try:
            fig.savefig(os.path.join(saveDir, label + extension))
        except Exception as exception:
            if verbose >= 1:
                print(exception)
    else:
        if verbose >= 1:
            print('Unsupported format. Figure not saved')

    if attemptSavingHTML:
        try:
            plot_mpl(fig, filename=os.path.join(saveDir, label + '.html'), auto_open=False)
        except Exception as exception:
            if verbose >= 1:
                print('Saving to iteractive HTML did not succeed')

    if close:
        try:
            plt.close(fig)
        except Exception as exception:
            if verbose >= 1:
                print(exception)
                print('Error while closing figure')

    return

def makeBarplotSingleDatasets(df_res_single, figsize=(10, 8), c=['green', 'navy', 'grey', 'crimson'], width=0.15, labelsAbove=False, saveName=None, dpi=300):

    oneax = None
    fig, ax = plt.subplots(figsize=figsize)
    for pos in np.arange(len(df_res_single)):
        df_temp = df_res_single.iloc[pos].unstack()
        bpos = np.array([-width*6/4, -width*2/4, width*2/4, width*6/4])
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
    ax.legend(oneax, ['ACDA', 'CDA', 'EN', 'EN-ACDA'], frameon=False, loc='upper right', fontsize=12)
    ax.set_ylim([-0.1, 1.095])
    ax.tick_params(axis='y', labelsize=16)

    fig.tight_layout()
    
    if not saveName is None:
        fig.savefig(saveName, facecolor='w', dpi=dpi)
    
    return fig

def makeBarplotCrossDatasets(df_res_cross, figsize=(12, 7), c=['green', 'navy', 'grey', 'crimson'], width=0.15, labelsAbove=False, saveName=None, dpi=300):

    oneax = None
    fig, ax = plt.subplots(figsize=figsize)
    for pos in np.arange(len(df_res_cross)):
        df_temp = df_res_cross.iloc[pos].unstack()
        bpos = np.array([-width*6/4, -width*2/4, width*2/4, width*6/4])
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
        ax.set_xticklabels([v[0].replace('_', ' ').capitalize() + '\n' + v[1].replace('_', ' ') for v in df_res_cross.index.values], rotation=90, va='top')
    else:
        ax.set_xticks([])

    ax.set_ylabel('Pearson corr. coef.', fontsize=16)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.legend(oneax, ['ACDA', 'CDA', 'EN', 'EN-ACDA'], frameon=False, loc='upper right', fontsize=12)
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
