import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.offline import plot as plot_offline
from plotly.offline import plot_mpl
from matplotlib import cm

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
