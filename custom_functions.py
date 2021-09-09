import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rcParams

# Function for ploting of pie charts
def func(pct, allvals):
        absolute = int(pct/100.*np.sum(allvals))
        return "{:.1f}%\n({:d})".format(pct, absolute)
  
# Framework for plotting and violin plots 
def plot_violin_marker(adata, markers, save=None, use_raw=True):
    for i in range(len(markers) // 2 + len(markers) % 2):
        if save is not None:
            sc.pl.violin(
                adata, 
                groupby='leiden', 
                keys=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                use_raw=use_raw, 
                rotation=90,
                size=3,
                save=save+"_"+str(i)+".png"
            )
        else:
            sc.pl.violin(
                adata, 
                groupby='leiden', 
                keys=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                use_raw=use_raw, 
                rotation=90, size=5,
            )

# Fframework for ploting and saving umap plots
def plot_umap_marker(adata, markers, color_map='RdGy_r', size=3, save=None, use_raw=True, frameon=True, title=None,
                     legend_loc='right margin', palette=None, vmin=None, vmax=None, show=None, ax=None):
    for i in range(len(markers) // 2 + len(markers) % 2):
        print(markers[(2*i):np.min([2*(i+1), len(markers)])])
        if save is not None:
            sc.pl.umap(
                adata, 
                color=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                size=size,
                use_raw=use_raw,
                color_map=color_map,
                palette=palette,
                frameon=frameon,
                title=title,
                legend_loc=legend_loc,
                vmin=vmin,
                vmax=vmax,
                show=None,
                ax=None,
                save=save+"_"+str(i)+".pdf"
            )
        else:
            sc.pl.umap(
                adata, 
                color=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                size=size,
                use_raw=use_raw,
                color_map=color_map,
                palette=palette,
                frameon=frameon,
                legend_loc=legend_loc,
                vmin=vmin,
                vmax=vmax,
                show=None,
                ax=None
            )

# Framework for plotting and saving the diffusion map plot
def plot_diffmap_marker(adata, markers, components='1,2', color_map='RdGy_r', palette=None, size=10, save=None, use_raw=True):
    for i in range(len(markers) // 2 + len(markers) % 2):
        print(markers[(2*i):np.min([2*(i+1), len(markers)])])
        if save is not None:
            sc.pl.diffmap(
                adata, 
                components=components,
                color=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                use_raw=use_raw,
                color_map=color_map,
                palette=palette,
                size=size,
                save=save+"_"+str(i)+".png"
            )
        else:
            sc.pl.diffmap(
                adata, 
                components=components,
                color=markers[(2*i):np.min([2*(i+1), len(markers)])], 
                use_raw=use_raw,
                color_map=color_map,
                palette=palette,
                size=size
            )

# Function of automatic calculation and ploting of percentage barplots over the conditions
def cell_percent(adata, cluster, condition, norm='index', xlabel='cell cluster', ylabel='cell count',
                 leg_loc='upper left', title=None, save=False, table=False):
    tmp = pd.crosstab(adata.obs[cluster],
                      adata.obs[condition], 
                      normalize=norm,
                      margins=True
    )
    tmp.plot.bar(stacked=True, color=['lime','blue','magenta','gainsboro',]).legend(loc=leg_loc)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if save is not None:
        plt.savefig(save+title+'.png')
    if table is True:
        if save is not None:
            tmp.to_csv(save+title+'.csv')
        return(tmp)

# Function for performing BH multiiplicity correction
def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr

# Functioin for performing Welch-ANOVA test for single-cell data
def welch_anova(adata):
        
    for i in range(adata.n_vars):
        gen = adata.var.index[i]

        aldh = ch[adata.obs['gfap_aldh']=='aldh_only',gen].X.toarray()
        gfap = ch[adata.obs['gfap_aldh']=='gfap_only',gen].X.toarray()
        both = adata[adata.obs['gfap_aldh']=='both',gen].X.toarray()
        args = (aldh, gfap, both)
        xt = np.concatenate(args)
        xt = xt.ravel()

        #create DataFrame
        df = pd.DataFrame({'score': xt,
                           'group': np.repeat(['aldh', 'gfap', 'both'],
                                              repeats=[len(aldh),len(gfap),len(both)])}) 

        #perform Welch's ANOVA
        oall = pg.welch_anova(dv='score',
                              between='group',
                              data=df)['p-unc']

        pergen = pg.pairwise_gameshowell(dv='score',
                                         between='group',
                                         data=df)['pval']

        gene_name.append(gen)
        pval.append(oall.loc[0])
        aldh_both_pval.append(pergen[0])
        aldh_gfap_pval.append(pergen[1])
        both_gfap_pval.append(pergen[2])

    d = pd.DataFrame()
    d['gene_name'] = gene_name
    d['pval'] = pval
    d['aldh_both_pval'] = aldh_both_pval
    d['aldh_gfap_pval'] = aldh_gfap_pval
    d['both_gfap_pval'] = both_gfap_pval

    # perform BH multipliciity correctioin
    d['pval_adj'] = cf.fdr(d['pval'])
    # sorrt by corrected pval
    dfin = d.sort_values(by=['pval_adj'])

    return dfin
    

def top_marker_as_csv(adata, key_rank_genes, key_groups, output_file):
    
    print(adata.uns[key_rank_genes].keys())
    dict_genes = adata.uns[key_rank_genes].copy()
    print(dict_genes['params'])
    tmp = list(set(adata.obs[key_groups]))
    if dict_genes['params']['reference'] != 'rest':
        tmp.remove(dict_genes['params']['reference'])
    print(tmp, type(tmp))
    del dict_genes['params']
    clusters = sorted(tmp)*len(dict_genes['scores'])
    for key, value in dict_genes.items():
        tmp_list = []
        for n in value:
            tmp_list = tmp_list + list(n)
        dict_genes[key] = tmp_list
    
    dict_genes['cluster'] = clusters
    dataframe = pd.DataFrame(dict_genes)
    dataframe.to_csv(output_file)
    return(dataframe)

def scale_data_5_75(data):
    mind = np.min(data)
    maxd = np.max(data)
    
    if maxd == mind:
        maxd=maxd+1
        mind=mind-1
        
    drange = maxd - mind
    return ((((data - mind)/drange*0.70)+0.05)*100)


def plot_enrich(data, n_terms=20, save=False):
    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError('Please input a Pandas Dataframe output by gprofiler.')
        
    if not np.all([term in data.columns for term in ['p_value', 'name', 'intersection_size']]):
        raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))
    
    data_to_plot = data.iloc[:n_terms,:].copy()
    data_to_plot['go.id'] = data_to_plot.index

    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap="cool", norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'bold'})

    sns.set(style="whitegrid")

    path = plt.scatter(x='recall', y="name", c='p_value', cmap='cool', 
                       norm=colors.LogNorm(min_pval, max_pval), 
                       data=data_to_plot, linewidth=1, edgecolor="grey", 
                       s=[(i+10)**1.5 for i in data_to_plot['scaled.overlap']])
    ax = plt.gca()
    ax.invert_yaxis()

    ax.set_ylabel('')
    ax.set_xlabel('Gene ratio', fontsize=14, fontweight='bold')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    #Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick)/6).astype(int)
    
    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick-1
    
    ticks_vals = [10**i for i in range(max_tick, min_tick-1, -tick_step)]
    ticks_labs = ['$10^{'+str(i)+'}$' for i in range(max_tick, min_tick-1, -tick_step)]

    #Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([0.8, 0.15, 0.03, 0.4])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0,0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs)
    cbar.set_label("Adjusted p-value", fontsize=14, fontweight='bold')

    #Size legend
    min_olap = data_to_plot['intersection_size'].min()
    max_olap = data_to_plot['intersection_size'].max()
    olap_range = max_olap - min_olap
    
    #Note: approximate scaled 5, 25, 50, 75 values are calculated
    #      and then rounded to nearest number divisible by 5
    size_leg_vals = [np.round(i/5)*5 for i in 
                          [min_olap, min_olap+(20/70)*olap_range, min_olap+(45/70)*olap_range, max_olap]]
    size_leg_scaled_vals = scale_data_5_75(size_leg_vals)

    
    l1 = plt.scatter([],[], s=(size_leg_scaled_vals[0]+10)**1.5, edgecolors='none', color='black')
    l2 = plt.scatter([],[], s=(size_leg_scaled_vals[1]+10)**1.5, edgecolors='none', color='black')
    l3 = plt.scatter([],[], s=(size_leg_scaled_vals[2]+10)**1.5, edgecolors='none', color='black')
    l4 = plt.scatter([],[], s=(size_leg_scaled_vals[3]+10)**1.5, edgecolors='none', color='black')

    labels = [str(int(i)) for i in size_leg_vals]

    leg = plt.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False, fontsize=12,
                     handlelength=1, loc = 'center left', borderpad = 1, labelspacing = 1.4,
                     handletextpad=2, title='Gene overlap', scatterpoints = 1,  bbox_to_anchor=(0.1, 1.5), 
                     facecolor='black')

    if save:
        plt.savefig(save, dpi=300, format='pdf', bbox_inches="tight")

    plt.show()
    
    
