import copy
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.patches import Patch
import matplotlib_venn
from itertools import chain
import upsetplot
from collections import Counter
import seaborn as sns
import colorcet as cc

"""Functions that only do the plotting part."""

# Define a few colour lists, shapes and hatches.
marker_shapes = ['o', 'D', 'P', 's', '*', 'X', '<', 'p', 'd', '^', 'v', '>', 'H', '$O$', '$D$', '$U$', '$Y$', '$N$']
hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']

tol_highcontrast = ['#004488', '#DDAA33', '#BB5566']

# Longer list of colours for categorical colours https://colorcet.holoviz.org/user_guide/Categorical.html
# Generated based on Glasbey, Chris; van der Heijden, Gerie & Toh, Vivian F. K. et al. (2007),
# “Colour displays for categorical images”, Color Research & Application 32.4: 304-309.
categ_colours = cc.glasbey_bw_minc_20
del categ_colours[19]  # The yellow colour is almost invisible.
glasbey_palettes = {'glasbey': categ_colours,
                    'glasbey_cool': cc.glasbey_cool,
                    'glasbey_warm': cc.glasbey_warm,
                    'glasbey_dark': cc.glasbey_dark,
                    'glasbey_light': cc.glasbey_light}


def venn_from_list(plot_list, label_list, plot_path, blob_colours=tol_highcontrast, title='',
                   scaled=True, linestyle='', number_size=11, xsize=5, ysize=5, formats=['pdf']):
    """
    Based on a list with the size of the sets and the respective labels, plot a non-scaled / scaled Venn diagram
    for up to three sets. If sets are given, the intersection will be done automatically.
    Choose non-scaled if the difference is too high.
    two sets: [a-b, b-a, a∩b]
    three sets: [a-b-c, b-a-c, a∩b-c, c-a-b, a∩c-b, b∩c-a, a∩b∩c]
     """
    if len(plot_list) > 3:
        print("ERROR, only making Venns for three or two sets")
        return
    if sum([type(x) == set for x in plot_list]) == len(plot_list):
        if len(plot_list) == 2:
            a, b = plot_list
            plot_list = [len(a - b), len(b - a), len(a & b)]
        elif len(plot_list) == 3:
            a, b, c = plot_list
            plot_list = [len(a - b - c), len(b - a - c), len(a & b - c), len(c-a-b), len(a & c - b), len(b & c - a), len(a & b & c)]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    if scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.5)
    elif not scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2_unweighted(subsets=plot_list, set_labels=label_list, ax=ax,
                                             set_colors=blob_colours, normalize_to=0.5)
    elif scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.5)
    elif not scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3_unweighted(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.5)
    if len(plot_list) == 3:
        matplotlib_venn.venn2_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=0.5)
    elif len(plot_list) == 7:
        matplotlib_venn.venn3_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=0.5)
    for text in v.set_labels:
        if text:
            text.set_fontsize(14)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(number_size)
    plt.title(title, size=16, y=1.15)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig((plot_path + title + '_Venn.'+form).replace(' ', ''), bbox_inches='tight', format=form)
    plt.close('All')


def basic_hist(plot_df, x_col, hue_col=None, hue_order=None, bin_num=None, title=None, output_path='', stat='count',
               cumulative=False, palette='tab10', binrange=None, xsize=12, ysize=8, colour='#2d63ad', font_s=14,
               ylabel=None, element='step', alpha=0.3, kde=False, legend_out=False, legend_title=True, fill=True,
               edgecolour=None, multiple='layer', shrink=1, hlines=[], vlines=[], discrete=False, grid=True,
               linewidth=None, formats=['pdf']):
    """
    Plots a basic layered histogram which allows for hue, whose order can be defined as well.
    If x_col is not a column in the df, it will be assumed that hue_col names all the columns which are supposed to be
    plotted.

    Args:
        stat: count: show the number of observations in each bin. frequency: show the number of observations divided
            by the bin width. probability or proportion: normalize such that bar heights sum to 1. percent: normalize
            such that bar heights sum to 100. density: normalize such that the total area of the histogram equals 1.
        element: {“bars”, “step”, “poly”}.
        multiple: {“layer”, “dodge”, “stack”, “fill”}
        discrete: If True, each data point gets their own bar with binwidth=1 and bin_num is ignored.
    """

    if x_col not in plot_df.columns:  # Reformat so that seaborn can interpret x_col as hue.
        plot_df = pd.DataFrame(list(chain(*[[[c, x] for x in plot_df[c].to_numpy()] for c in hue_order])),
                               columns=[hue_col, x_col])
    if palette and 'glasbey' in palette:
        palette = glasbey_palettes[palette][:len(set(plot_df[hue_col]))]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    hist = sns.histplot(data=plot_df, x=x_col, hue=hue_col, ax=ax, alpha=alpha, bins=bin_num if bin_num else 'auto',
                        linewidth=linewidth if linewidth is not None else (1 if fill else 3), color=colour, palette=palette if hue_col else None,
                        hue_order=hue_order, multiple=multiple, fill=fill, element=element, stat=stat,
                        discrete=discrete,
                        common_norm=False, cumulative=cumulative, binrange=binrange, kde=kde, shrink=shrink)
    for patch in ax.patches:
        if not edgecolour:
            clr = patch.get_facecolor()
            patch.set_edgecolor(clr)
        else:
            patch.set_edgecolor(edgecolour)
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(stat, fontsize=font_s+2)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=font_s+2)
    ax.set_xlabel(str(x_col), fontsize=font_s+2)
    if hue_col:
        plt.setp(hist.get_legend().get_texts(), fontsize=font_s)
        plt.setp(hist.get_legend().get_title(), fontsize=font_s+2)
        if not legend_title:
            hist.get_legend().set_title('')
        if legend_out:
            sns.move_legend(hist, prop={'size': font_s, 'weight': 'bold'}, loc='upper right',
                            bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            sns.move_legend(hist, prop={'size': font_s, 'weight': 'bold'}, loc='best')
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + str(x_col) + '_' + str(hue_col) + '_Hist.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def upset_plotter(inter_sets, max_groups=None, sort_by='cardinality', y_label='Intersection', title_tag='',
                  show_percent=False, plot_path='', min_degree=0, sort_categories_by='cardinality',
                  intersection_plot_elements=None, font_enhancer=0, element_size=None, formats=['pdf']):
    """
    Based on a dictionary with sets as values creates the intersection and an upsetplot.

    Args:
        max_groups: defines the maximum number of intersections plotted, sorted descending by size
        sort_categories_by: cardinality, degree or input.
        font_enhancer: what to add onto the default sizes, does not affect the title.
        intersection_plot_elements: height
        element_size: ~overall size and margins
    """
    # totals_plot_elements: ~size of the horizontal bars for total size
    if sort_categories_by == 'input':  # Reverse the order, to have it from top to bottom.
        inter_sets = {k: inter_sets[k] for k in list(inter_sets.keys())[::-1]}
        sort_categories_by = None
    intersection = upsetplot.from_contents(inter_sets)
    fig = plt.figure(figsize=(10 + int(len(inter_sets) / 2), 7 + int(len(inter_sets) / 2)))
    max_string = max([len(k) for k in inter_sets])
    if max_groups:
        upset = upsetplot.UpSet(intersection, show_counts=False, intersection_plot_elements=0, min_degree=min_degree,
                                with_lines=True, totals_plot_elements=max_groups - 2, element_size=60, sort_by=None,
                                sort_categories_by=sort_categories_by)
        cut_off = upset.intersections.sort_values(ascending=False).tolist()[max_groups]
        filtered_sets = set(upset.intersections[upset.intersections > cut_off].index.values)
        base_totals = upset.totals
        filtered_intersection = intersection[intersection.index.isin(filtered_sets)]
        filtered_upset = upsetplot.UpSet(filtered_intersection, show_counts=True, min_degree=min_degree, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets) if not intersection_plot_elements else intersection_plot_elements,
                                         facecolor="#010d4a",
                                         show_percentages=show_percent, totals_plot_elements=max(len(base_totals) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30 if not element_size else element_size, sort_by=sort_by,
                                         sort_categories_by=sort_categories_by)
        filtered_upset.totals = base_totals
    else:
        filtered_upset = upsetplot.UpSet(intersection, show_counts=True, min_degree=min_degree, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets) if not intersection_plot_elements else intersection_plot_elements,
                                         facecolor="#010d4a",
                                         show_percentages=show_percent, totals_plot_elements=max(len(inter_sets) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30 if not element_size else element_size, sort_by=sort_by,
                                         sort_categories_by=sort_categories_by)

    with plt.rc_context({'axes.titlesize': 20,
                         'axes.labelsize': 14+font_enhancer,
                         'xtick.labelsize': 11+font_enhancer,  # Ensures enough space between horizontal bars and UpSet.
                         'ytick.labelsize': 12+font_enhancer,
                         'font.size': 11+font_enhancer,  # For whatever reason font.size is only for the counts on top of the bars.
                         'legend.fontsize': 14+font_enhancer}):
        filtered_upset.plot(fig)

    ax = plt.gca()
    ax.set_ylabel(y_label)
    ax.set_title("UpSet " + title_tag + '\n#' + str(len(intersection)), fontsize=16, fontweight='bold', y=1.05)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig((plot_path + '_UpSet.'+form).replace(' ', ''), bbox_inches='tight', format=form)
    plt.close()


def stacked_bars(plot_df, x_col, y_cols, y_label='', title=None, output_path='', x_size=8, y_size=6,
                 rotation=None, palette=None, legend=True, fraction=False, numerate=False, sort_stacks=True,
                 legend_out=False, width=0.8, vertical=False, hatches=None, font_s=14, formats=['pdf']):
    """
    Plots a stacked barplot, with a stack for each y_col.

    Args:
        fraction: If True take all values as fraction of the row sum.
        hatches: If given assumes the colour list is meant for the x-axis.
    """
    plot_df = copy.deepcopy(plot_df)
    if vertical:  # To start the index at the top.
        plot_df = plot_df.iloc[::-1]
    count_dict = plot_df[y_cols].sum(axis=1).to_dict()
    if hatches and type(hatches) != list:
        hatches = hatches[:len(y_cols)]
        print(hatches)
    if fraction:
        y_label = 'Fraction of ' + y_label
        plot_df = plot_df.div(plot_df[y_cols].sum(axis=1).values, axis='rows')
    if sort_stacks:
        sorted_cols = list(plot_df.sum().sort_values(ascending=False).index)
        y_cols = sorted_cols
        plot_df.columns = pd.CategoricalIndex(plot_df.columns.values,
                                              ordered=True,
                                              categories=sorted_cols)
        plot_df = plot_df.sort_index(axis=1)
    if x_col not in plot_df.columns:  # Assumes the x_col is the index if the column doesn't exist.
        plot_df[x_col] = list(plot_df.index)
    if palette and 'glasbey' in palette:
        palette = glasbey_palettes[palette]
    if type(palette) == str:
        cmap = cm.get_cmap(palette, len(y_cols))
        palette = [cmap(i) for i in range(len(y_cols))]

    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    plot_df.plot(x=x_col, y=y_cols, kind='barh' if vertical else 'bar', stacked=True, ax=ax, alpha=1,
                 linewidth=2 if hatches else 0, color=palette, width=width, edgecolor='black' if hatches else None)
    if hatches:
        stacks = [thing for thing in ax.containers if isinstance(thing, mpl.container.BarContainer)]
        for s, stack in enumerate(stacks):  # For some reason the iteration is the stack layers.
            for p, patch in enumerate(stack):
                patch.set_hatch(hatches[s])
                patch.set_facecolor(palette[p])
    if legend_out:
        ax.legend(prop={'size': font_s, 'weight': 'bold'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
    else:
        ax.legend(prop={'size': font_s, 'weight': 'bold'})

    if hatches:  # TODO might not work when the y_cols are sorted in the function
        legend_elements = [Patch(facecolor='w', edgecolor='black', linewidth=2, label=c_label, hatch=hatches[p])
                           for p, c_label in enumerate(y_cols)]
        if legend_out:
            ax.legend(handles=legend_elements, prop={'size': 30}, loc='upper right',
                      bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            ax.legend(handles=legend_elements, prop={'size': 30})

    if not legend and ax.get_legend():
        ax.get_legend().remove()
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(y_label if not vertical else x_col, fontsize=font_s+2)
    ax.set_xlabel(x_col if not vertical else y_label, fontsize=font_s+2)
    if numerate:
        if not vertical:
            ax.set_xticklabels([x.get_text()+'\n(#' + str(int(count_dict[x.get_text()])) + ')' for x in ax.get_xmajorticklabels()],
                               ha='right' if rotation != 90 and rotation != 0 else 'center')
        else:
            ax.set_yticklabels([y.get_text()+'\n(#' + str(int(count_dict[y.get_text()])) + ')' for y in ax.get_ymajorticklabels()])
    if rotation is not None:
        ax.tick_params(axis='x', rotation=rotation)
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + '_'.join([str(x_col), str(y_label)]) + '_'+"Frac"*fraction+'StackedBars.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_violin(plot_df, y_col, x_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='',
                 numerate=False, ylim=None, palette=None, xsize=12, ysize=8, boxplot=False, boxplot_meanonly=False,
                 rotation=None, numerate_break=True, jitter=False, colour='#2d63ad', font_s=14, saturation=0.75,
                 jitter_colour='black', jitter_size=5, vertical_grid=False, legend_title=True, legend=True, grid=True,
                 formats=['pdf']):
    """
    Plots a basic violin plot which allows for hue, whose order can be defined as well.
    Use y_col=None and x_col=None for seaborn to interpret the columns as separate plots on the x-asis.

    Args:
        boxplot_meanonly: Remove all lines from the boxplot and show just the mean as horizontal line.
    """
    if palette and 'glasbey' in palette:
        palette = glasbey_palettes[palette]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if not boxplot:
        vio = sns.violinplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax,
                             color=colour if not palette else None, saturation=saturation,
                             palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    else:
        if boxplot_meanonly:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, hue_order=hue_order,
                              showfliers=False, showbox=False, showcaps=False, showmeans=True, meanline=True,
                              meanprops={'color': 'k', 'ls': '-', 'lw': 2}, medianprops={'visible': False},
                              whiskerprops={'visible': False}, zorder=1, saturation=saturation,
                              palette='tab10' if hue_col and not jitter_colour else jitter_colour)
        else:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, saturation=saturation,
                              color=colour if not palette else None, showfliers=False if jitter else True,
                              palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    if jitter:
        sns.stripplot(data=plot_df, x=x_col, y=y_col, jitter=True, ax=ax, hue=hue_col, hue_order=hue_order, zorder=10,
                      order=x_order, palette=jitter_colour, dodge=True, legend=False, edgecolor='black', linewidth=1,
                      size=jitter_size)
    ax.tick_params(axis='both', labelsize=font_s+4)
    ax.set_ylabel(y_col, fontsize=font_s+8)
    ax.set_xlabel(x_col, fontsize=font_s+8)
    if ylim:
        ax.set_ylim(ylim)
    if numerate:
        if not x_col:
            ax.set_xticklabels(['(#' + str((~plot_df[y_col].isna()).sum()) + ')' for x in ax.get_xmajorticklabels()])
        else:
            count_df = plot_df[[x_col, y_col]][~plot_df[y_col].isna()]
            x_counts = Counter(count_df[x_col].values)
            ax.set_xticklabels([x._text+'\n'*numerate_break+'(#'+str(x_counts[x._text])+')' for x in ax.get_xmajorticklabels()])
    if hue_col:
        plt.setp(vio.get_legend().get_texts(), fontsize=font_s)
        plt.setp(vio.get_legend().get_title(), fontsize=font_s+2)
        if not legend_title:
            vio.get_legend().set_title('')
        sns.move_legend(vio, prop={'size': 14, 'weight': 'bold'}, loc='best')
    if rotation:
        plt.xticks(rotation=rotation, ha='center')
    if vertical_grid:  # Fun part is minor ticks are always x5.
        for x in range(len(set(plot_df[x_col]))):
            plt.axvline(x+0.5, color='#f2f2f2', linewidth=1, zorder=0)
    if not legend:
        ax.get_legend().remove()
    plt.title(title, fontsize=22, fontweight='bold', y=1.02)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + str(x_col) + '_' + str(y_col) + '_' + str(hue_col) + '_Violin.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()



