#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Biota Technology.
# www.biota.com
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
"""
Parameters
----------
self.file string
    output directory name given from input
self.mpm dataframe
    mixing proportion result from gibbs
self.title string
    title
self.color string
    color scheme for heatmaps pulled from matplotlib
self.out_name string
    output name spaces are replaced with _
"""
import seaborn as sns
import matplotlib.pyplot as plt
import os


class ST_graphs:
    def __init__(self, mpm, output_dir,
                 title='Mixing Proportions', color='viridis'):
        self.file = output_dir
        self.mpm = mpm.T
        self.title = title
        self.color = color
        self.out_name = title.replace(" ", "_")

    def ST_heatmap(self, keep_unknowns=True, annot=True,
                   xlabel='Sources', ylabel='Sinks', vmax=1.0):
        """
        Default Plot for Gibbs method.
        Parameters
        ----------
        unknowns bool
            removes unknown column
        annot bool
            Defines visual proportions in plot
        xlabel string
            x label
        ylabel string
            y label
        vmax float
            determines the maximum color value of the plot
        returns
        --------
        none
        outputs a heatmap visualization in a SVG
        """
        prop = self.mpm
        if keep_unknowns:
            fp_suffix = "_heatmap.svg"
        else:
            fp_suffix = "_heatmap_nounknown.svg"
            prop = prop.drop(['Unknown'], axis=1)
            prop = prop.div(prop.sum(axis=1), axis=0)

        fig, ax = plt.subplots(figsize=((prop.shape[1] * 3 / 4)+4,
                                        (prop.shape[0] * 3 / 4)+4))
        sns.heatmap(prop, vmin=0.0, vmax=vmax, cmap=self.color,
                    annot=annot, linewidths=.5, ax=ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(self.title)
        plt.xticks(rotation=45, ha='right')

        plt.savefig(os.path.join(self.file, self.out_name + fp_suffix),
                    bbox_inches="tight")

    def ST_paired_heatmap(self, normalized=False, keep_unknowns=True,
                          transpose=False, annot=True, ylabel='Sinks',
                          heat_ratio=0.08, legend=True):
        """
        Parameters
        ----------
        normalized bool
            normalize each column to equal to 1 to represent the likelihood
        unknowns bool
            removes unknown column
        transpose bool
            transpose
        annot bool
            Defines visual proportions in plot
        xlabel string
            x label
        ylabel string
            y label
        heat ratio float
            ratio of example bar as compared to main columns.
            Should be much thinner than main columns
        returns
        --------
        none
        outputs a heatmap visualization in a svg defined by each
        individual column

        Any analysis should be done using a bin(n,x)

        Color examples:
        "viridis" "icefire" "vlag" "Spectral" "mako" "magma"
        "coolwarm" "rocket" "flare" "crest"
        "_r" reverses all of these

        vmax and min will show the maximum and minimum for the
        heat map settings. vmax=max is default and what we use here.
        vmin=min is also what we use here but the standard version
        will use vmin=0 in order to show the minimum possible
        and vmax=1.
        The reason I do not in this case is that these ranges
        are not particularly helpful to distinguishing the
        successful matches to each other.

        Paired heatmap should be used for assessing overlap in
        paired individual microbiomes or metabolomic overlap.
        Some examples of this may include organ donors and
        recipients, matching animal fecal samples to track
        movement patterns, overlap in mouse microbiomes
        caged together, and individual identification
        given a range of possible known samples.
        """
        prop = self.mpm
        if not keep_unknowns:
            prop = prop.drop(['Unknown'], axis=1)
            prop = prop.div(prop.sum(axis=1), axis=0)
        if normalized:
            prop = prop.div(prop.sum(axis=0), axis=1)
        tra = ""
        if transpose:
            prop = prop.T
            tra = "_Transposed"
        midpoint = len(prop.columns)/2
        midpoint = round(midpoint)
        ratios, g, axes = [], [], []
        for i in range(len(prop.columns)):
            ratios.append(1)
            axes.append("ax" + str(i))
            g.append("g" + str(i))
        ratios.append(heat_ratio)
        axes.append("axcb")
        fig, axes = plt.subplots(1, len(axes),
                                 gridspec_kw={'width_ratios': ratios},
                                 figsize=((prop.shape[1] * 3 / 4)+4,
                                          (prop.shape[0] * 3 / 4)+4))
        for i in range(len(prop.columns)):
            if i == 0:
                g[i] = sns.heatmap(prop.iloc[:, i:i + 1],
                                   vmin=0, cmap=self.color,
                                   cbar=False, annot=annot, ax=axes[i])
                g[i].set_xlabel("")
                g[i].set_ylabel(ylabel)
            elif i == midpoint:
                g[i] = sns.heatmap(prop.iloc[:, i:i+1], vmin=0,
                                   cmap=self.color, cbar=False,
                                   annot=annot, ax=axes[i])
                g[i].set_xlabel("Sources")
                g[i].set_ylabel("")
                g[i].set_yticks([])
                g[i].set_title(self.title)
            elif i == len(prop.columns) - 1:
                if legend:
                    g[i] = sns.heatmap(prop.iloc[:, i:i + 1], vmin=0,
                                       cmap=self.color, annot=annot,
                                       ax=axes[i],
                                       cbar_ax=axes[i + 1])
                else:
                    g[i] = sns.heatmap(prop.iloc[:, i:i + 1], vmin=0,
                                       cmap=self.color, annot=annot,
                                       ax=axes[i])
                g[i].set_xlabel("")
                g[i].set_ylabel("")
                g[i].set_yticks([])
            else:
                g[i] = sns.heatmap(prop.iloc[:, i:i + 1],
                                   vmin=0, cmap=self.color,
                                   cbar=False, annot=annot,
                                   ax=axes[i])
                g[i].set_xlabel("")
                g[i].set_ylabel("")
                g[i].set_yticks([])
        for ax in g:
            tl = ax.get_xticklabels()
            ax.set_xticklabels(tl, rotation=0)
            tly = ax.get_yticklabels()
            ax.set_yticklabels(tly, rotation=0)
        if normalized:
            if keep_unknowns:
                add_line = tra + "_pairedheatmap_normalized.svg"
            else:
                add_line = tra + "_pairedheatmap_nounknown_normalized.svg"
        else:
            if keep_unknowns:
                add_line = tra + "_pairedheatmap.svg"
            else:
                add_line = tra + "_pairedheatmap_nounknowns.svg"
        plt.savefig(os.path.join(self.file, self.out_name + add_line),
                    bbox_inches="tight")

    def ST_Stacked_bar(self, keep_unknowns=True, x_lab="Sink",
                       y_lab="Source Proportion", coloring=[], flipped=False):
        """
        Creates a Stacked bar plot for the user with direct svg save function
        Parameters
        ----------
        unknowns bool
            removes unknown column
        xlabel string
            x label
        ylabel string
            y label
        coloring string list
            string list of colors to encode the bar plot.
            must have and equal number of colors as sources
        flipped bool
            flips x and y axis
        returns
        --------
        none
        outputs a stacked bar visualization in a SVG

        color example list

        '#1f77b4'Blue, '#ff7f0e'Orange, '#2ca02c'Green, '#d62728'Red,
        '#9467bd'Purple, '#8c564b'Brown, '#e377c2'Pink, '#7f7f7f'Grey,
        '#bcbd22'Gold, '#17becf'Cyan
        make sure to use contrasting colors in order better illuminate
        your data above are some example codes to use
        """
        prop = self.mpm
        if not keep_unknowns:
            prop = prop.drop(['Unknown'], axis=1)
            prop = prop.div(prop.sum(axis=1), axis=0)
        if flipped:
            prop = prop.T
            y_lab_flip = x_lab
            x_lab_flip = y_lab
            y_lab = y_lab_flip
            x_lab = x_lab_flip
            prop = prop.div(prop.sum(axis=1), axis=0)
        prop = prop.reset_index()
        if len(coloring) != (prop.shape[1]-1):
            coloring = []
        if coloring == []:
            coloring = None
        prop.plot(kind='bar', x=prop.columns[0], stacked=True,
                  figsize=((prop.shape[1] * 3 / 4)+4,
                           (prop.shape[0] * 3 / 4)+4),
                  color=coloring)
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.title(self.title)
        plt.autoscale()
        plt.xticks(rotation=45, ha='right')
        fp_suffix = ".svg"
        if not keep_unknowns:
            fp_suffix = "_no_unknowns" + fp_suffix
        if flipped:
            fp_suffix = "_flipped" + fp_suffix
        fp_suffix = "_stacked_bar" + fp_suffix

        plt.savefig(os.path.join(self.file, self.out_name + fp_suffix),
                    bbox_inches="tight")
