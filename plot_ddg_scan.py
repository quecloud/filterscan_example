#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 17:05:16 2020

Plot output from filterscan.

@author: Quinton
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import re
from glob import glob

def scan_from_pdb(file, renumber=False, renumber_start=0, offset=0):
    temp={"WT" : [], "Resi" : [], "Mutation" : [], "ddG" : []}
    with open(file, "r") as f:
        for line in f:
            if line.startswith("DdGScan"):
                wt_resi_mut = line.split(" ")[2]
                ddG = float(line.split(" ")[-1].strip())
                resi = int(re.findall(r'\d+', wt_resi_mut)[0])
                og_resi = resi
                if renumber:
                    if resi >= renumber_start:
                        resi = resi - offset
                wt, mut = wt_resi_mut.split(str(og_resi))
                temp["WT"].append(wt)
                temp["Resi"].append(resi)
                temp["Mutation"].append(mut)
                temp["ddG"].append(ddG)
    result = pd.DataFrame.from_dict(temp, orient="columns")
    result["Source"] = file
    return result

def plot_means(x,y, color, data, group_name, y_name, median_width=0.5):
    #print(x,y, color, data)
    means = data.groupby([group_name])[y_name].mean()

    for n, mean in enumerate(means):
        plt.plot([n-median_width/2, n+median_width/2], [mean,mean], lw=2.5, color=color)

def MAD(dataset):
    if len(dataset) < 3:
        return 100000000
    else:
        med = np.median(dataset)
        return np.median(abs(dataset - med))

def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax[0].get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])

files = glob("/Users/quintondowling/sandbox/avitide_resin/filterscan/output/*/*.pdb")
data = []

for file in files:
    data.append(scan_from_pdb(file, renumber=True, renumber_start=704, offset=500))

df = pd.concat(data)
df["WT_Resi"] = df["WT"] + df["Resi"].astype(str)

i53_50_int_resis = [25,28,29,30,32,33,54,55,56,57,58,63,65,185,188,189,223,282,283,284,297,300,304,306,315]

df=df[~df["Resi"].isin(i53_50_int_resis)]

include_resis = [range(1,203), range(259,1000)]
for resi_set in include_resis:
    df2 = df[df["Resi"].isin(resi_set)].groupby(by=["WT_Resi", "Resi", "Mutation"]).mean().reset_index().pivot_table(columns='WT_Resi',index='Mutation',values='ddG')
    df2 = df2.dropna(axis=1, thresh=10).dropna(axis=0)
    
    #reorder columns by resi
    col_order = {"i" : [], "Resi" : []}
    for i, val in enumerate(df2.columns):
        col_order["Resi"].append(int(re.findall(r'\d+', val)[0]))
        col_order["i"].append(i)
    order_df = pd.DataFrame.from_dict(col_order, orient="columns")
    
    df2 = df2.iloc[:,order_df.sort_values("Resi")["i"]]
    
    
    plt.figure(figsize=(len(df2.columns)/4,5))
    sns.heatmap(df2, annot=False, cmap="hot", xticklabels=True, yticklabels=True, cbar_kws={'label': 'ddG (REU)'})
    
    line="sele interface, resi "
    for resi in df2.columns.str.findall(r'\d+'):
        line = line + str(resi[0]) + "+"

#SUPPRESS IDENT
for resi_set in include_resis:
    df2 = df[df["Resi"].isin(resi_set)].groupby(by=["WT_Resi", "Resi", "Mutation"]).mean().reset_index().pivot_table(columns='Resi',index='Mutation',values='ddG')
    df2 = df2.dropna(axis=1, thresh=10).dropna(axis=0)
    
    #reorder columns by resi
    col_order = {"i" : [], "Resi" : []}
    for i, val in enumerate(df2.columns):
        col_order["Resi"].append(val)
        col_order["i"].append(i)
    order_df = pd.DataFrame.from_dict(col_order, orient="columns")
    
    df2 = df2.iloc[:,order_df.sort_values("Resi")["i"]]
    
    
    plt.figure(figsize=(len(df2.columns)/4,5))
    sns.heatmap(df2, annot=False, cmap="hot", xticklabels=True, yticklabels=True, cbar_kws={'label': 'ddG (REU)'})
