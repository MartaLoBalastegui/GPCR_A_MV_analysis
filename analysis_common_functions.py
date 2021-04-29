# this file is used by the jupyter notebooks of the folder to build the figures
import os,sys,inspect
import re
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import urllib

import mdtraj as md
import itertools
from htmd.ui import *
from htmd.config import config
config(viewer='webgl')
import requests
import os
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator

import seaborn as sns
import matplotlib.pyplot as plt
import math

# this list is obtained from the jupyter notebook, all the names of class A receptors ordered alphabetically to be used as defoult parameter in some functions
all_gpcrs_axes=['5-hydroxytryptamine receptor 1A (HTR1A)', '5-hydroxytryptamine receptor 1B (HTR1B)', '5-hydroxytryptamine receptor 1D (HTR1D)', '5-hydroxytryptamine receptor 1E (HTR1E)', '5-hydroxytryptamine receptor 1F (HTR1F)', '5-hydroxytryptamine receptor 2A (HTR2A)', '5-hydroxytryptamine receptor 2B (HTR2B)', '5-hydroxytryptamine receptor 2C (HTR2C)', '5-hydroxytryptamine receptor 4 (HTR4)', '5-hydroxytryptamine receptor 5A (HTR5A)', '5-hydroxytryptamine receptor 6 (HTR6)', '5-hydroxytryptamine receptor 7 (HTR7)', 'adenosine A1 receptor (ADORA1)', 'adenosine A2a receptor (ADORA2A)', 'adenosine A2b receptor (ADORA2B)', 'adenosine A3 receptor (ADORA3)', 'adrenoceptor alpha 1A (ADRA1A)', 'adrenoceptor alpha 1B (ADRA1B)', 'adrenoceptor alpha 1D (ADRA1D)', 'adrenoceptor alpha 2A (ADRA2A)', 'adrenoceptor alpha 2B (ADRA2B)', 'adrenoceptor alpha 2C (ADRA2C)', 'adrenoceptor beta 1 (ADRB1)', 'adrenoceptor beta 2 (ADRB2)', 'adrenoceptor beta 3 (ADRB3)', 'angiotensin II receptor type 1 (AGTR1)', 'angiotensin II receptor type 2 (AGTR2)', 'apelin receptor (APLNR)', 'arginine vasopressin receptor 1A (AVPR1A)', 'arginine vasopressin receptor 1B (AVPR1B)', 'arginine vasopressin receptor 2 (AVPR2)', 'atypical chemokine receptor 1 (Duffy blood group) (ACKR1)', 'atypical chemokine receptor 2 (ACKR2)', 'atypical chemokine receptor 3 (ACKR3)', 'atypical chemokine receptor 4 (ACKR4)', 'bombesin receptor subtype 3 (BRS3)', 'bradykinin receptor B1 (BDKRB1)', 'bradykinin receptor B2 (BDKRB2)', 'C-C motif chemokine receptor 1 (CCR1)', 'C-C motif chemokine receptor 10 (CCR10)', 'C-C motif chemokine receptor 2 (CCR2)', 'C-C motif chemokine receptor 3 (CCR3)', 'C-C motif chemokine receptor 4 (CCR4)', 'C-C motif chemokine receptor 5 (gene/pseudogene) (CCR5)', 'C-C motif chemokine receptor 6 (CCR6)', 'C-C motif chemokine receptor 7 (CCR7)', 'C-C motif chemokine receptor 8 (CCR8)', 'C-C motif chemokine receptor 9 (CCR9)', 'C-C motif chemokine receptor like 2 (CCRL2)', 'C-X-C motif chemokine receptor 1 (CXCR1)', 'C-X-C motif chemokine receptor 2 (CXCR2)', 'C-X-C motif chemokine receptor 3 (CXCR3)', 'C-X-C motif chemokine receptor 4 (CXCR4)', 'C-X-C motif chemokine receptor 5 (CXCR5)', 'C-X-C motif chemokine receptor 6 (CXCR6)', 'C-X3-C motif chemokine receptor 1 (CX3CR1)', 'cannabinoid receptor 1 (CNR1)', 'cannabinoid receptor 2 (CNR2)', 'chemerin chemokine-like receptor 1 (CMKLR1)', 'cholecystokinin A receptor (CCKAR)', 'cholecystokinin B receptor (CCKBR)', 'cholinergic receptor muscarinic 1 (CHRM1)', 'cholinergic receptor muscarinic 2 (CHRM2)', 'cholinergic receptor muscarinic 3 (CHRM3)', 'cholinergic receptor muscarinic 4 (CHRM4)', 'cholinergic receptor muscarinic 5 (CHRM5)', 'coagulation factor II thrombin receptor (F2R)', 'coagulation factor II thrombin receptor like 2 (F2RL2)', 'complement C3a receptor 1 (C3AR1)', 'complement C5a receptor 1 (C5AR1)', 'complement component 5a receptor 2 (C5AR2)', 'cysteinyl leukotriene receptor 1 (CYSLTR1)', 'cysteinyl leukotriene receptor 2 (CYSLTR2)', 'dopamine receptor D1 (DRD1)', 'dopamine receptor D2 (DRD2)', 'dopamine receptor D3 (DRD3)', 'dopamine receptor D4 (DRD4)', 'dopamine receptor D5 (DRD5)', 'endothelin receptor type A (EDNRA)', 'endothelin receptor type B (EDNRB)', 'F2R like thrombin/trypsin receptor 3 (F2RL3)', 'F2R like trypsin receptor 1 (F2RL1)', 'follicle stimulating hormone receptor (FSHR)', 'formyl peptide receptor 1 (FPR1)', 'formyl peptide receptor 2 (FPR2)', 'formyl peptide receptor 3 (FPR3)', 'free fatty acid receptor 1 (FFAR1)', 'free fatty acid receptor 2 (FFAR2)', 'free fatty acid receptor 3 (FFAR3)', 'free fatty acid receptor 4 (FFAR4)', 'G protein-coupled bile acid receptor 1 (GPBAR1)', 'G protein-coupled estrogen receptor 1 (GPER1)', 'G protein-coupled receptor 1 (GPR1)', 'G protein-coupled receptor 101 (GPR101)', 'G protein-coupled receptor 12 (GPR12)', 'G protein-coupled receptor 132 (GPR132)', 'G protein-coupled receptor 135 (GPR135)', 'G protein-coupled receptor 139 (GPR139)', 'G protein-coupled receptor 141 (GPR141)', 'G protein-coupled receptor 142 (GPR142)', 'G protein-coupled receptor 146 (GPR146)', 'G protein-coupled receptor 148 (GPR148)', 'G protein-coupled receptor 149 (GPR149)', 'G protein-coupled receptor 15 (GPR15)', 'G protein-coupled receptor 150 (GPR150)', 'G protein-coupled receptor 151 (GPR151)', 'G protein-coupled receptor 152 (GPR152)', 'G protein-coupled receptor 153 (GPR153)', 'G protein-coupled receptor 160 (GPR160)', 'G protein-coupled receptor 161 (GPR161)', 'G protein-coupled receptor 162 (GPR162)', 'G protein-coupled receptor 17 (GPR17)', 'G protein-coupled receptor 171 (GPR171)', 'G protein-coupled receptor 173 (GPR173)', 'G protein-coupled receptor 174 (GPR174)', 'G protein-coupled receptor 176 (GPR176)', 'G protein-coupled receptor 182 (GPR182)', 'G protein-coupled receptor 183 (GPR183)', 'G protein-coupled receptor 19 (GPR19)', 'G protein-coupled receptor 20 (GPR20)', 'G protein-coupled receptor 21 (GPR21)', 'G protein-coupled receptor 22 (GPR22)', 'G protein-coupled receptor 25 (GPR25)', 'G protein-coupled receptor 26 (GPR26)', 'G protein-coupled receptor 27 (GPR27)', 'G protein-coupled receptor 3 (GPR3)', 'G protein-coupled receptor 31 (GPR31)', 'G protein-coupled receptor 32 (GPR32)', 'G protein-coupled receptor 33 (gene/pseudogene) (GPR33)', 'G protein-coupled receptor 34 (GPR34)', 'G protein-coupled receptor 35 (GPR35)', 'G protein-coupled receptor 37 (GPR37)', 'G protein-coupled receptor 37 like 1 (GPR37L1)', 'G protein-coupled receptor 39 (GPR39)', 'G protein-coupled receptor 4 (GPR4)', 'G protein-coupled receptor 42 (gene/pseudogene) (GPR42)', 'G protein-coupled receptor 45 (GPR45)', 'G protein-coupled receptor 50 (GPR50)', 'G protein-coupled receptor 52 (GPR52)', 'G protein-coupled receptor 6 (GPR6)', 'G protein-coupled receptor 61 (GPR61)', 'G protein-coupled receptor 62 (GPR62)', 'G protein-coupled receptor 63 (GPR63)', 'G protein-coupled receptor 65 (GPR65)', 'G protein-coupled receptor 68 (GPR68)', 'G protein-coupled receptor 75 (GPR75)', 'G protein-coupled receptor 78 (GPR78)', 'G protein-coupled receptor 82 (GPR82)', 'G protein-coupled receptor 83 (GPR83)', 'G protein-coupled receptor 84 (GPR84)', 'G protein-coupled receptor 85 (GPR85)', 'G protein-coupled receptor 87 (GPR87)', 'G protein-coupled receptor 88 (GPR88)', 'galanin receptor 1 (GALR1)', 'galanin receptor 2 (GALR2)', 'galanin receptor 3 (GALR3)', 'gastrin releasing peptide receptor (GRPR)', 'gonadotropin releasing hormone receptor (GNRHR)', 'GPR119 (GPR119)', 'GPR18 (GPR18)', 'GPR55 (GPR55)', 'growth hormone secretagogue receptor (GHSR)', 'histamine receptor H1 (HRH1)', 'histamine receptor H2 (HRH2)', 'histamine receptor H3 (HRH3)', 'histamine receptor H4 (HRH4)', 'hydroxycarboxylic acid receptor 1 (HCAR1)', 'hydroxycarboxylic acid receptor 2 (HCAR2)', 'hydroxycarboxylic acid receptor 3 (HCAR3)', 'hypocretin receptor 1 (HCRTR1)', 'hypocretin receptor 2 (HCRTR2)', 'KISS1 receptor (KISS1R)', 'leucine rich repeat containing G protein-coupled receptor 4 (LGR4)', 'leucine rich repeat containing G protein-coupled receptor 5 (LGR5)', 'leucine rich repeat containing G protein-coupled receptor 6 (LGR6)', 'leukotriene B4 receptor (LTB4R)', 'leukotriene B4 receptor 2 (LTB4R2)', 'Long-wave-sensitive opsin 1 (OPN1LW)', 'luteinizing hormone/choriogonadotropin receptor (LHCGR)', 'lysophosphatidic acid receptor 1 (LPAR1)', 'lysophosphatidic acid receptor 2 (LPAR2)', 'lysophosphatidic acid receptor 3 (LPAR3)', 'lysophosphatidic acid receptor 4 (LPAR4)', 'lysophosphatidic acid receptor 5 (LPAR5)', 'lysophosphatidic acid receptor 6 (LPAR6)', 'MAS related GPR family member D (MRGPRD)', 'MAS related GPR family member E (MRGPRE)', 'MAS related GPR family member F (MRGPRF)', 'MAS related GPR family member G (MRGPRG)', 'MAS related GPR family member X1 (MRGPRX1)', 'MAS related GPR family member X2 (MRGPRX2)', 'MAS related GPR family member X3 (MRGPRX3)', 'MAS related GPR family member X4 (MRGPRX4)', 'MAS1 proto-oncogene (MAS1)', 'MAS1 proto-oncogene like, G protein-coupled receptor (MAS1L)', 'Medium-wave-sensitive opsin 1 (OPN1MW)', 'melanin concentrating hormone receptor 1 (MCHR1)', 'melanin concentrating hormone receptor 2 (MCHR2)', 'melanocortin 1 receptor (MC1R)', 'melanocortin 2 receptor (MC2R)', 'melanocortin 3 receptor (MC3R)', 'melanocortin 4 receptor (MC4R)', 'melanocortin 5 receptor (MC5R)', 'Melanopsin (OPN4)', 'melatonin receptor 1A (MTNR1A)', 'melatonin receptor 1B (MTNR1B)', 'motilin receptor (MLNR)', 'neuromedin B receptor (NMBR)', 'neuromedin U receptor 1 (NMUR1)', 'neuromedin U receptor 2 (NMUR2)', 'neuropeptide FF receptor 1 (NPFFR1)', 'neuropeptide FF receptor 2 (NPFFR2)', 'neuropeptide S receptor 1 (NPSR1)', 'neuropeptide Y receptor Y1 (NPY1R)', 'neuropeptide Y receptor Y2 (NPY2R)', 'neuropeptide Y receptor Y4 (NPY4R)', 'neuropeptide Y receptor Y5 (NPY5R)', 'neuropeptides B and W receptor 1 (NPBWR1)', 'neuropeptides B and W receptor 2 (NPBWR2)', 'neurotensin receptor 1 (NTSR1)', 'neurotensin receptor 2 (NTSR2)', 'Olfactory receptor 1A1 (OR1A1)', 'Olfactory receptor 1G1 (OR1G1)', 'Olfactory receptor 2T11 (OR2T11)', 'opioid receptor delta 1 (OPRD1)', 'opioid receptor kappa 1 (OPRK1)', 'opioid receptor mu 1 (OPRM1)', 'opioid related nociceptin receptor 1 (OPRL1)', 'Opsin-3 (ECPN)', 'Opsin-5 (OPN5)', 'oxoeicosanoid receptor 1 (OXER1)', 'oxoglutarate receptor 1 (OXGR1)', 'oxytocin receptor (OXTR)', 'platelet activating factor receptor (PTAFR)', 'prokineticin receptor 1 (PROKR1)', 'prokineticin receptor 2 (PROKR2)', 'prolactin releasing hormone receptor (PRLHR)', 'prostaglandin D2 receptor (PTGDR)', 'prostaglandin D2 receptor 2 (PTGDR2)', 'prostaglandin E receptor 1 (PTGER1)', 'prostaglandin E receptor 2 (PTGER2)', 'prostaglandin E receptor 3 (PTGER3)', 'prostaglandin E receptor 4 (PTGER4)', 'prostaglandin F receptor (PTGFR)', 'prostaglandin I2 (prostacyclin) receptor (IP) (PTGIR)', 'purinergic receptor P2Y1 (P2RY1)', 'purinergic receptor P2Y10 (P2RY10)', 'purinergic receptor P2Y11 (P2RY11)', 'purinergic receptor P2Y12 (P2RY12)', 'purinergic receptor P2Y13 (P2RY13)', 'purinergic receptor P2Y14 (P2RY14)', 'purinergic receptor P2Y2 (P2RY2)', 'purinergic receptor P2Y8 (P2RY8)', 'pyrimidinergic receptor P2Y4 (P2RY4)', 'pyrimidinergic receptor P2Y6 (P2RY6)', 'pyroglutamylated RFamide peptide receptor (QRFPR)', 'relaxin/insulin like family peptide receptor 1 (RXFP1)', 'relaxin/insulin like family peptide receptor 2 (RXFP2)', 'relaxin/insulin like family peptide receptor 3 (RXFP3)', 'relaxin/insulin like family peptide receptor 4 (RXFP4)', 'Rhodopsin (OPN2)', 'Short-wave-sensitive opsin 1 (OPN1SW)', 'somatostatin receptor 1 (SSTR1)', 'somatostatin receptor 2 (SSTR2)', 'somatostatin receptor 3 (SSTR3)', 'somatostatin receptor 4 (SSTR4)', 'somatostatin receptor 5 (SSTR5)', 'sphingosine-1-phosphate receptor 1 (S1PR1)', 'sphingosine-1-phosphate receptor 2 (S1PR2)', 'sphingosine-1-phosphate receptor 3 (S1PR3)', 'sphingosine-1-phosphate receptor 4 (S1PR4)', 'sphingosine-1-phosphate receptor 5 (S1PR5)', 'succinate receptor 1 (SUCNR1)', 'tachykinin receptor 1 (TACR1)', 'tachykinin receptor 2 (TACR2)', 'tachykinin receptor 3 (TACR3)', 'thromboxane A2 receptor (TBXA2R)', 'thyroid stimulating hormone receptor (TSHR)', 'thyrotropin releasing hormone receptor (TRHR)', 'trace amine associated receptor 1 (TAAR1)', 'trace amine associated receptor 2 (gene/pseudogene) (TAAR2)', 'trace amine associated receptor 5 (TAAR5)', 'trace amine associated receptor 6 (TAAR6)', 'trace amine associated receptor 8 (TAAR8)', 'trace amine associated receptor 9 (gene/pseudogene) (TAAR9)', 'urotensin 2 receptor (UTS2R)', 'X-C motif chemokine receptor 1 (XCR1)']


def change_width_h(ax, new_value) :
    ''' Changes the height of the bars of 'ax' to 'new_value' '''
    for patch in ax.patches :
        current_width = patch.get_height()
        diff = current_width - new_value
        # we change the bar height
        patch.set_height(new_value)
        
        # we recenter the bar
        patch.set_y(patch.get_y() + diff * .5)

def change_width_v(ax, new_value) :
    '''Changes the width of the bars of 'ax' to 'new_value' '''
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value
        # we change the bar width
        patch.set_width(new_value)
        
        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
        
def sorted_var_info(gpcr_pos,gpcr_pos_count_df, list_binding_site):
    ''' Sort the data frame (names of the GPCRs) as descending number of variants. 
    the input parameters are, the original dataframe (0,1,2,3) the dataframe built
    with the counts of the variants (in order to make the horizontal sum) and the
    list with the GPCRnums that build each of the regions. Return all three elements
    sorted.'''
    gpcr_pos_count_df_sorted=gpcr_pos_count_df.sort_values(by="Count",ascending=False)
    gpcrs_sorted=list(gpcr_pos_count_df_sorted["GPCR"])
    
    gpcr_pos_sorted=pd.DataFrame(columns=list_binding_site,index=gpcrs_sorted).fillna(value=0)
    for gpcr in gpcr_pos_sorted.index:
        gpcr_pos_sorted.loc[gpcr]= gpcr_pos.loc[gpcr]

    return(gpcr_pos_sorted,gpcr_pos_count_df_sorted,gpcrs_sorted)
        
def defineRegion(listOfPositions):
    ''' Depending on the list of GPCRnum, set a diferent xlabel title'''
# bs=['7x37', '6x59', '5x461', '5x37', '23x49', '7x41', '3x29', '2x63', '2x53', '3x35', '6x44', '6x51', '7x35', '3x28', '7x38', '7x31', '6x48', '3x33', '6x52', '6x58', '7x39', '6x56', '3x40', '6x55', '45x50', '5x39', '45x51', '4x57', '5x36', '7x34', '4x56', '5x40', '7x30', '23x50', '2x64', '2x56', '3x37', '5x47', '6x54', '2x60', '2x59', '3x32', '3x36', '5x44', '45x52', '7x42', '5x43', '3x25']
# ibs=['1x59', '1x60', '12x48', '12x49', '12x50', '12x51', '2x34', '2x35', '2x36', '2x37', '2x38', '2x39', '2x40', '2x43', '3x49', '3x50', '3x53', '3x54', '3x55', '3x56', '34x50', '34x51', '34x52', '34x53', '34x54', '34x55', '34x56', '34x57', '4x39', '4x40', '5x61', '5x64', '5x65', '5x66', '5x67', '5x68', '5x69', '5x71', '5x72', '5x73', '5x74', '5x75', '6x24', '6x25', '6x26', '6x28', '6x29', '6x30', '6x32', '6x33', '6x36', '6x37', '6x38', '6x40', '7x53', '7x55', '7x56', '8x47', '8x48', '8x49']
# cm_naBS=['2x50', '3x39', '7x45', '7x46']
# cm_ionicLock=['3x50', '6x30']
# cm_rotTogg=['6x48']
    if listOfPositions[0]== '1x59':
        return("positions of the intracellular binding site ")
    elif listOfPositions[0]== '2x53':
        return("Positions of the binding site")
    elif listOfPositions[0]== '2x50': 
        return("Postions of conserved motifs")
    
def get_count_data(myval, dataFrame, yaxes=all_gpcrs_axes):
    ''' Build a dictionary and a data frame with the values under myval coming from the 
    orignal dataframe'''
    #build dict with values as 0
    gpcr_pos_count_sel=dict(zip(yaxes,[0]*len(yaxes)))
    for gpcridx in gpcr_pos_count_sel.keys():# iterate over the GPCR names
        # count how many times the cells has a value under or equal myval in the original df
        # higher than 0 (bening) lower than 3 (no-data)
        count=len([varval for varval in dataFrame.loc[gpcridx].values if varval >=myval and varval >0 and varval <3])
        gpcr_pos_count_sel[gpcridx]=count
        
    #build new data frame
    gpcr_pos_count_sel_df=pd.DataFrame(list(gpcr_pos_count_sel.items()),
                          columns=['GPCR','Count'])

    df_sel=pd.DataFrame.from_dict(gpcr_pos_count_sel,orient="index").T
    return (gpcr_pos_count_sel_df,df_sel)

def get_counts(dataFrame, varType, yaxes=all_gpcrs_axes):
    ''' Build a dictionary and a data frame with the values that are in the 'range'
    of the corresponding varType, coming from the orignal dataframe'''
    # select the correct scores for each varType
    scores={'benign':[0,1,2,4,5,6],'posDam':[1,2,4,5,6],'dam':[2,4,5,6], 'disPTMjoin':[4,5,6], 'PTMjoin':[5,6], 'join':[6]}
    myList=scores[varType]
    #build dict with values as 0
    gpcr_pos_count_sel=dict(zip(yaxes,[0]*len(yaxes)))
    for gpcridx in gpcr_pos_count_sel.keys():# iterate over the GPCR names
        # count how many times the cells has a valueof the correct list 
        count=len([varval for varval in dataFrame.loc[gpcridx].values if varval in myList])
        gpcr_pos_count_sel[gpcridx]=count
        
    #build new data frame
    gpcr_pos_count_sel_df=pd.DataFrame(list(gpcr_pos_count_sel.items()),
                          columns=['GPCR','Count'])

    df_sel=pd.DataFrame.from_dict(gpcr_pos_count_sel,orient="index").T
    return (gpcr_pos_count_sel_df,df_sel)

def get_count_dist(gpcr_ibspos_count_sel_df,df_sel):
    gpcr_ibspos_count_sel=df_sel.to_dict(orient="records")[0]
    abundance_dist_sel={}
    for gpcr,abundance in gpcr_ibspos_count_sel.items():
        if abundance in abundance_dist_sel:
            abundance_dist_sel[abundance]+=1
        else:
            abundance_dist_sel[abundance]=1

    minkey=min(abundance_dist_sel.keys())
    maxkey=max(abundance_dist_sel.keys())
    for mykey in range(minkey,maxkey):
        if mykey not in abundance_dist_sel:
            abundance_dist_sel[mykey]=0   
            
            
    abundance_dist_sel_df=pd.DataFrame(sorted(list(abundance_dist_sel.items()),key=lambda x: x[0]),
                          columns=['Count','GPCR num'])

    return abundance_dist_sel_df

def create_color_scale(fileInput):
    '''Use the colors from http://www.perbang.dk/rgbgradient/ copied
    into a file to build a list of colors'''
    colorscale=[]
    fileRead=open(fileInput,'r')
    for line in fileRead:
        if len(line)>4:
            new_col="#"+line[0:-1]
            colorscale.append(new_col)
    return(colorscale)


# Build general data frame for each region
def build_dataFrame(listGnum, offsetValue=3, yaxes=all_gpcrs_axes):
    '''Generate a dataframe filled with '3' (by defoult) rows as the GPCR names of the class and columns
    as the positions that define the region'''
    return pd.DataFrame(columns=listGnum,index=yaxes).fillna(value=offsetValue)

def fill_in_dataFrame(dataFrame, listGnum, df, yaxis=all_gpcrs_axes):
    '''Fill in the df with 0,1,2 depending on the impact score '''
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            value=df.loc[(df.Short == gpcr) & (df.GPCRdb == ibspos), 'Impact Prediction'].tolist()
            if value:#if there are matches
                valueOk=max(value)
                if not math.isnan(valueOk):
                    dataFrame[ibspos][gpcrname]=valueOk # we will take tha max in case there
                # are two rows matching
    return None


# heatmap and barplots graphs
def heatmap_only(dataFrame):
    '''Plot the data frame than contains the 0 (pale bening),1 (orange possibly damaging),
    2 (red damaging), 3(white not-data)'''
    fig, ax = plt.subplots(figsize=(40,40))
    ax2=sns.heatmap(dataFrame, cmap=["#DDC8A6","#F3752B","#A20021", "#EDEDF4"], linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

def heatmap_horizontal_sum(dataFrame, listGnum, score,  seeHeatmap=True, order_by_varnum=False, yaxis=all_gpcrs_axes):
    '''Plot the heatmap with the horizontal sum: sum of variants for each GPCR. Depending on
    the score value the sum can be for the damaging variants (2) or the possibly damaging (1).
    The input parameters are the original data frame (with 0,1,2,3) the list of the GPCRnum
    that for the region, score, T/F if you want to sort the GPCRs (yaxis) by descending number
    of variants. The last parameter is the list of names of the GPCRs (yaxis).
    Tune seeHeatmap to see it or not'''
    
    # build a dictionary with keys as the GPCR names values as zero
    gpcr_pos_count=dict(zip(all_gpcrs_axes,[0]*len(yaxis)))

    #score 1, possibly damaging, score 2 damaging
    for gpcridx in gpcr_pos_count.keys():# iterate over the keys (GPCR names)
        # sum the cells values of the dataframe depending on the score value
        count=len([varval for varval in dataFrame.loc[gpcridx].values if varval ==score])
        gpcr_pos_count[gpcridx]=count # set this new value to the corresponding key

    # new data frame, with the counts of each GPCR    
    gpcr_pos_count_df=pd.DataFrame(list(gpcr_pos_count.items()),
                          columns=['GPCR','Count'])

    # if True, order GPCRs by number of variants 
    if order_by_varnum:
        (dataFrame,gpcr_pos_count_df,yaxis)=sorted_var_info(dataFrame,gpcr_pos_count_df, listGnum)
    
    # represent the data frame gpcr_ibspos as barplot
    if seeHeatmap:
        fig, ax = plt.subplots(figsize=(40,40))
        axA=sns.heatmap(dataFrame, cmap=["#DDC8A6","#F3752B","#A20021", "#EDEDF4"], linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
        axA.xaxis.set_ticks_position('top')
        axA.set_xticklabels(labels=dataFrame.columns , rotation=90)
        regionTitle=defineRegion(listGnum)
        axA.xaxis.set_label_text(regionTitle)
        #axA.xaxis.set_label_position('top') 
        plt.show()

    #gpcr_ibspos_count_df 
    fig, ax = plt.subplots(figsize=(4,40))
    pal = sns.color_palette("Blues_d", len(gpcr_pos_count_df))
    rank = gpcr_pos_count_df["Count"].argsort().argsort() 
    axB = sns.barplot(data=gpcr_pos_count_df , x="Count",y="GPCR",palette=np.array(pal[::-1])[rank])
    # tune this if you want to see the GPCR names in the barplot
#     axB.get_yaxis().set_visible(False) # not visible
    axB.get_yaxis() # visible
    axB.xaxis.set_label_text("Number of variants")
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 

    plt.show()
    return None

def heatmap_horizontal_sum_mixed(dataFrame):
    '''Plot the heatmap with the horizontal sum: sum of variants for each GPCR, differenciating
    the damaging and the possibly damaging. Depending on the score value the sum can be for the
    damaging variants (2) or the possibly damaging (1). The input parameters is the original 
    data frame (with 0,1,2,3). Here the order will be the firts one, no by number of variants'''
    
    # build dictionaries and df with only damaging:
    (gpcr_ibspos_count_damaging_df,df_damaging)=get_count_data(2,dataFrame)#damaging
    # and only possibly damaging and damaging
    (gpcr_ibspos_count_all_df,df_all)=get_count_data(1,dataFrame)# possibly damaging and damaging

    abundance_dist_damaging_df=get_count_dist(gpcr_ibspos_count_damaging_df,df_damaging)
    abundance_dist_all_df=get_count_dist(gpcr_ibspos_count_all_df,df_all)
    
    #reprentation of both dataframes in a barplot
    grid = plt.GridSpec(1, 1)
    fig, ax = plt.subplots(figsize=(4,40))

    ax2=plt.subplot(grid[0, 0])
    axB = sns.barplot(data=gpcr_ibspos_count_all_df , x="Count",y="GPCR",color="#F3752B",label="Possibly-damaging variants")
    axB.get_yaxis()
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 
    # change_width_h(axB, 1)

    ax3=plt.subplot(grid[0, 0])
    axB = sns.barplot(data=gpcr_ibspos_count_damaging_df , x="Count",y="GPCR",color="#A20021",label="Damaging variants")
    axB.get_yaxis()
    axB.xaxis.set_label_text("Number of possibly/damaging variants")
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 
    # change_width_h(axB, 1)

    plt.legend()
    plt.show()
    return None

def heatmap_vert_sum_mixed(dataFrame, listGnum, downwards=False):
    gpcr_pos_count_bybs_all=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_damaging=dict(zip(listGnum,[0]*len(listGnum)))

    for pos in gpcr_pos_count_bybs_all.keys():
        count_all=len([varval for varval in dataFrame[pos].values if varval> 0 and varval<3 ])
        count_damaging=len([varval for varval in dataFrame[pos].values if varval> 1 and varval<3])
        gpcr_pos_count_bybs_all[pos]=count_all
        gpcr_pos_count_bybs_damaging[pos]=count_damaging
    
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_all_df=pd.DataFrame(list(gpcr_pos_count_bybs_all.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_damaging.items()),
                          columns=[axisText,'Count'])
    
    #representation downwards
    if downwards:
        grid = plt.GridSpec(1, 1)
        fig, ax = plt.subplots(figsize=(17,5))
        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax3)
    #     axB.get_xaxis()
        axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
        axB.yaxis.set_label_text("Number of variants")
        axB.invert_yaxis()
        change_width_v(axB, 0.9)

        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()

        plt.show()

    # representation upwards
    fig, ax = plt.subplots(figsize=(17,5))
    axB = sns.barplot(data=gpcr_pos_count_bybs_all_df, y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax)
    axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    #axB.invert_yaxis()
    change_width_v(axB, 0.9)

    axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    plt.legend()
    plt.show()
#     plt.savefig("/home/martalo/Documentos/TFM/GPCR_variants/Results/Intracel_Binding_site/ibs_count_positions.png")
    return None


# diseases
def fill_in_dataFrame_diseases(dataFrame, listGnum, df,yaxis=all_gpcrs_axes):
    '''Fill in the df with 0,1,2,.. depending on the number of related diseases.
    MAYBE tune with T/F to only the vars with impact score 2'''
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            value=df.loc[(df.Short == gpcr) & (df.GPCRdb == ibspos), 'related_diseases'].tolist()
            if value:#if there are matches               
                numValues=[x for x in value if x>=0]
                if numValues:
                    valueOk=max(numValues)
#                     print('Max value',valueOk)
                    dataFrame[ibspos][gpcrname]=valueOk # we will take tha max in case there
                # are two rows matching
            #  if value equal to NaN, te cell will be zero
    return None


# heatmap and barplots graphs
def heatmap_only_diseases(dataFrame):
    '''Plot the data frame than contains the color scale for the diseases'''
    color_map=create_color_scale("/home/martalo/Documentos/TFM/GPCR_variants/Data/colors_heatmap_diseases.txt")
    color_map.reverse()# white to 0, then from green to red
    fig, ax = plt.subplots(figsize=(40,40))
    
    ax2=sns.heatmap(dataFrame, cmap=color_map, linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

# ptm
def fill_in_dataFrame_ptm(dataFrame, listGnum, df,yaxis=all_gpcrs_axes):
    '''Fill in the df with 0,1, depending if the variant is a PTM. Return T/F depending
    if there is some value with ptm.'''
    p=re.compile(".*\((\w*)\)")
    somePTM=False
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            value=df.loc[(df.Short == gpcr) & (df.GPCRdb == ibspos), 'PTM'].tolist()             
#             print(value)
            if 'Yes' in value or 'yes' in value:
                somePTM=True
                dataFrame[ibspos][gpcrname]=1 # if it is a PTM position set 1
            #  if value equal to NaN or void, te cell will be zero
    return somePTM

def fill_in_dataFrame_ptm_type(dataFrame, listGnum, typePTM, df,yaxis=all_gpcrs_axes):
    '''Fill in the dataFrame with 0,7,8,9,10 depending of the PTM type of the variant. 
    0 - noPTM   7 - ModRes   8 - Glyco   9 - Lipid  10 - DB '''
    ptmScores={'Modified Residue':7, 'Glycosylation':8, 'Lipidation':9, 'Disulfide Bond':10}
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            value=df.loc[(df.Short == gpcr) & (df.GPCRdb == ibspos), typePTM].tolist()   
#             print(value)
            if len(value)>=1 and type(value[0])!=float:#check that the cell is not void or nan
#                 print('In the if:',gpcrname,ibspos,value)
#                 print('ok')
                dataFrame[ibspos][gpcrname]=ptmScores[typePTM] # if it is a PTM position set 
                # the corresponding score
            #  if value equal to NaN or void, te cell will be zero
    return None

def heatmap_only_ptm(dataFrame):
    '''Plot the data frame than contains the ptm data only, white (no PTM), pink (any PTM)'''
    fig, ax = plt.subplots(figsize=(40,40))
    # white to 0, pink to 1(PTM) (column PTM set to 'yes']
    ax2=sns.heatmap(dataFrame, cmap=["#EDEDF4","#F70F9A"], linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

def fill_in_dataFrame_ptm_all_types(dataFrame,dataFrameModRes, dataFrameGlyco, dataFrameLip,dataFrameDB, listGnum, yaxis=all_gpcrs_axes):
    '''Merge all dataframes and build the list of colors that will be used. We change the PTM
    type numbering in order to avoid extrange things with the color scale in the heatmap. We 
    put consequetive numbers 0-4'''
    colors=set()
    p=re.compile(".*\((\w*)\)")
    ptms=0
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            modRes=dataFrameModRes[ibspos][gpcrname] # 7 -> will be 1     
            glyco=dataFrameGlyco[ibspos][gpcrname] # 8  -> 2     
            lipid=dataFrameLip[ibspos][gpcrname]# 9   -> 3
            db=dataFrameDB[ibspos][gpcrname]# 10 -> 4
            listValues=[modRes,glyco,lipid,db]
#             listValues=[modRes,glyco,lipid]
            for element in listValues:
                if element!=0:
                    ptms+=1# if there is more than one type of PTM in that cell
#                     print(listValues)
            if ptms>=2:
                print('More than one PTM in',ibspos, gpcrname,':',listValues)
            # change the scoring to avoid scaling of the colormap
            if modRes !=0:
                dataFrame[ibspos][gpcrname]= 1 # soft pink
                colors.add("pink1")
            elif glyco !=0:
                dataFrame[ibspos][gpcrname]=2 # 
                colors.add("pink2")
            elif lipid !=0:
                dataFrame[ibspos][gpcrname]=3 # 
                colors.add("pink3")
            elif db !=0:
                dataFrame[ibspos][gpcrname]=4 # dark pink
                colors.add("pink4")
            else:
                dataFrame[ibspos][gpcrname]=0 # no ptm, white
                
            listValues=[]
            ptms=0
    return colors

def heatmap_only_ptm_all_types(dataFrame,colors,x=40,y=40):
    '''Plot the data frame than contains the info of the PTM types in different pinks'''    
    colormap=["#EDEDF4"]# zero values
    if 'pink1' in colors: #  were cells  with 7 -> 1 now
        colormap.append( "#fcacd0")#pink 1, soft
    if 'pink2' in colors:# there were cells with 8 -> 2
        colormap.append( "#f8589e")# pink 2
    if 'pink3' in colors: # there were cells with 9 -> 3
        colormap.append( "#d00d63")#pink 3 
    if 'pink4' in colors: # there were cells with 10 -> 4
        colormap.append( "#801f4a")#pink dark
    print(colormap)
    fig, ax = plt.subplots(figsize=(x,y))
    # no ptm  modRes  glyco   lipid   db
    # 0,         1,     2,      3,    4
    # white,  pink1, pink2,  pink3, pink4
    ax2=sns.heatmap(dataFrame, cmap=colormap, linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

def heatmap_vert_sum_ptms(dataFrame, listGnum, dbs=True, downwards=False):
    '''Horizontal barplot, sum each PTM type and plot them in different pinks'''
    gpcr_pos_count_bybs_all=dict(zip(listGnum,[0]*len(listGnum))) # DBs
    gpcr_pos_count_bybs_modRes=dict(zip(listGnum,[0]*len(listGnum))) # modRes
    gpcr_pos_count_bybs_glyco=dict(zip(listGnum,[0]*len(listGnum))) # glyco
    gpcr_pos_count_bybs_lipid=dict(zip(listGnum,[0]*len(listGnum))) # lipid

    for pos in gpcr_pos_count_bybs_all.keys():
        count_all=len([varval for varval in dataFrame[pos].values if varval==4 ])
        count_modRes=len([varval for varval in dataFrame[pos].values if varval==1])
        count_glyco=len([varval for varval in dataFrame[pos].values if varval==2])
        count_lipid=len([varval for varval in dataFrame[pos].values if varval==3])
        # store sums in the dictionaries position:sum of PTMs
        gpcr_pos_count_bybs_all[pos]=count_all
        gpcr_pos_count_bybs_modRes[pos]=count_modRes
        gpcr_pos_count_bybs_glyco[pos]=count_glyco
        gpcr_pos_count_bybs_lipid[pos]=count_lipid
    
    axisText=defineRegion(listGnum)
    # build dataframes
    gpcr_pos_count_bybs_all_df=pd.DataFrame(list(gpcr_pos_count_bybs_all.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_modRes_df=pd.DataFrame(list(gpcr_pos_count_bybs_modRes.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_glyco_df=pd.DataFrame(list(gpcr_pos_count_bybs_glyco.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_lipid_df=pd.DataFrame(list(gpcr_pos_count_bybs_lipid.items()),
                          columns=[axisText,'Count'])
    
    
    #representation downwards
    if downwards:
        grid = plt.GridSpec(1, 1)
        fig, ax = plt.subplots(figsize=(17,5))
        
        # DBs: they are the sum, as they are bigger in the case of the BS, we set this PTM as the max
        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#801f4a",label="Disulfide Bonds",ax=ax3)
    #     axB.get_xaxis()
        axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
        axB.yaxis.set_label_text("Number of variants")
        axB.invert_yaxis()
        change_width_v(axB, 0.9)
        
        # modRes
        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_modRes_df , y="Count",x=axisText,color="#fcacd0",label="Modified Residues",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()
        
        # Glyco
        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_glyco_df , y="Count",x=axisText,color="#f8589e",label="Glycosylation",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()
        
        # Lipid
        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_lipid_df , y="Count",x=axisText,color="#d00d63",label="Lipidation",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()

        plt.show()

    # representation upwards
    fig, ax = plt.subplots(figsize=(17,5))
    if dbs:
        axB = sns.barplot(data=gpcr_pos_count_bybs_all_df, y="Count",x=axisText,color="#801f4a",label="Disulfide Bond",ax=ax)
        axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
        #axB.get_xaxis().set_visible(False)
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)

    axB = sns.barplot(data=gpcr_pos_count_bybs_modRes_df , y="Count",x=axisText,color="#fcacd0",label="Modified Residues",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    axB = sns.barplot(data=gpcr_pos_count_bybs_glyco_df , y="Count",x=axisText,color="#f8589e",label="Glycosilation",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    axB = sns.barplot(data=gpcr_pos_count_bybs_lipid_df , y="Count",x=axisText,color="#d00d63",label="Lipidation",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    axB.yaxis.set_major_locator(MaxNLocator(integer=True)) # only integers in the y axis
    plt.legend()
    plt.show()
#     plt.savefig("/home/martalo/Documentos/TFM/GPCR_variants/Results/Intracel_Binding_site/ibs_count_positions.png")
    return None

# total: score, disease, ptm
def fill_in_dataFrame_total(dataFrame,dataFrameScore, dataFrameDis, dataFramePtm, listGnum, yaxis=all_gpcrs_axes):
    '''Change to 4,5,6 the cells with score 2 if they are related with a disease, ptm, or
    both.'''
    colors=set()
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            score=dataFrameScore[ibspos][gpcrname] # 0,1,2/3       
            dis=dataFrameDis[ibspos][gpcrname] # 0,inf         
            PTM=dataFramePtm[ibspos][gpcrname]#0,1          
            
            if score == 2: #red (2), blue (4), pink (5) or black (6)
                if dis !=0 and PTM ==1:
                    dataFrame[ibspos][gpcrname]=6 # black
                    colors.add("black")
                elif PTM ==1:
                    dataFrame[ibspos][gpcrname]=5 # pink
                    colors.add("pink")
                elif dis !=0:
                    dataFrame[ibspos][gpcrname]=4 # blue
                    colors.add("blue")
                else:
                    dataFrame[ibspos][gpcrname]=score # red
            else:
                dataFrame[ibspos][gpcrname]=score # 0,1/3 

            
    return colors

def heatmap_only_total(dataFrame,XaxisText, YaxisText,colors,x=40,y=40):
    '''Plot the data frame than contains the info of the score, diseases and PTM'''    
    colormap=["#DDC8A6","#F3752B","#A20021", "#EDEDF4"]
    if 'blue' in colors: # there are cells with 10
        colormap.append( "#09E9FB")
    if 'pink' in colors:# there are cells with 20
        colormap.append( "#F70F9A")
    if 'black' in colors: # there are cells with 30
        colormap.append( "#030303")
    print(colormap)
    fig, ax = plt.subplots(figsize=(x,y))
    # bening, poss, dam,nodata, dis,   ptm, dis-ptm
    # 0,      1,     2,    3,    4,   5,      6
    # pale, orange, red, white, blue, pink, black
    ax2=sns.heatmap(dataFrame, cmap=colormap, linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.set_label_text(XaxisText)
    ax2.yaxis.set_label_text(YaxisText)
    plt.xticks(rotation=90) 
    plt.show()
    return None

# filter aminergic
def filter_aminergic(dataFrame,aminDataFrame,aminList,listGnum,yaxis=all_gpcrs_axes):
    '''Build data frame with only the aminergic GPCRs'''
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        if gpcr in aminList:
#             print(gpcr)
            for ibspos in listGnum: # iterate over positions
                aminDataFrame[ibspos][gpcrname]=dataFrame[ibspos][gpcrname] # 0,1,2,3,4,5,6
    return aminDataFrame

def heatmap_vert_sum_mixed_total(dataFrame, listGnum, downwards=False):
    gpcr_pos_count_bybs_all=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_damaging=dict(zip(listGnum,[0]*len(listGnum)))
    
    gpcr_pos_count_bybs_dis=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_ptm=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_dis_ptm=dict(zip(listGnum,[0]*len(listGnum)))

    for pos in gpcr_pos_count_bybs_all.keys():
        # all: all variants with score 0,1 -> not counting any damaging kind
        count_all=len([varval for varval in dataFrame[pos].values if varval> 0 and varval<3 ])
        count_damaging=len([varval for varval in dataFrame[pos].values if varval==2])
        
        count_dis=len([varval for varval in dataFrame[pos].values if varval==4 ])
        count_ptm=len([varval for varval in dataFrame[pos].values if varval==5 ])
        count_dis_ptm=len([varval for varval in dataFrame[pos].values if varval==6 ])
        
        gpcr_pos_count_bybs_all[pos]=count_all
        gpcr_pos_count_bybs_damaging[pos]=count_damaging
        
        gpcr_pos_count_bybs_dis[pos]=count_dis
        gpcr_pos_count_bybs_ptm[pos]=count_ptm
        gpcr_pos_count_bybs_dis_ptm[pos]=count_dis_ptm
    
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_all_df=pd.DataFrame(list(gpcr_pos_count_bybs_all.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_damaging.items()),
                          columns=[axisText,'Count'])
    
    gpcr_pos_count_bybs_dis_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_ptm_df=pd.DataFrame(list(gpcr_pos_count_bybs_ptm.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_dis_ptm_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis_ptm.items()),
                          columns=[axisText,'Count'])
    
    #representation downwards, review colors
    if downwards:
        grid = plt.GridSpec(1, 1)
        fig, ax = plt.subplots(figsize=(17,5))

        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax3)
    #     axB.get_xaxis()
        axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
        axB.yaxis.set_label_text("Number of variants")
        axB.invert_yaxis()
        change_width_v(axB, 0.9)

        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()

        plt.show()

    # representation upwards
    # all-damaging(all kinds), score 0,1
    fig, ax = plt.subplots(figsize=(17,5))
    axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax)
    axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    #axB.invert_yaxis()
    change_width_v(axB, 0.9)

    # only damaging, score 2
    axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    #damaging + disease, score 4
    axB = sns.barplot(data=gpcr_pos_count_bybs_dis_df , y="Count",x=axisText,color="#09E9FB",label="Damaging variants related with a disease",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    #damaging + ptm, score 5
    axB = sns.barplot(data=gpcr_pos_count_bybs_ptm_df , y="Count",x=axisText,color="#F70F9A",label="Damaging variants and PTM",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    plt.legend()
    plt.show()
#     plt.savefig("/home/martalo/Documentos/TFM/GPCR_variants/Results/Intracel_Binding_site/ibs_count_positions.png")
    return None

def heatmap_vert_sum_mixed_total_allColors(dataFrame, listGnum, legendOut=False ,downwards=False, notBlue=False, notPink=False,notBlack=False, barWidth=0.9, xWidth=17):
    '''Plot the vertical sum of the variants (for each position) differenciating the
    benign, pos dam, only dam, dis, PTM and disPTM in different colors. 
    Legend can be outside the plot.
    The barplot can go downwards.'''
    gpcr_pos_count_bybs_all=dict(zip(listGnum,[0]*len(listGnum))) # 0,1,2,4,5,6 -> benign
    gpcr_pos_count_bybs_pos_dam=dict(zip(listGnum,[0]*len(listGnum))) # 1,2,4,5,6 -> poss dam    
    gpcr_pos_count_bybs_damaging=dict(zip(listGnum,[0]*len(listGnum))) # 2,4,5,6 -> dam 
    
    gpcr_pos_count_bybs_dis_join_ptm_join_dis_ptm=dict(zip(listGnum,[0]*len(listGnum))) # 4,5,6
    gpcr_pos_count_bybs_ptm_join_dis_ptm=dict(zip(listGnum,[0]*len(listGnum)))#5,6
    gpcr_pos_count_bybs_dis_ptm=dict(zip(listGnum,[0]*len(listGnum)))#6
    
    pos_dam=[1,2,4,5,6]
    dam=[2,4,5,6]
    damDisPTMjoin=[4,5,6]
    damPTMjoin=[5,6]
    for pos in gpcr_pos_count_bybs_all.keys():
        # all: all variants with score different than 3 (all variants with prediction)
        count_all=len([varval for varval in dataFrame[pos].values if varval!=3 ])# 0,1,2,4,5,6
        # bening: 
        count_pos_dam=len([varval for varval in dataFrame[pos].values if varval in pos_dam ])# 1,2,4,5,6
        # damaging 
        count_damaging=len([varval for varval in dataFrame[pos].values if varval in dam]) # 2,4,5,6
        # damaging + dis (4) and dam + PTM (5) and dam +PTM + dis (6)
        count_dis_join_ptm_join_dis_ptm=len([varval for varval in dataFrame[pos].values if varval in damDisPTMjoin ])
        # damaging + PTM (5) and dam + PTM + dis (6)
        count_join_dis_ptm=len([varval for varval in dataFrame[pos].values if varval in damPTMjoin ])
        # dam + PTM + dis (6)
        count_dis_ptm=len([varval for varval in dataFrame[pos].values if varval==6 ])
        
        gpcr_pos_count_bybs_all[pos]=count_all
        gpcr_pos_count_bybs_pos_dam[pos]=count_pos_dam
        gpcr_pos_count_bybs_damaging[pos]=count_damaging
        
        gpcr_pos_count_bybs_dis_join_ptm_join_dis_ptm[pos]=count_dis_join_ptm_join_dis_ptm
        gpcr_pos_count_bybs_ptm_join_dis_ptm[pos]=count_join_dis_ptm
        gpcr_pos_count_bybs_dis_ptm[pos]=count_dis_ptm
    
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_all_df=pd.DataFrame(list(gpcr_pos_count_bybs_all.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_pos_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_pos_dam.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_damaging.items()),
                          columns=[axisText,'Count'])
    
    gpcr_pos_count_bybs_dis_ptm_join_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis_join_ptm_join_dis_ptm.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_ptm_join_df=pd.DataFrame(list(gpcr_pos_count_bybs_ptm_join_dis_ptm.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_dis_ptm_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis_ptm.items()),
                          columns=[axisText,'Count'])
    

    # Benign -> 0
    fig, ax = plt.subplots(figsize=( xWidth,5))
    axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#DDC8A6",label="Benign variants",ax=ax)
    axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
    axB.yaxis.set_label_text("Number of variants")
    if downwards:
        axB.invert_yaxis()
        axB.get_xaxis().set_visible(False)
    change_width_v(axB, barWidth)
    
    # Poss damaging -> 1
    axB = sns.barplot(data=gpcr_pos_count_bybs_pos_damaging_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, barWidth)

    # Damaging (2,4,5,6) -> 'only damaging' 2 (red)
    axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, barWidth)
    
    # Damaging PTM Dis and join DisPTM (4,5,6) -> damaging + disease, score 4 (blue)
    if not notBlue:# there are not blue dots
        axB = sns.barplot(data=gpcr_pos_count_bybs_dis_ptm_join_df , y="Count",x=axisText,color="#09E9FB",label="Damaging variants related with a disease",ax=ax)
        #axB.get_xaxis().set_visible(False)
        axB.yaxis.set_label_text("Number of variants")
        change_width_v(axB, barWidth)
    
    # Damaging and join DisPTM (5,6)  -> damaging + PTM, score 5 (pink)
    if not notPink:
        axB = sns.barplot(data=gpcr_pos_count_bybs_ptm_join_df , y="Count",x=axisText,color="#F70F9A",label="Damaging variants at a PTM",ax=ax)
        #axB.get_xaxis().set_visible(False)
        axB.yaxis.set_label_text("Number of variants")
        change_width_v(axB, barWidth)
    
    #damaging + PTM + Dis, score 6 black
    if not notBlack:
        axB = sns.barplot(data=gpcr_pos_count_bybs_dis_ptm_df , y="Count",x=axisText,color="#030303",label="Damaging variants disease related and at PTM",ax=ax)
        #axB.get_xaxis().set_visible(False)
        axB.yaxis.set_label_text("Number of variants")
        change_width_v(axB, barWidth)
    
    axB.yaxis.set_major_locator(MaxNLocator(integer=True)) # only integers in the y axis

    if legendOut:
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        plt.legend()
    plt.show()
#     plt.savefig("/home/martalo/Documentos/TFM/GPCR_variants/Results/Intracel_Binding_site/ibs_count_positions.png")
    return None

# special funtions for CM region
def heatmap_only_total_CM(dataFrame,axisText,x=40,y=40):
    '''Plot the data frame than contains the info of the score, diseases and PTM, special for 
    CM aminergic: no 0 nor 4,5,6'''    
    colormap=["#F3752B","#A20021", "#EDEDF4"]
    fig, ax = plt.subplots(figsize=(x,y))
    # bening, poss, dam,nodata, dis,   ptm, dis-ptm
    # 0,      1,     2,    3,    4,   5,      6
    # pale, orange, red, white, blue, pink, black
    ax2=sns.heatmap(dataFrame, cmap=colormap, linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    ax.set(xlabel=axisText,ylabel='Class A GPCRs') 
    ax2.xaxis.set_label_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

# special functions BS region
def heatmap_vert_sum_mixed_total_BS(dataFrame, listGnum, downwards=False):
    gpcr_pos_count_bybs_all=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_damaging=dict(zip(listGnum,[0]*len(listGnum)))
    
    gpcr_pos_count_bybs_dis=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_ptm=dict(zip(listGnum,[0]*len(listGnum)))
    gpcr_pos_count_bybs_dis_ptm=dict(zip(listGnum,[0]*len(listGnum)))

    for pos in gpcr_pos_count_bybs_all.keys():
        # all: all variants with score 0,1 -> not counting any damaging kind
        count_all=len([varval for varval in dataFrame[pos].values if varval> 0 and varval<3 ])
        count_damaging=len([varval for varval in dataFrame[pos].values if varval==2])
        
        count_dis=len([varval for varval in dataFrame[pos].values if varval==4 ])
        count_ptm=len([varval for varval in dataFrame[pos].values if varval==5 ])
        count_dis_ptm=len([varval for varval in dataFrame[pos].values if varval==6 ])
        
        gpcr_pos_count_bybs_all[pos]=count_all
        gpcr_pos_count_bybs_damaging[pos]=count_damaging
        
        gpcr_pos_count_bybs_dis[pos]=count_dis
        gpcr_pos_count_bybs_ptm[pos]=count_ptm
        gpcr_pos_count_bybs_dis_ptm[pos]=count_dis_ptm
    
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_all_df=pd.DataFrame(list(gpcr_pos_count_bybs_all.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_damaging.items()),
                          columns=[axisText,'Count'])
    
    gpcr_pos_count_bybs_dis_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_ptm_df=pd.DataFrame(list(gpcr_pos_count_bybs_ptm.items()),
                          columns=[axisText,'Count'])
    gpcr_pos_count_bybs_dis_ptm_df=pd.DataFrame(list(gpcr_pos_count_bybs_dis_ptm.items()),
                          columns=[axisText,'Count'])
    
    #representation downwards, review colors
    if downwards:
        grid = plt.GridSpec(1, 1)
        fig, ax = plt.subplots(figsize=(17,5))

        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax3)
    #     axB.get_xaxis()
        axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
        axB.yaxis.set_label_text("Number of variants")
        axB.invert_yaxis()
        change_width_v(axB, 0.9)

        ax3=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax3)
        axB.get_xaxis()
        axB.yaxis.set_label_text("Number of variants")
        #axB.invert_yaxis()
        change_width_v(axB, 0.9)
        plt.tight_layout()

        plt.show()

    # representation upwards
    # all-damaging(all kinds), score 0,1
    fig, ax = plt.subplots(figsize=(17,5))
    axB = sns.barplot(data=gpcr_pos_count_bybs_all_df , y="Count",x=axisText,color="#F3752B",label="Possibly damaging variants",ax=ax)
    axB.set_xticklabels(labels=dataFrame.columns , rotation=90)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    #axB.invert_yaxis()
    change_width_v(axB, 0.9)

    # only damaging, score 2
    axB = sns.barplot(data=gpcr_pos_count_bybs_damaging_df , y="Count",x=axisText,color="#A20021",label="Damaging variants",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    #damaging + disease, score 4
    axB = sns.barplot(data=gpcr_pos_count_bybs_dis_df , y="Count",x=axisText,color="#09E9FB",label="Damaging variants related with a disease",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    #damaging + PTM, score 5
    axB = sns.barplot(data=gpcr_pos_count_bybs_ptm_df , y="Count",x=axisText,color="#F70F9A",label="Damaging variants PTM",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    #damaging + PTM + Dis, score 6
    axB = sns.barplot(data=gpcr_pos_count_bybs_dis_ptm_df , y="Count",x=axisText,color="#030303",label="Damaging, disease + PTM",ax=ax)
    #axB.get_xaxis().set_visible(False)
    axB.yaxis.set_label_text("Number of variants")
    change_width_v(axB, 0.9)
    
    plt.legend()
    plt.show()
#     plt.savefig("/home/martalo/Documentos/TFM/GPCR_variants/Results/Intracel_Binding_site/ibs_count_positions.png")
    return None

def heatmap_only_total_BS(dataFrame,colors,listGnum,amin=False,x=40,y=40):
    '''Plot the data frame than contains the info of the score, diseases and PTM.
    Add True to change the y axis label to aminergic'''    
    colormap=["#DDC8A6","#F3752B","#A20021", "#EDEDF4","#09E9FB","#F70F9A", '#030303']
    print("Benign: #DDC8A6 -> crudo")
    print("Pos Dam: #F3752B -> Orange")
    print("Dam: #A20021 -> Red")
    print("No data: #EDEDF4 -> white")
    print("Dam + Dis: #09E9FB -> blue")
    print("Dam + PTM: #F70F9A -> pink")
    print("Dam + Dis + PTM: #030303 -> black")
    
    print(colormap)
    fig, ax = plt.subplots(figsize=(x,y))
    # bening, poss, dam,nodata, dis,   ptm, dis-ptm
    # 0,      1,     2,    3,    4,   5,      6
    # pale, orange, red, white, blue, pink, black
    ax2=sns.heatmap(dataFrame, cmap=colormap, linewidths=0.1, annot=False, cbar=False,square=True,xticklabels=True, yticklabels=True,ax=ax)
    ax2.xaxis.set_ticks_position('top')
    axisText=defineRegion(listGnum)
    if amin:
        ax.set(xlabel=axisText, ylabel='Aminergic GPCRs')
    else: 
        ax.set(xlabel=axisText,ylabel='Class A GPCRs')
    
    ax2.xaxis.set_label_position('top')
    plt.xticks(rotation=90) 
    plt.show()
    return None

#### clustering 

def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))
    
def compute_serial_matrix(dist_mat,method="ward"):#[New]
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method)
    res_order = seriation(res_linkage, N, N + N-2)
    #seriated_dist = np.zeros((N,N))
    #a,b = np.triu_indices(N,k=1)
    #seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    #seriated_dist[b,a] = seriated_dist[a,b]
    
    return  res_order, res_linkage

def cluster_df(my_df):
    dist_mat_md = squareform(pdist(my_df))
    method="ward"
    res_order, res_linkage = compute_serial_matrix(dist_mat_md,method)

    df_order=[]
    for e in res_order:
        df_order.append(my_df.iloc[e].name)
    return my_df.loc[df_order]

def order_df_by_sorterlist(sorter,gpcr_bspos_count_df):
    sorterIndex = dict(zip(sorter,range(len(sorter))))

    gpcr_bspos_count_df['Rank'] = gpcr_bspos_count_df['GPCR'].map(sorterIndex)

    #gpcr_bspos_count_df.sort(['Player', 'Year', 'Tm_Rank'], ascending = [True, True, True], inplace = True)
    gpcr_bspos_count_df.sort_values(by="Rank",ascending=True,inplace = True)
    gpcr_bspos_count_df.drop('Rank', 1, inplace = True)
    return gpcr_bspos_count_df


def fill_in_dataFrame_total(dataFrame,dataFrameScore, dataFrameDis, dataFramePtm, listGnum, yaxis=all_gpcrs_axes):
    '''Change to 4,5,6 the cells with score 2 if they are related with a disease, ptm, or
    both.'''
    colors=set()
    p=re.compile(".*\((\w*)\)")
    for gpcrname in yaxis:# iteration over the GPCRs
        gpcr=re.search(p,gpcrname).group(1) # only the short name
        for ibspos in listGnum: # iterate over positions
            score=dataFrameScore[ibspos][gpcrname] # 0,1,2/3       
            dis=dataFrameDis[ibspos][gpcrname] # 0,inf         
            PTM=dataFramePtm[ibspos][gpcrname]#0,1          
            if score == 2: #red (2), blue (4), pink (5) or black (6)
                if dis !=0 and PTM ==1:
                    dataFrame[ibspos][gpcrname]=6 # black
                    colors.add("black")
                elif PTM ==1:
                    dataFrame[ibspos][gpcrname]=5 # pink
                    colors.add("pink")
                elif dis !=0:
                    dataFrame[ibspos][gpcrname]=4 # blue
                    colors.add("blue")
                else:
                    dataFrame[ibspos][gpcrname]=score # red
            else:
                dataFrame[ibspos][gpcrname]=score # 0,1/3 

            
    return colors


def get_count_data_damaging(myval, dataFrame, yaxes=all_gpcrs_axes):
    ''' Build a dictionary and a data frame with the damaging values equal to myval coming 
    from the orignal dataframe'''
    #build dict with values as 0
    gpcr_pos_count_sel=dict(zip(yaxes,[0]*len(yaxes)))
    for gpcridx in gpcr_pos_count_sel.keys():# iterate over the GPCR names
        # count how many times the cells has a value uequal myval in the original df
        count=len([varval for varval in dataFrame.loc[gpcridx].values if varval ==myval])
        gpcr_pos_count_sel[gpcridx]=count
        
    #build new data frame
    gpcr_pos_count_sel_df=pd.DataFrame(list(gpcr_pos_count_sel.items()),
                          columns=['GPCR','Count'])

    df_sel=pd.DataFrame.from_dict(gpcr_pos_count_sel,orient="index").T
    return (gpcr_pos_count_sel_df,df_sel)

def heatmap_total_horizontal_sum_mixed_cluster(dataFrame, small=False, dis=False, ptm=False, disPTM=False, legendOut=False, yText=False):
    ''' Build the horizontal barplot, sum for each GPCR, with the disease and PTM info.
    Add True as 2nd parameter if the heatmap will be for only aminergic (smaller df).
    Add True as 3rd parameter if there are ptm values to be displayed (pink).
    Add True as 4th parameter if there are Dis + ptm values to be displayed (black).
    Add True as 5th parameter if you want the legend out of the plot
    dd True as 6th parameter if you want to hide the ylabels'''
    # build dictionaries and df with only damaging:
    rows=list(dataFrame.index.values)
    (gpcr_ibspos_count_ben_df,df_ben)=get_counts(dataFrame,'benign',rows)#bening
    (gpcr_ibspos_count_damaging_df,df_damaging)=get_counts(dataFrame,'dam',rows)#damaging
    # and only possibly damaging and damaging
    (gpcr_ibspos_count_all_df,df_all)=get_counts(dataFrame,'posDam',rows)# possibly damaging and damaging
    
    # dam+dis
    (gpcr_ibspos_count_dis_df,df_dis)=get_counts(dataFrame,'disPTMjoin',rows)#damaging+disease
    # dam+ptm
    (gpcr_ibspos_count_ptm_df,df_ptm)=get_counts(dataFrame,'PTMjoin',rows)#damaging+ptm
    # dam+dis+ptm
    (gpcr_ibspos_count_dis_ptm_df,df_dis_ptm)=get_counts(dataFrame,'join',rows)#damaging+dis+ptm

    abundance_dist_damaging_df=get_count_dist(gpcr_ibspos_count_damaging_df,df_damaging)
    abundance_dist_all_df=get_count_dist(gpcr_ibspos_count_all_df,df_all)
    abundance_dist_bening_df=get_count_dist(gpcr_ibspos_count_ben_df,df_ben)
    
    abundance_dist_dis_df=get_count_dist(gpcr_ibspos_count_dis_df,df_dis)
    abundance_dist_ptm_df=get_count_dist(gpcr_ibspos_count_ptm_df,df_ptm)
    abundance_dist_dis_ptm_df=get_count_dist(gpcr_ibspos_count_dis_ptm_df,df_dis_ptm)
    
    #reprentation of both dataframes in a barplot
    grid = plt.GridSpec(1, 1)
    if not small:
        fig, ax = plt.subplots(figsize=(4,40))
    else: 
        fig, ax = plt.subplots(figsize=(6,11))
    #benign
    ax1=plt.subplot(grid[0, 0])
    axB = sns.barplot(data=gpcr_ibspos_count_ben_df , x="Count",y="GPCR",color="#DDC8A6",label="Benign variants")
    axB.get_yaxis()
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 
    if small: 
        change_width_h(axB, 0.9)
    if not yText:
        axB.get_yaxis().set_visible(False) # not visible
        
    # pos dam
    ax2=plt.subplot(grid[0, 0])
    axB = sns.barplot(data=gpcr_ibspos_count_all_df , x="Count",y="GPCR",color="#F3752B",label="Possibly-damaging variants")
    axB.get_yaxis()
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 
    if small: 
        change_width_h(axB, 0.9)
        
    # dam
    ax3=plt.subplot(grid[0, 0])
    axB = sns.barplot(data=gpcr_ibspos_count_damaging_df , x="Count",y="GPCR",color="#A20021",label="Damaging variants")
    axB.get_yaxis()
    axB.xaxis.set_label_text("Number of variants")
    axB.xaxis.set_ticks_position('top')
    #axB.xaxis.set_label_position('top') 
    if small: 
        change_width_h(axB, 0.9)
    # dis
    if dis:
        ax4=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_ibspos_count_dis_df , x="Count",y="GPCR",color="#09E9FB",label="Damaging variants related with diseases")
        axB.get_yaxis()
        axB.xaxis.set_label_text("Number of variants")
        axB.xaxis.set_ticks_position('top')
        #axB.xaxis.set_label_position('top') 
        if small: 
            change_width_h(axB, 0.9)
    
    # PTM
    if ptm:
        ax5=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_ibspos_count_ptm_df , x="Count",y="GPCR",color="#F70F9A",label="Damaging variants at PTMs")
        axB.get_yaxis()
        axB.xaxis.set_label_text("Number of variants")
        axB.xaxis.set_ticks_position('top')
        #axB.xaxis.set_label_position('top') 
        if small: 
            change_width_h(axB, 0.9)
    # PTM+Dis    
    if disPTM:
        ax6=plt.subplot(grid[0, 0])
        axB = sns.barplot(data=gpcr_ibspos_count_dis_ptm_df , x="Count",y="GPCR",color="#030303",label="Damaging variants disease related at PTMs")
        axB.get_yaxis()
        axB.xaxis.set_label_text("Number of variants")
        axB.xaxis.set_ticks_position('top')
        #axB.xaxis.set_label_position('top') 
        if small: 
            change_width_h(axB, 0.9)

    if legendOut:
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        plt.legend()
    plt.show()
    return None
# mapping
def gnum_all_to_one(num_type,gnum):
    (helix,posnum)=gnum.split(".")
    (bwnum,dbnum)=posnum.split("x")
    if num_type=="gpcrdb":
        gnum_ok=helix+"x"+dbnum
    else:
        gnum_ok=helix+"."+bwnum
    return(gnum_ok)

def open_json(json_filepath):
    json_file=open(json_filepath)
    json_str = json_file.read()
    json_dict=pd.io.json.loads(json_str)
    return json_dict

def gpcr_pdb_dict_fom_json(json_filepath,gnum_type="gpcrdb"):
    gpcr_pdbaa=open_json(json_filepath)
    gpcr_pdbaa_ok={}
    for gnum,v in gpcr_pdbaa.items():
        if "." in gnum or "x" in gnum:
            gnum=gnum_all_to_one(gnum_type,gnum)
        gpcr_pdbaa_ok[gnum]=v
    return gpcr_pdbaa_ok



def freq_score_dic(dataFrame, listGnum, score):
    '''Build a dictionary with keys as positions and values as the number of times the position
    appears with the corresponding score'''
    freq_dic={}
    freq_dic_aux={}
    
    gpcr_pos_count_bybs_damaging=dict(zip(listGnum,[0]*len(listGnum)))

    for pos in gpcr_pos_count_bybs_damaging.keys():
        count_damaging=len([varval for varval in dataFrame[pos].values if varval==score])
        gpcr_pos_count_bybs_damaging[pos]=count_damaging
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_damaging_df=pd.DataFrame(list(gpcr_pos_count_bybs_damaging.items()),
                          columns=[axisText,'Count'])
#     print(gpcr_pos_count_bybs_damaging_df)
    freq_dic_aux=gpcr_pos_count_bybs_damaging_df.set_index(axisText).T.to_dict('list')
    for key, element in freq_dic_aux.items():
        freq_dic[key]=int(freq_dic_aux[key][0])
    
    return freq_dic

def percentage_cells(gpcrs,freq_dict,score):
    '''Compute the % of the cells (gpcr x positions) that have the score'''
    den=len(gpcrs)*len(freq_dict.keys())
    num=0
    for key in freq_dict.keys():
        num+=freq_dict[key]
    return print('The',num/den*100,'% of the cells have score',score)

def freq_disease_dic(dataFrame, listGnum):
    '''Build a dictionary with keys as positions and values as the number of related diseases'''
    freq_dic={}
    freq_dic_aux={}
    
    gpcr_pos_count_bybs_diseases=dict(zip(listGnum,[0]*len(listGnum)))

    for pos in gpcr_pos_count_bybs_diseases.keys():
        count_diseases=dataFrame[pos].sum()# sum the values of the column,
        #(related diseases with variants in that pos)
        gpcr_pos_count_bybs_diseases[pos]=count_diseases
    axisText=defineRegion(listGnum)
    gpcr_pos_count_bybs_diseases_df=pd.DataFrame(list(gpcr_pos_count_bybs_diseases.items()),
                          columns=[axisText,'Count'])
#     print(gpcr_pos_count_bybs_diseases_df)
    freq_dic_aux=gpcr_pos_count_bybs_diseases_df.set_index(axisText).T.to_dict('list')
    for key, element in freq_dic_aux.items():
        freq_dic[key]=int(freq_dic_aux[key][0])
    
    return freq_dic


def create_reps_colormap(gpcr,pdbName,pdbPath,vars_display,freqs,norm=False, rep_chains=False):
    ''' Configuration of a view of a GPCR with its varints of the BS in VDW representantion and
    color depending on the number of times the variants appear in that position for all the GPCRs.
    norm is used to spread the freq values for all the 61 colors, otherwise not all the color scale
    may be reached'''
     
    pdb=pdbName
    affected_vars=[]
    unaffected_vars=[]
    unknown_vars=[]
    rep_chains=set()
    
    color_map=create_color_scale("/home/martalo/Documentos/TFM/GPCR_variants/Data/colors_heatmap_diseases.txt")
    color_map.reverse()# white to 0, then from green to red
    color_dic={}
    freq_set=set(freqs)
    freq_unique_sorted=list(freq_set)
    freq_unique_sorted.sort()
    norm_colors=[]
    print('Frequencies:',freq_unique_sorted)
    # less freq values than colors    
    if len(freq_unique_sorted)<=len(color_map):
        if not norm:# not normalized colors, maybe the top values are not reached
            for i in range(0,len(freq_unique_sorted)):
                color_dic[freq_unique_sorted[i]]=color_map[i]
#                 print('Freq:',freq_unique_sorted[i])
        else:#norm colors, the min and max will be taken
            jump=int(len(color_map)/len(freq_unique_sorted))
            for i in range(0,len(freq_unique_sorted)):
#                 print('Color position:',i*jump,'Color:',color_map[i*jump])
                norm_colors.append(color_map[i*jump])
#             print('Frequencies:')
            for i in range(0,len(freq_unique_sorted)):
#                 print(freq_unique_sorted[i])
                color_dic[freq_unique_sorted[i]]=norm_colors[i]
    else:
        print('You need more colors in the color scale')
    
    if rep_chains: 
        gpcr_sel="protein and (%s)" % " or ".join(rep_chains)
    else:
        gpcr_sel="protein"
    mypdb=pdb+".pdb"
    mypdbpath=os.path.join(pdbPath,mypdb)
    mol = Molecule(mypdbpath)
    mol.reps.add(sel=gpcr_sel,color="#f4f8ff",style='NewCartoon')

    for var in vars_display:
        mysel="(name CA and resid "+var["pdb_pos"]+" and chain "+var["chain"]+")"
        mol.reps.add(sel=mysel ,color=color_dic[var["value"]],style='VDW')

    mol.center()
    #mol.view()
    return mol