
import json
import pandas as pd
from htmd.ui import *
from htmd.config import config
config(viewer='webgl')
from create_csv import *


def create_reps_all(gpcr,var_info,path_to_project=False):
    
    name=var_info[gpcr]["name"]
    print("%s (%s)" % (name,gpcr))
    pdb=var_info[gpcr]["pdb"]
    is_model=var_info[gpcr]["is_model"]
    gpcr_vars=[]
    if len(list(var_info[gpcr].keys())) >1:
        if path_to_project:
            dictpath=os.path.join(path_to_project,"Results/gpcr_to_pdb/")
            gpcr_vars=var_to_display(ex_gpcr=gpcr,dictpath=dictpath)   
        else:
            gpcr_vars=var_to_display(gpcr)   
    affected_vars=[]
    unaffected_vars=[]
    unknown_vars=[]
    rep_chains=set()
    for var in gpcr_vars:
        affected_pred=var["affected_pred"]
        if affected_pred:
            if affected_pred=="?":
                unknown_vars.append("(resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
            else:
                affected_vars.append("(resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                #print("%s\t%s" % (var["gpcr_num"],", ".join(var["affected_items_pred"])))
        else:
            unaffected_vars.append("(resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
        rep_chains.add("chain %s" % var["pdb_pos"]["chain"])

    unknown_vars_s=" or ".join(unknown_vars)
    affected_vars_s=" or ".join(affected_vars)
    unaffected_vars_s=" or ".join(unaffected_vars)

    gpcr_sel="protein and (%s)" % " or ".join(rep_chains)

    #Create htmd molecule of pdb
    if path_to_project:
        mypdb=pdb+".pdb"
        pdb_path=os.path.join(path_to_project,"Data/pdb_files/")
        mypdbpath=os.path.join(pdb_path,mypdb)
    else:
        mypdbpath="/home/martalo/Documentos/TFM/Pharmacogenomics/Data/pdb_files/"+pdb+".pdb"
    mol = Molecule(mypdbpath)
    #mol.reps.remove()
    mol.reps.add(sel=gpcr_sel,color="#f4f8ff",style='NewCartoon')
    if unknown_vars_s:
        mol.reps.add(sel=unknown_vars_s ,color="#ffdd00",style='licorice')#yellow
    if affected_vars_s:
        mol.reps.add(sel=affected_vars_s ,color="#f9020a",style='licorice')#red
    if unaffected_vars_s:
        mol.reps.add(sel=unaffected_vars_s ,color="#00c631",style='licorice')#green
    mol.center()
    #mol.view()
    return mol
#{'interact', 'ptm', 'arrestinInt', 'sodiumP', 'microSwitch', 'DisGeNet_disease', 'activation', 'gprotInt', 'sift_poloP'}


# In[16]:


def create_reps_affected(gpcr,var_info,show_aff_items=False,restrict_vars=False,path_to_project=False):
    plot_only_del=True
    affected_items_d={'interact': {"color":"#3875B1","varsel":[]}, #dark blue
                      'ptm': {"color":"#F67D26","varsel":[]}, #orange
                      'arrestinInt+gprotInt': {"color":"#3EA233","varsel":[]}, #green 
                      'sodiumP': {"color":"#CD1E31","varsel":[]}, #red
                      'microSwitch': {"color":"#9161BB","varsel":[]}, #purple
                      'DisGeNet_disease': {"color":"#87554D","varsel":[]}, #brown
                      'activation':{"color":"#DD70C1","varsel":[]}, #pink
                      'sift_poloP':{"color":"#B8BF30","varsel":[]} #yellow-ish green
                     } #https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
    
    name=var_info[gpcr]["name"]
    pdb=var_info[gpcr]["pdb"]
    is_model=var_info[gpcr]["is_model"]
    gpcr_vars=[]
    if len(list(var_info[gpcr].keys())) >1:
        if path_to_project:
            dictpath=os.path.join(path_to_project,"Results/gpcr_to_pdb/")
            gpcr_vars=var_to_display(ex_gpcr=gpcr,dictpath=dictpath)   
        else:
            gpcr_vars=var_to_display(gpcr)   
    affected_vars=[]
    unaffected_vars=[]
    unknown_vars=[]
    rep_chains=set()
    for var in gpcr_vars:
        rep_chains.add("chain %s" % var["pdb_pos"]["chain"])
        if  restrict_vars and var["gpcr_num"] not in restrict_vars:
            continue
        affected_pred=var["affected_pred"]
        if affected_pred:
            if affected_pred=="?":
                unknown_vars.append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
            else:
                affected_vars.append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                if show_aff_items:
                    print("%s\t%s" % (var["gpcr_num"],", ".join(var["affected_items_pred"])))
                affected_items=var["affected_items_pred"] 
                if "sodiumP" in affected_items:
                    affected_items_d["sodiumP"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "microSwitch" in affected_items:
                    affected_items_d["microSwitch"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "activation" in affected_items:
                    affected_items_d["activation"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "ptm" in affected_items:
                    affected_items_d["ptm"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "arrestinInt" in affected_items or "gprotInt" in affected_items:
                    affected_items_d["arrestinInt+gprotInt"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "interact" in affected_items:
                    affected_items_d["interact"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "DisGeNet_disease" in affected_items:
                    affected_items_d["DisGeNet_disease"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                if "sift_poloP" in affected_items:
                    affected_items_d["sift_poloP"]["varsel"].append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")
                    continue
                    
        else:
            unaffected_vars.append("(name CA and resid "+var["pdb_pos"]["pos"]+" and chain "+var["pdb_pos"]["chain"]+")")

    unknown_vars_s=" or ".join(unknown_vars)
    affected_vars_s=" or ".join(affected_vars)
    unaffected_vars_s=" or ".join(unaffected_vars)
    gpcr_sel="protein and (%s)" % " or ".join(rep_chains)

    #Create htmd molecule of pdb
    if path_to_project:
        mypdb=pdb+".pdb"
        pdb_path=os.path.join(path_to_project,"Data/pdb_files/")
        mypdbpath=os.path.join(pdb_path,mypdb)
    else:
        mypdbpath="/home/martalo/Documentos/TFM/Pharmacogenomics/Data/pdb_files/"+pdb+".pdb"
    mol = Molecule(mypdbpath)
    #mol.reps.remove()
    mol.reps.add(sel=gpcr_sel,color="#f4f8ff",style='NewCartoon')
    if unaffected_vars_s:
        mol.reps.add(sel=unaffected_vars_s ,color="#f4f8ff",style='VDW')#green
    for aff_item in affected_items_d.keys():
        if affected_items_d[aff_item]["varsel"]:
            mysel=" or ".join(affected_items_d[aff_item]["varsel"])
            mol.reps.add(sel=mysel ,color=affected_items_d[aff_item]["color"],style='VDW')
    if unknown_vars_s:
        mol.reps.add(sel=unknown_vars_s ,color="#adadad",style='VDW')#grey
    mol.center()
    #mol.view()
    return mol
#{'interact', 'ptm', 'arrestinInt', 'sodiumP', 'microSwitch', 'DisGeNet_disease', 'activation', 'gprotInt', 'sift_poloP'}





def gpcrpdb_to_pdbgpcr(gpcr_pdbaa):
    pdb_gpcr={}
    for gnum,pdbpos in gpcr_pdbaa.items():
        if "." in gnum or "x" in gnum:
            gnum=gnum_all_to_one("gpcrdb",gnum)
            pdb_gpcr[int(pdbpos["pos"])]=gnum
    return pdb_gpcr
