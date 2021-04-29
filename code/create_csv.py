## this file is 'duplicated', the other one is in the folder 
## /Documentos/TFM/GPCR_variants/Code/Variants_general/create_csv.py
## and typically the one that is used is the other

import pandas as pd
import csv
from bioservices import UniProt
import requests
import urllib
from Bio.PDB.Polypeptide import *
import re
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
datapath=os.path.join(parentdir,"Data/studied_GPCR_vars")
resultspath=os.path.join(parentdir,"Results/studied_GPCR_vars")

def create_excel():
    with open(os.path.join(datapath,"proteins_and_families_raw.txt"),"r") as infile:
        with open(os.path.join(resultspath,'myprot_list.csv'), 'w') as outfile:    
            mywriter = csv.writer(outfile,delimiter=';')
            for line in infile:
                if len(line) - len(line.lstrip())> 8:
                    el_li=line.strip().split('","')
                    fam=el_li[2]
                    short_nm=el_li[10]
                    long_nm=el_li[11]
                    uniprot=el_li[15]
                    mywriter.writerow([fam,long_nm,short_nm,uniprot])
###########                
def obtain_struc_pdb(entry):
    struc=requests.get('http://gpcrdb.org/services/structure/protein/'+entry+'/representative/').json()
    if struc:
        if type(struc) is list:
            struccell_l=[]
            for struc_e in struc:
                struccell_l.append(struc_e["pdb_code"]+" ("+struc_e["state"]+")")
            struccell=", ".join(struccell_l)
        else:
            struccell=struc["pdb_code"]+" ("+struc["state"]+")"
        return struccell
    else:
        return False
                
def search_struc():   
    u = UniProt()             
    with open(os.path.join(resultspath,"myprot_list.csv"),"r") as infile:
        with open(os.path.join(resultspath,'myprot_list_struc.csv'), 'w') as outfile:   
            mywriter = csv.writer(outfile,delimiter=';') 
            myreader = csv.reader(infile, delimiter=';')
            for row in myreader:
               if row[0]=="Family":
                   mywriter.writerow(row+["Uniprot entry","Struc"])
                   continue
               struccell=""
               entry=""
               uprot=row[3]
               print("\n\n",row[2])
               if uprot:
                   data=u.quick_search("id:%s" % uprot)
                   if data:
                       entry=data[uprot]['Entry name'].lower()
                       struc_res=obtain_struc_pdb(entry)
                       if struc_res:
                           struccell=struc_res
                           print(struccell)
                       else:
                           temp=requests.get('http://gpcrdb.org/services/structure/template/'+entry).json()
                           if temp:
                               temp_res=obtain_struc_pdb(temp)
                               if temp_res:
                                   struccell="[Model]: "+temp_res
                                   print(struccell)
                               else:
                                   print("-----No struc for template")
                           else:
                               print("-----No template")
                   else:
                       print("-----Uprot ID not found")
               else:
                   print("-----No uprot ID")
            
               mywriter.writerow(row+[entry,struccell])

def uniprot_mapping(fromtype, totype, identifier):
    """Takes an identifier, and types of identifier 
    (to and from), and calls the UniProt mapping service
    Abbrebiations of Uniprot identifier types can be found here: https://www.uniprot.org/help/api_idmapping
    """
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
                'to':totype,
                'format':'tab',
                'query':identifier,
    }
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    #and grab the mapping
    response =  urllib.request.urlopen(url)
    #response.read() provides tab-delimited output of the mapping
    return (response.read())
    
def write_variant_row(myrow,entry, protT,mywriter):
    prot_var=protT.loc[entry]
    for index, var in prot_var.iterrows():
        newrow=myrow+[var["Ligandtype"],var["gene"], var["GPCRdb"], var["SequenceNumber"], var["Segment"], var["GPCRdbWT"], var["NMaa"], var["MutationType"], var["Allele Count"], var["Allele Frequency"], var["Allele Number"], var["Number of Homozygotes"], var["sift_word"], var["sift_score"], var["polyphen_word"], var["polyphen_score"], var["GProteinInteraction"], var["ArrestinInteraction"], var["ActivationPathway"], var["MicroSwitch"], var["SodiumPocket"], var["foldchangeMaxAbs"], var["diseaseId"], var["score"], var["diseaseName"], var["Type"], var["PTMsite"], var["LB_structure"], var["LB_fam"]]
        mywriter.writerow(newrow)                    

def get_gpcr_pos_info(entry):
    gpcr_pos_info={}
    pos_out=requests.get('http://gpcrdb.org/services/residues/extended/'+entry).json()
    for pos in pos_out:
        gpcr_pos_info[pos["sequence_number"]]={"protein_segment":pos["protein_segment"], "display_generic_number":pos["display_generic_number"]}
    return gpcr_pos_info
        

def check_aa_type_changed(fromAA,toAA):
    aa_groups={"A":"hydrophobic","C":"hydrophobic","F":"hydrophobic","I":"hydrophobic","L":"hydrophobic","M":"hydrophobic","V":"hydrophobic","W":"hydrophobic","Y":"hydrophobic",
     "F":"aromatic","H":"aromatic","W":"aromatic","Y":"aromatic",
     "S":"polar_uncharged","T":"polar_uncharged","N":"polar_uncharged","Q":"polar_uncharged",
     "P":"helix_breakers", "G":"helix_breakers",
     "D":"negative","E":"negative",
     "H":"positive","K":"positive","R":"positive"}
    MutationType="changed"
    if aa_groups[fromAA]==aa_groups[toAA]:
        MutationType="similar"
    return MutationType

                
                    
def find_variants():
    drugT = pd.read_csv(os.path.join(datapath,'AHauser_variant_data/GPCR_drug_targets.csv'),sep=";",index_col=0).dropna(how="all").fillna(value="")
    trialT = pd.read_csv(os.path.join(datapath,'AHauser_variant_data/GPCR_trial_targets.csv'),sep=";",index_col=0).dropna(how="all").fillna(value="")
    with open(os.path.join(resultspath,"myprot_list_struc.csv"),"r") as myprotfile:
        with open(os.path.join(resultspath,'myprot_list_strucNvar.csv'), 'w') as outfile:   
            myprotread = csv.reader(myprotfile, delimiter=';')
            mywriter = csv.writer(outfile,delimiter=';') 
            for myrow in myprotread:
                if myrow[0]=="Family":
                    newrow=myrow+['Ligandtype','gene', 'GPCRdb', 'SequenceNumber', 'Segment', 'GPCRdbWT', 'NMaa', 'MutationType', 'Allele Count', 'Allele Frequency', 'Allele Number', 'Number of Homozygotes', 'sift_word', 'sift_score', 'polyphen_word', 'polyphen_score', 'GProteinInteraction', 'ArrestinInteraction', 'ActivationPathway', 'MicroSwitch', 'SodiumPocket', 'foldchangeMaxAbs', 'diseaseId', 'score', 'diseaseName', 'Type', 'PTMsite', 'LB_structure', 'LB_fam']
                    mywriter.writerow(newrow)
                    continue
                entry=myrow[4]
                if not entry:
                    mywriter.writerow(myrow)
                    continue
                if entry in  drugT.index.unique():
                    write_variant_row(myrow,entry, drugT,mywriter)
                elif entry in  trialT.index.unique():
                    write_variant_row(myrow,entry, trialT,mywriter)
                else:
                    uprot_map=uniprot_mapping('ACC+ID', 'ENSEMBL_ID', myrow[3])
                    ens_id= uprot_map.decode().strip("\n").split("\t")[-1]
                    exac=requests.get("http://exac.hms.harvard.edu/rest/gene/variants_in_gene/"+ens_id).json()
                    gpcr_pos_info=get_gpcr_pos_info(entry)
                    found=False
                    for exac_var in exac:
                        if exac_var["category"] =="missense_variant" and exac_var["major_consequence"] == "missense_variant":
                            found=True
                            var_info=exac_var["HGVSp"]
                            mymatch=re.match("p\.([A-Za-z]*)(\d*)([A-Za-z]*)",var_info)
                            fromAA=three_to_one(mymatch.group(1).upper())
                            seqNum=mymatch.group(2)
                            toAA=three_to_one(mymatch.group(3).upper())
                            allele_count=exac_var["allele_count"]
                            allele_freq=exac_var["allele_freq"]
                            allele_num=exac_var["allele_num"]
                            hom_count=exac_var["hom_count"]
                            
                            Ligandtype="?"
                            gpcrNum=""
                            this_pos_info=gpcr_pos_info[int(seqNum)]
                            if this_pos_info["display_generic_number"]:
                                gpcrNum=this_pos_info["display_generic_number"]
                            Segment=this_pos_info["protein_segment"]
                            MutationType=check_aa_type_changed(fromAA,toAA)
                            sift_word="?"
                            sift_score="?"
                            polyphen_word="?"
                            polyphen_score="?"
                            GProteinInteraction="?"
                            ArrestinInteraction="?"
                            ActivationPathway="?"
                            MicroSwitch="?"
                            SodiumPocket="?"
                            foldchangeMaxAbs="?"
                            diseaseId="?"
                            score="?"
                            diseaseName="?"
                            Type="?"
                            PTMsite="?"
                            LB_structure="?"
                            LB_fam="?"
                            
                            mywriter.writerow(myrow + [Ligandtype, ens_id, gpcrNum, seqNum, Segment, fromAA, toAA, MutationType, allele_count, allele_freq, allele_num, hom_count, sift_word, sift_score, polyphen_word, polyphen_score, GProteinInteraction, ArrestinInteraction, ActivationPathway, MicroSwitch, SodiumPocket, foldchangeMaxAbs, diseaseId, score, diseaseName, Type, PTMsite, LB_structure, LB_fam])
                    if not found:
                        mywriter.writerow(myrow)
                        print("No variants found for %s" % entry)
                    
def download_strucs():
    with open(os.path.join(resultspath,"myprot_list_strucNvar.csv"),"r") as infile:
        myreader = csv.reader(infile, delimiter=';')
        last_name=False
        for row in myreader:
            if row[0]=="Family":
                continue
            name=row[1]
            if last_name != name:
                last_name=name
                struc_all=row[5]
                if not struc_all:
                    continue
                struc_li=struc_all.lstrip("[Model]: ").split(", ")
                for struc in struc_li:
                    if "Inactive" in struc:
                        struc_sel=struc.split()[0]
                    elif "Active" in struc:
                        struc_sel=struc.split()[0]
                    elif "Intermediate" in struc:
                        struc_sel=struc.split()[0]

def true_false_unknown(myvar,true_val):
    if myvar==true_val:
        myvar_res=True
    elif myvar=="?":
        myvar_res="?"
    else:
        myvar_res=False
    return myvar_res

def var_pos_in_pdb():
    '''Reads the file myprot_list_strucNvar.csv that contains the info of the variants detected in the GPCRs'''
    with open(os.path.join(resultspath,"myprot_list_strucNvar.csv"),"r") as infile:
        myreader = csv.reader(infile, delimiter=';')
        last_name=False
        var_info={}
        all_gpcrs=[]
        for row in myreader:
            if row[0]=="Family":
                continue
            name=row[1]
            if last_name != name:
                last_name=name
                sname=row[2]
                all_gpcrs.append(sname)
                var_info[sname]={}
                var_info[sname]["name"]=name
                struc_all=row[5]
                if not struc_all:
                    continue
                is_model=False
                if "[Model]" in struc_all:
                    is_model=True
                struc_li=struc_all.lstrip("[Model]: ").split(", ")
                for struc in struc_li:
                    if "Inactive" in struc:
                        struc_sel=struc.split()[0]
                    elif "Active" in struc:
                        struc_sel=struc.split()[0]
                    elif "Intermediate" in struc:
                        struc_sel=struc.split()[0]
                var_info[sname]["pdb"]=struc_sel
                var_info[sname]["is_model"]=is_model
                var_info[sname]["vars"]=[]
            this_var={"gpcrdb":row[8],"seqN":str(int(float(row[9]))),"var_aa":row[11] ,"wt_aa":row[12]}

            sift=row[18]
            poliP=row[20]
            sift_poloP_affected="?"
            if "deleterious" in sift or "damaging" in poliP:
                sift_poloP_affected=True
            elif "tolerated" in sift or "benign" in poliP:
                sift_poloP_affected=False
            this_var["sift_poloP"]=sift_poloP_affected
            
            DisGeNet_disease=row[30]
            DisGeNet_disease_val=True
            if DisGeNet_disease=="?":
                DisGeNet_disease_val="?"
            elif DisGeNet_disease=="":
                DisGeNet_disease_val=False
            this_var["DisGeNet_disease"]=DisGeNet_disease_val

            gprotInt=row[22]
            gprotInt_val=true_false_unknown(gprotInt,"putative")
            this_var["gprotInt"]=gprotInt_val

            arrestinInt=row[23]
            arrestinInt_val=true_false_unknown(arrestinInt,"putative")
            this_var["arrestinInt"]=arrestinInt_val

            activation=row[24]
            activation_val=true_false_unknown(activation,"yes")
            this_var["activation"]=activation_val

            microSwitch=row[25]
            microSwitch_val=true_false_unknown(microSwitch,"yes")
            this_var["microSwitch"]=microSwitch_val

            sodiumP=row[26]
            sodiumP_val=true_false_unknown(sodiumP,"yes")
            this_var["sodiumP"]=sodiumP_val

            ptm=row[32]
            ptm_val=true_false_unknown(ptm,"yes")
            this_var["ptm"]=ptm_val

            interact=row[33]
            interact_val=true_false_unknown(interact,"interacting")
            this_var["interact"]=interact_val

            var_info[sname]["vars"].append(this_var)

    return(var_info,all_gpcrs)

def gnum_all_to_one(num_type,gnum):
    '''It recives two stings, the 'num_type' and the  X.ZxY term. Depending on
    the num_type it returns XxY or XxZ'''
    (helix,posnum)=gnum.split(".") # separate X and ZxY
    (bwnum,dbnum)=posnum.split("x")# separate X and Y
    if num_type=="gpcrdb":
        gnum_ok=helix+"x"+dbnum # return XxY
    else:
        gnum_ok=helix+"."+bwnum # return XxZ
    return(gnum_ok)

def check_affected(pred_affect,pred_items,this_var,var_to_check):
    if this_var[var_to_check] and this_var[var_to_check]!="?":
        pred_items.append(var_to_check)
        pred_affect=True
    elif not this_var[var_to_check]:
        if pred_affect=="?":
            pred_affect=False
    return(pred_affect,pred_items)

    
def open_json(json_filepath):
    '''Opens and reads 'json_filepath' and parses it into a Python dictionary to be returned '''
    json_file=open(json_filepath)
    json_str = json_file.read()
    json_dict=pd.io.json.loads(json_str)# parse the json string into a Python dictionary
    return json_dict
    
def gpcr_pdb_dict_fom_json(json_filepath,gnum_type="gpcrdb"):
    '''Opens the json file and returns a dictionary, keys as the gnum 'helix x dbnum' and value as the items of the json file'''
    gpcr_pdbaa=open_json(json_filepath)
    gpcr_pdbaa_ok={}
    for gnum,v in gpcr_pdbaa.items():
        if "." in gnum or "x" in gnum:
            gnum=gnum_all_to_one(gnum_type,gnum)
        gpcr_pdbaa_ok[gnum]=v
    return gpcr_pdbaa_ok
    
def var_to_display(ex_gpcr="HRH2",define_pdb=False,gnum_type="gpcrdb",dictpath="/home/martalo/Documentos/TFM/GPCR_variants/Results/gpcr_to_pdb/"):#test if works
    '''Build a list of variants of the introduced GPCR in ex_gpcr, it also stores if the variant changes the aa type or not '''
    #ex_gpcr="CCR8"
    (var_info,all_gpcrs)=var_pos_in_pdb() # Reads the file myprot_list_strucNvar.csv that contains the info of the variants detected in the GPCRs
    
    if define_pdb:
        pdbname=define_pdb
    else:
        if "pdb" in var_info[ex_gpcr].keys():
            pdbname=var_info[ex_gpcr]["pdb"]
        else:
            print("No pdb found for %s"%ex_gpcr)
            return False
    json_filename=pdbname+"_conv.json"
    gpcr_pdbaa_ok=gpcr_pdb_dict_fom_json(dictpath+json_filename,gnum_type=gnum_type) # Opens the json file and returns a dictionary, keys as the gnum 'helix x dbnum' and value as the items of the json file


    #var_for_gpcr=[]
    var_for_gpcr_d={}
    for this_var in var_info[ex_gpcr]["vars"]:
        gnum=this_var["gpcrdb"]
        seqpos=this_var["seqN"]
        #print(seqpos,gnum1)
        if gnum:
            gnum_ok=gnum_all_to_one(gnum_type,gnum)
            if gnum_ok in gpcr_pdbaa_ok:
                pdb_pos=gpcr_pdbaa_ok[gnum_ok]
                pred_affect="?"
                pred_items=[]
                for var_to_check in ['gprotInt', 'arrestinInt', 'sodiumP', 'interact', 'microSwitch', 'DisGeNet_disease', 'activation', 'sift_poloP', 'ptm']:
                    (pred_affect,pred_items)=check_affected(pred_affect,pred_items,this_var,var_to_check)
                    
                #var_for_gpcr.append({"affected_pred":pred_affect, "affected_items_pred":pred_items , "pdb_pos":pdb_pos,"gpcr_num":gnum_ok})
            else:
                continue
        elif seqpos in gpcr_pdbaa_ok:
            gnum_ok=False
            pdb_pos=gpcr_pdbaa_ok[seqpos]
            pred_affect="?"
            pred_items=[]
            for var_to_check in ['gprotInt', 'arrestinInt', 'sodiumP', 'interact', 'microSwitch', 'DisGeNet_disease', 'activation', 'sift_poloP', 'ptm']:
                (pred_affect,pred_items)=check_affected(pred_affect,pred_items,this_var,var_to_check)
                
            #var_for_gpcr.append({"affected_pred":pred_affect, "affected_items_pred":pred_items , "pdb_pos":pdb_pos,"gpcr_num":gnum_ok})
        else:
            continue

        aa_type_changed=False
        if check_aa_type_changed(this_var['wt_aa'],this_var['var_aa'])=="changed":
            aa_type_changed=True
        if seqpos in var_for_gpcr_d:
            dict_pred_items=var_for_gpcr_d[seqpos]["affected_items_pred"]
            dict_pred_affect=var_for_gpcr_d[seqpos]["affected_pred"]
            dict_aa_type_changed=var_for_gpcr_d[seqpos]["aa_type_changed"]
            if pred_items != dict_pred_items:
                pred_items_def=list(set(pred_items)| set(dict_pred_items))
                var_for_gpcr_d[seqpos]["affected_items_pred"]=pred_items_def
            if pred_affect!= dict_pred_affect:
                if True in [pred_affect , dict_pred_affect]:
                    pred_affect_def=True
                elif False in [pred_affect , dict_pred_affect]:
                    pred_affect_def=False
                else:
                    pred_affect_def="?"
                var_for_gpcr_d[seqpos]["affected_pred"]=pred_affect_def
            if dict_aa_type_changed != aa_type_changed and aa_type_changed:
                var_for_gpcr_d[seqpos]["aa_type_changed"]=aa_type_changed
        else:
            var_for_gpcr_d[seqpos]={"affected_pred":pred_affect, "affected_items_pred":pred_items , "pdb_pos":pdb_pos,"gpcr_num":gnum_ok, "aa_type_changed":aa_type_changed}
            
    return list(var_for_gpcr_d.values())




#(var_info,all_gpcrs)=var_pos_in_pdb()
#out_name="var_info.json"
#with open("/home/martalo/Documentos/TFM/GPCR_variants/Results/gpcr_to_pdb/"+out_name, 'w') as outfile:
#    json.dump(var_info, outfile)
###########






                        
#create_excel()
#search_struc()
#find_variants()
#download_strucs()

