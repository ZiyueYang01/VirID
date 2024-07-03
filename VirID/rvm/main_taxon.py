import os
import sys
import logging
import subprocess
import pandas as pd
from os import system
from VirID.config.config import PROT_ACC2TAXID
from VirID.config.config import TAXONKIT_DB
from VirID.biolib.common import make_sure_path_exists

class main_taxon(object):
    def __init__(self,input_tsv,out_dir):
        """Instantiate the class."""
        self.input_tsv = input_tsv
        self.out_dir = out_dir
        self.out_taxon_tmp = os.path.join(self.out_dir,"step9_Taxon_temp")
        self.logger = logging.getLogger('timestamp')


    def accession2taxid(self):
        Tax_step1_accession_ID_data = pd.read_csv(self.input_tsv,sep='\t')
        make_sure_path_exists(self.out_taxon_tmp)
        Tax_step1_accession_ID_data2 = Tax_step1_accession_ID_data.drop_duplicates(subset = 'qseqid', keep = 'first')["NR_sseqid"]
        Tax_step1_accession_ID_data2.to_csv(self.out_taxon_tmp + "/diamond_nr_accID.txt", index = False, sep = "\t", header = False)

        system(f"grep -F -f {self.out_taxon_tmp}/diamond_nr_accID.txt {PROT_ACC2TAXID} \
             > {self.out_taxon_tmp}/diamond_nr_acc2taxid.tsv")

        system(f"cut -f3 {self.out_taxon_tmp}/diamond_nr_acc2taxid.tsv \
                |   sort -u \
                |   taxonkit --data-dir {TAXONKIT_DB} lineage \
                |   awk '$2>0' > {self.out_taxon_tmp}/diamond_nr_taxonomy_temp.tsv")
        
        args = ["taxonkit", "--data-dir",f"{TAXONKIT_DB}","reformat", "-f" ,"{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" ,"-F",
                f"{self.out_taxon_tmp}/diamond_nr_taxonomy_temp.tsv","-o",f"{self.out_taxon_tmp}/diamond_nr_taxonomy_temp_2.tsv"]
        env = os.environ.copy()
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)
        proc.communicate()

        taxid_table = os.path.join(self.out_taxon_tmp,"diamond_nr_taxonomy.tsv")
        if os.path.isfile(taxid_table):
            os.remove(taxid_table)


        system(f'cut -f1,3-  {self.out_taxon_tmp}/diamond_nr_taxonomy_temp_2.tsv \
            > {self.out_taxon_tmp}/diamond_nr_taxonomy.tsv')
        


    def merge_taxon(self):

        blastp_raw_table = pd.read_csv(self.input_tsv, sep='\t')
        
        acc2taxid_table =  pd.read_csv(f"{self.out_taxon_tmp}/diamond_nr_acc2taxid.tsv", sep='\t',encoding = "utf-8",names=['sseqid','NR_sseqid','taxid','gi'])

        taxonomy_table =  pd.read_csv(f"{self.out_taxon_tmp}/diamond_nr_taxonomy.tsv", sep='\t',encoding = "utf-8" ,
                                        names=['taxid','kindom','phylum','class','order','family','genus','species'])
    
        blastp_raw2taxid = pd.merge(left = blastp_raw_table , right = acc2taxid_table, on = 'NR_sseqid', how = "left",sort=True)
        blastp_taxonomy = pd.merge(left = blastp_raw2taxid , right = taxonomy_table, on = 'taxid', how = "left",sort=True)

        blastp_taxonomy = blastp_taxonomy.drop(columns=["sseqid","gi"])

        output_path = os.path.join(self.out_dir, "step9_Taxon")
        blastp_taxonomy.to_csv(output_path+".tsv", sep = '\t' , index = False)

        blastp_taxonomy_VirusesFlitered = blastp_taxonomy[blastp_taxonomy['kindom']=='Viruses']
        blastp_taxonomy_VirusesFlitered=blastp_taxonomy_VirusesFlitered[~blastp_taxonomy_VirusesFlitered['order'].isin(["Ortervirales"])]  
        blastp_taxonomy_VirusesFlitered=blastp_taxonomy_VirusesFlitered[~blastp_taxonomy_VirusesFlitered['genus'].isin(["retro"])]  
        blastp_taxonomy_VirusesFlitered=blastp_taxonomy_VirusesFlitered[~blastp_taxonomy_VirusesFlitered['species'].isin(["retro"])]  
        blastp_taxonomy_VirusesFlitered.to_csv(output_path+"_VirusesFlitered.tsv" , sep = '\t' , index = False)

        blastp_taxonomy_VirusesFlitered_qseqid = blastp_taxonomy_VirusesFlitered['qseqid']
        blastp_taxonomy_VirusesFlitered_qseqid.to_csv(output_path+"_VirusesFlitered_ID_list.txt" ,index = False,header = False)

        return output_path+"_VirusesFlitered_ID_list.txt"



    def run(self):

        self.accession2taxid()
        return self.merge_taxon()
        
