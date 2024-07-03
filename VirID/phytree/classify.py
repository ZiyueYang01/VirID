import os
from os import system
import logging
import pandas as pd
import numpy as np
from shutil import copy
from VirID.external.tree_build import Pplacer
from VirID.config.config import INTER_RdRP_CLUSTR_PATH
from VirID.phytree.tree_taxon import Tree_taxon
from VirID.biolib.common import make_sure_path_exists

class Classify(object):
    def __init__(self,classify_file,out_dir,clustr_name):
        self.classify_file = classify_file
        self.clustr_name = clustr_name
        self.out_dir = out_dir
        self.out_pair_path = os.path.abspath(os.path.join(self.out_dir, ".."))
        self.logger = logging.getLogger('timestamp')



    def draw_pic(self,lable_tree,label_csv,type):
        dir_path = os.path.dirname(__file__)
        make_sure_path_exists(f"{self.out_pair_path}/label_tree")
        system(f'Rscript {dir_path}/label_tree.r {lable_tree} {label_csv} \
             {self.clustr_name.replace("_RdRp","")} \
              {self.out_pair_path}/label_tree/{self.clustr_name.replace("_RdRp","")}_{type}_label_tree ')


    def _fast(self,ref_tsv):
        
        place_item = Pplacer('aa')
        ref_dir = os.path.join(INTER_RdRP_CLUSTR_PATH,self.clustr_name)
        json_out = os.path.join(ref_dir,self.clustr_name)+"_RAxML_info.txt"
        ref_tree = os.path.join(ref_dir,self.clustr_name)+"_RAxML_result.nwk"
        
        msa_file = self.classify_file

        pplacer_out = self.classify_file+"_pplacer.jplace" 
        pplacer_json_out = place_item.run(json_out,msa_file,ref_tree,pplacer_out)
        tree = place_item.tog(pplacer_json_out,self.classify_file+"_pplacer.nwk")


        query_taxon_fast_taxon,taxon_dic_csv = Tree_taxon(tree,ref_tsv,self.clustr_name.replace("_RdRp",""),
                                        'Taxon_set',self.classify_file+".fasta",self.out_dir).run()

        query_taxon_fast_taxon.columns=['qseqid','Virus_Taxon']
        query_taxon_fast_Human = Tree_taxon(tree,ref_tsv,'0','Human',self.classify_file+".fasta",self.out_dir).run()

        query_taxon_fast_Human.columns=['qseqid','Human']
        query_taxon_fast_Vertebrates = Tree_taxon(tree,ref_tsv,'0','Vertebrates',self.classify_file+".fasta",self.out_dir).run()

        query_taxon_fast_Vertebrates.columns=['qseqid','Vertebrates']
        query_taxon_fast_Plant = Tree_taxon(tree,ref_tsv,'0','Plant',self.classify_file+".fasta",self.out_dir).run()

        query_taxon_fast_Plant.columns=['qseqid','Plant']

        taxon_ill_merge_0 = pd.merge(left = query_taxon_fast_taxon , right = query_taxon_fast_Human, on = 'qseqid', how = "left",sort=True)
        taxon_ill_merge_1 = pd.merge(left = taxon_ill_merge_0 , right = query_taxon_fast_Vertebrates, on = 'qseqid', how = "left",sort=True)
        taxon_ill_merge_2 = pd.merge(left = taxon_ill_merge_1 , right = query_taxon_fast_Plant, on = 'qseqid', how = "left",sort=True)
        taxon_ill_merge = pd.merge(left = taxon_ill_merge_2 , right = taxon_dic_csv, on = 'qseqid', how = "left",sort=True)


        out_tsv = f"{self.out_dir}/{self.clustr_name}_query_Fast.csv"
        all_out_fast_tsv = f"{self.out_pair_path}/Contigs_classify_res_Fast.csv"
        copy(ref_tsv,out_tsv)
        with open(out_tsv,'a+', encoding='utf-8') as f:
            taxon_ill_merge.to_csv(f, mode='a+', index=False,  header=False)
        self.draw_pic(tree,out_tsv,"Fast")
        taxon_ill_merge.columns=['qseqid','Virus_Taxon','Human',"Vertebrates","Plant",'Realm','Kingdom','Phylum','Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species']
        aa_table = pd.read_csv(f"{self.out_dir}/{self.clustr_name}_query_prot.fasta_ORFlength.tsv",header=0,sep='\t',
                                encoding = "utf-8",names = ["qseqid","longest_aa_length"])
        query_taxon_aa_add = pd.merge(left = taxon_ill_merge, right = aa_table, on = 'qseqid',
                            how = 'left')[['qseqid','longest_aa_length','Virus_Taxon','Realm','Kingdom','Phylum',
                                           'Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species',
                                           'Human',"Vertebrates","Plant"]]
        output = query_taxon_aa_add.rename(columns={'longest_aa_length': 'Longest_aa_length'})
        output['RdRP_super_group'] = self.clustr_name.replace("_RdRp","")

        with open(all_out_fast_tsv,'a+', encoding='utf-8') as f:
            output.to_csv(f, mode='a+', index=False,  header=False)

    
    def run(self):
        self.logger.info('[phylogenetic_analysis] Contigs Taxon')
        ref_tsv = f"{INTER_RdRP_CLUSTR_PATH}/{self.clustr_name}/{self.clustr_name}.csv"
        self._fast(ref_tsv)


