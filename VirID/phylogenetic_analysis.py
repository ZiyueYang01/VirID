import logging
import os
import pandas as pd
from tqdm import tqdm
from shutil import copy,copytree,rmtree
import multiprocessing as mp
import numpy as np
from Bio import SeqIO
from os import system
from collections import defaultdict
from VirID.config.config import RdRP_DB_TABLE_PATH
from VirID.phytree.classify import Classify
from VirID.config.config import (AA_LENGTH_DICT,RdRP_ICTV_NR,Pathogenicity_VIRUS_CSV,INTER_RdRP_CLUSTR_PATH,Pathogenicity_VIRUS_GROUP,NR_DB_PATH,
                                    RdRP_DB_TABLE_PATH,RdRP_DB_PATH)
from VirID.external.fasta_process import Seqkit
from VirID.external.align_trim import Mafft,Trimal
from VirID.phytree.translate import Translate
from VirID.external.blast import Diamond
from VirID.biolib.common import make_sure_path_exists
from VirID.rvm.main_taxon import main_taxon

class Phylogenetic_analysis(object):
    def __init__(self,input_file,out_dir,options):
        self.input_file = input_file
        self.threads = options.threads
        self.translate_table = options.translate_table
        self.out_dir = out_dir
        self.out_pair_path = os.path.abspath(os.path.join(self.out_dir, ".."))
        self.group_list = options.group_list
        self.group_aa_length = options.group_aa_length
        self.logger = logging.getLogger('timestamp')


    def divide_clustr(self,tsv):
        clustr_dict = defaultdict(list)
        tsv_file = pd.read_csv(tsv,sep='\t',
                            names=["qseqid","qlen","RdRp_sseqid","RdRp_length","RdRp_pident","RdRp_matched_len","RdRp_evalue"])
        RdRp_db_tublr = pd.read_csv(RdRP_DB_TABLE_PATH, sep = ',' , usecols=['RdRp_sseqid','super_group'])
        super_group_outfile = pd.merge(tsv_file, RdRp_db_tublr, on = 'RdRp_sseqid',how = 'left')
        groups = super_group_outfile.groupby(super_group_outfile.super_group)
        clustr = super_group_outfile["super_group"].unique()

        for i in clustr:
            if i == "RdRp_sseqid":
                continue
            if len(self.group_list) > 0 :
                if (i in self.group_list) is False:
                    continue
            clustr_dir = os.path.join(self.out_dir,i)
            clustr_txt = f"{clustr_dir}/{i}.txt"
            make_sure_path_exists(clustr_dir)

            groups.get_group(i).to_csv(f"{clustr_dir}/{i}.tsv")

            with open(clustr_txt,'w+', newline='') as f:
                for a in groups.get_group(i)["qseqid"].tolist():
                    f.write(str(a)+'\n')
            clustr_dict[i] = clustr_txt
        return clustr_dict


    def Infectivity_aai(self,clustr_name,out_trimal_file):
        aai_test = os.path.join(self.out_dir,"ill_family_aai.csv")
        clustr_dir = os.path.join(self.out_dir,clustr_name)
        clustr_translate_file = f"{clustr_dir}/{clustr_name}_query_prot.fasta"
        query_list = [seq.id for seq in SeqIO.parse(clustr_translate_file, "fasta")]
        seq_trim_ID = [seq.id for seq in SeqIO.parse(out_trimal_file, "fasta")]
        seq_trim = SeqIO.to_dict(SeqIO.parse(out_trimal_file, "fasta"))

        q_ill = pd.DataFrame()

        Pathogenicity_table = pd.read_csv(Pathogenicity_VIRUS_CSV, sep = ',' , 
                                                 usecols=['Super-Group','Family','Taxon','ID','Human',
                                                          'Vertebrates','Plant'])
        Pathogenicity_superGroup = Pathogenicity_VIRUS_GROUP
        Pathogenicity_virus = Pathogenicity_table[["Taxon", "ID"]].groupby(Pathogenicity_table.Taxon)
        for q in [j for j in query_list if j in seq_trim_ID] :
            max_aai= 0
            max_aai_family = ""
            max_aai_one_identity = 0
            max_aai_id = ""
            query_str = str(seq_trim[q].seq)
            if (clustr_name in Pathogenicity_superGroup) :
                for family in Pathogenicity_superGroup[clustr_name]:
                    all_count = 0
                    max_ill_each = 0
                    max_ill_id =""
                    for ref_i in Pathogenicity_virus.get_group(family)["ID"].tolist():
                        if ref_i in seq_trim_ID:
                            ref_ill_str = str(seq_trim[ref_i].seq)
                            ill_count = 0
                            for i in range(len(query_str)):
                                if query_str[i] == ref_ill_str[i] and ref_ill_str[i]!="-":
                                    ill_count += 1
                            all_count+=ill_count/len(query_str)
                            if ill_count/len(query_str) > max_ill_each:
                                max_ill_each = ill_count/len(query_str)
                                max_ill_id = ref_i

                    if all_count/len(Pathogenicity_virus.get_group(family)["ID"].tolist()) > max_aai:
                        max_aai = all_count/len(Pathogenicity_virus.get_group(family)["ID"].tolist())
                        max_aai_family = family
                        max_aai_one_identity = max_ill_each
                        max_aai_id = max_ill_id
                q_ill = q_ill._append([[q,max_aai,max_aai_family,max_aai_one_identity,max_aai_id]])
        with open(aai_test,'a+', encoding='utf-8') as f:
            q_ill.to_csv(f, mode='a+',index=False,  header=False)


    def cal_aai(self,clustr_name,out_trimal_file):
        ref_fasta = f"{INTER_RdRP_CLUSTR_PATH}/{clustr_name}/{clustr_name}.fas"
        clustr_dir = os.path.join(self.out_dir,clustr_name)
        clustr_translate_file = f"{clustr_dir}/{clustr_name}_query_prot.fasta"
        query_list = [seq.id for seq in SeqIO.parse(clustr_translate_file, "fasta")]
        ref_list = [seq.id for seq in SeqIO.parse(ref_fasta, "fasta")]
        seq_trim_ID = [seq.id for seq in SeqIO.parse(out_trimal_file, "fasta")]
        seq_trim = SeqIO.to_dict(SeqIO.parse(out_trimal_file, "fasta"))
        del_list = []

        aai_list = pd.DataFrame()
        max_ref_aai=[]
        ref_trim_list = [j for j in ref_list if j in seq_trim_ID]

        for i in ref_trim_list:
            query_ref_str = str(seq_trim[i].seq)
            max_ref_aai_i=0
            for j in ref_trim_list:
                if i!=j:
                    aai_ref_str = str(seq_trim[j].seq)
                    count=0
                    for k in range(len(query_ref_str)):
                        if query_ref_str[k] == aai_ref_str[k] and aai_ref_str[k]!="-":
                            count+=1
                    max_ref_aai_i = count/len(query_ref_str) if count/len(query_ref_str) > max_ref_aai_i else max_ref_aai_i
            max_ref_aai.append(max_ref_aai_i)

        for q in [j for j in query_list if j in seq_trim_ID] :
            max = 0
            query_str = str(seq_trim[q].seq)
            for r in ref_list:
                ref_str = str(seq_trim[r].seq)
                count=0
                for i in range(len(query_str)):
                    if query_str[i] == ref_str[i] and ref_str[i]!="-":
                        count+=1
                max = count/len(query_str) if count/len(query_str) > max else max
            aai_list = aai_list._append([q,max])
            if (max < min(max_ref_aai)):
                del_list.append(q)

        aai_list.to_csv(os.path.join(self.out_dir,"aai.csv"), mode='a+',index=False,  header=False)
        return del_list
    

    def first_mafft(self,clustr_name,clustr_aa_length,clustr_query_file):
        clustr_dir = os.path.join(self.out_dir,clustr_name)
        clustr_all_fasta = f"{clustr_dir}/{clustr_name}_all.fasta"
        clustr_translate_file = f"{clustr_dir}/{clustr_name}_query_prot.fasta"
        ref_fasta = f"{INTER_RdRP_CLUSTR_PATH}/{clustr_name}/{clustr_name}.fas"

        translate_item = Translate(clustr_dir,True,clustr_aa_length,'build_tree',
                                   clustr_name=clustr_name,translate_table = self.translate_table)
        shortORF_ID,ORFlength_table = translate_item.run(clustr_query_file,clustr_translate_file)
        if os.path.getsize(clustr_translate_file) < 1:
            self.logger.warn(f'{clustr_name}:All contigs do not reach the threshold, end.')
            return 0
        if len(shortORF_ID) >0 :
            self.logger.warn(f'{clustr_name}:{len(shortORF_ID)} records contigs that do not reach the threshold')
        
        os.remove(clustr_translate_file+"_shortORF_exclude.txt")
        os.system(f'cat {clustr_translate_file} {ref_fasta} > {clustr_all_fasta}')
        mafft_item  = Mafft(self.threads)
        out_mafft_file = mafft_item.run(clustr_all_fasta)
        out_trimal_file =Trimal().run(out_mafft_file)
        return out_trimal_file


    def worker(self,clustr_name,clustr_txt):
        clustr_dir = os.path.join(self.out_dir,clustr_name)
        clustr_query_file = f"{clustr_dir}/{clustr_name}_query.fasta"
        clustr_all_fasta = f"{clustr_dir}/{clustr_name}_all.fasta"
        clustr_aa_length = AA_LENGTH_DICT[clustr_name]
        clustr_translate_file = f"{clustr_dir}/{clustr_name}_query_prot.fasta"

        if self.group_aa_length is not None:
            if clustr_name in self.group_aa_length:
                clustr_aa_length = int(self.group_aa_length[clustr_name])

        Seqkit().run(self.input_file,clustr_query_file,"grep",clustr_txt)
        os.remove(clustr_txt)
        out_trimal_file = self.first_mafft(clustr_name,clustr_aa_length,clustr_query_file)
        out_trimal_file_seq = [seq.seq for seq in SeqIO.parse(out_trimal_file, "fasta")]
        if len(out_trimal_file_seq[0]) < 50:
            out_trimal_file = self.first_mafft(clustr_name,300,clustr_query_file)

        query_list = [seq.id for seq in SeqIO.parse(clustr_translate_file, "fasta")]

        del_list = self.cal_aai(clustr_name,out_trimal_file)
        i=0

        seq_trim_ID = [seq.id for seq in SeqIO.parse(out_trimal_file, "fasta")]
        all_del_list = [j for j in query_list if j not in seq_trim_ID]

        while len(del_list) > 0:
            all_del_list = all_del_list + del_list
            if len(set(all_del_list)) ==  len(query_list):
                self.logger.warn(f'{clustr_name}:All contigs delete, end.')
                return 0
            self.logger.info(f'{clustr_name} delete {len(del_list)} seqs. ')
            with open(f"{clustr_dir}/{clustr_name}_del_{i}.txt", "w+") as opfile:
                    opfile.write("\n".join(all_del_list))
            clustr_all_fasta_del = f"{clustr_dir}/{clustr_name}_all_del_{i}.fasta"
            Seqkit().run(clustr_all_fasta,clustr_all_fasta_del,"grep-v",f"{clustr_dir}/{clustr_name}_del_{i}.txt")
            mafft_item  = Mafft(self.threads)
            out_mafft_file = mafft_item.run(clustr_all_fasta_del)
            out_trimal_file =Trimal().run(out_mafft_file)
            all_del_list = all_del_list + [j for j in query_list if j not in [seq.id for seq in SeqIO.parse(out_trimal_file, "fasta")]]
            del_list = self.cal_aai(clustr_name,out_trimal_file)
            i=i+1

        if i==0:
            out_trimal_file = f"{clustr_dir}/{clustr_name}_all.fasta_mafft.fasta.trim.fasta"
  
        else:
            out_trimal_file = f"{clustr_dir}/{clustr_name}_all_del_{i-1}.fasta_mafft.fasta.trim.fasta"
        
        self.Infectivity_aai(clustr_name,out_trimal_file)
        classify_item = Classify(out_trimal_file,clustr_dir,clustr_name)
        classify_item.run()


    def _novel_virus_identify(self,table_nvi):
        def getType(NR_matched_len,NR_pident,NR_length):
            if  ((NR_matched_len/NR_length >= 0.70) and NR_pident >=90) or NR_pident >=95:
                return "known"
            else:
                return "noval"

        def getPmatched(NR_matched_len,NR_length):
            return round((NR_matched_len/NR_length)*100,2)

        def getSpecies(qseqid,Virus_Type,NR_sseqid,Species):
            if Virus_Type == "noval" :
                return ""
            else:
                if NR_sseqid not in RdRP_ICTV_NR_table["NR_id"]:
                    if qseqid in NR_taxon.qseqid.values:
                        return NR_taxon[NR_taxon.qseqid == qseqid]["species"].values[0]
                    else:
                        return ""
                else:
                    return Species
                
        RdRP_ICTV_NR_table = pd.read_csv(RdRP_ICTV_NR,header=0,encoding = "utf-8",sep=',',
                        names =["Cluster_sseqid","ICTV_Realm","ICTV_Kingdom","ICTV_Phylum",
                                "ICTV_Subphylum","ICTV_Class","ICTV_Order",	"ICTV_Suborder","ICTV_Family",
                                "ICTV_Subfamily","ICTV_Genus","ICTV_Species","RdRP_id","Taxonid",
                                "NR_id","identity","NR_kindom","NR_phylum","NR_class","NR_order","NR_family","NR_genus","NR_species"] )
        NR_taxon = pd.read_csv(self.out_dir+"/step9_Taxon_VirusesFlitered.tsv",encoding = "utf-8",sep='\t')

        table_nvi['Virus_Type'] = table_nvi.apply(lambda x : getType(x['NR_matched_len'],x['NR_pident'],x['NR_length']),axis=1)
        table_nvi['pmatched'] = table_nvi.apply(lambda x : getPmatched(x['NR_matched_len'],x['NR_length']),axis=1)
        table_nvi['Species'] = table_nvi.apply(lambda x : getSpecies(x['qseqid'],x['Virus_Type'],x['NR_sseqid'],x['Species']),axis=1)
        return table_nvi



    def summary(self,input_csv):
        input_table = pd.read_csv(input_csv,header=None,encoding = "utf-8",
                                  names=['qseqid','Longest_aa_length','Virus_Taxon','Realm','Kingdom','Phylum',
                                           'Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species',
                                           'Human',"Vertebrates","Plant","RdRP_super_group"])

        seqkit = Seqkit()
        ID_txt = os.path.join(self.out_dir,"ID.txt")
        output_file_dir = os.path.join(self.out_pair_path,"results")
        with open(ID_txt,'w+', newline='') as f:
            for i in range(len(input_table)):
                f.write(str(input_table['qseqid'][i])+'\n')
        filename = os.path.split(input_csv)[1].split('.')[0]
        seqkit.run(self.input_file,f"{output_file_dir}/{filename}.fasta","grep",ID_txt)
        fieldnames = ["qseqid", "qlen", "sseqid", "stitle", "slen", "pident", "length", "evalue"]
        diamond = Diamond(NR_DB_PATH,self.threads,fieldnames,self.translate_table)
        if os.path.exists(os.path.join(self.out_pair_path,"assembly_and_basic_annotation")+"/step7_NR.tsv"):
            diamond_nr_out_tsv =os.path.join(self.out_pair_path,"assembly_and_basic_annotation")+"/step7_NR.tsv"
        else:
            diamond_nr_out_tsv = os.path.join(self.out_dir)+"/blast_nr.tsv"
            diamond.run(f"{output_file_dir}/{filename}.fasta",diamond_nr_out_tsv)
        nr_df = pd.read_csv(diamond_nr_out_tsv,sep='\t',names=["qseqid","qlen","NR_sseqid","NR_stitle","NR_length","NR_pident","NR_matched_len","NR_evalue"])
        nr_df.to_csv(os.path.join(self.out_dir)+"/blast_nr_name.tsv",sep = '\t',index = False)
        main_taxon(os.path.join(self.out_dir)+"/blast_nr_name.tsv",self.out_dir).run()

        if os.path.isfile(self.out_dir+"/ill_family_aai.csv"):
            ill_family_aai_table = pd.read_csv(self.out_dir+"/ill_family_aai.csv",header=None,encoding = "utf-8",sep=',',
                                    names =['qseqid','AAI','Pathogenicity_VIRUS_FAMILY','MAX_identity','MAX_ref_id'] )

            ill_family_aai_merge = pd.merge(left = input_table, right = ill_family_aai_table, on = 'qseqid',how = 'left')[['qseqid','AAI','Pathogenicity_VIRUS_FAMILY','MAX_identity','MAX_ref_id']]
        else:
            ill_family_aai_merge = pd.DataFrame()

        rdrp_res_table = pd.read_csv(self.out_dir+"/diamond_all.tsv",header=None,encoding = "utf-8",sep='\t',
                                 names = ["qseqid","qlen","RdRp_sseqid","RdRp_length","RdRp_pident","RdRp_matched_len","RdRp_evalue"] )
        
        nr_res_table = pd.read_csv(self.out_dir+"/blast_nr_name.tsv",header=0,encoding = "utf-8",sep='\t')
    
        RdRp_db_tublr = pd.read_csv(RdRP_DB_TABLE_PATH, sep = ',' , usecols=['RdRp_sseqid','super_group'])
        rdrp_table = pd.merge(left = rdrp_res_table, right = RdRp_db_tublr, on = 'RdRp_sseqid',how = 'left')[["qseqid","super_group"]]

        rdrp_nr_table = pd.merge(left = rdrp_table, right = nr_res_table, on = 'qseqid',how = 'right')

        input_table_merge = pd.merge(left = input_table, right = rdrp_nr_table, on = 'qseqid',how = 'left')

        nvi_nr_res_table = self._novel_virus_identify(input_table_merge)
        merge_table_illness = nvi_nr_res_table

        if not ill_family_aai_merge.empty:
            merge_table_illness_aai = pd.merge(merge_table_illness,ill_family_aai_merge,how='outer',on='qseqid')
            merge_table_illness_aai  = merge_table_illness_aai[['qseqid','qlen','Longest_aa_length','RdRP_super_group','Virus_Taxon',
                                                                'Realm','Kingdom','Phylum','Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species',
                                                                'Virus_Type','Human',"Vertebrates","Plant",'AAI','Pathogenicity_VIRUS_FAMILY','MAX_identity','MAX_ref_id']]
        else:
            merge_table_illness_aai = merge_table_illness[['qseqid','qlen','Longest_aa_length','RdRP_super_group','Virus_Taxon',
                                                           'Realm','Kingdom','Phylum','Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species',
                                                           'Virus_Type','Human',"Vertebrates","Plant"]]
        return merge_table_illness_aai


    def _compare(self):
        out_type = ["qseqid", "qlen", "sseqid", "slen", "pident", "length", "evalue"]

        if  os.path.exists(RdRP_DB_PATH+".dmnd") is False:
            Diamond(RdRP_DB_PATH,self.threads, 'a',self.translate_table).run(RdRP_DB_PATH+".fas",'a',model="makedb")
        diamond = Diamond(RdRP_DB_PATH,self.threads,out_type,self.translate_table)
        diamond_out_tsv = os.path.join(self.out_dir)+"/diamond_all.tsv"
        diamond.run(self.input_file,diamond_out_tsv)
        return diamond_out_tsv


    def draw_infectivity_distribution(self,label_csv,path):
        dir_path = os.path.dirname(__file__)
        system(f'Rscript {dir_path}/phytree/infectivity_distribution.r {label_csv} \
              {self.out_pair_path}/results/{path}')
        

    def run(self):
        self.logger.info('[phylogenetic_analysis] Classify known viruses contigs')
        if os.path.isfile(os.path.join(self.out_dir,"Contigs_classify_res_Fast.csv")):
            os.remove(os.path.join(self.out_dir,"Contigs_classify_res_Fast.csv"))
        if os.path.isfile(os.path.join(self.out_dir,"Contigs_classify_res_Accurate.csv")):
            os.remove(os.path.join(self.out_dir,"Contigs_classify_res_Accurate.csv"))
        if os.path.isfile(os.path.join(self.out_dir,"ill_family_aai.csv")):
            os.remove(os.path.join(self.out_dir,"ill_family_aai.csv"))

        diamond_out_tsv = self._compare()

        clustr_dict = self.divide_clustr(diamond_out_tsv)
        p_bar = tqdm(total=len(clustr_dict))
        p_bar.set_description('align')
        update = lambda *args: p_bar.update()
        pool = mp.Pool(len(clustr_dict))
        for clustr_name, clustr_txt in clustr_dict.items():
            self.logger.info(f'[phylogenetic_analysis] Translate-Align-Classify in clustr {clustr_name}')
            pool.apply_async(self.worker, (clustr_name,clustr_txt), callback=update)
        pool.close() 
        pool.join()

        res_path = os.path.join(self.out_pair_path,"results")
        res_Accurate_tsv = res_path+"/Contigs_classify_res_Accurate.csv"
        res_Fast_tsv = res_path+"/Contigs_classify_res_Fast.csv"


        if os.path.isfile(os.path.join(self.out_dir,"Contigs_classify_res_Fast.csv")):
            out = self.summary(os.path.join(self.out_dir,"Contigs_classify_res_Fast.csv"))
            with open(res_Fast_tsv,'w+', encoding='utf-8') as f:
                out.to_csv(f, index=False,  header=True)
        

        if os.path.isdir(os.path.join(res_path,"label_tree")):
            rmtree(os.path.join(res_path,"label_tree"))
        if os.path.isdir(os.path.join(self.out_dir,"label_tree")):
            copytree(os.path.join(self.out_dir,"label_tree"), os.path.join(res_path,"label_tree"))



