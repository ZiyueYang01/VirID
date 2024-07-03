import logging
import os
import pandas as pd
from shutil import copy
from VirID.config.config import RdRP_DB_TABLE_PATH
from VirID.biolib.common import make_sure_path_exists
from VirID.rvm.histogram import clustr_pic

class Summary(object):
    def __init__(self,options,input_tsv_dir,primary_tsv):
        self.options = options
        """Instantiate the class."""
        self.input_tsv_dir = input_tsv_dir
        self.primary_tsv = primary_tsv
        self.logger = logging.getLogger('timestamp')

    def _stitle_NR_standard(self,merged_table):
        for i in range(merged_table.shape[0]):
            merged_table.loc[i,'NR_stitle'] = merged_table.loc[i,'NR_stitle'].replace(merged_table.loc[i,'NR_sseqid']+' ','')
            merged_table.loc[i,'NR_stitle'] = merged_table.loc[i,'NR_stitle'][::-1].replace('[ ','@',1).replace(']','',1)[::-1]

        merged_table[['protein','virusesType']] = merged_table['NR_stitle'].str.split('@',expand = True)
        merged_table = merged_table[merged_table.columns[:3].tolist() + ['protein','virusesType'] + merged_table.columns[4:-2].tolist()]
        return merged_table
        

    def _novel_virus_identify(self,table_nvi):
        def getType(NR_matched_len,NR_pident,NR_length):
            if ((NR_matched_len/NR_length >= 0.70) and NR_pident >=90) or NR_pident >=95:
                return "known"
            else:
                return "noval"

        def getPmatched(NR_matched_len,NR_length):
            return round((NR_matched_len/NR_length)*100,2)

        table_nvi['virus_type'] = table_nvi.apply(lambda x : getType(x['NR_matched_len'],x['NR_pident'],x['NR_length']),axis=1)
        table_nvi['pmatched'] = table_nvi.apply(lambda x : getPmatched(x['NR_matched_len'],x['NR_length']),axis=1)
        return table_nvi


    def run(self):
        self.logger.info('Summary results')
        RdRp_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step6_RdRp.tsv"), 
                                 sep='\t', encoding = "utf-8",header=None,
                                 names = ["qseqid","RdRp_qlen","RdRp_sseqid","RdRp_length","RdRp_pident","RdRp_matched_len","RdRp_evalue"] )
        Taxed_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step9_Taxon_VirusesFlitered.tsv"), sep='\t',header=0,
                                  names=["qseqid","NR_qlen","NR_sseqid","NR_stitle","NR_length","NR_pident",
                                  "NR_matched_len",	"NR_evalue","longest_aa_length","taxid",
                                  "kindom","phylum","class","order","family","genus","species"])
        
        Taxed_table['NR_Virus'] = [x.split('[')[1].split(']')[0] for x in Taxed_table['NR_stitle']]


        if self.options.no_trim_contamination is False:
            ORFLength_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step10_blastn_trimed_orf.tsv"), 
                                          sep='\t',header=0,usecols=["qseqid",'longest_aa_length'])
        else:
            ORFLength_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step8_shortORF_exclude.tsv"),
                                          sep='\t',header=0,usecols=["qseqid",'longest_aa_length'])

        RdRp_db_tublr = pd.read_csv(RdRP_DB_TABLE_PATH, sep = ',' , usecols=['RdRp_sseqid','super_group'])



        merge_rdrp_with_taxon_table = pd.merge(left = Taxed_table , right = RdRp_table, on = 'qseqid', how = "left",sort=True)
        stitle_NR_standard_table = self._stitle_NR_standard(merge_rdrp_with_taxon_table)
        RdRP_info_merge_table = pd.merge(left = stitle_NR_standard_table, right = RdRp_db_tublr , on = 'RdRp_sseqid' , how = 'left' , sort=True)


        if self.options.no_trim_contamination is False:
            Blastn_trimed_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step10_blastn_trimed.tsv"), sep='\t',header=0,
                                            names = ["qseqid","Trimed_qlen"])
            trim_add = pd.merge(left = RdRP_info_merge_table, right = Blastn_trimed_table, on = 'qseqid',how = 'right')
            orflength_add = pd.merge(left = trim_add, right = ORFLength_table, on = 'qseqid',how = 'left')
            nvi = self._novel_virus_identify(orflength_add) 
            final_output = nvi[['qseqid','NR_qlen','Trimed_qlen','longest_aa_length_y','NR_sseqid','protein','NR_Virus',
                            'super_group','kindom','phylum','class','order','family','genus','species','virus_type']].rename(columns={'longest_aa_length_y':'longest_aa_length'})
        else:    
            orflength_add = pd.merge(left = RdRP_info_merge_table, right = ORFLength_table, on = 'qseqid',how = 'left')
            nvi = self._novel_virus_identify(orflength_add) 
            final_output = nvi[['qseqid','NR_qlen','longest_aa_length','NR_sseqid','protein','NR_Virus','super_group','kindom','phylum','class','order','family','genus','species','virus_type']]


        output = final_output.rename(columns={'longest_aa_length': 'Longest_aa_length','NR_qlen':'qlen',
                                        'protein':'NR_protein','super_group':'RdRP_super_group',
                                        'virus_type':'Virus_type'})

        if self.primary_tsv == "Primary_screen_res.tsv":
            rmdup_table = pd.read_csv(os.path.join(self.input_tsv_dir,"step11_rmdup_aniclust_ani80c40.tsv"), 
                                sep='\t',names = ["qseqid","other_id"])[['qseqid']]
            output_res = pd.merge(left = rmdup_table, right = output, on = 'qseqid',how = 'left')
            final_tsv = os.path.join(self.input_tsv_dir,'Primary_screen_res.tsv')
        else:
            output_res = output
            final_tsv = os.path.join(self.input_tsv_dir,'Primary_screen_keep_dup_res.tsv')


        output_res.to_csv(final_tsv, sep = '\t' , index = False)

        pair = os.path.abspath(os.path.join(self.input_tsv_dir, ".."))
        res_path = os.path.join(pair,"results")


        if self.primary_tsv == "Primary_screen_res.tsv":
            clustr_pic(output_res,os.path.join(res_path,"Contigs_clustr.svg"))
        else:
            clustr_pic(output_res,os.path.join(res_path,"Contigs_clustr_keep_dup.svg"))
        

        return output_res


