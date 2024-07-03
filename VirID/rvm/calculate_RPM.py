import logging
import os
from os import system
from VirID.external.reads_tool import Bowtie2,Samtools
from VirID.biolib.common import make_sure_path_exists
import pandas as pd
from shutil import copy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


class Calculate_RPM(object):
    def __init__(self,out_dir,threads,merged_table_tsv):
        """Instantiate the class."""
        self.out_dir = out_dir
        self.threads = threads
        self.merged_table_tsv = merged_table_tsv
        self.logger = logging.getLogger('timestamp')


    def calculate_mapping_counts(self,RPM_abundance_count_tsv):
        if os.path.exists(os.path.join(self.out_dir,"step4_rRNA_output.fq.1.gz")):
            norRNA_reads_list_1 = os.path.join(self.out_dir,"step4_rRNA_output.fq.1.gz")
            args = f'zcat {norRNA_reads_list_1} | grep -c "+"'
            norRNA_reads = int(os.popen(args).readlines()[0].replace("\n",""))*2
        else:
            norRNA_reads_list_1 = os.path.join(self.out_dir,"step4_rRNA_output.fq.gz")
            args = f'zcat {norRNA_reads_list_1} | grep -c "+"'
            norRNA_reads = int(os.popen(args).readlines()[0].replace("\n",""))

        merged_table = pd.read_csv(os.path.join(self.out_dir,self.merged_table_tsv), sep='\t',header=0)

        RPM_abundance_path = os.path.join(self.out_dir,"RPM_abundance")
        coverage_tabular = pd.read_csv(os.path.join(RPM_abundance_path,"coverage_tabular.tsv"), sep='\t', header=None,encoding = "utf-8",
                                        names=["qseqid","startpos","endpos","Numreads","Covbases","Coverage","Meandepth","meanbaseq","meanmapq"],
                                        usecols=['qseqid','Numreads','Covbases','Coverage','Meandepth'])

        coverage_table = pd.merge(left = merged_table, right = coverage_tabular, 
                              on = 'qseqid', how = "left")

        coverage_table.to_csv(os.path.join(self.out_dir,self.merged_table_tsv), sep = '\t' , index = False)
        pair = os.path.abspath(os.path.join(self.out_dir, ".."))
        res_path = os.path.join(pair,"results")
        make_sure_path_exists(res_path)
        res_tsv = res_path+"/"+self.merged_table_tsv
        copy(os.path.join(self.out_dir,self.merged_table_tsv),res_tsv)

        mapping_Table = pd.read_csv(RPM_abundance_count_tsv,sep='\t',encoding='utf-8',
                                    header=None, names=['qseqid','trimed_length','mapping_Counts','unmapping_Counts'])

        this_table = pd.merge(left = merged_table, right = mapping_Table, 
                              on = 'qseqid', how = "left", sort=True).sort_values(by = ['mapping_Counts'],ascending=False)
        sankey_path = res_path+"/RPM sankey.html"

        self.sankey(this_table,norRNA_reads,sankey_path)
        all_counts =  this_table['mapping_Counts'].sum()
        
        def getRPM(mapping_reads):
            return round(1000000 * (mapping_reads / norRNA_reads),2)

        mapping_counts = this_table.groupby(by = ['RdRP_super_group'])["mapping_Counts"].sum()

        mapping_counts = mapping_counts.to_frame()
        mapping_counts['RdRP_super_group'] = mapping_counts.index
        mapping_counts = mapping_counts.reset_index(drop = True)[['RdRP_super_group','mapping_Counts']]
        mapping_counts.loc[mapping_counts.shape[0]] = ['total_library',all_counts] 

        for i in range(mapping_counts.shape[0]):
            mapping_counts.loc[i,'Super_Group_RPM'] = getRPM(mapping_counts.loc[i,'mapping_Counts'])
        mapping_counts = mapping_counts.sort_values(by = ['Super_Group_RPM'],ascending=False)

        return mapping_counts 


    def sankey(self,this_table,norRNA_reads,path):
        tsv_input = this_table
        links_data = []
        tsv_input = tsv_input[['RdRP_super_group','phylum','class','order','family','genus','mapping_Counts']]
        list = ['RdRP_super_group','phylum','class','order','family','genus']

        for i in range(1,len(list)):
            tsv_input_sum = tsv_input.groupby([list[a] for a in range(i+1)],group_keys=False).sum().reset_index()
            for index,row in tsv_input_sum.iterrows():
                output = {"Node1": row[list[i-1]], "Node2": row[list[i]],  "Value": round(1000000 * (row['mapping_Counts'] / norRNA_reads),2),'group':i}
                links_data.append(output)
        
        pddata = pd.concat([pd.DataFrame({key: [item] for key, item in data.items()}) for data in links_data], ignore_index=True)
        pddata.to_csv(path+".csv",sep = ',',index = False)
        dir_path = os.path.dirname(__file__)
        system(f'Rscript {dir_path}/sankey.r {pddata} {path}')


    def RPM_auto_plot(self,RPM_table):
        RPM_abundance_path = os.path.join(self.out_dir,"RPM_abundance")
        data = RPM_table.sort_values(by=["Super_Group_RPM"],ascending=False).reset_index(drop=True)

        custom_params = {"axes.spines.right": False, "axes.spines.top": False}
        sns.set_theme(style = "ticks", rc = custom_params)
        fig , ax1 = plt.subplots(1,1,figsize=(8,5),sharex=True)
        RPM_plot = sns.barplot(data = data,x = 'Super_Group_RPM' , y = 'RdRP_super_group',ax=ax1, palette="muted", errorbar=('ci', 0),order=data['RdRP_super_group'])

        RPM_mean = data['Super_Group_RPM'].mean()//2.5
        for index,row in data.iterrows():
            RPM_plot.text(row.Super_Group_RPM+RPM_mean,row.name,row.Super_Group_RPM,ha="center",size=10)

        RPM_max = round(data['Super_Group_RPM'].max() * 1.05 / 1000) * 1000
        RPM_plot.set_xticks(np.linspace(0,RPM_max,10))
        RPM_plot.set_xticklabels([str(int(x/1000))+'k' for x in RPM_plot.get_xticks()])

        plt.xlabel("RPM")
        plt.ylabel("Super group")
        output_fig = RPM_plot.get_figure()
        RPM_abundance_count_png = os.path.join(RPM_abundance_path,"RPM_abundance_count.svg")
        output_fig.savefig(RPM_abundance_count_png,dpi=400,bbox_inches = 'tight')
        plt.close()
        pair = os.path.abspath(os.path.join(self.out_dir, ".."))
        res_path = os.path.join(pair,"results")
        if self.merged_table_tsv == "Primary_screen_keep_dup_res.tsv":
            res_svg = res_path+"/RPM_abundance_count_keep_dup.svg"
        else :
            res_svg = res_path+"/RPM_abundance_count.svg"
        copy(RPM_abundance_count_png,res_svg)


    def run(self,input_fasta,input_reads):
        self.logger.info('Calculate RPM')
        RPM_abundance_path = os.path.join(self.out_dir,"RPM_abundance")
        RPM_abundance_db_path = os.path.join(RPM_abundance_path,"db_bulid","candidate_viruses")
        make_sure_path_exists(RPM_abundance_path)
        make_sure_path_exists(os.path.join(RPM_abundance_path,"db_bulid"))
        RPM_abundance_sam = os.path.join(RPM_abundance_path,"RPM_abundance.sam")
        RPM_abundance_bam = os.path.join(RPM_abundance_path,"RPM_abundance.bam")
        RPM_abundance_sorted_bam = os.path.join(RPM_abundance_path,"RPM_abundance.sorted.bam")
        RPM_abundance_count_tsv = os.path.join(RPM_abundance_path,"RPM_abundance_count.tsv")
        coverage_tabular_tsv = os.path.join(RPM_abundance_path,"coverage_tabular.tsv")

        Bowtie2("build",self.threads).run(input_fasta,build_path=RPM_abundance_db_path)
        bowtie2_item = Bowtie2("local",self.threads)
        bowtie2_item.run(input_reads,out_file=RPM_abundance_sam,data_base=RPM_abundance_db_path)
        
        Samtools("view",self.threads).run(RPM_abundance_sam,RPM_abundance_bam)
        Samtools("sort",self.threads).run(RPM_abundance_bam,RPM_abundance_sorted_bam)
        Samtools("index",self.threads).run(RPM_abundance_sorted_bam)

        Samtools("coverage",self.threads).run(RPM_abundance_sorted_bam,coverage_tabular_tsv)
        Samtools("idxstats",self.threads).run(RPM_abundance_sorted_bam,RPM_abundance_count_tsv)

        RPM = self.calculate_mapping_counts(RPM_abundance_count_tsv)
        self.RPM_auto_plot(RPM)