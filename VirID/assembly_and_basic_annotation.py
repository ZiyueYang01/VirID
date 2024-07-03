import logging
import os
from os import system
import pandas as pd
import sys
from VirID.config.config import RdRP_DB_PATH, NR_DB_PATH, rRNA_DB_PATH,BBDUK_PATH
from VirID.external.fasta_process import Rmdup, Seqkit
from VirID.external.reads_tool import Bowtie2, Megahit
from VirID.external.blast import Diamond 
from VirID.phytree.translate import Translate
from VirID.rvm.main_taxon import main_taxon


class assembly_and_basic_annotation(object):
    def __init__(self,reads,out_dir,threads):
        """Instantiate the class."""
        self.reads = reads
        self.out_dir = out_dir
        self.threads = threads 
        self.logger = logging.getLogger('timestamp')

    def _reads_qc_single(self,reads):
        system(f'bbduk.sh \
            in={reads} \
            out={self.out_dir}/step1_QC.fq.gz \
            maq=10 qtrim=r trimq=10 ftl=5 minlen=90 > /dev/null 2>&1 ')
        system(f'bbduk.sh \
            in={self.out_dir}/step1_QC.fq.gz  \
            out={self.out_dir}/step2_QC.fq \
            entropy=0.5 entropywindow=50 entropyk=5 > /dev/null 2>&1 ')

        rmdup_item  = Rmdup()
        rmdup_item.run([f"{self.out_dir}/step2_QC.fq"], [f"{self.out_dir}/step3_QC_cdhit.fq"])
        return [str(f"{self.out_dir}/step3_QC_cdhit.fq")]


    def _reads_qc_double(self,reads_1,reads_2):
        system(f'bbduk.sh \
                in1={reads_1} in2={reads_2} \
                out1={self.out_dir}/step1_QC_1.fq.gz \
                out2={self.out_dir}/step1_QC_2.fq.gz \
                maq=10 qtrim=r trimq=10 ftl=5 minlen=90 > /dev/null 2>&1 ')
        system(f'bbduk.sh \
                in1={self.out_dir}/step1_QC_1.fq.gz \
                in2={self.out_dir}/step1_QC_2.fq.gz \
                out1={self.out_dir}/step2_QC_1.fq \
                out2={self.out_dir}/step2_QC_2.fq \
                entropy=0.5 entropywindow=50 entropyk=5 > /dev/null 2>&1')
        rmdup_item  = Rmdup()
        out_file = [str(f"{self.out_dir}/step3_QC_cdhit_1.fq"),str(f"{self.out_dir}/step3_QC_cdhit_2.fq")]
        rmdup_item.run(
            [str(f"{self.out_dir}/step2_QC_1.fq"),str(f"{self.out_dir}/step2_QC_2.fq")],
            out_file)
        return out_file
    

    def _project_name(self,str1,str2):
        ans = ''
        for i in zip(*[str1,str2]):	
            if len(set(i)) == 1:	         
                ans += i[0]
            else:					
                break
        return ans


    def _diamond_item(self,database,out_type,input_file,out_fasta,out_tsv,model):
        if os.path.exists(RdRP_DB_PATH+".dmnd") is False:
            Diamond(RdRP_DB_PATH,self.threads, 'a').run(RdRP_DB_PATH+".fas",'a',model="makedb")
        diamond_item= Diamond(database,self.threads, out_type)
        diamond_item.run(input_file,out_tsv,model)
        output_tsv_list = out_tsv+"_ID.txt"
        system(f'cut -f1 {out_tsv} |sort -u > {output_tsv_list}')
        seqkit_item = Seqkit()
        seqkit_item.run(input_file,out_fasta, "grep",output_tsv_list)
        os.remove(output_tsv_list)


    def qc_with_rm_rRNA(self):
        self.logger.info('[assembly_and_basic_annotation] Quality control of sequencing data')
        if len(self.reads) == 2:
            reads_1 = self.reads[0]
            reads_2 = self.reads[1]
            reads_1_name = os.path.basename(self.reads[0]).split('.')[0]
            reads_2_name = os.path.basename(self.reads[1]).split('.')[0]
            filename = self._project_name(reads_1_name, reads_2_name)
            after_qc_file = self._reads_qc_double(reads_1, reads_2)
            rm_rRNA_file = [os.path.join(self.out_dir, "step4_rRNA_output.fq.1.gz"),
                            os.path.join(self.out_dir, "step4_rRNA_output.fq.2.gz")]

        elif len(self.reads) == 1:
            after_qc_file = self._reads_qc_single(self.reads[0])
            filename = os.path.basename(self.reads[0 ]).split('.')[0]
            rm_rRNA_file = [os.path.join(self.out_dir, "step4_rRNA_output.fq.gz")]

        self.logger.info('[assembly_and_basic_annotation] Remove rRNA')
        out_rRNA_file = os.path.join(self.out_dir, "step4_rRNA_output.fq.gz")
        bowtie2_item = Bowtie2("fast-local",self.threads)
        bowtie2_item.run(after_qc_file,out_file=out_rRNA_file,data_base=rRNA_DB_PATH)
        return [rm_rRNA_file,filename]


    def run(self,options,rm_rRNA_file,filename):
        self.logger.info('[assembly_and_basic_annotation] Use megahit to splice reads into contigs')
        megahit_item = Megahit(self.threads)
        out_dir =  os.path.join(self.out_dir, "step5_Megahit_output")
        megahit_out_dir = megahit_item.run(rm_rRNA_file,300,out_dir)
        megahit_out_fasta = os.path.join(self.out_dir,"step5_Megahit_output.fasta")
        system(f'sed "s/ /_/g" {megahit_out_dir}/final.contigs.fa | sed "s/k/"{filename}"_k/g"| sed "s/=/_/g" > {megahit_out_fasta}')

        self.logger.info(f'[assembly_and_basic_annotation] Running diamond blastx to compare {RdRP_DB_PATH}')
        rdrp_out_type = ["qseqid", "qlen", "sseqid", "slen", "pident", "length", "evalue"]
        output_rdrp_tsv = os.path.join(self.out_dir,"step6_RdRp.tsv")
        rdrp_out_fasta = os.path.join(self.out_dir,"step6_RdRp.fasta")
        if options.ultra_sensitive is False:
            model = ""
        else:
            model = "ultra_sensitive"
        self._diamond_item(RdRP_DB_PATH,rdrp_out_type,megahit_out_fasta,rdrp_out_fasta,output_rdrp_tsv,model)
        if os.path.getsize(rdrp_out_fasta) < 1:
            self.logger.warn(f'[assembly_and_basic_annotation] No contigs are filtered out.')
            sys.exit()

        self.logger.info(f'[assembly_and_basic_annotation] Running diamond blastx to compare {NR_DB_PATH}')
        nr_out_type = ["qseqid", "qlen", "sseqid", "stitle", "slen", "pident", "length", "evalue"]
        output_nr_tsv = os.path.join(self.out_dir,"step7_NR.tsv")
        nr_out_fasta = os.path.join(self.out_dir,"step7_NR.fasta")
        self._diamond_item(NR_DB_PATH,nr_out_type,rdrp_out_fasta,nr_out_fasta,output_nr_tsv,model="")
        if os.path.getsize(nr_out_fasta) < 1:
            self.logger.warn(f'[assembly_and_basic_annotation] No contigs are filtered out.')
            sys.exit()

        self.logger.info(f'[assembly_and_basic_annotation] Remove contigs that cannot be translated into longer amino acid contigs.')
        translate_item = Translate(os.path.join(self.out_dir),False,200,'nt')
        translate_filter_file = os.path.join(self.out_dir,"step8_shortORF_exclude.fasta")
        translate_filter_tsv = os.path.join(self.out_dir,"step8_shortORF_exclude.tsv")
        shortORF_ID,ORFlength_tsv = translate_item.run(nr_out_fasta,translate_filter_file)
        if os.path.getsize(translate_filter_file) < 1:
            self.logger.warn(f'[assembly_and_basic_annotation] No contigs are filtered out.')
            sys.exit()
        fieldnames = ["qseqid","NR_qlen","NR_sseqid","NR_stitle","NR_length","NR_pident","NR_matched_len","NR_evalue"]
        translate_filter_df = pd.read_csv(output_nr_tsv,sep='\t',names=fieldnames)
        translate_filter_df.drop(translate_filter_df[translate_filter_df["qseqid"].isin(shortORF_ID)].index,
                                    inplace=True)
        df = pd.merge(left = translate_filter_df , right = ORFlength_tsv, on = 'qseqid', how = "left",sort=True)
        df.to_csv(translate_filter_tsv,sep = '\t',index = False)

        self.logger.info(f'[assembly_and_basic_annotation] Contigs annotation')
        accession_tax_VirusesFlitered_txt =  main_taxon(translate_filter_tsv,self.out_dir).run()
        accession_tax_VirusesFlitered_file = os.path.join(self.out_dir,"step9_Taxon_filter.fasta")
        seqkit_item = Seqkit()
        seqkit_item.run(translate_filter_file,accession_tax_VirusesFlitered_file,"grep",accession_tax_VirusesFlitered_txt)
        if os.path.getsize(accession_tax_VirusesFlitered_file) < 1:
            self.logger.warn(f'[assembly_and_basic_annotation] No contigs are filtered out.')
            sys.exit()

        return accession_tax_VirusesFlitered_file






        







