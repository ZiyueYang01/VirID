import pandas as pd
import numpy as np
import argparse
import os
import logging
from Bio import SeqIO
from VirID.phytree.translate import Translate
from VirID.config.config import NT_DB_PATH
from VirID.external.blast import Blastn

class Trim_contamination(object):
    def __init__(self,out_dir,threads,translate_table):
        """Instantiate the class."""
        self.out_dir = out_dir
        self.threads = threads 
        self.translate_table = translate_table
        self.logger = logging.getLogger('timestamp')
    
    def run(self,input_file):
        blastn_item = Blastn(self.threads)
        blastn_table_tsv = os.path.join(self.out_dir,"step10_blastn.tsv")
        blastn_table_trimed_tsv = os.path.join(self.out_dir,"step10_blastn_trimed.tsv")

        blastn_item.run(input_file,'nt',out_tsv=blastn_table_tsv,db_path=NT_DB_PATH)


        blastn_table = pd.read_csv(blastn_table_tsv, sep='\t', 
                                   names = ["qacc","qlen","NT_sseqid","NT_slen","NT_pident","NT_length",
                                            "qstart","qend","NT_sstart","NT_send","NT_evalue","NT_bitscore",'NT_qcovs'],encoding = "utf-8")


        ID_list = blastn_table['qacc'].to_list()
        dict_fasta = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
        output_nt_fasta = os.path.join(self.out_dir,"step10_blastn_trimed.fasta")

        len_dic = {}
        translate = Translate(self.out_dir,True,200,'nt',translate_table = self.translate_table)
        origin_orf_table = {}
        
        with open(output_nt_fasta, 'w') as fw:
            for seq_id in dict_fasta.keys():
                if seq_id in ID_list:

                    line = blastn_table.loc[blastn_table['qacc'] == seq_id]
                    qstart = int(line["qstart"])
                    qend = int(line["qend"])
                    origin_qlen = len(str(dict_fasta[seq_id].seq))
                    now_start = qstart if qstart < qend else qend
                    now_end = qend if qstart < qend else qstart
                    left = now_start-1
                    right = origin_qlen-(now_start-1) - int(line["NT_length"])

                    if left > right and left > 600:
                        left_seq = str(dict_fasta[seq_id].seq)[:left]
                        orf_left = len(translate.six_frame_translate(left_seq))

                        if int(orf_left)< 200:
                            continue
                        else:
                            origin_orf_table[seq_id]= orf_left
                        fw.write(">"+seq_id + '\n')
                        fw.write(left_seq + '\n')
                        len_dic[seq_id] = len(left_seq)

                    elif right > left and right > 600:
                        right_seq = str(dict_fasta[seq_id].seq)[now_end:]
                        orf_right = len(translate.six_frame_translate(right_seq))
                        if int(orf_right) < 200:
                            continue
                        else:
                            origin_orf_table[seq_id] = orf_right
                        fw.write(">"+seq_id + '\n')
                        fw.write(right_seq + '\n')
                        len_dic[seq_id] = len(right_seq)
                    else:
                        continue
                else:
                    fw.write(">"+seq_id + '\n')
                    fw.write(str(dict_fasta[seq_id].seq) + '\n')
                    len_dic[seq_id] = len(str(dict_fasta[seq_id].seq))
                    origin_orf_table[seq_id] = len(translate.six_frame_translate(str(dict_fasta[seq_id].seq)))

        origin_orf_table = pd.DataFrame.from_dict(origin_orf_table, orient='index',columns=['longest_aa_length'])
        origin_orf_table = origin_orf_table.reset_index().rename(columns = {'index':'qseqid'})
        origin_orf_table.to_csv(os.path.join(self.out_dir,"step10_blastn_trimed_orf.tsv"),sep='\t',index=False)


        len_df = pd.DataFrame.from_dict(len_dic, orient='index',columns=['Trimed_qlen'])
        len_df = len_df.reset_index().rename(columns = {'index':'qacc'})
        len_df.to_csv(blastn_table_trimed_tsv,index=False,sep='\t')

        return output_nt_fasta


