from Bio.Seq import Seq
from Bio import SeqIO
import regex as re
import os
import logging
import pandas as pd
import warnings
from VirID.external.blast import Blastp
from VirID.config.config import INTER_RdRP_CLUSTR_PATH
warnings.filterwarnings('ignore')


class Translate(object):
    def __init__(self, out_dir, outfile_aa,aa_length,mode,clustr_name="",translate_table=1):
        self.out_dir = out_dir
        self.outfile_aa = outfile_aa
        self.aa_length = aa_length
        self.mode = mode
        self.table = translate_table
        self.clustr_name = clustr_name
        self.logger = logging.getLogger('timestamp')
        
    def six_frame_translate(self, seq):
        amino_acids = []
        origin_seq=Seq(seq.strip())
        for seq in [origin_seq,origin_seq.reverse_complement()]:
            for frame in range(3):
                index = frame
                while index < len(str(seq)) - 6:
                    match = re.match('((ATG(?:\S{3})*?T(?:AG|AA|GA))|ATG(?:\S{3})*)', str(seq[index:]))
                    if match:
                        orf = Seq(match.group().strip())
                        prot = orf.translate(table=self.table)
                        if len(prot) > self.aa_length:
                            amino_acids.append(str(prot).replace('*',''))
                        index += len(orf)
                    else: index += 3
        if len(amino_acids)>0:
            if self.mode == "build_tree":
                res = self.blastp_res(amino_acids)
                if res!=[]:
                    return res
                else:
                    return []
            else:
                return max(amino_acids,key=len,default='')
        else:
            return []
    

    def blastp_res(self,amino_acids):
        aa_output_file = f'{self.out_dir}/aa.fas'
        out_tsv = f'{self.out_dir}/aa.tsv'
        db_path = f"{INTER_RdRP_CLUSTR_PATH}/{self.clustr_name}/db/{self.clustr_name}"
        db_fas = f"{INTER_RdRP_CLUSTR_PATH}/{self.clustr_name}/{self.clustr_name}.fas"
        with open(aa_output_file, 'w') as fw:
            count=1
            for seq in amino_acids:
                fw.write(f">aa_seq_{count}"+ '\n')
                fw.write(str(seq) + '\n')
                count+=1
        if os.path.isdir(f"{INTER_RdRP_CLUSTR_PATH}/{self.clustr_name}/db") is False:
            Blastp().run(db_fas,"makedb",
                        db_path=db_path)
        aa_blastp = Blastp().run(aa_output_file,"blastp",out_tsv,db_path)
        aa_table = pd.read_csv(aa_blastp,header=None,encoding = "utf-8",sep = '\t', 
                                names =['qaccver','saccver','pident','length',"mismatch","gapopen",
                                        'qstart','qend','sstart','send','evalue','bitscore','qlen','slen'])
        if os.path.getsize(aa_blastp):
            index = aa_table['pident'].idxmax()
            right_rdrp_aa = aa_table.loc[index,"qaccver"]
            all_aa =SeqIO.to_dict(SeqIO.parse(aa_output_file, "fasta"))
            return all_aa[right_rdrp_aa].seq
        else:
            return []


    def run(self, input_file, output_file):
        shortORF_file = output_file+"_shortORF_exclude.txt"
        ORFlength_tsv = output_file+"_ORFlength.tsv"
        ORFlength_table = pd.DataFrame()
        shortORF_ID = []
        with open(output_file, 'w') as fw:
            with open(shortORF_file, 'w') as fr:
                for seq_record in SeqIO.parse(input_file, "fasta"):
                    prot = self.six_frame_translate(str(seq_record.seq))
                    if prot:
                        fw.write(">"+seq_record.id + '\n')
                        line=pd.DataFrame({'qseqid':str(seq_record.id),
                                        'longest_aa_length':len(str(prot))},index=[0])
                        ORFlength_table=ORFlength_table._append(line,ignore_index=True)  

                        if self.outfile_aa is True:
                            fw.write(str(prot) + '\n')
                        else:
                            fw.write(str(seq_record.seq) + '\n')
                    else:
                        shortORF_ID.append(seq_record.id)
                        fr.write(seq_record.id + '\n')
        ORFlength_table.to_csv(ORFlength_tsv,sep = '\t' , index = False)
        return shortORF_ID,ORFlength_table