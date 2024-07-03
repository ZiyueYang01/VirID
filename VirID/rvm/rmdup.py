import os
from os import system
import logging
import pandas as pd
from VirID.external.blast import Blastn
from VirID.external.fasta_process import Seqkit
class Rmdup(object):
    def __init__(self,mode,out_dir,threads):
        """Instantiate the class."""
        self.mode = mode
        self.out_dir = out_dir
        self.threads = threads
        self.logger = logging.getLogger('timestamp')


    def run(self,input_file):
        dir_path = os.path.dirname(__file__)
        if self.mode!='phylogenetic_analysis':
            out_blastn_tsv = os.path.join(self.out_dir,"step11_rmdup_blastn_output.tsv")
            anicalc_tsv = os.path.join(self.out_dir,"step11_rmdup_anicalc.tsv")
            aniclust_ani80c40_tsv = os.path.join(self.out_dir,"step11_rmdup_aniclust_ani80c40.tsv")
            aniclust_id_ani80c40_txt = os.path.join(self.out_dir,"step11_rmdup_aniclust_id_ani80c40.txt")
            rmdup_res_fas = os.path.join(self.out_dir,"step11_rmdup_res.fas")
        else:
            out_blastn_tsv = os.path.join(self.out_dir,"rmdup_output.tsv")
            anicalc_tsv = os.path.join(self.out_dir,"rmdup_anicalc.tsv")
            aniclust_ani80c40_tsv = os.path.join(self.out_dir,"rmdup_aniclust_ani80c40.tsv")
            aniclust_id_ani80c40_txt = os.path.join(self.out_dir,"rmdup_aniclust_id_ani80c40.txt")
            rmdup_res_fas = os.path.join(self.out_dir,"rmdup_res.fas")

        blast_item = Blastn(self.threads)
        blast_item.run(input_file,"makedb",
                    db_path=os.path.join(self.out_dir,"rmdup","db"))
        blast_item.run(input_file,"blastn",
                    db_path=os.path.join(self.out_dir,"rmdup","db"),
                    out_tsv=out_blastn_tsv)


        system(f'python3 {dir_path}/anicalc.py -i {out_blastn_tsv} -o {anicalc_tsv}')
        system(f'python3 {dir_path}/aniclust.py --fna {input_file} --ani {anicalc_tsv} \
            --min_ani 80 --min_tcov 40 --out {aniclust_ani80c40_tsv}')
        system(f'cut -f1 {aniclust_ani80c40_tsv} > {aniclust_id_ani80c40_txt}')

        Seqkit().run(input_file,rmdup_res_fas,'grep',
                        ID_txt=aniclust_id_ani80c40_txt)
        return rmdup_res_fas

        

