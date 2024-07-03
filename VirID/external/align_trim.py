import logging
import os
import subprocess
from VirID.config.config import LOG_TASK



class Mafft(object):

    def __init__(self,threads):
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')
        self.threads = threads

    def run(self, input_file):
        env = os.environ.copy()
        # self.logger.log(LOG_TASK,'Using mafft to align query contigs and ref seqs..')
        msa_file = input_file+"_mafft.fasta"
        args = f' mafft --maxiterate 1000 --genafpair  --thread {self.threads} {input_file} > {msa_file}'
        mafft_log = input_file+"_mafft.log"
        with open(mafft_log, 'a+') as f_out_err:
            proc = subprocess.Popen(args, shell=True,stdout=f_out_err,
                                    stderr=f_out_err, env=env)
        proc.communicate()
        if proc.returncode != 0:
             self.logger.error(f'An error was encountered while running mafft in {input_file}.')
        return msa_file


class Trimal(object):

    def __init__(self):
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')

    def run(self, input_msa_file):

        env = os.environ.copy()
        output_msa_file = input_msa_file+".trim.fasta"
        args = ['trimal', '-in', input_msa_file, '-out', output_msa_file, '-automated1']
    
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,env=env)
        proc_out, proc_err = proc.communicate()
        if proc.returncode != 0:
             self.logger.error(
                f'An error was encountered while running trimal in {input_msa_file}.')

        
        return output_msa_file

