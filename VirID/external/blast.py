import logging
import os
import subprocess
from VirID.biolib.exceptions import ExternalException
from VirID.config.config import LOG_TASK

class Blastn(object):
    def __init__(self,threads):
        self.threads = threads
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')
    
    def run(self,input_file,model,out_tsv="",db_path=""):
        env = os.environ.copy()
        if model == "makedb":
            args = ['makeblastdb', '-in', input_file, 
                '-dbtype', 'nucl', '-out', db_path]
        elif model == "nt":
            args = ['blastn', '-query', input_file, '-db', db_path, 
                    '-out', out_tsv,'-evalue', "1E-10", "-max_hsps", str(1),
                    "-num_threads", str(self.threads), "-max_target_seqs", str(1),
                    "-outfmt", "6 qacc qlen sseqid slen pident length qstart qend sstart send evalue bitscore qcovs"]
        else:
            args = ['blastn', '-query', input_file, '-db', db_path, 
                    '-out', out_tsv,'-evalue', "1E-10", "-max_hsps", str(1),
                    "-num_threads", str(self.threads), 
                    "-outfmt", "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"]

        blastn_log = out_tsv+"_log.txt"
        with open(blastn_log, 'a+') as f_out_err:
            proc = subprocess.Popen(
                args, stdout=f_out_err, stderr=f_out_err, env=env)
        proc.communicate()
        if proc.returncode != 0:
            self.logger.error(
                f'An error was encountered while running blastn, please check {blastn_log}')
            raise ExternalException('blastn returned a non-zero exit code.')
        return out_tsv


class Blastp(object):
    def __init__(self):
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')
    
    def run(self,input_file,model,out_tsv="",db_path=""):
        env = os.environ.copy()
        if model == "makedb":
            args = ['makeblastdb', '-in', input_file, 
                '-dbtype', 'prot', '-out', db_path]
        else:
            args = ['blastp', '-query', input_file, '-db', db_path, 
                    '-out', out_tsv,'-evalue', "1E-10", "-max_hsps", str(1),
                    "-outfmt", "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"]

        blastp_log = out_tsv+"_log.txt"
        with open(blastp_log, 'a+') as f_out_err:
            proc = subprocess.Popen(
                args, stdout=f_out_err, stderr=f_out_err, env=env)
        proc.communicate()
        if proc.returncode != 0:
            self.logger.error(
                f'An error was encountered while running blastp, please check {blastp_log}')
            raise ExternalException('blastp returned a non-zero exit code.')
        return out_tsv


class Diamond(object):

    def __init__(self, database_path,threads,out_type,translate_table):
        self.threads = threads
        self.database_path = database_path
        self.out_type = out_type
        self.translate_table = translate_table
        self.logger = logging.getLogger('timestamp')

    def run(self, origin_file, output_tsv,model=""):
        env = os.environ.copy()

        if model=='ultra_sensitive':
            args = ['diamond','blastx', '-q', origin_file, '-d', self.database_path, '-o', 
                 output_tsv, '-e', '1E-4', '--query-gencode',str(self.translate_table),'-k', str(1), '-p', str(self.threads),'--ultra-sensitive', '-f',str(6)]
        elif model == "makedb":
            args = ['diamond','makedb', '--in', origin_file, '-d', self.database_path]
        else:
            args = ['diamond','blastx', '-q', origin_file, '-d', self.database_path, '-o', 
                 output_tsv, '-e', '1E-4', '--query-gencode',str(self.translate_table), '-k', str(1), '-p', str(self.threads),'-f',str(6)]
        
        if model != "makedb":
            for a in self.out_type:
                args.append(a)
        
        diamond_log = output_tsv+"_log.txt"
        with open(diamond_log, 'a+') as f_out_err:
            proc = subprocess.Popen(
                args, stdout=f_out_err, stderr=f_out_err, env=env)
        proc.communicate()

        if proc.returncode != 0:
            self.logger.error(
                'An error was encountered while running Diamond.')
            raise ExternalException('Diamond returned a non-zero exit code.')
        if model != "makedb" and not os.path.isfile(output_tsv):
            self.logger.error(
                'An error was encountered while running Diamond.')
            raise ExternalException(
                'Diamond output file is missing: {}'.format(output_tsv))
