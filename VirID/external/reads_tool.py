import logging
import os
import subprocess
import shutil
from VirID.biolib.exceptions import ExternalException
from VirID.config.config import LOG_TASK


class Bowtie2(object):
    def __init__(self,model,threads):

        self.threads = threads 
        self.model = model
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')


    def run(self,input_file,out_file="",build_path="",data_base="",co1_file=""):
        env = os.environ.copy()
        if self.model == "fast-local":
            if len(input_file) == 2:
                    args = ['bowtie2', '--fast-local', '--threads', str(self.threads),
                        '--un-conc-gz', out_file,'-x', data_base, '-1', 
                        input_file[0], '-2', input_file[1]]
            elif len(input_file) == 1:
                    args = ['bowtie2', '--fast-local', '--threads', str(self.threads),
                        '--un-gz', out_file,'-x', data_base,'-U',input_file[0]]
        elif self.model == "build":
            args = ["bowtie2-build",input_file,build_path]
        elif self.model == "fast":
            if len(input_file) == 2:
                args = ['bowtie2', '--fast','--threads', str(self.threads),
                '-x',data_base,'-1', input_file[0], '-2', input_file[1],'--no-unal',
                '-S',out_file,'--al-conc',co1_file]
            elif len(input_file) == 1:
                args = ['bowtie2', '--fast','--threads', str(self.threads),
                '-x',data_base,'-U', input_file[0],'--no-unal',
                '-S',out_file,'--al',co1_file]
        elif self.model == "local":
            if len(input_file) == 2:
                args = ['bowtie2', '--local','--threads', str(self.threads),
                '-x',data_base,'-1', input_file[0], '-2', input_file[1],
                '-S',out_file]
            elif len(input_file) == 1:
                args = ['bowtie2', '--local','--threads', str(self.threads),
                '-x',data_base,'-U',input_file[0],'-S',out_file]
        proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise ExternalException('An error was encountered while running bowtie2.')



class Megahit(object):
    def __init__(self, threads):
        self.threads = threads 
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')

    def run(self,input_file,cut_len,out_dir):
        env = os.environ.copy()

        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir) 

        if len(input_file) == 2:
                args = ['megahit', '-1', input_file[0],'-2', input_file[1],
                        '--num-cpu-threads', str(self.threads), '--memory', str(1),
                        '--min-contig-len',str(cut_len),
                        '-o', out_dir]
        elif len(input_file) == 1:
                args = ['megahit', '-r', input_file[0],
                        '--num-cpu-threads', str(self.threads), '--memory', str(1),
                        '--min-contig-len',str(cut_len),
                        '-o', out_dir]

        proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)

        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise ExternalException('An error was encountered while '
                                    f'running megahit, please check {stderr}')
        return out_dir


class Samtools(object):
    def __init__(self, model,threads):
        self.model = model
        self.threads = threads 
        self.logger = logging.getLogger('timestamp')

    def run(self,input_file,output_file=""):
        env = os.environ.copy()
        # self.logger.info(f'Running Samtools  with model {self.model}')
        if self.model == "view":
            args = ["samtools", "view","-bSF4","-@",str(self.threads),"-o",output_file,input_file]
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)
        elif self.model == "sort":
            args = ["samtools","sort","-@",str(self.threads),"-o",output_file,input_file]
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)
        elif self.model == "index":
            args = ["samtools", "index","-@",str(self.threads),input_file]
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=env)
        elif self.model == "idxstats":
            args = f'samtools idxstats {input_file} > {output_file}'
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True, env=env)
        elif self.model == "coverage":
            args = f'samtools coverage -H  {input_file} -o  {output_file}'
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True, env=env)

        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise ExternalException('An error was encountered while '
                                    f'running Samtools, please check {stderr}')
