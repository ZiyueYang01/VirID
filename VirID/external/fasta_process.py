import logging
import os
import subprocess
from VirID.biolib.exceptions import ExternalException
class Seqkit(object):

    def __init__(self):
        self.logger = logging.getLogger('timestamp')

    def run(self, input_file, output_file, model, ID_txt={}):
        env = os.environ.copy()
        args = []
        if model == 'grep':
            args = ["seqkit", "grep", '-f', ID_txt, input_file, '-o', output_file]
        if model == 'grep-v':
            args = ["seqkit", "grep", '-f', ID_txt, input_file,'-v', '-o', output_file]
        elif model == 'rmdup':
            args = ["seqkit", "rmdup", '-s', '-o', output_file, input_file]
        elif model == 'fq2fa':
            args = ["seqkit", "fq2fa", '-o', output_file, input_file]

        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,env=env)
        proc.communicate()
        if proc.returncode != 0:
            raise ExternalException('An error was encountered while '
                                    'running seqkit, check the log')
        if not os.path.isfile(output_file):
            self.logger.error('seqkit returned a zero exit code but no output '
                              'file was generated.')


class Rmdup(object):
    def __init__(self):
        self.logger = logging.getLogger('timestamp')
    
    def run(self,input_file,out_file):
        env = os.environ.copy()

        if len(input_file) == 2:
            args = ['cd-hit-est', '-i', input_file[0], '-i2', input_file[1], 
                    '-o', out_file[0],'-o2', out_file[1]]
        elif len(input_file) == 1:
            args = ['cd-hit-est', '-i', input_file[0],'-o', out_file[0]]

        proc = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            raise ExternalException(f'Error deleting duplicate reads,{stderr}')
