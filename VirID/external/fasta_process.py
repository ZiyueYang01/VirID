import logging
import os
import subprocess
from os import system
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
