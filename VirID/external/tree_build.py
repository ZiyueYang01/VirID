import logging
import os
import subprocess
from Bio import SeqIO
from VirID.biolib.exceptions import ExternalException


class Pplacer(object):
    def __init__(self,classify_data_type):
        self.classify_data_type = classify_data_type
        self.logger = logging.getLogger('timestamp')

    def run(self, json_out, msa_file, ref_tree, pplacer_out):
        env = os.environ.copy()
        dir_path = os.path.dirname(__file__)
        package_path = os.path.join(dir_path,"package")

        if self.classify_data_type == 'nt':
            model = "GTR"
        else:
            model = "LG"
        args = ['pplacer', '-m', model, '-s', json_out,'-t', ref_tree, '-o', pplacer_out, msa_file]
        placer_log = pplacer_out+"_placer_log.txt"

        with open(placer_log, 'a+') as f_out_err:
            proc = subprocess.Popen(args, stdout=f_out_err,
                                    stderr=f_out_err, env=env)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise ExternalException('An error was encountered while '
                                    'running pplacer, check the log '
                                    'file: {}'.format(placer_log))
        if not os.path.isfile(json_out):
            self.logger.error('pplacer returned a zero exit code but no output '
                              'file was generated.')
            raise ExternalException
        return pplacer_out


    def tog(self, pplacer_json_out, tree_file):
        dir_path = os.path.dirname(__file__)
        package_path = os.path.join(dir_path,"package")
        """ Convert the pplacer json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            TogException: If a non-zero exit code is returned, or the tree file
                          isn't output.
        """

        args = ['guppy', 'tog', '-o', tree_file, pplacer_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise ExternalException(proc_err)

        if not os.path.isfile(pplacer_json_out):
            self.logger.error('tog returned a zero exit code but no output '
                              'file was generated.')
            raise ExternalException

        return tree_file
