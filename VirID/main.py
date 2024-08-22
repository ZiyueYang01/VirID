import logging
import os
import sys
import subprocess
from multiprocessing import Process
from tqdm import tqdm
from shutil import copy
from VirID.config.config import LOG_TASK
from VirID.rvm.trim_contamination import Trim_contamination
from VirID.rvm.summary import Summary
from VirID.rvm.rmdup import Rmdup
from VirID.rvm.calculate_RPM import Calculate_RPM
from .phylogenetic_analysis import Phylogenetic_analysis
from VirID.biolib.common import (check_dependencies,
                                  check_file_exists,
                                  make_sure_path_exists)
from VirID.biolib.exceptions import VirIDException
from .assembly_and_basic_annotation import assembly_and_basic_annotation

class OptionsParser(object):

    def __init__(self, version):
        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.version = version
        self.genomes_to_process = None


    def assembly_and_basic_annotation(self,options):

        self.logger.log(LOG_TASK,'START Primary Screen PART...')

        assembly_and_basic_annotation_path = os.path.join(options.out_dir,'assembly_and_basic_annotation')
        res_path = os.path.join(options.out_dir,"results")
        assembly_and_basic_annotation_fasta = f"{res_path}/Primary_screen_res.fasta"
        assembly_and_basic_annotation_fasta_keep_dup = f"{res_path}/Primary_screen_keep_dup_res.fasta"
        check_file_exists(options.i)
        make_sure_path_exists(assembly_and_basic_annotation_path)
        make_sure_path_exists(os.path.join(options.out_dir,'results'))

        if options.i2 is not None:
            check_file_exists(options.i2)
            reads_list = [options.i,options.i2]
        else:
            reads_list = [options.i]

        assembly_and_basic_annotation_item = assembly_and_basic_annotation(reads_list,assembly_and_basic_annotation_path,
                                                                           options.threads,options.translate_table)
        rm_rRNA_file,filename = assembly_and_basic_annotation_item.qc_with_rm_rRNA()

        accession_tax_VirusesFlitered_file = assembly_and_basic_annotation_item.run(options,rm_rRNA_file,filename)

        if len(reads_list) == 2:
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step1_QC_1.fq.gz"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step1_QC_2.fq.gz"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step2_QC_1.fq"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step2_QC_2.fq"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step3_QC_cdhit_1.fq"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step3_QC_cdhit_1.fq.clstr"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step3_QC_cdhit_1.fq2.clstr"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step3_QC_cdhit_2.fq"))
        else:
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step1_QC.fq.gz"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step2_QC.fq"))
            os.remove(os.path.join(assembly_and_basic_annotation_path,"step3_QC_cdhit.fq"))

        if not os.path.isfile(accession_tax_VirusesFlitered_file):
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()
        elif os.path.getsize(accession_tax_VirusesFlitered_file) < 1:
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()

        if options.no_trim_contamination is False:
            self.logger.info('[assembly_and_basic_annotation] Cut the sequence contamination at both ends of contigs')
            output_nt_fasta = Trim_contamination(assembly_and_basic_annotation_path,
                                                 options.threads,options.translate_table).run(accession_tax_VirusesFlitered_file)
            res_file = output_nt_fasta
        else:
            res_file = accession_tax_VirusesFlitered_file

        if not os.path.isfile(res_file):
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()
        elif os.path.getsize(res_file) < 1:
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()

        rmdup_res_fas = Rmdup(options.subparser_name,assembly_and_basic_annotation_path,options.threads).run(res_file)
        copy(os.path.join(assembly_and_basic_annotation_path,'step11_rmdup_res.fas'),assembly_and_basic_annotation_fasta)
        rmdup_res_fas_keep_dup = res_file
        copy(rmdup_res_fas_keep_dup,assembly_and_basic_annotation_fasta_keep_dup)
        
        Summary(options,assembly_and_basic_annotation_path,"Primary_screen_res.tsv").run()
        Summary(options,assembly_and_basic_annotation_path,"Primary_screen_keep_dup_res.tsv").run()

        Calculate_RPM(assembly_and_basic_annotation_path,
                    options.threads,"Primary_screen_res.tsv").run(rmdup_res_fas,rm_rRNA_file)
        Calculate_RPM(assembly_and_basic_annotation_path,
            options.threads,"Primary_screen_keep_dup_res.tsv").run(rmdup_res_fas_keep_dup,rm_rRNA_file)
        
        self.logger.log(LOG_TASK,'END Primary Screen PART...')

        return rmdup_res_fas



    def phylogenetic_analysis(self,options):
        self.logger.log(LOG_TASK,'START Contigs Classify  PART...')
        if os.path.getsize(options.classify_i) < 1:
            self.logger.warn(f'There is no sequence in the file.')
            return 0
        contigs_classify_path = os.path.join(options.out_dir,'phylogenetic_analysis')
        make_sure_path_exists(contigs_classify_path)
        make_sure_path_exists(os.path.join(options.out_dir,'results'))
        check_file_exists(options.classify_i)

        if options.keep_dup is False:
            rmdup_res_fas = Rmdup(options.subparser_name,
                        contigs_classify_path,options.threads).run(options.classify_i)
        else:
            rmdup_res_fas = options.classify_i

        classify_item = Phylogenetic_analysis(rmdup_res_fas,
                                         contigs_classify_path,options)
        classify_item.run()
        self.logger.log(LOG_TASK,'END Contigs Classify  PART...')


    def parse_options(self, options):
        if sys.version_info.major < 3:
            raise VirIDException('Python 2 is no longer supported.')

        if hasattr(options, 'threads') and options.threads < 1:
            self.logger.warning(
                'You cannot use less than 1 CPU, defaulting to 1.')
            options.cpus = 1

        if options.subparser_name == 'end_to_end':
            check_dependencies(['bbduk.sh',
                    'diamond','seqkit','bowtie2','taxonkit',
                    'megahit','makeblastdb','blastn',
                    'blastp','mafft','trimal','pplacer','guppy'])
            
            contigs_file = self.assembly_and_basic_annotation(options)
            options.classify_i = contigs_file
            self.phylogenetic_analysis(options)

        elif options.subparser_name == 'assembly_and_basic_annotation':
            check_dependencies(['bbduk.sh',
                    'diamond','seqkit','bowtie2','taxonkit',
                    'megahit','makeblastdb','blastn',
                    'blastp'])

            self.assembly_and_basic_annotation(options)
        elif options.subparser_name == 'phylogenetic_analysis':
            check_dependencies(['bbduk.sh',
                    'diamond','seqkit','makeblastdb','blastn','taxonkit',
                    'blastp','mafft','trimal','pplacer','guppy'])

            self.phylogenetic_analysis(options)

        elif options.subparser_name == 'install':
            parent_path = os.path.dirname(os.path.realpath(__file__))
            args = f'pip install -r {parent_path}/requirements.txt'
            subprocess.call(args)
        else:
            self.logger.error('Unknown VirID command: "' +
                              options.subparser_name + '"\n')
            sys.exit(1)

        return 0

