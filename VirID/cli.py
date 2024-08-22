import argparse
from contextlib import contextmanager
from VirID.biolib.custom_help_formatter import CustomHelpFormatter

@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc,
                            formatter_class=CustomHelpFormatter)


@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive '
                                      f'{"required" if required else "optional"} '
                                      f'arguments')
    yield group.add_mutually_exclusive_group(required=required)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)

def __reads(group,required):
    group.add_argument(
        '-i', help="Single-end sequencing data only #i,paired-end with files in #i2")

def __reads_2(group):
    group.add_argument(
        '-i2', help="Paired-end sequencing data with files in #i")

def __classify_i(group,required):
    group.add_argument(
        '-classify_i', help="Contigs files to be classified")

def __out_dir(group, required):
    group.add_argument('-out_dir', type=str, default=None, required=required,
                       help="Directory to output files")

    
def __threads(group):
    group.add_argument('--threads', type=int,default=1,
                       help='Threads num')


def __translate_table(group):
    group.add_argument('--translate_table', type=int,default=1,
                       help='Genetic code')
    
def __group_list(group):
    group.add_argument('--group_list', nargs='*', default='',
                       help='Multiple specific Super-Group with commas spacing,eg: Astro-Poty,Hepe-Virga')

def __no_trim_contamination(group):
    group.add_argument('--no_trim_contamination',  action='store_true', default=False,
    help='It is not necessary to compare nt library to cut the sequence')

def __keep_dup(group):
    group.add_argument('--keep_dup',  action='store_true', default=False,
    help='Preserve redundant sequences in the evolutionary analysis module')
    
def __ultra_sensitive(group):
    group.add_argument('--ultra_sensitive',  action='store_true', default=False,
    help='Select ultra sensitive mode when comparing RdRP libraries')

def __group_aa_length(group):
    group.add_argument('--group_aa_length', type=str, nargs='+',default=None,
    help='Choose different length thresholds for each supercluster')

def __help(group):
    group.add_argument('-h', '--help', action="help", help="show help message")

def get_main_parser():
    main_parser = argparse.ArgumentParser(
        prog='VirID', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')
    with subparser(sub_parsers, 'end_to_end', '') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __reads(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __reads_2(grp)
            __threads(grp)
            __translate_table(grp)
            __no_trim_contamination(grp)
            __keep_dup(grp)
            __group_list(grp)
            __group_aa_length(grp)
            __ultra_sensitive(grp)
    with subparser(sub_parsers, 'assembly_and_basic_annotation', '') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __reads(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __reads_2(grp)
            __translate_table(grp)
            __threads(grp)
            __no_trim_contamination(grp)
            __ultra_sensitive(grp)
    with subparser(sub_parsers, 'phylogenetic_analysis', '') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __classify_i(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __threads(grp)
            __translate_table(grp)
            __keep_dup(grp)
            __group_list(grp)
            __group_aa_length(grp)

    return main_parser


