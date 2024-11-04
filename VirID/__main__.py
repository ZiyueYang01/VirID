import logging
import sys
import traceback
from VirID import  __version__
from VirID.cli import get_main_parser
from VirID.biolib.logger import logger_setup
from VirID.biolib.exceptions import *
from VirID.main import OptionsParser


def print_help():
    print('''\

              ...::: VirID :::...

  An Integrated Platform for Comprehensive RNA Virus Characterization.

  VirID  <medthod>  [options]
  
  Methods:
  
    end_to_end                           -> Complete workflow to discover potential RNA virus sequences and perform further phylogenetic characterisation.
                                            ( end_to_end  = assembly_and_basic_annotation + phylogenetic_analysis )
    assembly_and_basic_annotation        -> Identification of potential putative RNA virus sequences.
    phylogenetic_analysis                -> Phylogenetic analysis and pathogenicity estimation on putative (or other input) sequences.
    
  Tools:
  
    install -> install python package
    download -> data download
    
  Detail: VirID <Method> -h for command specific help
    ''')

def main():
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"VirID: version {__version__} ")
        sys.exit(0)

    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()
        if ('group_aa_length' in vars(args).keys()) and (args.group_aa_length is not None):
            i = iter(args.group_aa_length)
            args.group_aa_length = dict(zip(i, i))

    logger_setup(args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None,
                 "VirID.log", "VirID", __version__, False,
                 hasattr(args, 'debug') and args.debug)
    logger = logging.getLogger('timestamp')

    try:
        gt_parser = OptionsParser(__version__)
        gt_parser.parse_options(args)
    except SystemExit:
        logger.error('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except VirIDExit as e:
        if len(str(e)) > 0:
            logger.error('{}'.format(e))
        logger.error('Controlled exit resulting from an unrecoverable error or warning.')
        sys.exit(1)
    except (VirIDException, ExternalException) as e:
        msg = 'Controlled exit resulting from an unrecoverable error or warning.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)
    except Exception as e:
        msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)


if __name__ == '__main__':
    a = main()
    
