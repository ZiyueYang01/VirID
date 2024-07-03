import logging
import os
from .exceptions import FileNotFound, VirIDException
import sys

def check_file_exists(input_file):
    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exist: ' + input_file)
        raise FileNotFound('Input file does not exist: ' + input_file)
    return True

def make_sure_path_exists(path):
    if not path:
        return True
    elif os.path.isdir(path):
        return True
    try:
        os.makedirs(path)
        return True
    except OSError:
        logger = logging.getLogger('timestamp')
        logger.error('Specified path could not be created: ' + path)
        raise VirIDException('Specified path could not be created: ' + path)


def is_executable(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, _fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_executable(exe_file):
                return exe_file
    return None

def check_on_path(program, exit_on_fail=True):
    if which(program):
        return True

    if exit_on_fail:
        print('%s is not on the system path.' % program)
        sys.exit(1)
    return False


def check_dependencies(programs, exit_on_fail=True):
    for program in programs:
        if not check_on_path(program, exit_on_fail):
            return False
    return True
