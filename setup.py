from setuptools import setup, find_packages
import os

def package_files(dir_list):
    paths = []
    for dir in dir_list:
        for (path, directories, filenames) in os.walk(dir):
            for filename in filenames:
                filepath = os.path.join('..', path, filename)
                paths.append(filepath)
                full_path = os.path.join(path, filename)
    return paths



dir_path = os.path.dirname(__file__)
extra_files = package_files([dir_path+'/VirID/'])



setup(
    name='VirID',
    version='2.0.0',
    author="YangZiyue&ShanYongtao",
    author_email="yangzy58@mail2.sysu.edu.cn",
    packages=find_packages(),
    package_data={'': extra_files},

    install_requires=[
        "Bio",
        "biopython",
        "DendroPy",
        "matplotlib",
        "numpy",
        "pandas",
        "regex",
        "seaborn",
        "tqdm",
    ],
    entry_points={
        'console_scripts': [
            'VirID = VirID.__main__:main',
        ],
    },
)
