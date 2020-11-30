from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'scidat/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='scidat',
      version=read_version(),
      description='Download-Annotate-TCGA: Facilitates the download of data and annotation with metadata from TCGA',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Ariane Mora',
      author_email='ariane.n.mora@gmail.com',
      url='https://github.com/ArianeMora/scidat',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/ArianeMora/scidat/issues",
          "Documentation": "https://github.com/ArianeMora/scidat",
          "Source Code": "https://github.com/ArianeMora/scidat",
      },
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='annotation',
      packages=['scidat'],
      entry_points={
          'console_scripts': [
              'scidat = scidat.__main__:main'
          ]
      },
      install_requires=['sciutil', 'pandas', 'numpy', 'scibiomart', 'sciviso'],
      python_requires='>=3.6',
      data_files=[("", ['LICENSE', 'tests/data/d3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv',
                        'tests/data/jhu-usc.edu_KIRC.HumanMethylation450.6.lvl-3.TCGA-CZ-5989-01A-11D-1670-05.gdc_hg38.txt',
                        'tests/data/clinical.txt', 'tests/data/sample_sheet.txt',
                        'tests/data/mutation_TCGA-KN-8422_20200410.tsv', 'tests/data/manifest.tsv',
                        'tests/data/gdc-client', 'tests/data/mutation_TCGA-A3-3308_20200410.tsv'])]
      )