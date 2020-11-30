###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import argparse
import sys
import os

from scidat import __version__, API
from sciutil import SciUtil


def print_help():
    lines = [f'scidat {__version__}',
             'Use --h for help information.',
             'usage: ',
             '--download ',
             '--mf [manifest file downloaded from TCGA ]',
             '--dd [directory to save the sub manifests and download files to]',
             '--gdc [directory with gdc_client file]',
             '--cf [clinical file]',
             '--sf [sample file]',
             '--af [annotation file]',
             '--o [output directory]',
             '--rna [t or f depending on if your manifest has RNAseq data]',
             '--meth [t or f depending on if your manifest has DNA Methylation data]',
             '--r [list of features (optional)]',
             '--c [list of clinical features (optional)]',
             '--maxp [num threads for download]',
             '--download [t or f whether to download the data from TCGA]',
             '--mutations [t or f whether to download mutation data for the samples as well]'
             ]
    print('\n'.join(lines))


def run(args):
    # This currently only supports the API not download and annotate separately
    api = API(args.mf, args.gdc, args.cf, args.sf, args.dd, args.md, args.af, args.r, args.c, args.maxp)
    # Check if they want to download the data
    if args.download == 't':
        api.u.warn_p(['Downloading data... WARN this may take some time.'])
        api.download_data_from_manifest()
    # Now we want to first build the annotation
    api.build_annotation()
    # Check if they want to also download mutation data
    if args.mutations == 't':
        api.download_mutation_data()
    # Build the rna df if we have RNA data
    if args.rna == 't':
        api.build_rna_df()
        df = api.get_rna_df()
        df = api.add_gene_metadata_to_df(df)
        # Save this to the output directory
        df.to_csv(f'{args.o}TCGA-scidat_RNA.csv', index=False)

    # Check if there were methylation data in the input
    if args.meth == 't':
        api.minify_meth_files(args.dd, args.o)
        api.build_meth_df(args.o)
        meth_df = api.get_meth_df()
        meth_df.to_csv(f'{args.o}TCGA-scidat_DNA-methylation.csv', index=False)

    if args.meth == 't' and args.rna == 't':
        merged_df = api.merge_rna_meth_values(meth_df, df, index_col='external_gene_name',
                                              merge_col='external_gene_name')
        merged_df.to_csv(f'{args.o}TCGA-scidat_DNA-methylation_merged_RNA.csv', index=False)

    api.u.dp(['FINISHED.'])


def gen_parser():
    parser = argparse.ArgumentParser(description='scidat')
    parser.add_argument('--md', type=str, help='Meta data directory (where all the required files are located e.g. '
                                               'the gdc file and the annotation files.)')
    parser.add_argument('--mf', type=str, help='File name of manifest file e.g.(gdc_manifest_20200416_055430.txt)')
    parser.add_argument('--gdc', type=str, help='File name of GDC client (for unix you need dir/./gdc-client')
    parser.add_argument('--cf', type=str, help='File name of clinical file (e.g. clinical.tsv)')
    parser.add_argument('--sf', type=str, help='File name of the sample file (e.g. gdc_sample_sheet.2020-04-16.tsv)')
    parser.add_argument('--af', type=str, help='Annotation file name (e.g. tcga_hsapiens_gene_ensembl-GRCh38.p13.csv).')
    parser.add_argument('--dd', type=str, help='Download directory')
    parser.add_argument('--o', type=str, help='Output directory')
    parser.add_argument('--rna', type=str, help='Either "t" or "f" for whether there are RNAseq datasets in '
                                                'the download manifest.')
    parser.add_argument('--meth', type=str, help='Either "t" or "f" for whether there are DNA Methylation datasets in '
                                                'the download manifest.')
    parser.add_argument('--r', type=str, default=None, help='Comma separated list of required features (defaults to None, '
                                              'format: "x,y,z")')
    parser.add_argument('--c', type=str, default=None, help='Comma separated list of clinical columns e.g "x,y,z"')
    parser.add_argument('--maxp', type=str, default=1, help='Maximum number of manifests to download in one file, the smaller '
                                                 'this number is the faster it will go, but it will start to use all '
                                                 'your processes see documentation for more info on how to set this param.')
    parser.add_argument('--download', type=str, default='f', help='Either "t" or "f" for whether to download or not.')
    parser.add_argument('--mutations', type=str, default='f', help='Either "t" or "f" for whether to also download mutation data.')

    return parser

def check_args(args):
    u = SciUtil()
    # Validate the input arguments.
    if not args.af or not os.path.isfile(args.af):
        u.err_p([f'The annotation file could not be located, file passed: {args.af}'])
        sys.exit(1)
    if not args.sf or not os.path.isfile(args.sf):
        u.err_p([f'The sample file could not be located, file passed: {args.sf}'])
        sys.exit(1)
    if not args.cf or not os.path.isfile(args.cf):
        u.err_p([f'The clinical file could not be located, file passed: {args.cf}'])
        sys.exit(1)
    if not args.gdc or not os.path.isfile(args.gdc):
        u.err_p([f'The GDC installer could not be located, file passed: {args.gdc}'])
        sys.exit(1)
    if not args.mf or not os.path.isfile(args.mf):
        u.err_p([f'The manifest file could not be located, file passed: {args.mf}'])
        sys.exit(1)
    if not args.md or not os.path.isdir(args.md):
        u.err_p([f'The meta data directory could not be located, passed: {args.md}'])
        sys.exit(1)
    if not args.dd or not os.path.isdir(args.dd):
        u.err_p([f'The download directory could not be located, passed: {args.dd}'])
        sys.exit(1)
    if not args.o or not os.path.isdir(args.o):
        u.err_p([f'The output directory could not be located, passed: {args.o}'])
        sys.exit(1)
    if args.download != 't' and args.download != 'f':
        u.err_p([f'The download argument passed is not supported: {args.download}, '
                 f'filetype must be "t" for downloading or "f" for not downloading the data.'])
        sys.exit(1)
    if args.download != 't' and args.download != 'f':
        u.err_p([f'The download argument passed is not supported: {args.download}, '
                 f'filetype must be "t" for downloading or "f" for not downloading the data.'])
        sys.exit(1)
    if args.meth != 't' and args.meth != 'f':
        u.err_p([f'The meth argument passed is not supported: {args.meth}, '
                 f'meth must be "t" for including methylation data or "f" for not.'])
        sys.exit(1)
    if args.rna != 't' and args.rna != 'f':
        u.err_p([f'The rna argument passed is not supported: {args.rna}, '
                 f'rna must be "t" for including RNAseq data or "f" for not.'])
        sys.exit(1)

def main(args=None):
    parser = gen_parser()
    u = SciUtil()
    if args:
        sys.argv = args
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'scidat v{__version__}')
        sys.exit(0)
    else:
        print(f'scidat v{__version__}')
        args = parser.parse_args(args)
        # Check the args
        check_args(args)
        # Otherwise we have need successful so we can run the program
        u.dp(['Running scidat...'
              '\nSaving to output directory: ', args.o,
              '\nInput data located in meta directory: ', args.md,
              '\nAnnotation file: ', args.af,
              '\nManifest file:', args.mf,
              '\nClinical file: ', args.cf,
              '\nSample file: ', args.sf,
              '\nGDC file:', args.gdc,
              '\nRNA data:', args.rna,
              '\nDNA Methylation data:', args.meth,
              '\nDownloading data:', args.download,
              '\nJoining mutation data:', args.mutations,
              '\nRunning download with num. threads:', args.maxp,
              '\nDownload directory:', args.dd,
              '\nRequired features:', args.r,
              '\nClinical features:', args.c
              ])
        # RUN!
        run(args)
    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()
    # ----------- Example below -----------------------

    # root_dir = '../tests/data/'
    # gdc_client = root_dir + './gdc-client'
    # clinical_file = root_dir + 'clinical.txt'
    # sample_file = root_dir + 'sample_sheet.txt'
    # manifest_file = root_dir + 'manifest.tsv'
    # annotation_file = root_dir + 'tcga_hsapiens_gene_ensembl-GRCh38.p13.csv'
    # main(['--md', root_dir,
    #       '--mf', manifest_file,
    #       '--gdc', gdc_client,
    #       '--sf', sample_file,
    #       '--af', annotation_file,
    #       '--cf', clinical_file,
    #       '--dd', './',
    #       '--o', './',
    #       '--rna', 't',
    #       '--meth', 't',
    #       '--download', 'f'])