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
import os
import sys

from scidat import __version__

 # = '/Users/ariane/Documents/data/tcga/kidney/'
 # = f'{base_dir}./gdc-client'
 # = f'{base_dir}downloads/'
# manifest = f'{base_dir}gdc_manifest_20200402_003828.txt'

def print_help():
    lines = [f'scidat {__version__}',
             'usage: ',
             '--download ',
             '      [path to manifest file downloaded from TCGA ]',
             '      [directory to save the sub manifests]',
             '      [directory with gdc_client file]',
             '      [directory to download the files to]',
             '--annotate ',
             '      [directory with TCGA downloaded data]',
             '      [list of columns to annotate each file with (or blank for the default)]']
    print('\n'.join(lines))


def main(args=None):
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('newick', type=str, help='path to the newick file')

    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'phylodm v{__version__}')
        sys.exit(0)
    else:
        print(f'phylodm v{__version__}')
        args = parser.parse_args(args)

        # Validate the input arguments.
        if not os.path.isfile(args.newick):
            print(f'The input file cannot be found: {args.newick}')
            sys.exit(1)
        if args.method not in {'pd', 'node'}:
            print('Invalid choice for method, select either: pd or node')
            sys.exit(1)
        out_dir = os.path.dirname(args.output)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()