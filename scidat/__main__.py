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

from scidat import __version__


def print_help():
    lines = [f'scidat {__version__}',
             'WARNING: scidat is currently not for CLI. See notebook examples for use.',
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

    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'scidat v{__version__}')
        sys.exit(0)

    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()