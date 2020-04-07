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

import shutil
import tempfile
import unittest

from scidat import Annotate


class TestAnnotate(unittest.TestCase):

    def setUp(self):
        # Setup temp dir
        self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')
        self.meta_dir = '/Users/ariane/Documents/code/tcga-format/tests/data/'
        self.clinical_file = self.meta_dir + 'clinical.txt'
        self.sample_file = self.meta_dir + 'sample_sheet.txt'
        self.manifest_file = self.meta_dir + 'manifest.tsv'
        self.mutations_file = None

        # (self, output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
        #                  sep='_', mutations_file=None, clin_cols=None
        self.annotator = Annotate(self.tmp_dir, self.clinical_file, self.sample_file, self.manifest_file,
                                  ['count', 'm450'])

    def tearDown(self):
        # Delete temp dr
        shutil.rmtree(self.tmp_dir)

    def test_annotate(self):
        self.annotator.build_annotation()
        # Now we want to check that the dataframe
        self.annotator.save_annotation(self.meta_dir, 'out')
