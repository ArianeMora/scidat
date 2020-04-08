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

import os
import shutil
import tempfile
import unittest

from scidat import Annotate


class TestAnnotate(unittest.TestCase):

    def setUp(self):
        # Setup temp dir
        self.local = True
        if not self.local:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')
        else:
            self.tmp_dir = '../tests/data/tmp/'
        meta_dir = '../tests/data/'
        clinical_file = meta_dir + 'clinical.txt'
        sample_file = meta_dir + 'sample_sheet.txt'
        manifest_file = meta_dir + 'manifest.tsv'
        mutations_file = None

        # (self, output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
        #                  sep='_', mutations_file=None, clin_cols=None
        self.annotator = Annotate(self.tmp_dir, clinical_file, sample_file, manifest_file,
                                  ['count', 'm450'])

    def tearDown(self):
        if not self.local:
            # Delete temp dr
            shutil.rmtree(self.tmp_dir)

    def test_annotate(self):

        self.annotator.build_annotation()

        # Now we want to check that the dataframe
        self.annotator.save_annotation(self.tmp_dir, 'out')

        # Check file exists
        self.assertEqual(os.path.exists(self.annotator.u.generate_label([self.tmp_dir, 'out'], '.csv')), True)

        # Check our annotation has what we expect
        print(self.annotator.annotated_file_dict)
