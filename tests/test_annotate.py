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
        self.dir_tmp = tempfile.mkdtemp(prefix='scidatannotate_tmp_')
        self.meta_dir = '/Users/ariane/Documents/vae_paper/data/rcc/raw/meta_tcga_kirp-kirc/'
        self.processed_dir = '/Users/ariane/Documents/vae_paper/data/rcc/processed/tcga_kirp-kirc/'
        self.sci_dir = '/Users/ariane/Documents/vae_paper/data/rcc/processed/scigacrux/'

        self.clinical_file = self.meta_dir + 'clinical.cart.2020-03-11/clinical.tsv'
        self.sample_file = self.meta_dir + 'gdc_sample_sheet.2020-03-11.tsv'
        self.manifest_file = self.meta_dir + 'gdc_manifest_20200311_043417.txt'
        self.mutations_file = self.meta_dir + 'mutations_20200311.tsv'

        # (self, output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
        #                  sep='_', mutations_file=None, clin_cols=None
        self.annotator = Annotate(self.dir_tmp, self.clinical_file, self.sample_file, self.manifest_file,
                                  ['count', 'm450'])

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_annotate(self):
        self.annotator.build_annotation()
        # Now we want to check that the dataframe
