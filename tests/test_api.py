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

from scidat.api import API

import os
import shutil
import tempfile
import unittest


class TestAPI(unittest.TestCase):
    def setUp(self):
        # Flag to set data to be local so we don't have to download them repeatedly. ToDo: Remove when publishing.
        self.local = True
        if self.local:
            self.tmp_dir = '../tests/data/tmp/'
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')

        # def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
        #                  max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_'):

        self.data_dir = '../tests/data/'
        manifest_file = self.data_dir + 'manifest.tsv'
        gdc_client = self.data_dir + './gdc-client'
        meta_dir = '../tests/data/'
        clinical_file = meta_dir + 'clinical.txt'
        sample_file = meta_dir + 'sample_sheet.txt'
        manifest_file = meta_dir + 'manifest.tsv'

        self.download = API(manifest_file, gdc_client, clinical_file, sample_file, download_dir=self.tmp_dir,
                            max_cnt=1)

    def tearDown(self):
        if not self.local:
            # Delete temp dir
            shutil.rmtree(self.tmp_dir)

