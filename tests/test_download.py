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

from scidat.download import Download

import os
import shutil
import tempfile
import unittest


class TestDownload(unittest.TestCase):

    def setUp(self):
        # Flag to set data to be local so we don't have to download them repeatedly. ToDo: Remove when publishing.
        self.local = True

        self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')
        self.data_dir = '../tests/data/'
        manifest_file = self.data_dir + 'manifest.tsv'
        gdc_client = self.data_dir + './gdc-client'
        if self.local:
            # manifest_file, split_manifest_dir, download_dir, gdc_client, max_cnt=100)
            self.download = Download(manifest_file, self.data_dir, self.data_dir, gdc_client, max_cnt=1)
        else:
            self.download = Download(manifest_file, self.tmp_dir, self.tmp_dir, gdc_client, max_cnt=1)
            self.data_dir = self.tmp_dir

    def tearDown(self):
        # Delete temp dir
        shutil.rmtree(self.tmp_dir)

    def test_download(self):
        files_pre = os.listdir(self.tmp_dir)
        self.download.download()

        # Now there should be an extra file in the downloads dir
        files_post = os.listdir(self.tmp_dir)
        self.assertEqual(len(files_pre), len(files_post) - 1)

        # Check file name
        self.assertEqual(files_post[0], '001ae925-102c-4818-8eb0-c8d2e5726e7c')

    def test_check_downloads(self):
        if not self.local:
            self.download.download()

        download_status = self.download.check_downloads(self.data_dir + 'download_status.csv')
        download_status.sort()

        # Check the download status was correctly assigned
        self.assertEqual(download_status[0][0], '001ae925-102c-4818-8eb0-c8d2e5726e7c')
        self.assertEqual(download_status[0][5], 'True')
        self.assertEqual(download_status[-1][0], '06b98f02-9853-4c2f-87fa-1e0a454ef4bc')
        self.assertEqual(download_status[-1][5], 'False')

        # Check the file was written with the download status
        self.assertEqual(os.path.exists(self.data_dir + 'download_status.csv'), True)

    def test_copy_downloads_to_new_dir(self):
        files_pre = os.listdir(self.data_dir)

        if not self.local:
            self.download.download()

        self.download.copy_downloads_to_new_dir(self.data_dir)

        files_post = os.listdir(self.data_dir)

        self.assertEqual(len(files_pre), len(files_post) - 2)

        assert 'jhu-usc.edu_KIRP.HumanMethylation450.15.lvl-3.TCGA-UZ-A9PS-05A-11D-A42K-05.gdc_hg38.txt' in files_post
        assert 'd3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv' in files_post
        assert 'd3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv' not in files_pre