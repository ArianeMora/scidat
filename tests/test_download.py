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
        self.local = False
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        if self.local:

            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_') + '/'
        self.data_dir = os.path.join(THIS_DIR, 'data/')
        manifest_file = self.data_dir + 'manifest.tsv'
        gdc_client = self.data_dir + './gdc-client'
        self.download = Download(manifest_file, self.tmp_dir, self.tmp_dir, gdc_client, max_cnt=1)

    def tearDown(self):
        # Delete temp dir
        shutil.rmtree(self.tmp_dir)

    def test_download(self):
        self.download.download()

        # Now there should be an extra file in the downloads dir
        files_post = os.listdir(self.tmp_dir)

        # Check file name
        self.assertEqual('001ae925-102c-4818-8eb0-c8d2e5726e7c' in files_post, True)

    def test_check_downloads(self):
        # Save us from having to download them again
        if self.local:
            os.system(f'cp -r {self.data_dir}001ae925-102c-4818-8eb0-c8d2e5726e7c {self.tmp_dir}')
            os.system(f'cp -r {self.data_dir}19601351-3c26-4293-b87d-97222cd64a19 {self.tmp_dir}')
        else:
            self.download.download()

        self.download.copy_downloads_to_new_dir(self.tmp_dir)
        # Run the download check
        download_status = self.download.check_downloads(self.data_dir + 'download_status.csv')
        download_status.sort()

        # Check the download status was correctly assigned
        self.assertEqual(download_status[0][0], '001ae925-102c-4818-8eb0-c8d2e5726e7c')
        self.assertEqual(download_status[0][5], 'True')
        self.assertEqual(download_status[-1][0], '19601351-3c26-4293-b87d-97222cd64a19')
        self.assertEqual(download_status[-1][5], 'True')

        # Check the file was written with the download status
        self.assertEqual(os.path.exists(self.data_dir + 'download_status.csv'), True)

    def test_copy_downloads_to_new_dir(self):
        # Copy the files to the tmp dir
        if self.local:
            os.system(f'cp -r {self.data_dir}001ae925-102c-4818-8eb0-c8d2e5726e7c {self.tmp_dir}')
            os.system(f'cp -r {self.data_dir}19601351-3c26-4293-b87d-97222cd64a19 {self.tmp_dir}')
        else:
            self.download.download()

        files_pre = os.listdir(self.tmp_dir)

        self.download.copy_downloads_to_new_dir(self.tmp_dir)

        files_post = os.listdir(self.tmp_dir)

        self.assertEqual('jhu-usc.edu_KIRC.HumanMethylation450.6.lvl-3.TCGA-CZ-5989-01A-11D-1670-05.gdc_hg38.txt' in
                         files_post, True)
        self.assertEqual('d3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv' in files_post, True)
        self.assertEqual('d3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv' in files_pre, False)

    def test_download_data_using_api(self):
        # Here we want to download the mutation files
        case_ids = ['C3N-00310', 'TCGA-KN-8422', 'TCGA-CZ-5989-11A', 'TCGA-A4-8312-01A']
        self.download.download_data_using_api(case_ids, 'mutation')

        labels = []
        for c in case_ids:
            labels.append(self.download.u.generate_label(['mutation', c], '.tsv'))
        # Let's check if the files were downloaded
        files = os.listdir(self.tmp_dir)
        overlap = 0
        for f in files:
            if f in labels:
                overlap += 1
        # There should only be two files with mutation data
        self.assertEqual(overlap, len(labels) - 2)

