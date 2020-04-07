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
import gzip
import os
import shutil
from subprocess import Popen

from sciutil import *


class SciException(Exception):

    def __init__(self, message=''):
        Exception.__init__(self, message)


class DownloadException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class Download:

    def __init__(self, manifest_file, split_manifest_dir, download_dir, gdc_client, max_cnt=100, sciutil=None):
        self.u = SciUtil() if not sciutil else sciutil
        self.download_dir = download_dir
        self.split_manifest_dir = split_manifest_dir
        self.gdc_client = gdc_client
        self.manifest_file = manifest_file
        self.max_cnt = max_cnt
        self.downloaded_files = []

    def check_args(self) -> None:
        """
        Check that the directories supplied by the user are locations on their computer
        Returns
        -------

        """

        if not os.path.isdir(self.download_dir):
            raise DownloadException(f'ARGPARSE ERR: download directory was not a directory: {self.download_dir}')

        if not os.path.isfile(self.manifest_file):
            raise DownloadException(f'ARGPARSE ERR: manifest file does not exist: {self.manifest_file}')

        if not os.path.isfile(self.gdc_client):
            raise DownloadException(f'ARGPARSE ERR: gdc-client file does not exist: {self.gdc_client}')

        if not os.path.isdir(self.split_manifest_dir):
            raise DownloadException(f'ARGPARSE ERR: the directory you selected to save your '
                                    f'manifest files to was not a directory: {self.split_manifest_dir}')

        # Now check they all have trailing slash for directories
        self.download_dir = self.u.check_dir_format(self.download_dir)
        self.split_manifest_dir = self.u.check_dir_format(self.split_manifest_dir)

    @staticmethod
    def write_file(filepath: str, lines: list) -> None:
        with open(filepath, 'w+') as f:
            for line in lines:
                f.write(line)

    @staticmethod
    def run_cmds(cmds: list) -> None:
        processes = [Popen(cmd, shell=True) for cmd in cmds]
        for p in processes:
            p.wait()

    def copy_downloads_to_new_dir(self, processed_dir):
        files = os.listdir(self.download_dir)
        count = 0
        for f in files:
            if os.path.isdir(self.download_dir + f):
                try:
                    # Get the downloaded folders
                    folders = os.listdir(self.download_dir + f)
                    file_idx = 0
                    for df in folders:
                        # TCGA also downloads log files we don't want to copy those across
                        if 'log' not in df:
                            file_idx += 1
                            count += 1
                            dir_str = self.u.check_dir_format(self.download_dir + f)
                            if 'gz' in df:
                                with gzip.open(dir_str + df, 'rb') as f_in:
                                    with open(processed_dir + df[:-3] + '.tsv', 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                            else:
                                # If it is a methylation file it isn't zipped so we just copy it across
                                os.system('cp ' + dir_str + df + ' ' + processed_dir + df)

                except Exception as e:
                    self.u.warn_p(["Unable to process", f])

    def check_downloads(self, ouput_file=None) -> list:
        """
        Checks if all the files in the manifest were downloaded.
        Prints out message of number of successful downloads and failed downloads.
        Optionally writes this to a file.
        
        Returns: defaultdict -> keys:str file name, value:boolean success state of download
        -------

        Parameters
        ----------
        ouput_file : str -> the path of the output file.

        """
        downloaded_files = os.listdir(self.download_dir)
        download_status = []
        success = []
        fail = []

        with open(self.manifest_file, 'r+') as f:
            first = True
            for line in f:
                # Remove trailing new line character before splitting
                line = line[:-1].split('\t')
                if first:
                    download_status.append(line + ['download_status'])
                    first = False
                else:
                    file_id = line[0]
                    if file_id not in downloaded_files:
                        download_status.append(line + ['False'])
                        fail.append(line[1])
                    else:
                        download_status.append(line + ['True'])
                        success.append(line[1])
        if len(fail) == 0:
            self.dp(["Successfully downloaded all files: no.", len(success)])
        else:
            self.u.err_p(["\tSuccessfully downloaded: ", len(success), '\n', "Failed: ", len(fail), '\n\nFailed Files:\n' + '\n'.join(fail)])

        if ouput_file:
            with open(ouput_file, 'w+') as f:
                for d in download_status:
                    f.write('\t'.join(d) + '\n')

        # Remove the header
        return download_status[1:]

    def download(self) -> None:
        cmds = []
        manifest_filename = self.manifest_file.split(self.u.dir_sep)[-1]
        with open(self.manifest_file, 'r+') as f:
            hdr = None
            cnt = 0
            to_write = []
            file_cnt = 0

            for line in f:
                if hdr is None:
                    hdr = line
                    to_write.append(hdr)
                else:
                    if cnt < self.max_cnt:
                        to_write.append(line)
                        cnt += 1
                    else:
                        new_manifest = f'{self.split_manifest_dir}{file_cnt}_{manifest_filename}'
                        self.write_file(new_manifest, to_write)
                        cmds.append(f'{self.gdc_client} download -m {new_manifest} -d {self.download_dir}')
                        file_cnt += 1
                        cnt = 0
                        to_write = [hdr]
        self.run_cmds(cmds)
