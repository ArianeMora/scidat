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

    @staticmethod
    def check_dir_trailing_slash(dir_str: str) -> str:
        if dir_str[-1] != '/':
            return dir_str + '/'
        return dir_str

    @staticmethod
    def write_file(filepath: str, lines: list) -> None:
        with open(filepath, 'w+') as f:
            for line in lines:
                f.write(line)

    @staticmethod
    def run_cmds(cmds: list) -> None:
        processes = [Popen(cmd, shell=True) for cmd in cmds]
        for p in processes: p.wait()

    def copy_downloads_to_new_dir(self, new_dir):
        return

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
                line = line.split('\t')
                if first:
                    download_status.append(line + ['Download Status'])
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
        manifest_filename = self.manifest_file.split('/')[-1]
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
