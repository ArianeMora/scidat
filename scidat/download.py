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

from subprocess import Popen


class SciException(Exception):

    def __init__(self, message=''):
        Exception.__init__(self, message)


class DownloadException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class Download:

    def __init__(self, manifest_file, manifest_dir, download_dir, gdc_client, max_cnt=100):
        self.download_dir = download_dir
        self.manifest_dir = manifest_dir
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
    def download(cmds: list) -> None:
        processes = [Popen(cmd, shell=True) for cmd in cmds]
        for p in processes: p.wait()

    def run(self) -> None:
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
                        new_manifest = f'{self.manifest_dir}{file_cnt}_{manifest_filename}'
                        self.write_file(new_manifest, to_write)
                        cmds.append(f'{self.gdc_client} download -m {new_manifest} -d {self.download_dir}')
                        file_cnt += 1
                        cnt = 0
                        to_write = [hdr]
        self.download(cmds)