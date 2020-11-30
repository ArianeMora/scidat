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
import pandas as pd
import numpy as np
import shutil
from subprocess import Popen

from sciutil import *


class DownloadException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class Download:
    """
    Class to enable the download of data from TCGA.
    """
    def __init__(self, manifest_file, split_manifest_dir, download_dir, gdc_client, max_cnt=100, sciutil=None):
        self.u = SciUtil() if not sciutil else sciutil
        self.download_dir = download_dir
        self.split_manifest_dir = split_manifest_dir
        self.gdc_client = gdc_client
        self.manifest_file = manifest_file
        self.max_cnt = 6000
        self.downloaded_files = []
        self.file_type = None
        self.check_args()

    def add_windows_slash(self, change_str):
        if change_str and '/' in change_str and (self.file_type == 'windows' or 'C:' in change_str):
            return change_str.replace('/', '\\')
        return change_str

    def check_args(self) -> None:
        """
        Check that the directories supplied by the user are locations on their computer
        Returns
        -------

        """
        if self.gdc_client and 'exe' in self.gdc_client:
            self.file_type = 'windows'
        self.download_dir = self.add_windows_slash(self.download_dir)
        self.manifest_file = self.add_windows_slash(self.manifest_file)
        self.split_manifest_dir = self.add_windows_slash(self.split_manifest_dir)
        self.gdc_client = self.add_windows_slash(self.gdc_client)

        if not os.path.isdir(self.download_dir):
            raise DownloadException(f'ARGPARSE ERR: download directory was not a directory: {self.download_dir}')

        if not os.path.isfile(self.manifest_file):
            raise DownloadException(f'ARGPARSE ERR: manifest file does not exist: {self.manifest_file}')

        gdc_client = self.gdc_client.replace('/./', '/')
        if not os.path.isfile(gdc_client):
            raise DownloadException(f'ARGPARSE ERR: gdc-client file does not exist: {self.gdc_client}')
        if not os.access(gdc_client, os.X_OK):
            raise DownloadException(f'ARGPARSE ERR: gdc-client file is not executable: {self.gdc_client}')
        if not os.path.isdir(self.split_manifest_dir):
            raise DownloadException(f'ARGPARSE ERR: the directory you selected to save your '
                                    f'manifest files to was not a directory: {self.split_manifest_dir}')

        # Now check they all have trailing slash for directories
        # self.download_dir = self.u.check_dir_format(self.download_dir)
        # self.split_manifest_dir = self.u.check_dir_format(self.split_manifest_dir)

    def download_data_using_api(self, case_ids: list, data_type: str) -> None:
        """
        Download data using the API from TCGA: refer to \
        https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/#simple-somatic-mutation-endpoint-examples

        Parameters
        ----------
        case_ids
        data_type

        Returns
        -------

        """
        cmds = []
        files = []

        for case_id in case_ids:
            if data_type == 'mutation':
                file_name = self.u.generate_label([self.download_dir, data_type, case_id], '.tsv')
                cmds.append(self.gen_mutation_api_str(
                    case_id, file_name))
                files.append(file_name)

        self.run_cmds(cmds)

        if data_type == 'mutation':
            self.format_mutation_files(files, case_ids)

    def add_mutation_values(self, label: str, df: pd.DataFrame):
        """
        Helper function for the elow function
        Parameters
        ----------
        label
        df

        Returns
        -------

        """
        cols = []
        for c in df.columns:
            if label in c:
                cols.append(c)

        values = df[cols].values
        non_nans = []
        non_nans_idx = []
        for val_lst in values:
            i = 0
            non_nan_idxs = []
            label_vals = []
            for v in val_lst:
                if isinstance(v, str):
                    non_nan_idxs.append(str(i))
                    label_vals.append(v)
                i += 1
            non_nans.append(','.join(label_vals))
            non_nans_idx.append(','.join(non_nan_idxs))
        return non_nans, non_nans_idx

    def format_mutation_files(self, files: list, case_ids: list) -> None:
        """
        Document that we only keep amino acid changes when there is a ssm consequence documented.

        Parameters
        ----------
        files

        Returns
        -------

        """
        # Since the files are very ugly by default, we want to only keep some cols
        cols_to_keep = ['ssm.ssm_id', 'ssm.consequence.0.transcript.gene.symbol',
                        'ssm.consequence.0.transcript.gene.gene_id', 'ssm.genomic_dna_change', 'ssm.mutation_subtype']
        transcript_cols = ['consequence_type', 'aa_change']

        # We also want to provide the user with a summary of the cases that had mutations and which didn't
        cases_without_mutations = []

        # Since we potentially have multiple transcripts we need to combine the last two columns into a sting/
        # separable list so that we can use it for downstream analyses.
        file_idx = 0
        for f in files:
            try:
                df = pd.read_csv(f, sep='\t')
                formatted_df = pd.DataFrame()
                formatted_df['case_id'] = np.array([case_ids[file_idx] for _ in range(len(df))])
                for c in cols_to_keep:
                    formatted_df[c] = df[c].values
                # If it is non empty we will keep it otherwise delete the file
                for col in transcript_cols:
                    formatted_df[col], formatted_df[col + '_transcript_idx'] = self.add_mutation_values(col, df)

                self.u.save_df(formatted_df, f, sep='\t')

            except Exception as e:
                # Delete the file ToDo: Generalise the split
                cases_without_mutations.append(f.split('_')[-2])
                os.remove(f)
            file_idx += 1
        self.u.dp(["Num cases with mutations: ", len(files) - len(cases_without_mutations),
                   "\nNum cases without mutations: ", len(cases_without_mutations), "\n\nCases: \n",
                   '\n'.join(cases_without_mutations)])

    @staticmethod
    def gen_mutation_api_str(case_id: str, output_file: str) -> str:
        """
        Copied the string from the TCGA website and updated the case_id
        See simple somatic mutations
        https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/#simple-somatic-mutation-endpoint-examples
        Parameters
        ----------
        case_id
        output_file

        Returns
        -------

        """
        return f'curl "https://api.gdc.cancer.gov/ssm_occurrences?format=tsv&fields=ssm.ssm_id,ssm.genomic_dna_change,' \
               f'ssm.mutation_subtype,ssm.consequence.transcript.gene.gene_id,ssm.consequence.transcript.gene.symbol,' \
               f'ssm.consequence.transcript.aa_change,ssm.consequence.transcript.consequence_type&' \
               f'size=5000&filters=%7B%0D%0A%22op%22%3A%22in%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22' \
               f'case.submitter_id%22%2C%0D%0A%22value%22%3A%5B%0D%0A%22{case_id}%22%0D%0A%5D%0D%0A%7D%0D%0A%7D" ' \
               f'> {output_file}'

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
        processed_dir = self.add_windows_slash(processed_dir)
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
                                if self.file_type == 'windows' or 'C:' in dir_str:
                                    print('xcopy ' + dir_str + df + ' ' + processed_dir + df)
                                    os.system('xcopy ' + dir_str + df + ' ' + processed_dir + df)
                                else:
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
            self.u.dp(["Successfully downloaded all files: no.", len(success)])
        else:
            self.u.err_p(["\tSuccessfully downloaded: ", len(success), '\n', "Failed: ", len(fail), '\n\nFailed Files:\n' + '\n'.join(fail)])

        if ouput_file:
            with open(ouput_file, 'w+') as f:
                for d in download_status:
                    f.write('\t'.join(d) + '\n')

        # Remove the header
        return download_status[1:]

    def download(self, manifest_filename=None, manifest_path=None) -> None:
        cmds = []
        manifest_filename = manifest_filename if manifest_filename is not None else self.manifest_file.split(self.u.dir_sep)[-1]
        manifest_path = manifest_path + manifest_filename if manifest_path is not None else self.manifest_file

        with open(manifest_path, 'r+') as f:
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
            # We want to download the last part of the files
            new_manifest = f'{self.split_manifest_dir}{file_cnt}_{manifest_filename}'
            self.write_file(new_manifest, to_write)
            cmds.append(f'{self.gdc_client} download -m {new_manifest} -d {self.download_dir}')
        self.run_cmds(cmds)
