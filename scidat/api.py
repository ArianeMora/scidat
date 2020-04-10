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
from sciutil import SciUtil
from scidat import Download
from scidat import Annotate

import time
import pandas as pd
import numpy as np

class API:

    def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
                 max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_'):
        self.u = SciUtil() if not sciutil else sciutil
        self.download_dir = download_dir
        self.manifest_file = manifest_file
        self.gdc_client = gdc_client
        self.clinical_file = clinical_file
        self.sample_file = sample_file
        self.requires_lst = [] if not requires_lst else requires_lst
        self.max_cnt = max_cnt
        self.split_manifest_dir = split_manifest_dir
        self.meta_dir = meta_dir
        self.sep = sep
        self.clin_cols = clin_cols
        self.download = Download(self.manifest_file, self.split_manifest_dir, self.download_dir, self.gdc_client,
                                 self.max_cnt, self.u)
        self.annotate = Annotate(self.meta_dir, self.clinical_file, self.sample_file, self.manifest_file,
                                 self.requires_lst, self.sep, self.clin_cols, self.u)

    def download_mutation_data(self, case_ids=None) -> None:
        # If there are no case ids then we'll download the mutation data for all cases in the annotation data.
        if not case_ids:
            case_ids = case_ids if case_ids is not None else list(self.annotate.get_cases())
            self.u.warn_p(
                ["Warning: you didn't provide any case ids. Are you sure you want to download all the cases?\n"
                 "Total number: ", len(case_ids), "\nTerminate the program now if you wish to stop. Waiting 2 seconds."])
            time.sleep(2)
            self.u.warn_p(["Continuing..."])

        self.download.download_data_using_api(case_ids, 'mutation')

    def build_annotation(self) -> None:
        self.annotate.build_annotation()

    def save_annotation(self) -> None:
        """
        Wrapper for annotation saver. Can use the detailed one as well this just provides easy access.
        Returns
        -------

        """
        # save_annotation(self, dir_str: str, filename_str: str, sep = '\t', list_sep = ',') -> None:
        self.annotate.save_annotation(self.meta_dir, 'annotation')

    def get_genes_with_mutations(self, case_ids=None, mutation_file_dir=None) -> np.array:
        """

        Returns
        -------
        A matrix with genes as the first column and case ids for the other columns.
        a cell i (gene) x j (case) = 1 if a mutation exists else 0.
        """
        case_ids = case_ids if case_ids is not None else list(self.annotate.get_cases())
        mutation_file_dir = mutation_file_dir if mutation_file_dir is not None else self.download_dir
        if self.mutation_df is None:
            self.annotate.build_mutation_df()
        return []

    def get_cases_with_mutations(self, mutations=None) -> list:
        """

        Parameters
        ----------
        mutations

        Returns
        -------

        """

    def get_methylation_beta_values(self, case_ids=None) -> np.array:
        """

        Parameters
        ----------
        case_ids

        Returns
        -------

        """

    def get_rnaseq_count_values(self, case_ids=None) -> np.array:
        """

        Parameters
        ----------
        case_ids

        Returns
        -------

        """

    def get_cases_with_meta(self, meta:dict) -> list:
        """

        Parameters
        ----------
        meta

        Returns
        -------

        """

    def get_merged_beta_count_df(self, case_ids=None) -> pd.DataFrame:
        """

        Parameters
        ----------
        case_ids

        Returns
        -------

        """