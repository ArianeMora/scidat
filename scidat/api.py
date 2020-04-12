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

from functools import reduce
import time
import pandas as pd
import numpy as np
import os


class API:

    def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
                 max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_'):
        self.mutation_file_dir = None
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
        self.rna_df = None
        self.meth_df = None
        self.rna_meth_df = None
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

    def build_rna_df(self, rna_seq_dir=None) -> np.array:
        """

        Parameters
        ----------
        case_ids
        rna_seq_dir

        Returns
        -------

        """
        rna_seq_dir = rna_seq_dir if rna_seq_dir is not None else self.download_dir
        files = os.listdir(rna_seq_dir)
        count_files = []
        for f in files:
            if 'counts' in f and '.DS' not in f and '.txt' in f:
                count_files.append(rna_seq_dir + f)

        first_file = count_files[0][:-3] + 'gz'
        df = pd.read_csv(count_files[0], sep='\t', header=None)
        df.columns = ['id', self.annotate.annotated_file_dict[first_file]['label']]

        # Add all others.
        for f in count_files[1:]:
            try:
                df_tmp = pd.read_csv(f, sep='\t', header=None)
                name = self.annotate.annotated_file_dict[f[:-3] + 'gz']['label']
                df_tmp.columns = ['id', name]
                df = df.join(df_tmp.set_index('id'), on='id')
            except:
                self.u.warn_p(["Unable to process: ", f])
        self.rna_df = df

    def minify_meth_files(self, processed_dir=None, cpg_min_dir=None):
        cpg_min_dir = cpg_min_dir if cpg_min_dir is not None else self.download_dir
        processed_dir = processed_dir if processed_dir is not None else self.download_dir

        files = os.listdir(processed_dir)
        count_files = []
        cpg_files = []

        for f in files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_files.append(f)
            elif 'min' not in f and '.DS' not in f:
                count_files.append(f)

        self.u.dp(["Minified ", len(cpg_files), "\n Number of count files to process: ", len(count_files)])

        cpg_to_gene = {}

        for f in cpg_files:
            print(f)
            try:
                dups = []
                idx_dict = {}
                df_tmp = pd.read_csv(processed_dir + f, sep='\t')
                i = 0
                beta_values = df_tmp['Beta_value'].values
                cpg_ids = df_tmp['Composite Element REF'].values
                for t in df_tmp['Gene_Symbol'].values:
                    transcripts = t.split(';')
                    beta_value = beta_values[i]
                    cpg_id = cpg_ids[i]
                    for transcript in transcripts:
                        # We want to keep track of the CPG ids so that we can map between them
                        if cpg_to_gene.get(cpg_id):
                            if transcript not in cpg_to_gene.get(cpg_id):
                                cpg_to_gene[cpg_id].append(transcript)
                        else:
                            cpg_to_gene[cpg_id] = [transcript]
                        if idx_dict.get(transcript):
                            if beta_value not in idx_dict[transcript]:
                                dups.append(transcript)
                            else:
                                idx_dict[transcript].append(beta_value)
                        else:
                            idx_dict[transcript] = [beta_value]

                    i += 1

                # We need to go through our dictionary and we'll only keep the maximum methylation value
                transcript_rows = []
                beta_rows = []
                for transcript, beta_values in idx_dict.items():
                    max_beta = 0
                    for b in beta_values:
                        if b > max_beta:
                            max_beta = b
                    transcript_rows.append(transcript)
                    beta_rows.append(max_beta)

                meth_df = pd.DataFrame()
                meth_df['id'] = transcript_rows
                meth_df[f] = beta_rows
                meth_df.to_csv(cpg_min_dir + f, index=False)
            except Exception as e:
                self.u.err_p([f, e])

    def build_meth_df(self, cpg_min_dir: str, df=None, drop_empty_rows=False) -> None:
        """
        Here we build a dataframe based on the methylation data, we can add it to an existing dataframe
        and just join on the columns or create a new dataframe. Here we also choose whether to keep empty rows or only
        keep genes that have full data.

        Parameters
        ----------
        cpg_min_dir
        df
        drop_empty_rows

        Returns
        -------

        """
        files = os.listdir(cpg_min_dir)
        df = pd.DataFrame() if df is None else df
        for f in files:
            try:
                df_tmp = pd.read_csv(cpg_min_dir + f)
                name = self.annotate.annotated_file_dict[f]['label']
                df_tmp.columns = ['id', name]

                df = df.join(df_tmp.set_index('id'), on='gene_id')
            except Exception as e:
                print(e)
        if drop_empty_rows:
            df = df.dropna()
        return df

    def merge_rna_meth_values(self, meth_df: pd.DataFrame, rna_df: pd.DataFrame, merge_col='gene_id') -> pd.DataFrame:
        self.rna_meth_df = meth_df.join(rna_df, on=merge_col)
        return self.rna_meth_df

    def save_annotation(self) -> None:
        """
        Wrapper for annotation saver. Can use the detailed one as well this just provides easy access.
        Returns
        -------

        """
        # save_annotation(self, dir_str: str, filename_str: str, sep = '\t', list_sep = ',') -> None:
        self.annotate.save_annotation(self.meta_dir, 'annotation')

    def get_genes_with_mutations(self, case_ids=None, id_type='symbol') -> np.array:
        """

        Returns
        -------
        A matrix with genes as the first column and case ids for the other columns.
        a cell i (gene) x j (case) = 1 if a mutation exists else 0.
        """

        column = 'ssm.consequence.0.transcript.gene.gene_id'
        if id_type == 'symbol':
            column = 'ssm.consequence.0.transcript.gene.symbol'
        return self.get_mutation_values_on_filter(column, case_ids, 'case_id')

    def get_mutation_values_on_filter(self, column: str, filter_values: list, filter_column: str) -> list:
        mutation_df = self.annotate.get_mutation_df()
        if mutation_df is None:
            self.u.warn_p(["You have not yet built the mutation dataframe. Please run: "
                           "\nbuild_mutation_df(self, mutation_dir: str, "
                           "output_file=None, mutation_prefix='mutation', sep='\t', case_ids=None) -> None"])
            return
        if not filter_values or not filter_column:
            # Now we want to query the gene ID column:
            return list(set(mutation_df[column].values))
        # Otherwise we need to filter by gene ids
        idxs = []
        for value in mutation_df[column].values:
            if value in filter_values:
                idxs.append(True)
            else:
                idxs.append(False)
        return list(set(mutation_df[column].values[np.where(idxs == True)]))

    def get_cases_with_mutations(self, gene_list=None, id_type='symbol') -> list:
        """

        Parameters
        ----------
        gene_list

        Returns
        -------

        """
        filter_column = 'ssm.consequence.0.transcript.gene.gene_id'
        if id_type == 'symbol':
            filter_column = 'ssm.consequence.0.transcript.gene.symbol'

        return self.get_mutation_values_on_filter('case_id', gene_list, filter_column)

    def get_values_from_files(self, case_ids=None, file_type='counts') -> np.array:
        """

        Parameters
        ----------
        file_type
        case_ids

        Returns
        -------

        """

    def get_rna_df(self):
        if self.rna_df is not None:
            return self.rna_df
        else:
            self.build_rna_df()
            return self.rna_df

    def get_cases_with_meta(self, meta: dict, method="all") -> list:
        """
        Here we want to allow the user to specify some dictionary of metadata and retrieve all the cases with that data
        for example: gender: ["male"], race: ["white", "asian"], status:["dead"]
        Parameters
        ----------
        meta
        method: str being "any" i.e. that case meets any of the criteria in the meta_dict or "all" has to meet all of
                the criteria.
        Returns
        -------

        """
        cases_lst = []
        for key, values in meta.items():
            if not isinstance(values, list):
                self.u.err_p([self.u.msg.msg_data_type("get_cases_with_meta", values, "list")])
                return
            cases_lst.append(list(set(self.annotate.clin_df['submitter_id'][self.annotate.clin_df[key].isin(values)])))

        if method == "all":
            return list(reduce(set.intersection, [set(item) for item in cases_lst]))

        elif method == "any":
            # Otherwise just return all the cases that had any of those metadata
            return list(set([item for sublist in cases_lst for item in sublist]))

    def get_merged_rna_meth_df(self, case_ids=None) -> pd.DataFrame:
        """

        Parameters
        ----------
        case_ids

        Returns
        -------

        """
        if self.rna_meth_df is None:
            self.u.err_p([self.u.msg.msg_data_gen("get_merged_rna_meth_df",
                "rna_meth_df", ["minify_meth_files", "build_meth_df", "build_rna_df", "merge_rna_meth_df"])])
            return

        if case_ids is not None:
            return self.rna_meth_df[self.rna_meth_df['case_id'].isin(case_ids)]

        return self.rna_meth_df