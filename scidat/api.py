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
from sciutil import SciUtil, SciException
from scidat import Download
from scidat import Annotate

from functools import reduce
import time
from typing import Tuple
import pandas as pd
import numpy as np
import os


class APIException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class API:
    """
    API provides an interface to the download and annotate classes. They can also be used independently. Finally  it adds
    helpful interfaces that query between the data-structures.
    """

    def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, download_dir, meta_dir, requires_lst=None, clin_cols=None,
                 max_cnt=100, sciutil=None, split_manifest_dir=None, sep='_'):
        self.mutation_file_dir = None
        self.u = SciUtil() if not sciutil else sciutil
        self.download_dir = download_dir
        self.manifest_file = manifest_file
        self.gdc_client = gdc_client
        self.clinical_file = clinical_file
        self.sample_file = sample_file
        self.requires_lst = [] if not requires_lst else requires_lst
        self.max_cnt = max_cnt
        self.split_manifest_dir = split_manifest_dir if split_manifest_dir is not None else self.download_dir
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

    def download_data_from_manifest(self):
        self.download.download()
        # Copy downloads
        self.download.copy_downloads_to_new_dir(self.download_dir)

    def download_mutation_data(self, case_ids=None) -> None:
        # If there are no case ids then we'll download the mutation data for all cases in the annotation data.
        if not case_ids:
            if self.annotate.annotated_file_dict is None:
                msg = self.u.msg.msg_data_gen("download_mutation_data", "annotate.case_ids", ["build_annotation"])
                self.u.err_p([msg])
                raise APIException(msg)
            case_ids = case_ids if case_ids is not None else list(self.annotate.get_cases())
            self.u.warn_p(
                ["Warning: you didn't provide any case ids. Are you sure you want to download all the cases?\n"
                 "Total number: ", len(case_ids), "\nTerminate the program now if you wish to stop. Waiting 2 seconds."])
            time.sleep(2)
            self.u.warn_p(["Continuing..."])

        self.download.download_data_using_api(case_ids, 'mutation')

    def build_mutation_df(self, mutation_dir=None):
        mutation_dir = self.download_dir if not mutation_dir else mutation_dir
        self.annotate.build_mutation_df(mutation_dir)

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
        # first check if they have performed the annotation step
        if self.annotate.annotated_file_dict is None:
            msg = self.u.msg.msg_data_gen("build_rna_df", "annotate.annotated_file_dict", ["build_annotation"])
            self.u.err_p([msg])
            raise APIException(msg)
        rna_seq_dir = rna_seq_dir if rna_seq_dir is not None else self.download_dir
        files = os.listdir(rna_seq_dir)
        count_files = []
        for f in files:
            if 'counts' in f and '.DS' not in f and '.tsv' in f:
                count_files.append(f)

        first_file = count_files[0][:-3] + 'gz'
        df = pd.read_csv(rna_seq_dir + count_files[0], sep='\t', header=None)
        df.columns = ['id', self.annotate.annotated_file_dict[first_file]['label']]

        # Add all others.
        for f in count_files[1:]:
            try:
                df_tmp = pd.read_csv(rna_seq_dir + f, sep='\t', header=None)
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
        cpg_files = []

        for f in files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_files.append(f)

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

    def build_meth_df(self, cpg_min_dir: str, df=None, drop_empty_rows=False, join_id='gene_id') -> None:
        """
        Here we build a dataframe based on the methylation data, we can add it to an existing dataframe
        and just join on the columns or create a new dataframe. Here we also choose whether to keep empty rows or only
        keep genes that have full data.

        Parameters
        ----------
        cpg_min_dir
        df
        drop_empty_rows
        join_id: Tells us which column we want to join on, this must be a column in DF if we're joining to that
        Returns
        -------

        """
        if self.annotate.annotated_file_dict is None:
            msg = self.u.msg.msg_data_gen("build_meth_df", "annotate.annotated_file_dict", ["build_annotation"])
            self.u.err_p([msg])
            raise APIException(msg)
        files = os.listdir(cpg_min_dir)
        cpg_files = []

        for f in files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_files.append(f)

        start_idx = 0
        if df is None and len(cpg_files) > 0:
            df = pd.read_csv(cpg_min_dir + cpg_files[0])
            name = self.annotate.annotated_file_dict[cpg_files[0]]['label']
            df.columns = ['id', name]
            start_idx = 1

        for f in cpg_files[start_idx:]:
            try:
                df_tmp = pd.read_csv(cpg_min_dir + f)
                name = self.annotate.annotated_file_dict[f]['label']
                df_tmp.columns = ['id', name]

                df = df.join(df_tmp.set_index('id'), on=join_id)
            except Exception as e:
                print(e)
        if drop_empty_rows:
            df = df.dropna()
        self.meth_df = df

    def merge_rna_meth_values(self, meth_df: pd.DataFrame, rna_df: pd.DataFrame, index_col='id', merge_col='gene_id') -> pd.DataFrame:
        self.rna_meth_df = rna_df.join(meth_df.set_index(index_col), on=merge_col)
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
        # Throws an annotation API exception if the mutation df hasn't been created.
        mutation_df = self.annotate.get_mutation_df()

        if not filter_values or not filter_column:
            # Now we want to query the gene ID column:
            return list(set(mutation_df[column].values))
        # Otherwise we need to filter by gene ids
        idxs = []
        for value in mutation_df[filter_column].values:
            if value in filter_values:
                idxs.append(True)
            else:
                idxs.append(False)
        return list(set(mutation_df[column].values[idxs]))

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

    def get_values_from_df(self, df: pd.DataFrame, gene_id_column: str, case_ids=None, gene_ids=None,
                           column_name_includes=None, column_name_method="all") -> Tuple[np.array, list, pd.DataFrame]:
        """

        Parameters
        ----------
        df
        case_ids
        gene_ids
        gene_id_column

        Returns
        -------

        """
        if case_ids is None and gene_ids is None:
            return df.values, list(df.columns), df
        # Otherwise we only want to get the columns with the case_id in it
        columns = [] if case_ids is not None else list(df.columns)
        for c in df.columns:
            for case in case_ids:
                # ToDo: WARN this won't be stable if we let users choose their own format for the filenames
                if case == c.split(self.sep)[-1]:
                    columns.append(c)
                    break
        # Lets check if we also have a column requirement
        if column_name_includes is not None:

            new_columns = []
            for c in columns:
                i = 0
                for req in column_name_includes:
                    if req in c:
                        i += 1
                if column_name_method == 'all':
                    if i == len(column_name_includes):
                        new_columns.append(c)
                elif len(column_name_includes) > 0:
                    new_columns.append(c)
            columns = new_columns
        # Lets check if they are filtering on the gene id as well
        if gene_ids is not None:
            idxs = []
            i = 0
            for gene_id in df[gene_id_column].values:
                if gene_id in gene_ids:
                    idxs.append(i)
                i += 1

            new_df = pd.DataFrame(df[[gene_id_column] + columns].values[idxs], columns=[gene_id_column] + columns)
            return new_df.values, list(new_df.columns), new_df

        return df[[gene_id_column] + columns].values, [gene_id_column] + columns, df[[gene_id_column] + columns]

    def get_rna_df(self):
        if self.rna_df is not None:
            return self.rna_df
        else:
            self.build_rna_df()
            return self.rna_df

    def get_meth_df(self):
        if self.meth_df is not None:
            return self.meth_df
        msg = self.u.msg.msg_data_gen("get_meth_df", "meth_df", ["minify_meth_files", "build_meth_df"])
        self.u.err_p([msg])
        raise APIException(msg)

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
        if self.annotate.annotated_file_dict is None:
            msg = self.u.msg.msg_data_gen("get_cases_with_meta", "annotate.clinical_df", ["build_annotation"])
            self.u.err_p([msg])
            raise APIException(msg)

        for key, values in meta.items():
            if not isinstance(values, list):
                msg = self.u.msg.msg_data_type("get_cases_with_meta", values, "list")
                self.u.err_p([msg])
                raise APIException(msg)

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