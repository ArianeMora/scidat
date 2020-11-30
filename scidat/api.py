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

    def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, download_dir, meta_dir,
                 annotation_file: str, requires_lst=None, clin_cols=None,
                 max_cnt=100, sciutil=None, split_manifest_dir=None, sep='_'):
        if gdc_client and '.exe' in gdc_client:
            self.file_type = 'windows'
        else:
            self.file_type = 'unix'
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
        self.annotation_file = annotation_file
        self.setup_paths()
        self.annotation_df = pd.read_csv(self.annotation_file)
        self.annotation_df.set_index('ensembl_gene_id', inplace=True)
        self.download = Download(self.manifest_file, self.split_manifest_dir, self.download_dir, self.gdc_client,
                                 self.max_cnt, self.u)
        self.annotate = Annotate(self.meta_dir, self.clinical_file, self.sample_file, self.manifest_file,
                                 self.requires_lst, self.sep, self.clin_cols, self.u)

    def add_windows_slash(self, change_str):
        if change_str and '/' in change_str and (self.file_type == 'windows' or 'C:' in change_str):
            return change_str.replace('/', '\\')
        return change_str

    def setup_paths(self):
        self.meta_dir = self.add_windows_slash(self.meta_dir)
        self.download_dir = self.add_windows_slash(self.download_dir)
        self.split_manifest_dir = self.add_windows_slash(self.split_manifest_dir)
        self.gdc_client = self.add_windows_slash(self.gdc_client)
        self.clinical_file = self.add_windows_slash(self.clinical_file)
        self.sample_file = self.add_windows_slash(self.sample_file)
        self.manifest_file = self.add_windows_slash(self.manifest_file)
        self.annotation_file = self.add_windows_slash(self.annotation_file)

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

    def minify_meth_files(self, processed_dir=None, cpg_min_dir=None, beta_cutoff=0.5):
        cpg_min_dir = cpg_min_dir if cpg_min_dir is not None else self.download_dir
        processed_dir = processed_dir if processed_dir is not None else self.download_dir

        files = os.listdir(processed_dir)
        cpg_files = []

        for f in files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_files.append(f)

        for f in cpg_files:
            print(f)
            try:
                df_tmp = pd.read_csv(processed_dir + f, sep='\t')
                i = 0
                gene_symbols = []
                for g in df_tmp['Gene_Symbol'].values:
                    gene_symbols.append(g.split(';')[0])

                meth_df = pd.DataFrame()
                meth_df['Composite Element REF'] = df_tmp['Composite Element REF'].values
                meth_df['external_gene_name'] = gene_symbols
                meth_df['beta_value'] = df_tmp['Beta_value'].values
                meth_df.to_csv(cpg_min_dir + f, index=False)
            except Exception as e:
                self.u.err_p([f, e])

    def build_meth_df(self, cpg_min_dir: str, df=None, drop_empty_rows=False, join_id='external_gene_name',
                      meth_id='external_gene_name') -> None:
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
        files_parsed = []
        cpg_to_genes = {}
        cpg_to_beta = {}
        for f in cpg_files[start_idx:]:
            try:
                df_tmp = pd.read_csv(cpg_min_dir + f)
                name = self.annotate.annotated_file_dict[f]['label']
                df_tmp.columns = ['Composite Element REF', meth_id, name]
                gene_names = df_tmp[meth_id].values
                beta_values = df_tmp[name].values
                # Build dict from each one
                for i, cpg in enumerate(df_tmp['Composite Element REF'].values):
                    # Don't keep cpgs that don't have gene names
                    if gene_names[i] != '.':
                        if cpg_to_genes.get(cpg) == None:
                            cpg_to_genes[cpg] = gene_names[i]
                            cpg_to_beta[cpg] = {}
                        if cpg_to_genes.get(cpg) != gene_names[i]:
                            self.u.warn_p(['WARN mismatch between gene names for same CpG: ', cpg, gene_names[i],
                                           cpg_to_genes.get(cpg)])
                        cpg_to_beta[cpg][name] = beta_values[i]
                files_parsed.append(name)
            except Exception as e:
                self.u.warn_p(['WARNING: minify meth files, unable to parse: ', f,
                               '\nIs this file missisng from your sample or clinical file?'])
        df = pd.DataFrame()
        # Now we want to build the df with each of our columns for the beta values and the gene names should be
        # always the same
        f_betas = {}
        gene_names = []
        cpgs = []
        for cpg, beta_vals in cpg_to_beta.items():
            for f in files_parsed:
                if not f_betas.get(f):
                    f_betas[f] = []
                f_betas[f].append(beta_vals.get(f))
            gene_names.append(cpg_to_genes.get(cpg))
            cpgs.append(cpg)
        df['cpg-id'] = cpgs
        df[join_id] = gene_names
        for f, vals in f_betas.items():
            df[f] = vals

        if drop_empty_rows:
            df = df.dropna()
        self.meth_df = df

    def merge_rna_meth_values(self, meth_df: pd.DataFrame, rna_df: pd.DataFrame, index_col='external_gene_name',
                              merge_col='external_gene_name') -> pd.DataFrame:
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

    def get_mutations_for_specific_gene(self, column: str, filter_values: list, filter_column: str, gene: str,
                                        id_type='symbol')\
            -> list:
        gene_column = 'ssm.consequence.0.transcript.gene.gene_id'
        if id_type == 'symbol':
            gene_column = 'ssm.consequence.0.transcript.gene.symbol'
        # Throws an annotation API exception if the mutation df hasn't been created.
        mutation_df = self.annotate.get_mutation_df()
        genes = mutation_df[gene_column].values
        # Otherwise we need to filter by gene ids
        idxs = []
        i = 0
        for value in mutation_df[filter_column].values:
            found = False
            for f in filter_values:
                if isinstance(value, str) and f in value:
                    found = True
                    break

            if genes[i] == gene and found:
                idxs.append(i)
            i += 1
        if column is None:
            return idxs

        return list(set(mutation_df[column].values[idxs]))

    @staticmethod
    def get_columns_without_cases(df: pd.DataFrame, cases: list) -> list:
        """
        Here we want to return the columns without the cases in it, this enables us to have a "control" group.

        Parameters
        ----------
        df
        cases

        Returns
        -------

        """
        columns = []
        for c in df.columns:
            found = False
            for case in cases:
                if case in c:
                    found = True
                    break
            if not found:
                columns.append(c)
        return columns

    def filter_columns_on_gene_value(self, df: pd.DataFrame, gene_id: str, id_column: str,
                                   cutoff_value: float, method: str, ensembl_id=True) -> list:
        """
        Returns the cases that have the gene value
        Parameters
        ----------
        gene_id
        cutoff_value
        method

        Returns
        -------

        """
        if method not in ['greater', 'less']:
            self.u.err_p([f'ARG ERROR: filter_cases_on_gene_value, you passed: {method}, but value must be one of: ',
                          ','.join(['greater', 'less'])])
            raise APIException(f'ARG ERROR: filter_cases_on_gene_value, you passed: {method}, but value must'
                               f' be one of:{",".join(["greater", "less"])}')
        gene_idx = None
        i = 0
        for g in df[id_column].values:
            if ensembl_id:
                # Removes the last tag.
                g = g.split('.')[0]
            if g == gene_id:
                gene_idx = i
                break
            i += 1
        if gene_idx is None:
            self.u.err_p([f'Gene {gene_id} not found. Please check your gene id column {id_column} is correct. If you'
                          f' have an ensembl ID make sure you select ensembl_id=True if you dont want the .X, '
                          f'set ensembl_id=Flase if you want to check versions as well.'])
            return
        values = df.values[gene_idx]
        c_i = 0
        columns = []
        if method == 'less':
            for c in df.columns:
                if c != gene_id:
                    if values[c_i] < cutoff_value:
                        columns.append(c)
                        break
                c_i += 1
        elif method == 'greater':
            for c in df.columns:
                if c != gene_id:
                    if values[c_i] < cutoff_value:
                        columns.append(c)
                        break
                c_i += 1
        return columns

    @staticmethod
    def add_meta_to_cases(df: pd.DataFrame, cases: list, meta_str: str, alt_str=None, sep='_') -> pd.DataFrame:
        """
        Adds a meta tag to the end of the column
        Parameters
        ----------
        df
        cases
        meta_str
        sep

        Returns
        -------

        """
        new_df = pd.DataFrame()
        for c in df.columns:
            found = False
            for case in cases:
                if case in c:
                    new_df[f'{c}{sep}{meta_str}'] = df[c].values
                    found = True
                    break
            if not found:
                if alt_str is not None:
                    new_df[f'{c}{sep}{alt_str}'] = df[c].values
                else:
                    new_df[c] = df[c].values
        return new_df

    def get_mutation_values_on_filter(self, column: str, filter_values: list, filter_column: str) -> list:
        # Throws an annotation API exception if the mutation df hasn't been created.
        mutation_df = self.annotate.get_mutation_df()

        if not filter_values or not filter_column:
            # Now we want to query the gene ID column:
            return list(set(mutation_df[column].values))

        # Otherwise we need to filter by gene ids
        idxs = []
        for value in mutation_df[filter_column].values:
            found = False
            for f in filter_values:
                if isinstance(value, str) and f in value:
                    found = True
                    break
            idxs.append(found)

        if column is None:
            return idxs

        return list(set(mutation_df[column].values[idxs]))

    @staticmethod
    def create_annotation_for_columns(df: pd.DataFrame, includes: list, new_annotation: list) -> list:
        """
        Makes it easy to create annotations for the columns based on annotations stored in the column name
        Parameters
        ----------
        df
        includes
        new_annotation

        Returns
        -------

        """
        annotations = []
        for c in df.columns:
            i = 0
            found = False
            for value in includes:
                if '_' in value and value in c:
                    annotations.append(new_annotation[i])
                    found = True
                elif f'_{value}' in c:
                    annotations.append(new_annotation[i])
                    found = False
                i += 1
            if not found:
                annotations.append('None')

        return annotations

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
                if case in c: # == c.split(self.sep)[-1]:
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

    def get_files_with_meta(self, meta: dict, method="all") -> list:
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
        files_lst = []
        if self.annotate.annotated_file_dict is None:
            msg = self.u.msg.msg_data_gen("get_cases_with_meta", "annotate.clinical_df", ["build_annotation"])
            self.u.err_p([msg])
            raise APIException(msg)

        for filename, case in self.annotate.annotated_file_dict.items():
            count_meeting_req = 0
            for key, values in meta.items():
                if not isinstance(values, list):
                    msg = self.u.msg.msg_data_type("get_cases_with_meta", values, "list")
                    self.u.err_p([msg])
                    raise APIException(msg)
                if case[key] in values:
                    if method == 'any':
                        files_lst.append(case['label'])
                        break
                    else:
                        count_meeting_req += 1
            if method == 'all' and count_meeting_req == len(meta):
                files_lst.append(case['label'])

        return list(set(files_lst))

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

    def add_gene_metadata_to_df(self, df:pd.DataFrame, left_id='id', right_id='ensembl_gene_id'):
        """ Adds gene name information to the RNAseq dataframe """
        # We need to remove the '.' from the ensembl ids
        ensembl_ids = []
        for g in df[left_id].values:
            ensembl_ids.append(g.split('.')[0])

        df[right_id] = ensembl_ids
        df.set_index(right_id, inplace=True)

        joined_df = df.join(self.annotation_df, how='left', on=right_id)
        return joined_df

    def set_gene_annotation_file(self, filepath):
        """ Set the annotation file to override the default. """
        self.annotation_df = pd.read_csv(filepath)