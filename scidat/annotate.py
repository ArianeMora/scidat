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

# https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/
from collections import defaultdict
import os
import pandas as pd
from pybiomart import Server
from typing import Tuple

from sciutil import *


class SciException(Exception):

    def __init__(self, message=''):
        Exception.__init__(self, message)


class AnnotateException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class Annotate:

    def __init__(self, output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
                 sep='_', clin_cols=None, sciutil=None):

        self.u = SciUtil() if not sciutil else sciutil
        self.output_dir = output_dir
        self.clinical_file = clinical_file
        self.sample_file = sample_file
        self.manifest_file = manifest_file
        self.file_types = file_types
        if clin_cols is None:
            clin_cols = ['submitter_id', 'project_id', 'age_at_index', 'gender', 'race', 'vital_status', 'tumor_stage']
        self.annotation_cols = ['case_files', 'tumor_stage_num', 'race', 'gender', 'project_id']
        self.clin_cols = clin_cols
        self.clin_df = None
        self.mutation_df = None
        self.annotated_file_dict = None
        self.sample_df = None
        self.sep = sep
        self.case_to_mutation = None
        self.case_to_file = None
        self.gene_annot_dict = None
        self.ens_to_name = None

    def run_setup(self) -> None:
        self.clin_df = self.setup_clin_df()
        self.sample_df = self.setup_sample_df()
        self.annotated_file_dict, self.case_to_file = self.setup_file_dict()
        self.update_clin_data()

    def setup_clin_df(self) -> pd.DataFrame:
        clin_df = pd.read_csv(self.clinical_file, sep='\t')
        clin_df = clin_df[self.clin_cols]
        clin_df = clin_df.drop_duplicates()
        self.u.dp(["Clinical dataframe"])
        self.u.dp([clin_df.head()])

        return clin_df

    def setup_sample_df(self) -> pd.DataFrame:
        # Now we want to add to this dataframe our other information (mutation) and the file ids
        sample_df = pd.read_csv(self.sample_file, sep='\t')
        self.u.dp(["Sample dataframe"])
        self.u.dp([sample_df.head()])

        return sample_df

    def setup_file_dict(self) -> Tuple[dict, defaultdict]:
        # Make a dictionary that connects our cases and our files
        file_dict = {}
        case_to_file = defaultdict(list)

        sample_types = self.sample_df['Sample Type'].values
        case_ids = self.sample_df['Case ID'].values
        files = self.sample_df['File Name'].values
        project_ids = self.sample_df['Project ID'].values
        i = 0
        for f in files:
            # At this point we could filter but it is assumed that the data was pre-filtered by the user
            # prior to downloading the manifest
            case_id = project_ids[i] + '_' + case_ids[i]
            file_dict[f] = {'case_id': case_id, 'sample_type': sample_types[i]}
            case_to_file[case_id].append(f)

            i += 1
        return file_dict, case_to_file

    def update_clin_data(self) -> None:
        normals, tumors, case_files = [], [], []
        project_ids = self.clin_df['project_id'].values

        # Now lets iterate through out clinical data
        i = 0
        for c in self.clin_df['submitter_id'].values:
            # Lets now add the new information
            tumor = 0
            normal = 0
            c_files = []
            case = project_ids[i] + '_' + c
            for f in self.case_to_file[case]:
                if self.annotated_file_dict[f]['sample_type'] == 'Solid Tissue Normal':
                    normal += 1
                elif self.annotated_file_dict[f]['sample_type'] == 'Primary Tumor':
                    tumor += 1
                c_files.append(f)
            normals.append(normal)
            tumors.append(tumor)
            case_files.append(c_files)
            i += 1
        self.clin_df['normal_samples'] = normals
        self.clin_df['tumor_samples'] = tumors
        self.clin_df['case_files'] = case_files

        # Lastly we want to convert our tumor stage into a numeric value (i.e. 1-4 instead of a string)
        stages = []
        for stage in self.clin_df['tumor_stage'].values:
            if 'iv' in stage:
                stages.append(4)
            elif 'iii' in stage:
                stages.append(3)
            elif 'ii' in stage:
                stages.append(2)
            else:
                stages.append(1)
        self.clin_df['tumor_stage_num'] = stages

    def build_annotation(self, annotation_cols=None) -> None:
        # Setup the necessary files.
        self.run_setup()
        # Assign each of our files with a meaningful name based on the meta data.
        # Lets now go through the cases and assign these files with meaningful names
        # The naming format is: tumor-type_gender_race__stage
        if not annotation_cols:
            # use default columns
            annotation_cols = self.annotation_cols
        else:
            self.annotation_cols = annotation_cols
        for value in self.clin_df[annotation_cols].values:
            for f in value[0]:
                i = 0
                for c in annotation_cols:
                    self.annotated_file_dict[f][c] = value[i]
                    i += 1
        # Now lets build our new column dictionary
        for filename, meta in self.annotated_file_dict.items():
            file_type = None
            for ft in self.file_types:
                if ft in filename:
                    file_type = ft
            # ToDo: Update this --> won't work with user specified annotation columns
            label = f'{meta["project_id"]}{self.sep}{meta["sample_type"]}{self.sep}{meta["gender"]}{self.sep}' \
                    f'{meta["race"]}{self.sep}{meta["tumor_stage_num"]}{self.sep}{file_type}' \
                    f'{self.sep}{meta["case_id"]}'
            # Remove any spaces from the label
            self.annotated_file_dict[filename]['label'] = label.replace(' ', '')

    def build_mutation_df(self, mutation_dir: str, mutation_prefix='mutation', sep='\t',
                          case_ids=None) -> None:
        # Create a dataframe of all the mutation
        files = os.listdir(mutation_dir)
        mutation_files = []
        for f in files:
            if mutation_prefix in f:
                if case_ids:
                    for c in case_ids:
                        if c in f:
                            mutation_files.append(mutation_dir + f)
                else:
                    mutation_files.append(mutation_dir + f)
        # Check there were any files
        if len(mutation_files) < 1:
            if case_ids:
                self.u.err_p(["Error in build_mutation_df!\nYou had no files in: ", mutation_dir, "\nwith the prefix: "
                             , mutation_prefix, "\nAnd case ids in: ", ','.join(case_ids),
                              "\nPlease check this directory exists and you have downloaded the mutation data."])
            else:
                self.u.err_p(["Error in build_mutation_df!\nYou had no files in: ", mutation_dir, "\nwith the prefix: "
                              , mutation_prefix, "\nPlease check this directory exists and you have downloaded "
                              "the mutation data."])
            return
        # Now we want to go through each of the files and add it to a dataframe
        dfs = []

        for f in mutation_files[1:]:
            dfs.append(pd.read_csv(f, sep=sep))
        # Now we just concatenate them together
        self.mutation_df = pd.concat(dfs)

    def build_gene_info_file(self, organism_lbl: str, gene_ids: list, output_dir: str) -> Tuple[str, pd.DataFrame]:
        """
        https://jrderuiter.github.io/pybiomart/

        organism_lbl must be one of the datasets available in ensembl i.e. hsapiens_gene_ensembl,
        ToDo: allow the user to have choice in their query
        Parameters
        ----------
        organism_lbl

        Returns
        -------

        """
        server = Server(host='http://www.ensembl.org')

        dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets[organism_lbl])

        gene_info_df = dataset.query(attributes=["ensembl_gene_id", "external_gene_name",  "percentage_gene_gc_content",
                                  "chromosome_name", "start_position","end_position", "strand", "go_id",
                                  "entrezgene_id"],
                      filters={'ensembl_gene_id': gene_ids})
        output_path = self.u.generate_label([output_dir, "gene_info"], ".tsv")
        # Save this to a tsv file
        self.u.save_df(gene_info_df, output_path, sep='\t')
        return output_path, gene_info_df

    def build_gene_annot_dict(self, gene_info_file):
        if not os.path.exists(gene_info_file):
            self.u.err_p(["Your gene info file doesn't exist, please run function: "
                          "build_gene_info_file. \nERR FILENAME: ", gene_info_file])
            return
        # Load csv just created, we want a dictionary on gene id, with values for start, end and go terms (as a list)
        gene_dict = {}
        ens_to_name = {}
        ensembl_id, id_idx, gc_idx, chr_idx, start_idx, end_idx, strand_idx, go_idx, ncbi_idx = 0, 1, 2, 3, 4, 5, 6, 7, 8
        end_err, chr_err, start_err, strand_err = 0, 0, 0, 0

        with open(gene_info_file, 'r+') as fp:
            hdr = True
            for line in fp:
                line = line.split('\t')
                try:
                    if not hdr:
                        g_id = line[id_idx].strip()
                        ens_id = line[ensembl_id].strip()
                        # Here if we have multiple go terms they wil be separated by a comma which we change to a pipe
                        # so we can save it as a tab
                        go_term = line[go_idx].strip().replace(',', '|')
                        chrom = 'chr' + line[chr_idx].strip()
                        # i.e. we're only keeping normal chrs.
                        if chrom:
                            gc_content = float(line[gc_idx].strip())
                            start = int(line[start_idx].strip())
                            end = int(line[end_idx].strip())
                            strand = int(line[strand_idx].strip())
                            ncbi = line[ncbi_idx].strip()
                            # Check if we've already added it
                            gene = gene_dict.get(g_id)
                            if gene:
                                # Check if we have the same chrom
                                if gene['chr'] == chrom:
                                    if gene['start'] == start:
                                        if gene['end'] == end:
                                            if gene['direction'] == strand:
                                                if len(go_term) > 1 and len(go_term.split(':')) > 1 and 'GO' in go_term:
                                                    if go_term not in gene['go_terms']:
                                                        gene['go_terms'].append(
                                                            int(go_term.split(':')[1]))  # go_term.split('/')[-1])
                                            else:
                                                strand_err += 1
                                        else:
                                            end_err += 1
                                    else:
                                        start_err += 1
                                else:
                                    chr_err += 1

                            else:
                                ens_to_name[ens_id] = g_id
                                gene_dict[g_id] = {
                                    'id': g_id,
                                    'chr': chrom,
                                    'gc': gc_content,
                                    'start': start,
                                    'end': end,
                                    'direction': strand,
                                    'go_terms': [],
                                    'ncbi': ncbi # Used for mapping to kegg pathways.
                                }
                                if len(go_term.split(':')) > 1 and 'GO' in go_term:
                                    gene_dict[g_id]['go_terms'].append(
                                        int(go_term.split(':')[1]))
                    hdr = False
                except Exception as e:
                    print(line, e)
        # Print out our summary
        self.u.warn_p(
            ["Sorted our gene information, gene dictionary length: ", len(gene_dict), "\n Error log: \t chr: ", chr_err,
             "\t end: ", end_err, "\t start: ", start_err])

        self.gene_annot_dict = gene_dict
        self.ens_to_name = ens_to_name

        return gene_dict, ens_to_name

    def build_roi(self):
        """
        Here we make regions of interest for our peaks to map to.

        We sort this in the same way that our peak files are sorted, otherwise we would not be going
        through it efficiently.
        """
        chr_dict = {}
        start_err = 0

        for gene_id, values in self.gene_annot_dict.items():
            chrom = values['chr']
            start = values['start']
            # Lets make this have a "fake" start based on the TSS
            if values['direction'] < 0:
                start = values['end']

            # Check if we already have the chrom.
            if chr_dict.get(chrom):
                # Check if we have that start already
                if chr_dict[chrom].get(start):
                    chr_dict[chrom][start].append(gene_id)
                else:
                    chr_dict[chrom][start] = [gene_id]
            else:
                chr_dict[chrom] = {}
                chr_dict[chrom][start] = [gene_id]
        # Again, lets use the ordering of the index keys
        chrs_sorted = list(chr_dict.keys())
        chrs_sorted.sort()
        # Now lets make a list of the gene regions of interest
        gene_rois = []
        for c in chrs_sorted:
            starts_labels = list(chr_dict[c].keys())
            starts_labels.sort()
            for i in starts_labels:
                genes = chr_dict[c][i]
                for g_id in genes:
                    self.gene_annot_dict[g_id]['chr'] = self.gene_annot_dict[g_id]['chr']
                    gene_rois.append(self.gene_annot_dict[g_id])

        return gene_rois

    def add_gene_metadata_to_df(self, df: pd.DataFrame, drop_empty_rows=False) -> pd.DataFrame:
        names = []
        terms = []
        gc_content = []
        ncbi = []
        for g_id in df['id'].values:
            g = g_id.split('.')[0]
            if self.ens_to_name.get(g) is None:
                names.append(None)
                terms.append(None)
                gc_content.append(None)
                ncbi.append(None)
            else:
                gene_info = self.gene_annot_dict[self.ens_to_name.get(g)]
                names.append(self.ens_to_name.get(g))
                terms.append(gene_info['go_terms'])
                gc_content.append(gene_info['gc'])
                ncbi_saved = gene_info['ncbi']
                if ncbi_saved != 'NA':
                    ncbi.append(ncbi_saved)
                else:
                    ncbi.append(None)

        df['gene_id'] = names
        df['gc'] = gc_content
        df['go_terms'] = terms
        df['ncbi'] = ncbi
        if drop_empty_rows:
            df = df.dropna()
        return df

    # Getter functions
    def get_cases(self):
        return list(self.case_to_file.keys())

    def get_mutation_df(self) -> pd.DataFrame:
        if self.mutation_df is None:
            self.u.warn_p(["No mutation dataframe built. Please run: \nbuild_mutation_df(self, mutation_dir: str, "
                           "output_file=None, mutation_prefix='mutation', sep='\t', case_ids=None) -> None"])
            return None
        return self.mutation_df

    # Saving functions
    def save_mutation_df(self, file_path: str, sep='\t') -> None:
        self.u.save_df(file_path, sep)

    def save_annotated_clinical_df(self, dir_str: str, filename: str) -> None:
        if '.csv' not in filename:
            filename += '.csv'
        self.u.save_df(self.clin_df, self.u.generate_label([dir_str, filename]))

    def save_annotation(self, dir_str: str, filename_str: str, sep='\t', list_sep=',') -> None:
        hdr = ['filename']
        first = True
        rows = []

        for filename, meta in self.annotated_file_dict.items():
            row = [filename]
            for m, v in meta.items():
                if first:
                    hdr.append(str(m))
                if isinstance(v, list):
                    row.append(list_sep.join([str(i) for i in v]))
                else:
                    row.append(str(v))
            rows.append(row)
            first = False

        with open(self.u.generate_label([dir_str, filename_str], '.csv'), 'w+') as f:
            f.write(sep.join(hdr) + '\n')
            for row in rows:
                f.write(sep.join(row) + '\n')

        self.u.dp(["Saved annotation file to: ", filename_str])