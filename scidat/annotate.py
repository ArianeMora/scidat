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
from typing import Tuple

from sciutil import SciException, SciUtil


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
        self.file_types = file_types if file_types is not None else ['450', 'counts']
        if clin_cols is None:
            clin_cols = ['project_id', 'age_at_index', 'gender', 'race', 'vital_status',
                         'tumor_stage', 'days_to_death']
        self.annotation_cols = ['case_files', 'tumor_stage_num', 'race', 'gender', 'project_id', 'days_to_death']
        self.clin_cols = clin_cols
        self.case_submitter_id = 'submitter_id'
        self.clin_df = None
        self.mutation_df = None
        self.annotated_file_dict = None
        self.sample_df = None
        self.sep = sep
        self.case_to_mutation = None
        self.case_to_file = None
        self.file_type = None
        self.check_args()

    def set_case_submitter_id(self, case_submitter_id):
        self.case_submitter_id = case_submitter_id

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
        self.output_dir = self.add_windows_slash(self.output_dir)
        self.clinical_file = self.add_windows_slash(self.clinical_file)
        self.sample_file = self.add_windows_slash(self.sample_file)
        self.manifest_file = self.add_windows_slash(self.manifest_file)

        if not os.path.isdir(self.output_dir):
            raise AnnotateException(f'ARGPARSE ERR: output directory in annotate was not a directory: {self.output_dir}')

        if not os.path.isfile(self.manifest_file):
            raise AnnotateException(f'ARGPARSE ERR: manifest file does not exist: {self.manifest_file}')

        if not os.path.isfile(self.sample_file):
            raise AnnotateException(f'ARGPARSE ERR: sample sheet file does not exist: {self.sample_file}')

        if not os.path.isfile(self.clinical_file):
            raise AnnotateException(f'ARGPARSE ERR: clinical file does not exist: {self.clinical_file}')

    def run_setup(self) -> None:
        self.clin_df = self.setup_clin_df()
        self.sample_df = self.setup_sample_df()
        self.annotated_file_dict, self.case_to_file = self.setup_file_dict()
        self.update_clin_data()

    def setup_clin_df(self) -> pd.DataFrame:
        clin_df = pd.read_csv(self.clinical_file, sep='\t')
        # set the submitter id column
        new_clin_cols = []
        # First check if the normal plain submitter_id column isn't in our columns
        if self.case_submitter_id not in list(clin_df.columns):
            self.u.warn_p([f'Warning: the column you selected for submitter ID {self.case_submitter_id} '
                           f'was not in the colums. Resetting submitter_id Available columns: ',
                           ', '.join(list(clin_df.columns)), '\n Run: annotate.set_case_submitter_id() to setup. \n'
                                                             'Continuing with automatic selection.'])
            self.case_submitter_id = None
        else:
            new_clin_cols.append(self.case_submitter_id)
        for c in clin_df:
            if 'submitter_id' in c and self.case_submitter_id is None:
                self.case_submitter_id = c
                self.u.dp(["Submitter ID set as: ", c])
                new_clin_cols.append(c)
        new_clin_cols += self.clin_cols
        self.clin_cols = new_clin_cols
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
        file_ids = self.sample_df['File ID'].values
        i = 0
        for f in files:
            # At this point we could filter but it is assumed that the data was pre-filtered by the user
            # prior to downloading the manifest
            case_id = project_ids[i] + '_' + case_ids[i]
            file_dict[f] = {'case_id': case_id, 'sample_type': sample_types[i], 'submitter_id': case_ids[i],
                            'file_id': file_ids[i], 'project_id': project_ids[i]}
            case_to_file[case_id].append(f)
            i += 1
        return file_dict, case_to_file

    def update_clin_data(self) -> None:
        normals, tumors, case_files = [], [], []
        project_ids = self.clin_df['project_id'].values

        # Now lets iterate through out clinical data
        i = 0
        for c in self.clin_df[self.case_submitter_id].values:
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
            if 'htseq.counts.gz' in filename:
                file_type = 'htseq.counts'
            elif 'HumanMethylation450' in filename:
                file_type = 'HumanMethylation450'
            else:
                print("FILETYPE NOT FOUND", filename)
            # ToDo: Update this --> won't work with user specified annotation columns
            # ToDo: WARN: Must have the case_id as the last element
            label = f'{meta["project_id"]}{self.sep}{meta["sample_type"]}{self.sep}{meta["gender"]}{self.sep}' \
                    f'{meta["race"]}{self.sep}{meta["tumor_stage_num"]}{self.sep}{file_type}{self.sep}' \
                    f'{meta["days_to_death"]}{self.sep}{meta["case_id"]}{self.sep}{meta["file_id"]}'

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
            raise AnnotateException("Build mutation exception")
        # Now we want to go through each of the files and add it to a dataframe
        dfs = []

        for f in mutation_files:
            dfs.append(pd.read_csv(f, sep=sep))
        # Now we just concatenate them together and drop any rows that are all nans
        self.mutation_df = pd.concat(dfs).dropna(how='all')

    # Getter functions
    def get_cases(self):
        cases = list(self.case_to_file.keys())
        case_ids = []
        for c in cases:
            case_ids.append(c.split('_')[-1])
        return case_ids

    def get_mutation_df(self) -> pd.DataFrame:
        if self.mutation_df is None:
            err_msg = self.u.msg.msg_data_gen("get_mutation_df", "self.mutation_df",
                                              ["download_data_using_api()", "build_mutation_df()"])
            self.u.err_p([err_msg])
            raise AnnotateException(err_msg)

        return self.mutation_df

    # Saving functions
    def save_mutation_df(self, dir_str: str, filename: str, sep='\t') -> None:
        self.u.save_df(self.mutation_df, self.u.generate_label([dir_str, filename], '.tsv'), sep)

    def save_annotated_clinical_df(self, dir_str: str, filename: str) -> None:
        self.u.save_df(self.clin_df, self.u.generate_label([dir_str, filename], '.csv'))

    def save_annotation(self, dir_str: str, filename: str, sep='\t', list_sep=',') -> None:
        hdr = ['filename']
        first = True
        rows = []

        for fname, meta in self.annotated_file_dict.items():
            row = [fname]
            for m, v in meta.items():
                if first:
                    hdr.append(str(m))
                if isinstance(v, list):
                    row.append(list_sep.join([str(i) for i in v]))
                else:
                    row.append(str(v))
            rows.append(row)
            first = False

        with open(self.u.generate_label([dir_str, filename], '.csv'), 'w+') as f:
            f.write(sep.join(hdr) + '\n')
            for row in rows:
                f.write(sep.join(row) + '\n')

        self.u.dp(["Saved annotation file to: ", filename])