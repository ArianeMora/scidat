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

from collections import defaultdict
import os
import pandas as pd
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
                 sep='_', mutations_file=None, clin_cols=None):
        self.u = SciUtil()
        self.output_dir = output_dir
        self.clinical_file = clinical_file
        self.sample_file = sample_file
        self.manifest_file = manifest_file
        self.mutations_file = mutations_file
        self.file_types = file_types
        if clin_cols is None:
            clin_cols = ['submitter_id', 'project_id', 'age_at_index', 'gender', 'race', 'vital_status', 'tumor_stage']
        self.annotation_cols = ['case_files', 'tumor_stage_num', 'mutations', 'race', 'gender', 'project_id']
        self.clin_cols = clin_cols
        self.clin_df = None
        self.mutation_df = None
        self.file_dict = None
        self.sample_df = None
        self.sep = sep

    def run_setup(self) -> None:
        self.clin_df = self.setup_clin_df()
        self.sample_df = self.setup_sample_df()
        self.file_dict, self.case_to_file = self.setup_file_dict()
        self.mutation_df, self.case_to_mutation = self.setup_mutation_df()
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

    def setup_mutation_df(self) -> Tuple[pd.DataFrame, defaultdict]:
        # Check if we also have mutation data
        if self.mutations_file:
            mutation_df = pd.read_csv(self.mutations_file, sep='\t')
            self.u.dp(["Muttaion dataframe"])
            self.u.dp([mutation_df.head()])

            case_to_mutation = defaultdict(list)
            mutations = mutation_df['Consequences'].values
            i = 0
            for cases in mutation_df['Cases'].values:
                cases = cases.split(',')
                mutation = mutations[i]
                for c in cases:
                    case_to_mutation[c].append(mutation)
                i += 1
            return mutation_df, case_to_mutation

        # Otherwise return an empty dataframe
        return pd.DataFrame(), {}

    def update_clin_data(self) -> None:
        normals, tumors, mutations, case_files = []
        project_ids = self.clin_df['project_id'].values

        # Now lets iterate through out clinical data
        i = 0
        for c in self.clin_df['submitter_id'].values:
            # Lets now add the new information
            mutations.append(self.case_to_mutation.get(c))
            tumor = 0
            normal = 0
            c_files = []
            case = project_ids[i] + '_' + c
            for f in self.case_to_file[case]:
                if self.file_dict[f]['sample_type'] == 'Solid Tissue Normal':
                    normal += 1
                elif self.file_dict[f]['sample_type'] == 'Primary Tumor':
                    tumor += 1
                c_files.append(f)
            normals.append(normal)
            tumors.append(tumor)
            case_files.append(c_files)
            i += 1
        self.clin_df['normal_samples'] = normals
        self.clin_df['tumor_samples'] = tumors
        self.clin_df['case_files'] = case_files
        # Check if we even have mutations
        if len(self.mutation_df.values) > 1:
            self.clin_df['mutations'] = mutations
            # Update our mutations to only keep the first element (otherwise it is too specific)
            # maybe later on we'll change this
            # ToDo: Potential user selection here
            mutations = []
            for m in self.clin_df['mutations']:
                if m is not None:
                    mutations.append(m[0])
                else:
                    mutations.append(None)

            self.clin_df['mutations'] = mutations
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
        # The naming format is: tumor-type_gender_race_mutation_stage
        if not annotation_cols:
            # use default columns
            annotation_cols = self.annotation_cols
        else:
            self.annotation_cols = annotation_cols
        for value in self.clin_df[annotation_cols].values:
            for f in value[0]:
                i = 0
                for c in annotation_cols:
                    self.file_dict[f][c] = value[i]
                    i += 1
        # Now lets build our new column dictionary
        for filename, meta in self.file_dict.items():
            file_type = None
            for ft in self.file_types:
                if ft in filename:
                    file_type = ft
            label = f'{meta["project_id"]}{self.sep}{meta["sample_type"]}{self.sep}{meta["gender"]}{self.sep}' \
                    f'{meta["race"]}{self.sep}{meta["stage"]}{self.sep}{meta["mutation"]}{self.sep}{file_type}' \
                    f'{self.sep}{meta["case_id"]}'
            # Remove any spaces from the label
            self.file_dict[filename]['label'] = label.replace(' ', '')

    def save_annotated_clinical_df(self, dir_str: str, filename: str) -> None:
        if '.csv' not in filename:
            filename += '.csv'
        self.u.save_df(self.clin_df, self.u.generate_label([dir_str, filename]))

    def save_annotation(self, dir_str: str, filename: str, sep=',') -> None:
        hdr = ['filename']
        first = True
        rows = []

        for filename, meta in self.file_dict.items():
            row = [filename]
            for m, v in meta.items():
                if first:
                    hdr.append(str(m))
                row.append(str(v))
            rows.append(row)

        with open(self.u.generate_label([dir_str, filename]), 'w+') as f:
            f.write(sep.join(hdr))
            for row in rows:
                f.write(sep.join(row))

        self.u.dp(["Saved annotation file to: ", filename])
