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

from scidat.api import API, APIException
from scidat.annotate import AnnotateException

import shutil
import tempfile
import unittest

from scidat import API
from sciutil import Biomart

import pandas as pd
import numpy as np

import os

meta_dir = '/Users/ariane/Documents/data/tcga/kidney/'
output_folder = '/Users/ariane/Documents/code/scivae_paper/data/rcc_datasets/'
manifest_file = f'{meta_dir}gdc_manifest_20200416_055430.txt'
gdc_client = f'{meta_dir}./gdc-client'
download_dir = f'{meta_dir}downloads/'

clinical_file = f'{meta_dir}clinical.tsv'
sample_file = f'{meta_dir}gdc_sample_sheet.2020-04-07.tsv'
processed_dir = f'{meta_dir}processed/'


class TestAPI(unittest.TestCase):

    def setUp(self):
        # Flag to set data to be local so we don't have to download them repeatedly. ToDo: Remove when publishing.
        self.local = False
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        if self.local:

            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_') + '/'

        # def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
        #                  max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_'):

        self.data_dir = os.path.join(THIS_DIR, 'data/')
        gdc_client = self.data_dir + './gdc-client'
        clinical_file = self.data_dir + 'clinical.txt'
        sample_file = self.data_dir + 'sample_sheet.txt'
        manifest_file = self.data_dir + 'manifest.tsv'

        self.api = API(manifest_file, gdc_client, clinical_file, sample_file, self.tmp_dir, self.tmp_dir,
                            max_cnt=1, requires_lst=['counts', 'm450'])

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_download(self):

        self.api.download_data_from_manifest()
        # Now there should be an extra file in the downloads dir
        files_post = os.listdir(self.tmp_dir)

        # Check file name
        self.assertEqual(files_post[0], 'd3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv')
        # Run the download check
        download_status = self.api.download.check_downloads(self.tmp_dir + 'download_status.csv')
        download_status.sort()

        # Check the download status was correctly assigned
        self.assertEqual(download_status[0][0], '001ae925-102c-4818-8eb0-c8d2e5726e7c')
        self.assertEqual(download_status[0][5], 'True')
        self.assertEqual(download_status[-1][0], '19601351-3c26-4293-b87d-97222cd64a19')
        self.assertEqual(download_status[-1][5], 'True')

        # Check the file was written with the download status
        self.assertEqual(os.path.exists(self.tmp_dir + 'download_status.csv'), True)

    def test_download_mutation_data(self):
        # Here we want to download the mutation files
        case_ids = ['TCGA-A3-3308', 'TCGA-KN-8422', 'TCGA-CZ-5989-11A', 'TCGA-A4-8312-01A']
        with self.assertRaises(APIException):
            self.api.download_mutation_data()

        # Now we want to first build the annotation
        self.api.build_annotation()
        self.api.download_mutation_data(case_ids)

        labels = []
        for c in case_ids:
            labels.append(self.api.u.generate_label(['mutation', c], '.tsv'))
        # Let's check if the files were downloaded
        files = os.listdir(self.tmp_dir)
        overlap = 0
        print(labels, files)
        for f in files:
            if f in labels:
                overlap += 1
        # There should only be two files with mutation data
        self.assertEqual(overlap, len(labels) - 2)

        # self.api.download_mutation_data()

    def test_build_rna_df(self):
        with self.assertRaises(APIException):
            self.api.build_rna_df()

        # Now we want to first build the annotation
        self.api.build_annotation()

        # Now lets copy over the data to the download dir (i.e. the temp dir)
        os.system(f'cp {self.data_dir}d3f73c0f-d518-4e91-b038-a4360495ee27.htseq.counts.tsv {self.tmp_dir}')

        # Lets now build it and run our tests
        self.api.build_rna_df()
        df = self.api.get_rna_df()
        self.assertEqual("TCGA-KIRC_PrimaryTumor_male_asian_3_counts_TCGA-KIRC_TCGA-A3-3308", df.columns[1])
        self.assertEqual(df['TCGA-KIRC_PrimaryTumor_male_asian_3_counts_TCGA-KIRC_TCGA-A3-3308'].values[3], 753)

    def test_minify_meth_files(self):
        self.api.minify_meth_files(self.data_dir, self.tmp_dir)
        files = os.listdir(self.data_dir)
        output_files = os.listdir(self.tmp_dir)
        cpg_in = []
        cpg_out = []
        for f in files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_in.append(f)
        for f in output_files:
            if 'HumanMethylation450' in f and 'min' not in f and '.DS' not in f:
                cpg_out.append(f)

        self.assertEqual(cpg_in[0], cpg_out[0])

    def test_build_meth_df(self):
        self.api.minify_meth_files(self.data_dir, self.tmp_dir)
        # Check if we haven't annotated our files we have an exception
        with self.assertRaises(APIException):
            self.api.build_meth_df(self.tmp_dir, join_id='id')

        # Build annotation and then build the df again
        self.api.build_annotation()
        self.api.build_meth_df(self.tmp_dir, join_id='id')
        meth_df = self.api.get_meth_df()

        # build rnaseq df first then join the meth df
        self.api.build_rna_df(self.data_dir)
        df = self.api.get_rna_df()

        # Lets get annotation information from biomart, first we'll get the gene ids from the rnaseq df
        bm = Biomart()
        df = bm.build_add_metadata(df, 'hsapiens_gene_ensembl', 'id', self.tmp_dir)
        # Now we want to add our methylation data to our dataframe, here we're being very strict and dropping any null
        # rows
        self.api.build_meth_df(self.tmp_dir, df, drop_empty_rows=True)
        meth_df = self.api.get_meth_df()

        # Let's now run some checks
        self.assertEqual(len(meth_df), 5)
        self.assertEqual(meth_df.values[0][2], 'TSPAN6')
        self.assertEqual(meth_df.values[0][1], 2881)
        self.assertEqual(meth_df.values[0][-1], 0.0617698530665544)

    def test_merge_rna_meth_values(self):
        # Run almost the same as above
        self.api.minify_meth_files(self.data_dir, self.tmp_dir)
        # Build annotation and then build the df again
        self.api.build_annotation()
        self.api.build_meth_df(self.tmp_dir, join_id='id')
        meth_df = self.api.get_meth_df()

        # build rnaseq df first then join the meth df
        self.api.build_rna_df(self.data_dir)
        df = self.api.get_rna_df()

        # Lets get annotation information from biomart, first we'll get the gene ids from the rnaseq df
        bm = Biomart()
        df = bm.build_add_metadata(df, 'hsapiens_gene_ensembl', 'id', self.tmp_dir)
        merged_df = self.api.merge_rna_meth_values(meth_df, df, index_col='id', merge_col='gene_id')

        self.assertEqual(merged_df.values[1][2], 'TNMD')
        self.assertEqual(merged_df.values[1][1], 6)

        self.assertEqual(merged_df.values[-1][2], 'CYP51A1')
        self.assertEqual(merged_df.values[-1][1], 132)
        self.assertEqual(merged_df.values[-1][-1], 0.935303719030208)

        # Note we have much more rows in this merge because we don't drop the NA cols as we do in the prev test
        self.assertEqual(len(merged_df), 18)

    def test_get_mutation_values_on_filter(self):
        filter_col = 'ssm.consequence.0.transcript.gene.symbol'
        gene_ids = ['WASL', 'CDCP2', 'BHLHE40']
        # First check if the exception is raised
        with self.assertRaises(AnnotateException):
            self.api.get_mutation_values_on_filter('ssm.consequence.0.transcript.gene.gene_id', gene_ids, filter_col)

        # This will tell us we need to first download and create our mutation data frame
        self.api.build_mutation_df(self.data_dir)

        gene_ids = self.api.get_mutation_values_on_filter('ssm.consequence.0.transcript.gene.gene_id', gene_ids, filter_col)
        gene_ids.sort()
        #['ENSG00000173114', 'ENSG00000115414', 'ENSG00000157184'] --> these are the IDs that should be in the mutation df
        self.assertEqual(gene_ids[0], 'ENSG00000115414')

        changes = ['Small deletion']
        gene_changes = self.api.get_mutation_values_on_filter(filter_col, changes, 'ssm.mutation_subtype')
        self.assertEqual(gene_changes[0], 'FOLH1B')
        self.assertEqual(len(gene_changes), 1)

    def test_get_genes_with_mutations(self):
        # This will tell us we need to first download and create our mutation data frame
        self.api.build_mutation_df(self.data_dir)
        genes = self.api.get_genes_with_mutations()

        # Now lets get them just for one case
        genes_case = self.api.get_genes_with_mutations(['TCGA-A3-3308'])

        self.assertEqual(len(genes), 11)
        self.assertEqual(len(genes_case), 6)

    def test_get_values_from_df(self):
        # build rnaseq df first then join the meth df
        self.api.build_annotation()
        self.api.build_rna_df(self.data_dir)
        df = self.api.get_rna_df()

        values, columns, df = self.api.get_values_from_df(df, 'id', ['TCGA-A3-3308'])
        self.assertEqual(len(columns), 2)
        self.assertEqual(values[0][1], 2881)

        # Lets try just getting some gene IDs of interest
        values, columns, df = self.api.get_values_from_df(df, 'id', ['TCGA-A3-3308'], ['ENSG00000000938.11', 'ENSG00000000005.5'])
        print(df)
        self.assertEqual(len(df), 2)

        # Lets also check if we do a oclumn must include
        values, columns, df = self.api.get_values_from_df(df, 'id', ['TCGA-A3-3308'], None, column_name_includes=['Solid'])
        print(df)
        self.assertEqual(len(df), 2)

    def test_get_cases_with_meta(self):
        # Now we need to run the annotation building (thats what the previous API exception was for
        metadata = {'race': 'white'}
        with self.assertRaises(APIException):
            self.api.get_files_with_meta(metadata, "any")
        self.api.build_annotation()

        # Test first with badly formated metadata
        with self.assertRaises(APIException):
            self.api.get_files_with_meta(metadata, "any")

        # Now we want to check the functionality we expect
        metadata = {'race': ['asian']}

        cases = self.api.get_files_with_meta(metadata, "any")
        cases.sort()
        # self.assertEqual(cases[-1], 'TCGA-CZ-5989')
        # self.assertEqual(cases[0], 'TCGA-A3-3308')
        # self.assertEqual(len(cases), 2)

        female_stage1 = self.api.get_files_with_meta({'gender': ['female'], 'tumor_stage_num': [1, 2],
                                                      'project_id': ['TCGA-KIRP']})
        self.assertEqual(len(female_stage1), 0)
        male_cases = self.api.get_files_with_meta({'gender': ['male'], 'tumor_stage_num': [1, 2, 3, 4],
                                                      'project_id': ['TCGA-KIRP']})
        self.assertEqual(len(male_cases), 1)
