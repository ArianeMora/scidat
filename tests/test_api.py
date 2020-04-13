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
from sciutil import Biomart

import os
import shutil
import tempfile
import unittest


class TestAPI(unittest.TestCase):
    def setUp(self):
        # Flag to set data to be local so we don't have to download them repeatedly. ToDo: Remove when publishing.
        self.local = True
        if self.local:
            self.tmp_dir = '../tests/data/tmp/'
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')

        # def __init__(self, manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
        #                  max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_'):

        self.data_dir = '../tests/data/'
        gdc_client = self.data_dir + './gdc-client'
        meta_dir = '../tests/data/'
        clinical_file = meta_dir + 'clinical.txt'
        sample_file = meta_dir + 'sample_sheet.txt'
        manifest_file = meta_dir + 'manifest.tsv'

        self.api = API(manifest_file, gdc_client, clinical_file, sample_file, download_dir=self.tmp_dir,
                            max_cnt=1, requires_lst=['counts', 'm450'])

    def tearDown(self):
        if not self.local:
            # Delete temp dir
            shutil.rmtree(self.tmp_dir)

    def test_download_mutation_data(self):
        # Here we want to download the mutation files
        case_ids = ['TCGA-A3-3308', 'TCGA-KN-8422', 'TCGA-CZ-5989-11A', 'TCGA-A4-8312-01A']
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

    def test_build_rna_df(self):
        with self.assertRaises(APIException):
            self.api.build_rna_df()

        # Now we want to first build the annotation
        self.api.build_annotation()
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
        #['ENSG00000173114', 'ENSG00000115414', 'ENSG00000157184'] --> these are the IDs that should be in the mutation df
        self.assertEqual(gene_ids[0], 'ENSG00000157184')

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
        values, columns = self.api.get_values_from_df(df, ['TCGA-A3-3308'])
        self.assertEqual(len(columns), 1)
        self.assertEqual(values[0][0], 2881)

    def test_get_cases_with_meta(self):
        # Now we need to run the annotation building (thats what the previous API exception was for
        metadata = {'race': 'white'}
        with self.assertRaises(APIException):
            self.api.get_cases_with_meta(metadata, "any")
        self.api.build_annotation()

        # Test first with badly formated metadata
        with self.assertRaises(APIException):
            self.api.get_cases_with_meta(metadata, "any")

        # Now we want to check the functionality we expect
        metadata = {'race': ['asian']}

        cases = self.api.get_cases_with_meta(metadata, "any")
        cases.sort()
        self.assertEqual(cases[-1], 'TCGA-CZ-5989')
        self.assertEqual(cases[0], 'TCGA-A3-3308')
        self.assertEqual(len(cases), 2)