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

import os
import pandas as pd
import shutil
import tempfile
import unittest

from scidat.annotate import AnnotateException
from scidat import Annotate


class TestAnnotate(unittest.TestCase):

    def setUp(self):
        # Setup temp dir
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.local = False
        if self.local:
            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')

        meta_dir = os.path.join(THIS_DIR, 'data/')
        clinical_file = meta_dir + 'clinical.txt'
        sample_file = meta_dir + 'sample_sheet.txt'
        manifest_file = meta_dir + 'manifest.tsv'
        self.meta_dir = meta_dir
        # (self, output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
        #                  sep='_', mutations_file=None, clin_cols=None
        self.annotator = Annotate(self.tmp_dir, clinical_file, sample_file, manifest_file,
                                  ['count', 'm450'])

    def tearDown(self):
        if not self.local:
            # Delete temp dr
            shutil.rmtree(self.tmp_dir)

    def test_annotate(self):

        self.annotator.build_annotation()

        # Now we want to check that the dataframe
        self.annotator.save_annotation(self.tmp_dir, 'out')

        # Check file exists
        self.assertEqual(os.path.exists(self.annotator.u.generate_label([self.tmp_dir, 'out'], '.csv')), True)

        # Check our annotation has what we expect
        print(self.annotator.annotated_file_dict)

    def test_build_mutation_df(self):
        mutation_dir = self.tmp_dir
        with self.assertRaises(AnnotateException):
            self.annotator.build_mutation_df(mutation_dir, 'cakes')

        mutation_dir = self.meta_dir
        # Now use the proper dir.

        self.annotator.build_mutation_df(mutation_dir)
        mutation_df = self.annotator.get_mutation_df()
        self.assertListEqual(['TCGA-KN-8422', 1.0, 'SIGLEC14', 'ENSG00000163125', 'chr3:g.33544725G>A', 'Single base substitution',
                'intron_variant,upstream_gene_variant,missense_variant,downstream_gene_variant', '0,1,3', 'Q379P',
                '2,3,4,5,8,9'], list(mutation_df.values[0]))

        self.assertEqual(len(mutation_df), 11)
        self.assertListEqual(['TCGA-A3-3308', 6.0, 'CDCP2', 'ENSG00000115414', 'chr1:g.160173678C>G', 'Single base substitution',
                         'splice_region_variant,missense_variant,missense_variant,missense_variant,missense_variant',
                        '0,1,12,23,29,30', 'W134Lfs*8', '1,2,3,10,15,17,18,19,20'], list(mutation_df.values[-1]))

    def test_get_mutation_df(self):
        # Check we get a raise error
        with self.assertRaises(AnnotateException):
            self.annotator.get_mutation_df()

        mutation_dir = self.meta_dir
        # Now use the proper dir.
        self.annotator.build_mutation_df(mutation_dir)
        mutation_df = self.annotator.get_mutation_df()
        self.assertListEqual(
            ['TCGA-KN-8422', 1.0, 'SIGLEC14', 'ENSG00000163125', 'chr3:g.33544725G>A', 'Single base substitution',
             'intron_variant,upstream_gene_variant,missense_variant,downstream_gene_variant', '0,1,3', 'Q379P',
             '2,3,4,5,8,9'], list(mutation_df.values[0]))

    def test_save_mutation_df(self):
        mutation_dir = self.meta_dir
        # Now use the proper dir.
        self.annotator.build_mutation_df(mutation_dir)
        file_path = self.annotator.u.generate_label([self.tmp_dir, "mut_df"], '.tsv')
        self.annotator.save_mutation_df(self.tmp_dir, "mut_df")
        self.assertEqual(os.path.exists(file_path), True)

    def test_save_annotation(self):
        self.annotator.build_annotation()
        file_path = self.annotator.u.generate_label([self.tmp_dir, "annot_df"], '.csv')
        self.annotator.save_annotation(self.tmp_dir, "annot_df")
        self.assertEqual(os.path.exists(file_path), True)

    def test_save_annotated_clinical_df(self):
        self.annotator.build_annotation()
        file_path = self.annotator.u.generate_label([self.tmp_dir, "clin_df"], '.csv')
        self.annotator.save_annotated_clinical_df(self.tmp_dir, "clin_df")
        self.assertEqual(os.path.exists(file_path), True)

    def test_get_cases(self):

        self.annotator.build_annotation()

        cases = self.annotator.get_cases()
        cases.sort()
        self.assertEqual(cases[0], "TCGA-A3-3308")
        self.assertEqual(cases[-1], "TCGA-KN-8422")
        self.assertEqual(len(cases), 4)

