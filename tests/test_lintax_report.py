import importlib.util
from pathlib import Path
import tempfile
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'bin' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lin = load('report_lin', 'report-lin.py')
report = load('global_report', 'global_report.py')


class LINtaxTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.groups = self.write('groups.tsv',
            'LINgroup_Name\tLINgroup_prefix\nA\t0,1\nB\t0,10\nAbsent\t0,9\n')
        # Put 0,10 first to expose the old raw startswith("0,1") bug.
        self.taxonomy = self.write('data.txt',
            'LIN\ttaxid_LIN\n0,10,0\t1,31,32\n0,1,0\t1,21,22\n')
        self.kreport = self.write('sample.kreport',
            '20\t1\t1\tU\t0\tunclassified\n'
            '80\t4\t1\tR\t1\troot\n'
            '40\t2\t1\tS\t21\t  A\n'
            '20\t1\t1\tS1\t22\t    strain A\n'
            '20\t1\t0\tS\t31\t  B\n'
            '20\t1\t1\tS1\t32\t    strain B\n')
        self.koutput = self.write('sample.koutput',
            'C\tr1\t1\t5|5\t1:1\n'
            'C\tr2\t21\t10|20\t21:1\n'
            'C\tr3\t22\t20|30\t22:1\n'
            'C\tr4\t32\t15|25\t32:1\n'
            'U\tr5\t0\t10|10\t0:1\n')

    def write(self, name, text):
        path = self.work / name
        path.write_text(text)
        return path

    def make(self):
        return lin.make_report(self.groups, self.taxonomy, self.kreport, self.koutput).set_index('LINgroup_Name')

    def test_component_prefix_matching_and_descendant_lengths(self):
        frame = self.make()
        self.assertEqual(frame.loc['A', 'Assigned_reads'], 2)
        self.assertEqual(frame.loc['A', 'Total_reads_length'], 80)
        self.assertEqual(frame.loc['B', 'Assigned_reads'], 1)
        self.assertEqual(frame.loc['B', 'Total_reads_length'], 40)
        self.assertEqual(frame.loc['A', 'Percentage_assigned_reads'], 50)

    def test_direct_counts_are_distinct_from_clade_counts(self):
        frame = self.make()
        self.assertEqual(frame.loc['A', 'Unique_Assigned_reads'], 1)
        self.assertEqual(frame.loc['B', 'Unique_Assigned_reads'], 0)
        self.assertEqual(frame.loc['B', 'Assigned_reads'], 1)

    def test_totals_and_missing_groups_are_filled(self):
        frame = self.make()
        total = frame.loc['Total_reads']
        self.assertEqual(total['Assigned_reads'], 4)
        self.assertEqual(total['Unique_Assigned_reads'], 1)
        self.assertEqual(total['Total_reads_length'], 130)
        self.assertFalse(total[['Assigned_reads', 'Unique_Assigned_reads', 'Total_reads_length']].isna().any())
        self.assertEqual(frame.loc['Absent', 'Assigned_reads'], 0)
        self.assertEqual(frame.loc['Absent', 'LINtax_note'], 'group_prefix_missing_from_taxonomy')

    def test_total_row_never_counts_as_a_hit(self):
        path = self.work / 'LINreport.txt'
        self.make().reset_index().to_csv(path, index=False)
        hits = report.read_lintax(str(path))
        self.assertEqual(set(hits.Raw_hit), {'A', 'B'})

    def test_empty_kraken_files(self):
        self.kreport.write_text('')
        self.koutput.write_text('')
        frame = self.make()
        self.assertTrue((frame.Assigned_reads == 0).all())
        self.assertTrue((frame.Total_reads_length == 0).all())

    def test_bracketed_taxids(self):
        self.assertEqual(lin.taxid_tokens('[1, 21, 22]'), (1, 21, 22))
        self.assertEqual(lin.taxid_tokens('1,21,22'), (1, 21, 22))

    def test_inconsistent_taxonomy_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Inconsistent taxonomy'):
            lin.build_taxonomy(pd.DataFrame({'LIN': ['0,1', '0,1'],
                                            'taxid_LIN': ['1,21', '1,31']}))


if __name__ == '__main__':
    unittest.main()
