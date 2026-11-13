import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('global_report', ROOT / 'bin/global_report.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def write(self, name, text):
        path = self.work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def test_classified_counts_do_not_double_count_clades(self):
        path = self.write('firstpass.kreport',
            '20\t20\t20\tU\t0\tunclassified\n'
            '80\t100\t5\tR\t1\troot\n'
            '76\t95\t15\tS\t1280\t  S. aureus\n'
            '64\t80\t80\tS1\t9001\t    strain\n')
        self.assertEqual(report.parse_kraken_report(str(path)),
                         {'Classified_reads': 100, 'Kraken_note': ''})

    def test_custom_root_no_target_taxid_and_zero_classified(self):
        path = self.write('custom.kreport',
            '100\t100\t10\tR\t999\tcustom root\n'
            '90\t90\t90\tS\t9001\t  strain\n')
        self.assertEqual(report.parse_kraken_report(str(path))['Classified_reads'], 100)
        path.write_text('100\t100\t100\tU\t0\tunclassified\n')
        self.assertEqual(report.parse_kraken_report(str(path)),
                         {'Classified_reads': 0, 'Kraken_note': ''})

    def test_invalid_report_has_diagnostic(self):
        path = self.write('invalid.kreport', 'invalid\n')
        self.assertEqual(report.parse_kraken_report(str(path))['Kraken_note'],
                         'kraken_report_no_valid_rows')

    def test_placeholders_are_not_strain_hits(self):
        for contents in ('strain\tgcov\n', 'strain\tgcov\nNA\t0\n'):
            path = self.write('placeholder.tsv', contents)
            self.assertTrue(report.read_strainge(str(path)).empty)

    def test_strainscan_depth_column_variants(self):
        for column in ('Predicted_Depth', 'Predicted_Depth (Enet)'):
            self.write('strainscan_output/final_report.txt',
                       f'Strain_Name\t{column}\nDRX014091_S99\t11.31\n')
            rows = report.read_strainscan(str(self.work / 'strainscan_output'))
            self.assertEqual(rows.iloc[0]['Evidence_value'], 11.31)

    def mock_tool(self, name, source):
        path = self.write('mock_bin/' + name, f'#!{sys.executable}\n' + source)
        path.chmod(0o755)

    def run_wrapper(self, script, *args, **env):
        return subprocess.run(['bash', str(ROOT / 'bin' / script), *map(str, args)],
                              cwd=self.work, env=dict(os.environ,
                                  PATH=str(self.work / 'mock_bin') + ':' + os.environ['PATH'],
                                  **env), capture_output=True, text=True, check=True)

    def test_strainscan_report_and_error_statuses(self):
        self.mock_tool('strainscan', '''import os, sys
from pathlib import Path
args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
os.chdir('/')  # Emulate a launcher that changes directory.
assert all(Path(args[k]).is_absolute() for k in ('-i', '-j', '-d', '-o'))
assert Path(args['-i']).read_text() == 'reads'
assert Path(args['-j']).read_text() == 'reads'
if os.environ['MOCK_REPORT'] == 'yes':
    Path(args['-o'], 'final_report.txt').write_text('Strain_Name\\tPredicted_Depth\\nDRX014091_S99\\t11.31\\n')
print('diagnostic stdout')
print('diagnostic stderr', file=sys.stderr)
sys.exit(int(os.environ['MOCK_EXIT']))
''')
        self.write('R1.fastq', 'reads')
        self.write('R2.fastq', 'reads')
        for exit_code, has_report, expected in (
                ('0', 'yes', 'success'), ('2', 'yes', 'report_with_error'),
                ('2', 'no', 'failed'), ('0', 'no', 'no_report')):
            with self.subTest(expected=expected):
                old_report = self.work / 'strainscan_output/final_report.txt'
                if old_report.exists():
                    old_report.unlink()
                self.run_wrapper('run_strainscan.sh', 'sample.S99.5', 'R1.fastq',
                                 'R2.fastq', self.work / 'db',
                                 MOCK_EXIT=exit_code, MOCK_REPORT=has_report)
                with (self.work / 'strainscan_output/strainscan_status.tsv').open() as f:
                    status = next(csv.DictReader(f, delimiter='\t'))
                self.assertEqual(status['status'], expected)
                self.assertIn('diagnostic stderr',
                              (self.work / 'strainscan_output/strainscan.stderr.log').read_text())
                if has_report == 'yes':
                    self.assertEqual(len(report.read_strainscan(str(old_report.parent))), 1)

    def test_strainge_dotted_sample_and_failure_placeholders(self):
        self.mock_tool('straingst', '''import os, sys
from pathlib import Path
if os.environ['MOCK_EXIT'] != '0':
    print('database mismatch', file=sys.stderr)
    sys.exit(int(os.environ['MOCK_EXIT']))
prefix = Path(sys.argv[sys.argv.index('-o') + 1])
prefix.with_suffix('.strains.tsv').write_text('strain\\tgcov\\nDRX014091_S99\\t0.99\\n')
prefix.with_suffix('.stats.tsv').write_text('sample\\tpkmers\\noriginal\\t100\\n')
''')
        self.write('sample.hdf5', 'hdf5')
        sample = 'DRX014091_S99.5_1_ST5476_CCNA'
        self.run_wrapper('run_strainge.sh', sample, 'sample.hdf5', 'db.hdf5', MOCK_EXIT='0')
        output = self.work / (sample + '.strains.tsv')
        self.assertEqual(report.read_strainge(str(output)).iloc[0]['Raw_hit'], 'DRX014091_S99')
        self.assertTrue((self.work / (sample + '.stats.tsv')).exists())
        self.run_wrapper('run_strainge.sh', sample, 'sample.hdf5', 'db.hdf5', MOCK_EXIT='1')
        self.assertTrue(report.read_strainge(str(output)).empty)
        self.assertIn('straingst_failed_exit_1',
                      (self.work / (sample + '.stats.tsv')).read_text())

    def test_summary_schema_and_real_tool_counts(self):
        kraken = self.write('firstpass.kreport', '100\t100\t100\tS\t9001\tstrain\n')
        strainge = self.write('sample.strains.tsv', 'strain\tgcov\nNA\t0\n')
        self.write('strainscan_output/final_report.txt',
                   'Strain_Name\tPredicted_Depth\nDRX014091_S99\t11.31\n')
        summary, evidence = self.work / 'summary.csv', self.work / 'evidence.tsv'
        result = subprocess.run([sys.executable, str(ROOT / 'bin/global_report.py'),
                        '-kr', str(kraken), '-st', str(strainge),
                        '-sc', str(self.work / 'strainscan_output'),
                        '-s', 'sample.S99.5', '-o', str(summary), '-e', str(evidence)],
                       capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with summary.open() as f:
            reader = csv.DictReader(f)
            self.assertEqual(reader.fieldnames, [
                'Sample', 'Classified_reads', 'Kraken_note', 'Tools_reporting_n',
                'Tools_reporting', 'Tool_hit_counts', 'Total_tool_hits', 'Notes'])
            row = next(reader)
        self.assertEqual(row['Classified_reads'], '100')
        self.assertEqual(row['Tools_reporting'], 'StrainScan')
        self.assertEqual(row['Total_tool_hits'], '1')
        self.assertEqual(row['Kraken_note'], '-')


if __name__ == '__main__':
    unittest.main()
