"""Regression checks for taxonomy parsing; no external dependencies required."""
import importlib.util
from pathlib import Path
import unittest


spec = importlib.util.spec_from_file_location(
    'extract_kraken_reads',
    Path(__file__).resolve().parents[1] / 'bin' / 'extract_kraken_reads.py',
)
extract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extract)


def row(taxid, depth, rank='S'):
    return f'100.00\t10\t5\t{rank}\t{taxid}\t{"  " * depth}taxon\n'


def descendants(node):
    return {child.taxid for child in node.children}.union(
        *(descendants(child) for child in node.children)
    )


class ReportTreeTests(unittest.TestCase):
    def test_custom_root_without_taxid_one(self):
        nodes = extract.parse_report_tree(
            [row(0, 0, 'U'), row(1280, 0), row(9001, 1, '-')], [1280]
        )
        self.assertIsNone(nodes[1280].parent)
        self.assertEqual(descendants(nodes[1280]), {9001})

    def test_standard_root_and_sibling_branch(self):
        nodes = extract.parse_report_tree(
            [row(1, 0, 'R'), row(2, 1, 'D'), row(1280, 2),
             row(9001, 3, '-'), row(9002, 4, '-'), row(9999, 2)],
            [1, 1280, 9999],
        )
        self.assertEqual(descendants(nodes[1280]), {9001, 9002})
        self.assertEqual(nodes[1280].parent.taxid, 2)
        self.assertEqual(nodes[9999].parent.taxid, 2)
        self.assertEqual(descendants(nodes[1]), {2, 1280, 9001, 9002, 9999})

    def test_independent_roots_and_skipped_depths(self):
        nodes = extract.parse_report_tree(
            [row(1280, 0), row(9001, 3, '-'), row(9999, 0), row(9998, 1)],
            [1280, 9999],
        )
        self.assertEqual(descendants(nodes[1280]), {9001})
        self.assertIsNone(nodes[9999].parent)
        self.assertEqual(descendants(nodes[9999]), {9998})

    def test_absent_target_or_unclassified_only(self):
        for lines in ([], [row(0, 0, 'U')], [row(9999, 0)]):
            with self.subTest(lines=lines):
                self.assertEqual(extract.parse_report_tree(lines, [1280]), {})

    def test_krakenuniq_report(self):
        nodes = extract.parse_report_tree(
            ['100\t10\t5\t100\t20\t1280\tspecies\tS. aureus\n',
             '100\t10\t5\t100\t20\t9001\tstrain\t  strain A\n'],
            [1280],
        )
        self.assertEqual(descendants(nodes[1280]), {9001})


if __name__ == '__main__':
    unittest.main()
