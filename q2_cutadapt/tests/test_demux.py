# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import os
import pathlib
import shutil
import tempfile
import unittest

import pandas as pd

from q2_cutadapt._demux import (_build_demux_command, _rename_files,
                                _write_barcode_fasta,
                                _write_empty_fastq_to_mux_barcode_in_seq_fmt)
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequenceDirFmt,
    MultiplexedPairedEndBarcodeInSequenceDirFmt)
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat)
from qiime2 import Artifact, CategoricalMetadataColumn
from qiime2.util import redirected_stdio
from qiime2.plugin.testing import TestPluginBase


class TestDemuxSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def assert_demux_results(self, exp_samples_and_barcodes, exp_results,
                             obs_demuxed_art):
        obs_demuxed = obs_demuxed_art.view(
            SingleLanePerSampleSingleEndFastqDirFmt)
        obs_demuxed_seqs = obs_demuxed.sequences.iter_views(FastqGzFormat)
        zipped = itertools.zip_longest(exp_samples_and_barcodes.iteritems(),
                                       exp_results, obs_demuxed_seqs)
        for (sample_id, barcode), exp, (filename, fmt) in zipped:
            filename = str(filename)
            self.assertTrue(sample_id in filename)
            self.assertTrue(barcode in filename)
            with gzip.open(str(fmt), 'rt') as fh:
                obs = ''.join(fh.readlines())
            self.assertEqual(exp, obs)

    def assert_untrimmed_results(self, exp, obs_untrimmed_art):
        obs_untrimmed = obs_untrimmed_art.view(
            MultiplexedSingleEndBarcodeInSequenceDirFmt)
        obs_untrimmed = obs_untrimmed.file.view(FastqGzFormat)
        with gzip.open(str(obs_untrimmed), 'rt') as fh:
            obs = ''.join(fh.readlines())
        self.assertEqual(exp, obs)

    def setUp(self):
        super().setUp()
        self.demux_single_fn = self.plugin.methods['demux_single']

        muxed_sequences_fp = self.get_data_path('forward.fastq.gz')
        self.muxed_sequences = Artifact.import_data(
            'MultiplexedSingleEndBarcodeInSequence', muxed_sequences_fp)

    def test_typical(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                                      obs_untrimmed_art)

    def test_all_matched(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'GGGG'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c'],
                                     name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample c
            '@id6\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('', obs_untrimmed_art)

    # NOTE: this test used to check for an exception because it was
    # possible to generate a completely empty output dir with no fastq.gz
    # files in it. As of cutadapt 3 this is no longer possible, because output
    # files are generated for every sample (and we must specify at least one
    # sample in order for the barcodes to be valid QIIME 2 Metadata). Rather
    # than remove the test, we will retool it here.
    def test_none_matched(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['TTTT'], name='Barcode',
                      index=pd.Index(['sample_d'], name='id')))

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), [''], obs_demuxed_art)
        self.assert_untrimmed_results('@id1\nAAAAACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id2\nCCCCACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id3\nAAAAACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id4\nCCCCACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id5\nCCCCACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                                      obs_untrimmed_art)

    def test_error_tolerance_filtering(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAG', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a has no reads (bc we misspelled the barcode)
            '',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp,
                                  obs_demuxed_art)
        self.assert_untrimmed_results('@id1\nAAAAACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id3\nAAAAACGTACGT\n+\nzzzzzzzzzzzz\n'
                                      '@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                                      obs_untrimmed_art)

    def test_error_tolerance_high_enough_to_prevent_filtering(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAG', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata,
                                     error_rate=0.25)

        # This test should yield the same results as test_typical, above
        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                                      obs_untrimmed_art)

    def test_extra_barcode_in_metadata(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'GGGG', 'TTTT'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c',
                                      'sample_d'], name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample c
            '@id6\nACGTACGT\n+\nzzzzzzzz\n',
            # sample d is empty bc no reads matched the barcode TTTT
            '', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        exp_samples_and_barcodes = pd.Series(['AAAA', 'CCCC', 'GGGG', 'TTTT'],
                                             index=['sample_a', 'sample_b',
                                                    'sample_c', 'sample_d'])
        self.assert_demux_results(exp_samples_and_barcodes, exp,
                                  obs_demuxed_art)
        self.assert_untrimmed_results('', obs_untrimmed_art)

    def test_variable_length_barcodes(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAAA', 'CCCCCC', 'GGGG'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c'],
                                     name='id')))
        muxed_sequences_fp = self.get_data_path('variable_length.fastq.gz')
        muxed_sequences = Artifact.import_data(
            'MultiplexedSingleEndBarcodeInSequence', muxed_sequences_fp)
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample c
            '@id6\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('', obs_untrimmed_art)

    def test_batch_size(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata,
                                     batch_size=1)

        # This test should yield the same results as test_typical, above,
        # the fact that we are batching shouldn't impact the final results
        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                                      obs_untrimmed_art)

    def test_invalid_batch_size(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))

        with self.assertRaisesRegex(ValueError, '5.*cannot be greater.*2'):
            self.demux_single_fn(self.muxed_sequences, metadata, batch_size=5)

    def test_batch_size_odd_number_of_samples(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'GGGG'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c'],
                      name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample c
            '@id6\nACGTACGT\n+\nzzzzzzzz\n', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata,
                                     batch_size=2)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('', obs_untrimmed_art)

    def test_min_length(self):
        metadata = CategoricalMetadataColumn(
            # The third barcode is meant to completely remove the only GGGG
            # coded sequence
            pd.Series(['AAAA', 'CCCC', 'GGGGACGTACGT'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c'],
                      name='id')))
        exp = [
            # sample a
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample c is empty because the barcode matched the entire
            # read, which removed everything.
            '', ]

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_single_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results('', obs_untrimmed_art)


class TestDemuxPaired(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def assert_demux_results(self, exp_samples_and_barcodes, exp_results,
                             obs_demuxed_art):
        obs_demuxed = obs_demuxed_art.view(
            SingleLanePerSamplePairedEndFastqDirFmt)
        obs_demuxed_seqs = obs_demuxed.sequences.iter_views(FastqGzFormat)
        # Since we are working with fwd/rev reads, duplicate each list elem
        exp = [x for x in exp_samples_and_barcodes.iteritems() for _ in (0, 1)]
        zipped = itertools.zip_longest(exp, exp_results, obs_demuxed_seqs)
        for (sample_id, barcode), exp, (filename, fmt) in zipped:
            filename = str(filename)
            self.assertTrue(sample_id in filename)
            self.assertTrue(barcode in filename)
            with gzip.open(str(fmt), 'rt') as fh:
                obs = ''.join(fh.readlines())
            self.assertEqual(exp, obs)

    def assert_untrimmed_results(self, exp, obs_untrimmed_art):
        obs_untrimmed = obs_untrimmed_art.view(
            MultiplexedPairedEndBarcodeInSequenceDirFmt)

        # first check the fwd reads
        obs_untrimmed_f = obs_untrimmed.forward_sequences.view(FastqGzFormat)
        with gzip.open(str(obs_untrimmed_f), 'rt') as fh:
            obs_untrimmed_f = ''.join(fh.readlines())
        self.assertEqual(exp[0], obs_untrimmed_f)

        # next check the rev reads
        obs_untrimmed_r = obs_untrimmed.reverse_sequences.view(FastqGzFormat)
        with gzip.open(str(obs_untrimmed_r), 'rt') as fh:
            obs_untrimmed_r = ''.join(fh.readlines())
        self.assertEqual(exp[1], obs_untrimmed_r)

    def setUp(self):
        super().setUp()
        self.demux_paired_fn = self.plugin.methods['demux_paired']

        muxed_sequences_f_fp = self.get_data_path('forward.fastq.gz')
        muxed_sequences_r_fp = self.get_data_path('reverse.fastq.gz')

        with tempfile.TemporaryDirectory() as temp:
            shutil.copy(muxed_sequences_f_fp, temp)
            shutil.copy(muxed_sequences_r_fp, temp)
            self.muxed_sequences = Artifact.import_data(
                'MultiplexedPairedEndBarcodeInSequence', temp)

    # Just one proof-of-concept test here - the single-end test suite
    # covers the edge cases.
    def test_typical(self):
        metadata = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='Barcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a, fwd
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample a, rev
            '@id1\nGGGGTGCATGCA\n+\nzzzzzzzzzzzz\n'
            '@id3\nGGGGTGCATGCA\n+\nzzzzzzzzzzzz\n',
            # sample b, fwd
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample b, fwd
            '@id2\nTTTTTGCATGCA\n+\nzzzzzzzzzzzz\n'
            '@id4\nTTTTTGCATGCA\n+\nzzzzzzzzzzzz\n'
            '@id5\nTTTTTGCATGCA\n+\nzzzzzzzzzzzz\n', ]
        exp_untrimmed = ['@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                         '@id6\nTTTTTGCATGCA\n+\nzzzzzzzzzzzz\n']

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_paired_fn(self.muxed_sequences, metadata)

        self.assert_demux_results(metadata.to_series(), exp, obs_demuxed_art)
        self.assert_untrimmed_results(exp_untrimmed, obs_untrimmed_art)

    def test_di_typical(self):
        forward_barcodes = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='ForwardBarcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        reverse_barcodes = CategoricalMetadataColumn(
            pd.Series(['GGGG', 'TTTT'], name='ReverseBarcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        exp = [
            # sample a, fwd
            '@id1\nACGTACGT\n+\nzzzzzzzz\n'
            '@id3\nACGTACGT\n+\nzzzzzzzz\n',
            # sample a, rev
            '@id1\nTGCATGCA\n+\nzzzzzzzz\n'
            '@id3\nTGCATGCA\n+\nzzzzzzzz\n',
            # sample b, fwd
            '@id2\nACGTACGT\n+\nzzzzzzzz\n'
            '@id4\nACGTACGT\n+\nzzzzzzzz\n'
            '@id5\nACGTACGT\n+\nzzzzzzzz\n',
            # sample a, fwd
            '@id2\nTGCATGCA\n+\nzzzzzzzz\n'
            '@id4\nTGCATGCA\n+\nzzzzzzzz\n'
            '@id5\nTGCATGCA\n+\nzzzzzzzz\n', ]
        exp_untrimmed = ['@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                         '@id6\nTTTTTGCATGCA\n+\nzzzzzzzzzzzz\n']

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_paired_fn(self.muxed_sequences,
                                     forward_barcodes=forward_barcodes,
                                     reverse_barcodes=reverse_barcodes)

        self.assert_demux_results(forward_barcodes.to_series(), exp,
                                  obs_demuxed_art)
        self.assert_untrimmed_results(exp_untrimmed, obs_untrimmed_art)

    def test_mixed_orientation_success(self):
        # sample_a and sample_b have reads in both fwd and rev directions.
        # sample_c only has reads in the fwd direction.
        # sample_d only has reads in the rev direction.
        forward_barcodes = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'GGGG', 'TTTT'], name='ForwardBarcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c',
                                      'sample_d'], name='id')))
        mixed_orientation_sequences_f_fp = self.get_data_path(
            'mixed-orientation/forward.fastq.gz')
        mixed_orientation_sequences_r_fp = self.get_data_path(
            'mixed-orientation/reverse.fastq.gz')
        with tempfile.TemporaryDirectory() as temp:
            shutil.copy(mixed_orientation_sequences_f_fp, temp)
            shutil.copy(mixed_orientation_sequences_r_fp, temp)
            mixed_orientation_sequences = Artifact.import_data(
                'MultiplexedPairedEndBarcodeInSequence', temp)

        with redirected_stdio(stderr=os.devnull):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_paired_fn(mixed_orientation_sequences,
                                     forward_barcodes=forward_barcodes,
                                     mixed_orientation=True)
        exp = [
            # sample_a fwd
            '@id1\nACGTACGT\n+\nyyyyyyyy\n' \
            '@id3\nACGTACGT\n+\nyyyyyyyy\n',
            # sample_a rev
            '@id1\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n' \
            '@id3\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n',
            # sample_b fwd
            '@id4\nACGTACGT\n+\nyyyyyyyy\n' \
            '@id2\nACGTACGT\n+\nyyyyyyyy\n',
            # sample_b rev
            '@id4\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n' \
            '@id2\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n',
            # sample_c fwd
            '@id5\nACGTACGT\n+\nyyyyyyyy\n',
            # sample_c rev
            '@id5\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n',
            # sample_d fwd
            '@id6\nACGTACGT\n+\nyyyyyyyy\n',
            # sample_d rev
            '@id6\nTGCATGCATGCA\n+\nzzzzzzzzzzzz\n', ]

        # We want to be sure that the validation is 100%, not just `min`,
        obs_demuxed_art.validate(level='max')
        # checkpoint assertion for the above `validate` - nothing should fail
        self.assertTrue(True)

        self.assert_demux_results(forward_barcodes.to_series(), exp,
                                  obs_demuxed_art)

        # Everything should match, so untrimmed should be empty
        self.assert_untrimmed_results(['', ''], obs_untrimmed_art)

    def test_di_mismatched_barcodes(self):
        forward_barcodes = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'ACGT'], name='ForwardBarcode',
                      index=pd.Index(['sample_a', 'sample_b', 'sample_c'],
                                     name='id')))
        reverse_barcodes = CategoricalMetadataColumn(
            pd.Series(['GGGG', 'TTTT'], name='ReverseBarcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))

        with self.assertRaisesRegex(ValueError, 'do not have.*sample_c'):
            self.demux_paired_fn(self.muxed_sequences,
                                 forward_barcodes=forward_barcodes,
                                 reverse_barcodes=reverse_barcodes)

    def test_di_duplicate_barcode_pairs(self):
        forward_barcodes = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC', 'AAAA', 'CCCC'], name='ForwardBarcode',
                      index=pd.Index(
                            ['sample_a', 'sample_b', 'sample_d', 'sample_c'],
                            name='id')))
        reverse_barcodes = CategoricalMetadataColumn(
            pd.Series(['GGGG', 'TTTT', 'GGGG', 'TTTT'], name='ReverseBarcode',
                      index=pd.Index(
                            ['sample_a', 'sample_b', 'sample_d', 'sample_c'],
                            name='id')))

        with self.assertRaisesRegex(
                ValueError, 'duplicate barcode.*sample_c.*sample_d'):
            self.demux_paired_fn(self.muxed_sequences,
                                 forward_barcodes=forward_barcodes,
                                 reverse_barcodes=reverse_barcodes)

    def test_multiple_orientations_dual_indices(self):
        forward_barcodes = CategoricalMetadataColumn(
            pd.Series(['AAAA', 'CCCC'], name='ForwardBarcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))
        reverse_barcodes = CategoricalMetadataColumn(
            pd.Series(['GGGG', 'TTTT'], name='ReverseBarcode',
                      index=pd.Index(['sample_a', 'sample_b'], name='id')))

        mixed_orientation_sequences_f_fp = self.get_data_path(
            'mixed-orientation/forward.fastq.gz')
        mixed_orientation_sequences_r_fp = self.get_data_path(
            'mixed-orientation/reverse.fastq.gz')

        # These files have forward and reverse reads mixed together in the same
        # file
        with tempfile.TemporaryDirectory() as temp:
            shutil.copy(mixed_orientation_sequences_f_fp, temp)
            shutil.copy(mixed_orientation_sequences_r_fp, temp)
            mixed_orientation_sequences = Artifact.import_data(
                'MultiplexedPairedEndBarcodeInSequence', temp)

        with self.assertRaisesRegex(ValueError,
                                    'Dual-indexed barcodes for mixed '
                                    'orientation reads are not supported.'):
            obs_demuxed_art, obs_untrimmed_art = \
                self.demux_paired_fn(mixed_orientation_sequences,
                                     forward_barcodes=forward_barcodes,
                                     reverse_barcodes=reverse_barcodes,
                                     mixed_orientation=True)


class TestDemuxUtilsSingleEnd(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.fastq_fp = \
            self.get_data_path('forward.fastq.gz')
        self.fastq = FastqGzFormat(self.fastq_fp, mode='r')
        self.seqs_dir_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt()
        self.seqs_dir_fmt.file.write_data(self.fastq, FastqGzFormat)
        self.barcode_series = pd.Series(['A', 'G'],
                                        index=['sample_a', 'sample_b'])
        self.per_sample_dir_fmt = SingleLanePerSampleSingleEndFastqDirFmt()
        self.untrimmed_dir_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt()

    def test_build_demux_command(self):
        with tempfile.NamedTemporaryFile() as barcode_fasta:
            obs = _build_demux_command(self.seqs_dir_fmt,
                                       {'fwd': barcode_fasta, 'rev': None},
                                       self.per_sample_dir_fmt,
                                       self.untrimmed_dir_fmt,
                                       0.1,
                                       2)
            self.assertTrue(barcode_fasta.name in obs[2])
        self.assertTrue('0.1' in obs[4])
        self.assertTrue('2' in obs[6])
        self.assertTrue(str(self.per_sample_dir_fmt) in obs[8])
        self.assertTrue(str(self.untrimmed_dir_fmt) in obs[10])
        self.assertEqual(str(self.seqs_dir_fmt.file.view(FastqGzFormat)),
                         obs[11])

    def test_rename_files_single(self):
        for fn in ['sample_a.1.fastq.gz', 'sample_b.1.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.seqs_dir_fmt, self.per_sample_dir_fmt,
                      self.barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        counter = 0
        for fn, (sample_id, barcode) in zip(seqs,
                                            self.barcode_series.iteritems()):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))
            counter += 1
        self.assertEqual(counter, 2)

    def test_rename_files_extra_samples_in_barcode_map(self):
        barcode_series = pd.Series(['A', 'G', 'C'],
                                   index=['sample_a', 'sample_b', 'sample_c'])

        for fn in ['sample_a.1.fastq.gz', 'sample_b.1.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.seqs_dir_fmt, self.per_sample_dir_fmt,
                      barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        counter = 0
        for fn, (sample_id, barcode) in zip(seqs, barcode_series.iteritems()):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))
            counter += 1
        self.assertEqual(counter, 2)

    def test_write_empty_fastq_to_mux_barcode_in_seq_fmt(self):
        _write_empty_fastq_to_mux_barcode_in_seq_fmt(self.untrimmed_dir_fmt)
        forward = self.untrimmed_dir_fmt.file.view(FastqGzFormat)
        self.assertTrue(forward.path.exists())
        self.assertTrue(forward.path.is_file())

    # Not really a single-end or paired-end specific test, but this is as
    # good a place as any to stash it
    def test_write_barcode_fasta(self):
        with tempfile.NamedTemporaryFile() as fh:
            _write_barcode_fasta(self.barcode_series, fh)
            fasta = open(fh.name).read()
            for (sample_id, barcode) in self.barcode_series.iteritems():
                self.assertTrue(sample_id in fasta)
                self.assertTrue(barcode in fasta)


class TestDemuxUtilsPairedEnd(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.fastq_fp = \
            self.get_data_path('forward.fastq.gz')
        self.fastq = FastqGzFormat(self.fastq_fp, mode='r')
        self.seqs_dir_fmt = MultiplexedPairedEndBarcodeInSequenceDirFmt()
        self.seqs_dir_fmt.forward_sequences.write_data(self.fastq,
                                                       FastqGzFormat)
        self.seqs_dir_fmt.reverse_sequences.write_data(self.fastq,
                                                       FastqGzFormat)
        self.barcode_series = pd.Series(['A', 'G'],
                                        index=['sample_a', 'sample_b'])
        self.per_sample_dir_fmt = SingleLanePerSampleSingleEndFastqDirFmt()
        self.untrimmed_dir_fmt = MultiplexedPairedEndBarcodeInSequenceDirFmt()

    def test_build_demux_command(self):
        with tempfile.NamedTemporaryFile() as barcode_fasta:
            obs = _build_demux_command(self.seqs_dir_fmt,
                                       {'fwd': barcode_fasta, 'rev': None},
                                       self.per_sample_dir_fmt,
                                       self.untrimmed_dir_fmt,
                                       0.1,
                                       2)
            self.assertTrue(barcode_fasta.name in obs[2])
        self.assertTrue('0.1' in obs[4])
        self.assertTrue('2' in obs[6])
        self.assertTrue(str(self.per_sample_dir_fmt) in obs[8])  # fwd
        self.assertTrue(str(self.per_sample_dir_fmt) in obs[12])  # rev
        self.assertTrue(str(self.untrimmed_dir_fmt) in obs[10])  # fwd
        self.assertTrue(str(self.untrimmed_dir_fmt) in obs[14])  # rev
        exp_f = str(self.seqs_dir_fmt.forward_sequences.view(FastqGzFormat))
        self.assertEqual(exp_f, obs[15])
        exp_r = str(self.seqs_dir_fmt.reverse_sequences.view(FastqGzFormat))
        self.assertEqual(exp_r, obs[16])

    def test_build_di_demux_command(self):
        with tempfile.NamedTemporaryFile() as barcode_fasta_f:
            with tempfile.NamedTemporaryFile() as barcode_fasta_r:
                obs = _build_demux_command(self.seqs_dir_fmt,
                                           {'fwd': barcode_fasta_f,
                                            'rev': barcode_fasta_r},
                                           self.per_sample_dir_fmt,
                                           self.untrimmed_dir_fmt,
                                           0.1,
                                           2)
                self.assertTrue(barcode_fasta_f.name in obs[2])
                self.assertTrue('--pair-adapters' == obs[11])
                self.assertTrue(barcode_fasta_r.name in obs[13])

    def test_rename_files(self):
        for fn in ['sample_a.1.fastq.gz', 'sample_a.2.fastq.gz',
                   'sample_b.1.fastq.gz', 'sample_b.2.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.seqs_dir_fmt, self.per_sample_dir_fmt,
                      self.barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        exp = [('sample_a', 'A'), ('sample_a', 'A'),
               ('sample_b', 'G'), ('sample_b', 'G')]
        counter = 0
        for fn, (sample_id, barcode) in zip(seqs, exp):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))
            counter += 1
        self.assertEqual(counter, 4)

    def test_rename_files_extra_samples_in_barcode_map(self):
        barcode_series = pd.Series(['A', 'G', 'C'],
                                   index=['sample_a', 'sample_b', 'sample_c'])

        for fn in ['sample_a.1.fastq.gz', 'sample_a.2.fastq.gz',
                   'sample_b.1.fastq.gz', 'sample_b.2.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.seqs_dir_fmt, self.per_sample_dir_fmt,
                      barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        exp = [('sample_a', 'A'), ('sample_a', 'A'),
               ('sample_b', 'G'), ('sample_b', 'G')]
        counter = 0
        for fn, (sample_id, barcode) in zip(seqs, exp):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))
            counter += 1
        self.assertEqual(counter, 4)

    def test_write_empty_fastq_to_mux_barcode_in_seq_fmt(self):
        _write_empty_fastq_to_mux_barcode_in_seq_fmt(self.untrimmed_dir_fmt)
        forward = self.untrimmed_dir_fmt.forward_sequences.view(FastqGzFormat)
        self.assertTrue(forward.path.exists())
        self.assertTrue(forward.path.is_file())

        reverse = self.untrimmed_dir_fmt.reverse_sequences.view(FastqGzFormat)
        self.assertTrue(reverse.path.exists())
        self.assertTrue(reverse.path.is_file())


if __name__ == '__main__':
    unittest.main()
