# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def cutadapt_demux_single(use):
    sequence_url = ('https://qiime2-data.s3.us-west-2.amazonaws.com/'
                    'usage-examples/cutadapt/single_multiplexed.qza')
    metadata_url = ('https://qiime2-data.s3.us-west-2.amazonaws.com/'
                    'usage-examples/cutadapt/barcodes')

    seqs = use.init_artifact_from_url('seqs', sequence_url)
    md = use.init_metadata_from_url('md', metadata_url)
    barcodes = use.get_metadata_column('barcodes', 'BarcodeSequence', md)

    per_sample_sequences, untrimmed_sequences = use.action(

            use.UsageAction(plugin_id='cutadapt', action_id='demux_single'),
            use.UsageInputs(seqs=seqs, barcodes=barcodes),
            use.UsageOutputNames(per_sample_sequences='per_sample_sequences',
                                 untrimmed_sequences='untrimmed_sequences')
            )
    per_sample_sequences.assert_output_type('SampleData[SequencesWithQuality]')
    untrimmed_sequences.assert_output_type(
        'MultiplexedSingleEndBarcodeInSequence'
        )
