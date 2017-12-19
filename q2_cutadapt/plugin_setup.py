# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Plugin,
    MetadataCategory,
    Float,
    Range,
    Int,
    List,
    Str,
    Bool,
)
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequence,
    MultiplexedPairedEndBarcodeInSequence,
)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
)

import q2_cutadapt
import q2_cutadapt._demux
import q2_cutadapt._trim


plugin = Plugin(
    name='cutadapt',
    version=q2_cutadapt.__version__,
    website='https://github.com/qiime2/q2-cutadapt',
    package='q2_cutadapt',
    description='This QIIME 2 plugin uses cutadapt to work with '
                'adapters (e.g. barcodes, primers) in sequence data.',
    short_description='Plugin for removing adapter sequences, primers, and '
                      'other unwanted sequence from sequence data.',
    citation_text='Martin, M. (2011). Cutadapt removes adapter sequences from '
                  'high-throughput sequencing reads. EMBnet.Journal, 17(1), '
                  'pp. 10-12.\ndoi:http://dx.doi.org/10.14806/ej.17.1.200',
)

plugin.methods.register_function(
    function=q2_cutadapt._trim.trim_single,
    inputs={
        'demultiplexed_sequences': SampleData[SequencesWithQuality],
    },
    parameters={
        'cores': Int % Range(1, None),
        'adapter': List[Str],
        'front': List[Str],
        'anywhere': List[Str],
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
        'indels': Bool,
        'times': Int % Range(1, None),
        'overlap': Int % Range(1, None),
        'match_read_wildcards': Bool,
        'match_adapter_wildcards': Bool,
    },
    outputs=[
        ('trimmed_sequences', SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        'demultiplexed_sequences': 'The single-end sequences to be trimmed.',
    },
    parameter_descriptions={
        'cores': 'Number of CPU cores to use.',
        'adapter': 'Sequence of an adapter ligated to the 3\' end. The '
                   'adapter and any subsequent bases are trimmed. If a `$` '
                   'is appended, the adapter is only found if it is a suffix '
                   'of the read.',
        'front': 'Sequence of an adapter ligated to the 5\' end. The adapter '
                 'and any preceding bases are trimmed. Partial matches at the '
                 '5\' end are allowed. If a `^` character is prepended, the '
                 'adapter is only found if it is a prefix of the read.',
        'anywhere': 'Sequence of an adapter that may be ligated to the 5\' or '
                    '3\' end. Both types of matches as described under '
                    '`adapter` and `front` are allowed. If the first base of '
                    'the read is part of the match, the behavior is as with '
                    '`front`, otherwise as with `adapter`. This option is '
                    'mostly for rescuing failed library preparations - do '
                    'not use if you know which end your adapter was ligated '
                    'to.',
        'error_rate': 'Maximum allowed error rate.',
        'indels': 'Allow indels in alignments.',
        'times': 'Remove up to `times` adapters from each read.',
        'overlap': 'Require `overlap` overlap between read and adapter for '
                   'an adapter to be found.',
        'match_read_wildcards': 'Interpret IUPAC wildcards in reads.',
        'match_adapter_wildcards': 'Interpret IUPAC wildcards in adapters.',
    },
    output_descriptions={
        'trimmed_sequences': 'The resulting trimmed sequences.',
    },
    name='Find and remove adapters in demultiplexed single-end sequences.',
    description='Search demultiplexed single-end sequences for adapters and '
                'remove them. The parameter descriptions in this method are '
                'adapted from the official cutadapt docs - please see those '
                'docs at https://cutadapt.readthedocs.io for complete '
                'details.',
)

plugin.methods.register_function(
    function=q2_cutadapt._trim.trim_paired,
    inputs={
        'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        'cores': Int % Range(1, None),
        'adapter_f': List[Str],
        'front_f': List[Str],
        'anywhere_f': List[Str],
        'adapter_r': List[Str],
        'front_r': List[Str],
        'anywhere_r': List[Str],
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
        'indels': Bool,
        'times': Int % Range(1, None),
        'overlap': Int % Range(1, None),
        'match_read_wildcards': Bool,
        'match_adapter_wildcards': Bool,
    },
    outputs=[
        ('trimmed_sequences', SampleData[PairedEndSequencesWithQuality]),
    ],
    input_descriptions={
        'demultiplexed_sequences': 'The paired-end sequences to be trimmed.',
    },
    parameter_descriptions={
        'cores': 'Number of CPU cores to use.',

        'adapter_f': 'Sequence of an adapter ligated to the 3\' end. The '
                     'adapter and any subsequent bases are trimmed. If a `$` '
                     'is appended, the adapter is only found if it is a '
                     'suffix of the read. Search in forward read.',
        'front_f': 'Sequence of an adapter ligated to the 5\' end. The '
                   'adapter and any preceding bases are trimmed. Partial '
                   'matches at the 5\' end are allowed. If a `^` character '
                   'is prepended, the adapter is only found if it is a '
                   'prefix of the read. Search in forward read.',
        'anywhere_f': 'Sequence of an adapter that may be ligated to the 5\' '
                      'or 3\' end. Both types of matches as described under '
                      '`adapter` and `front` are allowed. If the first base '
                      'of the read is part of the match, the behavior is as '
                      'with `front`, otherwise as with `adapter`. This option '
                      'is mostly for rescuing failed library preparations - '
                      'do not use if you know which end your adapter was '
                      'ligated to. Search in forward read.',
        'adapter_r': 'Sequence of an adapter ligated to the 3\' end. The '
                     'adapter and any subsequent bases are trimmed. If a `$` '
                     'is appended, the adapter is only found if it is a '
                     'suffix of the read. Search in reverse read.',
        'front_r': 'Sequence of an adapter ligated to the 5\' end. The '
                   'adapter and any preceding bases are trimmed. Partial '
                   'matches at the 5\' end are allowed. If a `^` character '
                   'is prepended, the adapter is only found if it is a '
                   'prefix of the read. Search in reverse read.',
        'anywhere_r': 'Sequence of an adapter that may be ligated to the 5\' '
                      'or 3\' end. Both types of matches as described under '
                      '`adapter` and `front` are allowed. If the first base '
                      'of the read is part of the match, the behavior is as '
                      'with `front`, otherwise as with `adapter`. This '
                      'option is mostly for rescuing failed library '
                      'preparations - do not use if you know which end your '
                      'adapter was ligated to. Search in reverse read.',
        'error_rate': 'Maximum allowed error rate.',
        'indels': 'Allow indels in alignments.',
        'times': 'Remove up to `times` adapters from each read.',
        'overlap': 'Require `overlap` overlap between read and adapter for '
                   'an adapter to be found.',
        'match_read_wildcards': 'Interpret IUPAC wildcards in reads.',
        'match_adapter_wildcards': 'Interpret IUPAC wildcards in adapters.',
    },
    output_descriptions={
        'trimmed_sequences': 'The resulting trimmed sequences.',
    },
    name='Find and remove adapters in demultiplexed paired-end sequences.',
    description='Search demultiplexed paired-end sequences for adapters and '
                'remove them. The parameter descriptions in this method are '
                'adapted from the official cutadapt docs - please see those '
                'docs at https://cutadapt.readthedocs.io for complete '
                'details.',
)

plugin.methods.register_function(
    function=q2_cutadapt._demux.demux_single,
    inputs={
        'seqs': MultiplexedSingleEndBarcodeInSequence,
    },
    parameters={
        'barcodes': MetadataCategory,
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
    },
    outputs=[
        ('per_sample_sequences', SampleData[SequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedSingleEndBarcodeInSequence),
    ],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata category listing the per-sample '
                    'barcodes.',
        'error_rate': 'The level of error tolerance, specified as the maximum '
                      'allowable error rate. The default value specified by '
                      'cutadapt is 0.1 (=10%), which is greater than '
                      '`demux emp-*`, which is 0.0 (=0%).',
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'untrimmed_sequences': 'The sequences that were unmatched to '
                               'barcodes.',
    },
    name='Demultiplex single-end sequence data with barcodes in-sequence.',
    description='Demultiplex sequence data (i.e., map barcode reads to '
                'sample ids). Barcodes are expected to be located within the '
                'sequence data (versus the header, or a separate barcode '
                'file).',
)

plugin.methods.register_function(
    function=q2_cutadapt._demux.demux_paired,
    inputs={
        'seqs': MultiplexedPairedEndBarcodeInSequence,
    },
    parameters={
        'forward_barcodes': MetadataCategory,
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
    },
    outputs=[
        ('per_sample_sequences', SampleData[PairedEndSequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedPairedEndBarcodeInSequence),
    ],
    input_descriptions={
        'seqs': 'The paired-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'forward_barcodes': 'The sample metadata category listing the '
                            'per-sample barcodes for the forward reads.',
        'error_rate': 'The level of error tolerance, specified as the maximum '
                      'allowable error rate.',
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'untrimmed_sequences': 'The sequences that were unmatched to '
                               'barcodes.',
    },
    name='Demultiplex paired-end sequence data with barcodes in-sequence.',
    description='Demultiplex sequence data (i.e., map barcode reads to '
                'sample ids). Barcodes are expected to be located within the '
                'sequence data (versus the header, or a separate barcode '
                'file).',
)
