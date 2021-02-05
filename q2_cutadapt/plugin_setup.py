# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Plugin,
    Citations,
    MetadataColumn,
    Categorical,
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
    citations=Citations.load('citations.bib', package='q2_cutadapt')
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
        'minimum_length': Int % Range(1, None),
        'discard_untrimmed': Bool,
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
                   'is appended, the adapter is only found if it is at the '
                   'end of the read. If your sequence of interest is "framed" '
                   'by a 5\' and a 3\' adapter, use this parameter to define '
                   'a "linked" primer - see https://cutadapt.readthedocs.io '
                   'for complete details.',
        'front': 'Sequence of an adapter ligated to the 5\' end. The adapter '
                 'and any preceding bases are trimmed. Partial matches at the '
                 '5\' end are allowed. If a `^` character is prepended, the '
                 'adapter is only found if it is at the beginning of the '
                 'read.',
        'anywhere': 'Sequence of an adapter that may be ligated to the 5\' or '
                    '3\' end. Both types of matches as described under '
                    '`adapter` and `front` are allowed. If the first base of '
                    'the read is part of the match, the behavior is as with '
                    '`front`, otherwise as with `adapter`. This option is '
                    'mostly for rescuing failed library preparations - do '
                    'not use if you know which end your adapter was ligated '
                    'to.',
        'error_rate': 'Maximum allowed error rate.',
        'indels': 'Allow insertions or deletions of bases when matching '
                  'adapters.',
        'times': 'Remove multiple occurrences of an adapter if it is '
                 'repeated, up to `times` times.',
        'overlap': 'Require at least `overlap` bases of overlap between read '
                   'and adapter for an adapter to be found.',
        'match_read_wildcards': 'Interpret IUPAC wildcards (e.g., N) in '
                                'reads.',
        'match_adapter_wildcards': 'Interpret IUPAC wildcards (e.g., N) in '
                                   'adapters.',
        'minimum_length': 'Discard reads shorter than specified value. Note, '
                          'the cutadapt default of 0 has been overridden, '
                          'because that value produces empty sequence '
                          'records.',
        'discard_untrimmed': 'Discard reads in which no adapter was found.',
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
        'minimum_length': Int % Range(1, None),
        'discard_untrimmed': Bool,
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
                     'is appended, the adapter is only found if it is at the '
                     'end of the read. Search in forward read. If your '
                     'sequence of interest is "framed" by a 5\' and a 3\' '
                     'adapter, use this parameter to define a "linked" primer '
                     '- see https://cutadapt.readthedocs.io for complete '
                     'details.',
        'front_f': 'Sequence of an adapter ligated to the 5\' end. The '
                   'adapter and any preceding bases are trimmed. Partial '
                   'matches at the 5\' end are allowed. If a `^` character '
                   'is prepended, the adapter is only found if it is at the '
                   'beginning of the read. Search in forward read.',
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
                     'is appended, the adapter is only found if it is at the '
                     'end of the read. Search in reverse read. If your '
                     'sequence of interest is "framed" by a 5\' and a 3\' '
                     'adapter, use this parameter to define a "linked" primer '
                     '- see https://cutadapt.readthedocs.io for complete '
                     'details.',
        'front_r': 'Sequence of an adapter ligated to the 5\' end. The '
                   'adapter and any preceding bases are trimmed. Partial '
                   'matches at the 5\' end are allowed. If a `^` character '
                   'is prepended, the adapter is only found if it is at the '
                   'beginning of the read. Search in reverse read.',
        'anywhere_r': 'Sequence of an adapter that may be ligated to the 5\' '
                      'or 3\' end. Both types of matches as described under '
                      '`adapter` and `front` are allowed. If the first base '
                      'of the read is part of the match, the behavior is as '
                      'with `front`, otherwise as with `adapter`. This '
                      'option is mostly for rescuing failed library '
                      'preparations - do not use if you know which end your '
                      'adapter was ligated to. Search in reverse read.',
        'error_rate': 'Maximum allowed error rate.',
        'indels': 'Allow insertions or deletions of bases when matching '
                  'adapters.',
        'times': 'Remove multiple occurrences of an adapter if it is '
                 'repeated, up to `times` times.',
        'overlap': 'Require at least `overlap` bases of overlap between read '
                   'and adapter for an adapter to be found.',
        'match_read_wildcards': 'Interpret IUPAC wildcards (e.g., N) in '
                                'reads.',
        'match_adapter_wildcards': 'Interpret IUPAC wildcards (e.g., N) in '
                                   'adapters.',
        'minimum_length': 'Discard reads shorter than specified value. Note, '
                          'the cutadapt default of 0 has been overridden, '
                          'because that value produces empty sequence '
                          'records.',
        'discard_untrimmed': 'Discard reads in which no adapter was found.',
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
        'barcodes': MetadataColumn[Categorical],
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
        'batch_size': Int % Range(0, None),
        'minimum_length': Int % Range(1, None),
    },
    outputs=[
        ('per_sample_sequences', SampleData[SequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedSingleEndBarcodeInSequence),
    ],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata column listing the per-sample '
                    'barcodes.',
        'error_rate': 'The level of error tolerance, specified as the maximum '
                      'allowable error rate. The default value specified by '
                      'cutadapt is 0.1 (=10%), which is greater than '
                      '`demux emp-*`, which is 0.0 (=0%).',
        'batch_size': 'The number of samples cutadapt demultiplexes '
                      'concurrently. Demultiplexing in smaller batches will '
                      'yield the same result with marginal speed loss, and '
                      'may solve "too many files" errors related to sample '
                      'quantity. Set to "0" to process all samples at once.',
        'minimum_length': 'Discard reads shorter than specified value. Note, '
                          'the cutadapt default of 0 has been overridden, '
                          'because that value produces empty sequence '
                          'records.',
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
        'forward_barcodes': MetadataColumn[Categorical],
        'reverse_barcodes': MetadataColumn[Categorical],
        'error_rate': Float % Range(0, 1, inclusive_start=True,
                                    inclusive_end=True),
        'batch_size': Int % Range(0, None),
        'minimum_length': Int % Range(1, None),
        'mixed_orientation': Bool,
    },
    outputs=[
        ('per_sample_sequences', SampleData[PairedEndSequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedPairedEndBarcodeInSequence),
    ],
    input_descriptions={
        'seqs': 'The paired-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'forward_barcodes': 'The sample metadata column listing the '
                            'per-sample barcodes for the forward reads.',
        'reverse_barcodes': 'The sample metadata column listing the '
                            'per-sample barcodes for the reverse reads.',
        'error_rate': 'The level of error tolerance, specified as the maximum '
                      'allowable error rate.',
        'batch_size': 'The number of samples cutadapt demultiplexes '
                      'concurrently. Demultiplexing in smaller batches will '
                      'yield the same result with marginal speed loss, and '
                      'may solve "too many files" errors related to sample '
                      'quantity. Set to "0" to process all samples at once.',
        'minimum_length': 'Discard reads shorter than specified value. Note, '
                          'the cutadapt default of 0 has been overridden, '
                          'because that value produces empty sequence '
                          'records.',
        'mixed_orientation': 'Handle demultiplexing of mixed orientation '
                             'reads (i.e. when forward and reverse reads '
                             'coexist in the same file).'
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
