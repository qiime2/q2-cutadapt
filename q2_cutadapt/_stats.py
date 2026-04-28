# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from os import PathLike
from typing import Any, Mapping

import pandas as pd
import qiime2


JsonReportPath = str | PathLike[str]


def _adapter_sequence(adapter: dict[str, Any]) -> str:
    sequences: list[str] = []
    for end in ('five_prime_end', 'three_prime_end'):
        end_report = adapter.get(end)
        if end_report is not None:
            sequence = end_report.get('sequence')
            if sequence is not None and sequence not in sequences:
                sequences.append(sequence)

    if len(sequences) == 0:
        sequence = adapter.get('sequence', adapter.get('name'))
        if sequence is None:
            raise ValueError('Cutadapt adapter report has no sequence.')
        return str(sequence)
    elif len(sequences) == 1:
        return sequences[0]
    else:
        # Linked adapters contain a 5-prime and 3-prime sequence in one block.
        return '...'.join(sequences)


def _add_adapter_percentages(
    row: dict[str, int | float],
    report: dict[str, Any],
    read_key: str,
    read_label: str,
    adapter_sequences: list[str],
    read_denominator: int,
    total_denominator: int,
) -> None:
    adapter_reports = report.get(read_key)
    if adapter_reports is None:
        return

    for adapter in adapter_reports:
        sequence = _adapter_sequence(adapter)
        if sequence not in adapter_sequences:
            adapter_sequences.append(sequence)
        count = adapter.get('total_matches', 0)
        column = f'{read_label}-{sequence}'
        row[column] = row.get(column, 0) + _percent(
            count, read_denominator)
        total_column = f'total-{sequence}'
        row[total_column] = row.get(total_column, 0) + _percent(
            count, total_denominator)


def _percent(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0
    else:
        return numerator / denominator * 100


def summarize_cutadapt_json_reports(
    json_reports: Mapping[str, JsonReportPath],
) -> qiime2.Metadata:
    rows: dict[str, dict[str, int | float]] = {}
    base_columns = [
        'reads-before',
        'percent-reads-after',
        'bases-before',
        'percent-bases-after',
    ]
    adapter_sequences: list[str] = []

    for sample_id, json_report in json_reports.items():
        with open(json_report) as fh:
            report = json.load(fh)

        read_counts = report['read_counts']
        basepair_counts = report['basepair_counts']

        reads_before = read_counts['input']
        reads_after = read_counts['output']
        bases_before = basepair_counts['input']
        bases_after = basepair_counts['output']
        total_reads_before = reads_before
        if basepair_counts.get('input_read2') is not None:
            total_reads_before *= 2

        row = {
            'reads-before': reads_before,
            'percent-reads-after': _percent(reads_after, reads_before),
            'bases-before': bases_before,
            'percent-bases-after': _percent(bases_after, bases_before),
        }

        _add_adapter_percentages(row, report, 'adapters_read1', 'read1',
                                 adapter_sequences, reads_before,
                                 total_reads_before)
        _add_adapter_percentages(row, report, 'adapters_read2', 'read2',
                                 adapter_sequences, reads_before,
                                 total_reads_before)

        rows[sample_id] = row

    result = pd.DataFrame.from_dict(rows, orient='index').fillna(0)
    adapter_columns = [
        f'{read}-{sequence}'
        for sequence in adapter_sequences
        for read in ('read1', 'read2', 'total')
    ]
    result = result.reindex(columns=base_columns + adapter_columns,
                            fill_value=0)
    result.index.name = 'sample-id'
    return qiime2.Metadata(result)
