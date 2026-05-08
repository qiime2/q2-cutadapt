# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from os import PathLike
from typing import Mapping

import pandas as pd
import qiime2


JsonReportPath = str | PathLike[str]

_READ_LABELS = {'adapters_read1': 'R1', 'adapters_read2': 'R2'}


def _percent(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0
    return numerator / denominator * 100


def _adapter_column(read_label: str, adapter: dict) -> str:
    five = adapter.get('five_prime_end')
    three = adapter.get('three_prime_end')
    if five and three:
        end = "5' or 3'"
        sequence = five['sequence']
    elif three:
        end = "3'"
        sequence = three['sequence']
    else:
        end = "5'"
        sequence = five['sequence']
    return f"{end} {read_label} {sequence}"


def _summarize_cutadapt_json_reports(
    json_reports: Mapping[str, JsonReportPath],
) -> qiime2.Metadata:
    rows: dict[str, dict[str, int | float]] = {}
    columns: list[str] = [
        'reads-before',
        'percent-reads-after',
        'bases-before',
        'percent-bases-after',
        'percent-bases-quality-trimmed',
    ]

    for sample_id, json_report_fp in json_reports.items():
        with open(json_report_fp) as fh:
            report = json.load(fh)

        read_counts = report['read_counts']
        basepair_counts = report['basepair_counts']

        reads_before = read_counts['input']
        reads_after = read_counts['output']
        bases_before = basepair_counts['input']
        bases_after = basepair_counts['output']
        bases_quality_trimmed = basepair_counts.get('quality_trimmed') or 0

        row: dict[str, int | float] = {
            'reads-before': reads_before,
            'percent-reads-after': _percent(reads_after, reads_before),
            'bases-before': bases_before,
            'percent-bases-after': _percent(bases_after, bases_before),
            'percent-bases-quality-trimmed': _percent(
                bases_quality_trimmed, bases_before),
        }

        for column, count in (
            ('percent-r1-with-adapter', read_counts.get('read1_with_adapter')),
            ('percent-r2-with-adapter', read_counts.get('read2_with_adapter')),
        ):
            if count is not None:
                row[column] = _percent(count, reads_before)
                if column not in columns:
                    columns.append(column)

        for read_key, read_label in _READ_LABELS.items():
            for adapter in report.get(read_key) or []:
                column = _adapter_column(read_label, adapter)
                row[column] = _percent(
                    adapter.get('total_matches', 0), reads_before)
                if column not in columns:
                    columns.append(column)

        rows[sample_id] = row

    result = pd.DataFrame.from_dict(rows, orient='index').fillna(0)
    result = result.reindex(columns=columns, fill_value=0)
    result.index.name = 'sample-id'
    return qiime2.Metadata(result)
