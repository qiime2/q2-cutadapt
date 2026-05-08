# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from os import PathLike
from typing import Any, Mapping, NamedTuple

import pandas as pd
import qiime2


JsonReportPath = str | PathLike[str]


class AdapterSpec(NamedTuple):
    read_key: str
    column: str


def _add_adapter_percentages(
    row: dict[str, int | float],
    report: dict[str, Any],
    read_key: str,
    adapter_specs: list[AdapterSpec],
    read_denominator: int,
) -> None:
    adapter_reports = report.get(read_key)
    if adapter_reports is None:
        return

    for adapter_spec, adapter in zip(adapter_specs, adapter_reports):
        count = adapter.get('total_matches', 0)
        row[adapter_spec.column] = _percent(count, read_denominator)


def _percent(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0
    else:
        return numerator / denominator * 100


def _extend_adapter_specs(
    adapter_specs: list[AdapterSpec],
    adapters: list[str] | None,
    end_label: str,
    read_label: str,
    read_key: str,
) -> None:
    for adapter in adapters or []:
        adapter_specs.append(
            AdapterSpec(read_key, f"{end_label} {read_label} {adapter}")
        )


def make_adapter_specs(
    adapter_f: list[str] | None = None,
    front_f: list[str] | None = None,
    anywhere_f: list[str] | None = None,
    adapter_r: list[str] | None = None,
    front_r: list[str] | None = None,
    anywhere_r: list[str] | None = None,
) -> list[AdapterSpec]:
    adapter_specs: list[AdapterSpec] = []
    _extend_adapter_specs(adapter_specs, adapter_f, "3'", 'R1',
                          'adapters_read1')
    _extend_adapter_specs(adapter_specs, front_f, "5'", 'R1',
                          'adapters_read1')
    _extend_adapter_specs(adapter_specs, anywhere_f, "5' or 3'", 'R1',
                          'adapters_read1')
    _extend_adapter_specs(adapter_specs, adapter_r, "3'", 'R2',
                          'adapters_read2')
    _extend_adapter_specs(adapter_specs, front_r, "5'", 'R2',
                          'adapters_read2')
    _extend_adapter_specs(adapter_specs, anywhere_r, "5' or 3'", 'R2',
                          'adapters_read2')
    return adapter_specs


def summarize_cutadapt_json_reports(
    json_reports: Mapping[str, JsonReportPath],
    adapter_specs: list[AdapterSpec],
) -> qiime2.Metadata:
    rows: dict[str, dict[str, int | float]] = {}
    base_columns = [
        'reads-before',
        'percent-reads-after',
        'bases-before',
        'percent-bases-after',
    ]
    read1_adapter_specs = [
        spec for spec in adapter_specs if spec.read_key == 'adapters_read1'
    ]
    read2_adapter_specs = [
        spec for spec in adapter_specs if spec.read_key == 'adapters_read2'
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

        row = {
            'reads-before': reads_before,
            'percent-reads-after': _percent(reads_after, reads_before),
            'bases-before': bases_before,
            'percent-bases-after': _percent(bases_after, bases_before),
        }

        _add_adapter_percentages(row, report, 'adapters_read1',
                                 read1_adapter_specs, reads_before)
        _add_adapter_percentages(row, report, 'adapters_read2',
                                 read2_adapter_specs, reads_before)

        rows[sample_id] = row

    result = pd.DataFrame.from_dict(rows, orient='index').fillna(0)
    adapter_columns = [spec.column for spec in adapter_specs]
    result = result.reindex(columns=base_columns + adapter_columns,
                            fill_value=0)
    result.index.name = 'sample-id'
    return qiime2.Metadata(result)
