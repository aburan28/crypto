#!/usr/bin/env python3
"""Freeze the two root11 factor-base finalists and exact IC candidates."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
IDS = {
    2816: 'IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b',
    3072: 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a',
}


def manifest_path(k):
    folder = ('n83_intermediate_factor_base_holdout_20261017' if k == 2816
              else 'n83_smallk_challenge_20261017')
    return ROOT / 'research' / folder / f'candidate_{IDS[k]}.json'


def choose():
    for k, candidate_id in IDS.items():
        manifest = json.loads(manifest_path(k).read_text())
        if (manifest['candidate_id'] != candidate_id or
                manifest['factor_base']['actual_usable_points'] != 166 * k or
                manifest['factor_base']['effective_columns'] != k or
                manifest['implementation']['source_sha256']['ic'] != SOURCE_SHA256 or
                manifest['implementation']['algorithm_flags']['target_threads'] != 1):
            raise ValueError(f'K{k} exact candidate manifest differs')
    arms = tuple((f'k{k}root11', k, 'root11_dynamic', candidate_id)
                 for k, candidate_id in IDS.items())
    return arms, {
        'schema': 'n83.root11-finalists-choice/v1',
        'policy': 'fresh_pairwise_K2816_vs_K3072_after_intermediate_screen',
        'ordered_candidate_ids': [arm[3] for arm in arms],
    }
