#!/usr/bin/env python3
"""Resolve the sealed twelve-point holdout's candidate list from two panels."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
K992_ID = 'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094'
K1107_ID = 'IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa'
K1024_ID = 'IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289'
K1280_ID = 'IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


five_run = module('n83_final_arms_five_run',
                  ROOT / 'research/n83_smallk_challenge_20261017/run.py')
five_select = module('n83_final_arms_five_select',
                     ROOT / 'research/n83_smallk_challenge_20261017/select.py')
k992_run = module('n83_final_arms_k992_run',
                  ROOT / 'research/n83_k992_cold_followup_20261017/run.py')
k992_select = module('n83_final_arms_k992_select',
                     ROOT / 'research/n83_k992_cold_followup_20261017/select.py')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel(validated_path, selection_path, selector, schema, targets_sha):
    validated = json.loads(validated_path.read_text())
    selection = json.loads(selection_path.read_text())
    if (validated.get('schema') != schema or
            validated.get('source_sha256') != SOURCE_SHA256 or
            validated.get('targets_sha256') != targets_sha or
            selection != selector.select(validated)):
        raise ValueError('upstream panel evidence or frozen selection differs')
    selected = (selection['online']['selected_candidate_id'] if
                selection.get('admission_status') == 'PASS' and
                selection['online']['status'] == 'SELECTED' else None)
    if selected is not None and selected not in validated['candidate_ids']:
        raise ValueError('upstream selected candidate was absent')
    return validated, selection, selected, {
        'run_id': validated['run_id'],
        'head_sha': validated['head_sha'],
        'validated_sha256': digest(validated_path),
        'selection_sha256': digest(selection_path),
        'admission_status': selection['admission_status'],
        'online_status': selection['online']['status'],
        'selected_candidate_id': selected,
    }


def unique_names(chosen):
    """Keep candidate IDs intact when two exact configurations share a K label."""
    used = set()
    result = []
    for name, k, layout, candidate_id in chosen:
        if name in used:
            name = f'{name}_h{candidate_id[-12:]}'
        if name in used:
            raise ValueError('candidate run labels still collide')
        used.add(name)
        result.append((name, k, layout, candidate_id))
    return tuple(result)


def choose(five_validated, five_selection, k992_validated, k992_selection):
    f, _, f_selected, f_record = panel(
        five_validated, five_selection, five_select,
        'n83.smallk-challenge-validated/v1', five_run.TARGETS_SHA256)
    k, _, k_selected, k_record = panel(
        k992_validated, k992_selection, k992_select,
        'n83.k992-cold-followup-validated/v1', k992_run.TARGETS_SHA256)
    if (k['candidate_ids'] != [arm[3] for arm in k992_run.ARMS] or
            f['candidate_ids'][:4] != [arm[3] for arm in five_run.FIXED_ARMS] or
            f['candidate_ids'][4] != f['choice']['comparison_candidate_id']):
        raise ValueError('upstream candidate manifests or fifth arm differ')
    f_arm_id = f_selected or K1280_ID
    k_arm_id = k_selected or K1024_ID
    all_arms = (*k992_run.ARMS, *five_run.FIXED_ARMS,
                *five_run.finalist.ARMS)
    by_id = {arm[3]: arm for arm in all_arms}
    if len(by_id) != len(all_arms):
        raise ValueError('candidate registry has a duplicate ID')
    requested = (K992_ID, K1107_ID, f_arm_id, k_arm_id)
    if any(item not in by_id for item in requested):
        raise ValueError('chosen candidate ID has no exact arm configuration')
    arms = unique_names(tuple(by_id[item] for item in dict.fromkeys(requested)))
    if len(arms) < 3 or len(arms) > 4:
        raise ValueError('deduplicated holdout arm count differs')
    return arms, {
        'schema': 'n83.final-online-holdout-choice/v1',
        'source_sha256': SOURCE_SHA256,
        'five_arm_panel': f_record,
        'k992_panel': k_record,
        'five_arm_fallback_used': f_selected is None,
        'k992_fallback_used': k_selected is None,
        'ordered_candidate_ids': [arm[3] for arm in arms],
    }
