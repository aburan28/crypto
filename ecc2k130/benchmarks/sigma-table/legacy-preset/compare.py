"""Control-versus-candidate comparison of ECC_PACKED_SIGMA_TABLE on one RTX PRO 6000.

    /path/to/modal run compare.py

Builds the audited RTX PRO 6000 preset (control, byte-identical to production
apart from the unused table knob) and three candidates from this checkout:

    control    batch 32, 256 threads, two blocks/SM, networks (the preset)
    table      batch 32, 512 threads, one block/SM, shared-memory nibble tables
    table64    batch 64, 512 threads, one block/SM, tables
    control64  batch 64, 256 threads, two blocks/SM, networks

The table needs 68 KB of shared memory per block, so it runs one 512-thread
block per SM, which keeps the register budget (128) and resident warps (16)
of the preset. Each build also compiles the GPU arithmetic test, and the
image build disassembles every client so the walk kernel's instruction mix is
recorded. On the GPU: every arithmetic test must pass, the client integration
test must pass for the control and table clients, then complete walk
benchmarks alternate control and candidates after one excluded warm-up, each
at the preset grid of 192,512 workers (so a batch-64 run performs twice the
updates of a batch-32 run). Writes result.json.
"""
import json
import re
import statistics
import subprocess
from pathlib import Path

import modal

HERE = Path(__file__).parent
ROOT = HERE.parent.parent          # the ecc2k130 checkout
REMOTE = '/root/ecc2k130'
GENCODE = '-gencode arch=compute_120,code=sm_120'
COMMON = ('STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0 '
          'PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 '
          'PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 '
          'PACKED_DIRECT_REDUCE=1 PACKED_ADD_COMBINE=0')
COMMON_DEFS = ('-DECC_STREAM_KARAT=0 -DECC_SMEM_SPILL=0 '
               '-DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_CACHE_DENOM=1 -DECC_PACKED_BY_VALUE=1 '
               '-DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_POLY_CHAIN=1 -DECC_PACKED_UNROLL_INV=1 '
               '-DECC_PACKED_PAIR_PRODUCTS=1 -DECC_PACKED_POLY_STATE=1 -DECC_PACKED_DIRECT_REDUCE=1 '
               '-DECC_PACKED_ADD_COMBINE=0')
MODES = {
    'control':   dict(batch=32, threads=256, minBlocks=2, table=0),
    'table':     dict(batch=32, threads=512, minBlocks=1, table=1),
    'table64':   dict(batch=64, threads=512, minBlocks=1, table=1),
    'control64': dict(batch=64, threads=256, minBlocks=2, table=0),
}
WORKERS = 192512
STEPS, LAUNCHES = 1024, 32


def knobs(m):
    return (f'BATCH={m["batch"]} THREADS={m["threads"]} MINBLOCKS={m["minBlocks"]} '
            f'PACKED_SIGMA_TABLE={m["table"]} {COMMON}')


def defs(m):
    return (f'-DECC_BATCH={m["batch"]} -DECC_THREADS={m["threads"]} -DECC_MINBLOCKS={m["minBlocks"]} '
            f'-DECC_PACKED_SIGMA_TABLE={m["table"]} {COMMON_DEFS}')


buildCommands = ['nvcc --version']
for name, m in MODES.items():
    buildCommands += [
        f'cd {REMOTE} && (make -B gpu ARCH="{GENCODE}" {knobs(m)} > /root/build-{name}.log 2>&1 || (cat /root/build-{name}.log; exit 1))',
        f'cp {REMOTE}/ecc2k130 /root/client-{name}',
        f'cd {REMOTE} && nvcc -O3 -std=c++17 {GENCODE} {defs(m)} src/testpackedcuda.cu -o /root/test-{name}',
        f'cuobjdump -sass /root/client-{name} > /root/sass-{name}.txt',
    ]

app = modal.App('ecc2k130-sigma-table')
image = (modal.Image.from_registry('nvidia/cuda:13.0.0-devel-ubuntu24.04', add_python='3.12')
         .entrypoint([])
         .apt_install('build-essential')
         .add_local_dir(ROOT, remote_path=REMOTE, copy=True,
                        ignore=['ecc2k130-cpu', 'ecc2k130', 'build/*', '__pycache__', '*.pyc', 'aws/*',
                                'benchmarks/*', '.claude/*'])
         .run_commands(*buildCommands))


def sh(cmd, timeout):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    return dict(command=cmd, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)


def walkMix(path):
    """Opcode histogram of the walk kernel (helpers are embedded in its section)."""
    text = Path(path).read_text()
    sections = text.split('Function : ')
    section = next((s for s in sections[1:] if 'walk' in s.split('\n', 1)[0]), sections[-1] if len(sections) > 1 else text)
    ops = {}
    for line in section.splitlines():
        m = re.match(r'\s*/\*[0-9a-f]+\*/\s+(?:@!?U?P\w+\s+)?([A-Z][A-Z0-9_.]*)', line)
        if m:
            op = m.group(1)
            ops[op] = ops.get(op, 0) + 1
    fam = {}
    for op, n in ops.items():
        key = ('LOP3' if op.startswith('LOP3') else 'IADD3' if op.startswith('IADD3') else
               'IMAD.IADD' if op.startswith('IMAD.IADD') else 'IMAD.WIDE' if op.startswith('IMAD.WIDE') else
               'IMAD.SHL' if op.startswith('IMAD.SHL') else 'IMAD' if op.startswith('IMAD') else
               'SHF' if op.startswith('SHF') else 'IADD' if op.startswith('IADD') else
               'LDS' if op.startswith('LDS') else 'STS' if op.startswith('STS') else op.split('.')[0])
        fam[key] = fam.get(key, 0) + n
    return dict(instructions=sum(ops.values()), families=dict(sorted(fam.items(), key=lambda kv: -kv[1])))


def resources(name):
    log = Path(f'/root/build-{name}.log').read_text()
    walk = re.search(r"Compiling entry function '(\S*walk\S*)'.*?Used (\d+) registers.*?(?=Compiling entry|\Z)", log, re.S)
    spill = re.findall(r'(\d+) bytes spill stores, (\d+) bytes spill loads', log)
    stack = re.findall(r'(\d+) bytes stack frame', log)
    return dict(log=log[-4000:], walkRegisters=int(walk.group(2)) if walk else None,
                spillPairs=spill, stackFrames=stack)


@app.function(image=image, gpu='RTX-PRO-6000', timeout=5400)
def run():
    result = dict(kind='ECC_PACKED_SIGMA_TABLE control/candidate comparison; complete scalar walk updates',
                  gpu=sh('nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit --format=csv', 60),
                  nvcc=sh('nvcc --version', 60), modes=MODES, workers=WORKERS, steps=STEPS, launches=LAUNCHES,
                  mix={m: walkMix(f'/root/sass-{m}.txt') for m in MODES},
                  resources={m: resources(m) for m in MODES},
                  arithmetic={}, integration={}, walks=[])
    for m in MODES:
        t = sh(f'/root/test-{m}', 600)
        result['arithmetic'][m] = t
        print('arithmetic', m, t['returncode'], t['stdout'][-700:], t['stderr'][-600:], flush=True)
        if t['returncode'] != 0 or 'PASS' not in t['stdout']:
            result['valid'] = False
            return result
    for m in ('control', 'table'):
        t = sh(f'cd {REMOTE} && python3 codegen/testpackedclient.py /root/client-{m}', 1200)
        result['integration'][m] = t
        print('integration', m, t['returncode'], t['stdout'][-700:], t['stderr'][-900:], flush=True)
        if t['returncode'] != 0:
            result['valid'] = False
            return result
    order = [('control', True), ('control', False), ('table', False), ('control64', False), ('table64', False),
             ('control', False), ('table', False), ('control64', False), ('table64', False),
             ('control', False), ('table', False), ('control64', False), ('table64', False)]
    for m, warmup in order:
        cmd = f'/root/client-{m} --curve 131 --packed --bench --steps {STEPS} --launches {LAUNCHES} --verify 0 --threads {WORKERS}'
        p = sh(cmd, 1800)
        out = p['stdout']
        rate = re.findall(r'^\s*finished:\s+(\S+) M it/s,', out, re.M)
        iters = re.findall(r'M it/s\s+(\d+) iterations', out)
        identity = re.findall(r'packed sigma table: (\d)', out)
        resident = re.findall(r'(\d+) block\(s\) of (\d+) packed threads resident per SM', out)
        expected = WORKERS * MODES[m]['batch'] * STEPS * LAUNCHES
        row = dict(mode=m, warmup=warmup, returncode=p['returncode'],
                   rateMPerSecond=float(rate[0]) if len(rate) == 1 else None,
                   completedUpdates=int(iters[-1]) if iters else None, expectedUpdates=expected,
                   tableIdentity=identity[0] if identity else None, resident=resident[0] if resident else None,
                   stdout=out[-2500:], stderr=p['stderr'][-1000:])
        row['valid'] = (p['returncode'] == 0 and row['rateMPerSecond'] is not None and row['completedUpdates'] == expected
                        and row['tableIdentity'] == str(MODES[m]['table']) and 'MISMATCH' not in out and 'stopping:' not in out)
        result['walks'].append(row)
        print('walk', m, 'warmup' if warmup else 'timed', row['rateMPerSecond'], row['completedUpdates'], row['resident'], row['valid'], flush=True)
    result['valid'] = all(r['valid'] for r in result['walks'])
    return result


@app.local_entrypoint()
def main():
    r = run.remote()
    out = HERE / 'result.json'
    out.write_text(json.dumps(r, indent=1) + '\n')
    print('valid', r.get('valid'))
    for m in MODES:
        mix = r['mix'][m]
        print('mode', m, 'walk instructions', mix['instructions'], 'registers', r['resources'][m]['walkRegisters'],
              'spills', r['resources'][m]['spillPairs'][:2])
        print('   families', {k: v for k, v in list(mix['families'].items())[:10]})
    timed = [w for w in r.get('walks', []) if not w['warmup']]
    medians = {}
    for m in MODES:
        rates = [w['rateMPerSecond'] for w in timed if w['mode'] == m and w['valid']]
        medians[m] = statistics.median(rates) if rates else None
        print('walk mode', m, 'rates', rates, 'median', medians[m])
    c = [w['rateMPerSecond'] for w in timed if w['mode'] == 'control' and w['valid']]
    for m in ('table', 'table64', 'control64'):
        a = [w['rateMPerSecond'] for w in timed if w['mode'] == m and w['valid']]
        if c and a and len(c) == len(a):
            print(m, 'paired differences vs control in %:', [round(100 * (y - x) / x, 3) for x, y in zip(c, a)],
                  'median change: %.3f%%' % (100 * (statistics.median(a) - statistics.median(c)) / statistics.median(c)))
    print('wrote', out)
    if not r.get('valid'):
        raise SystemExit('comparison invalid')
