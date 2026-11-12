#!/usr/bin/env python3
"""Check electrostatic output against existing membrane/solid regression systems.

Double-precision MPI SPDYN and NumPy are required. Example from GENESIS root:
  python3 tests/regression_test/test_ele_existing.py --spdyn src/spdyn/spdyn
Optional --baseline /path/to/unchanged/spdyn adds an independent control run.
"""

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np

from test_ele_output import dcd_numeric_records, read_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spdyn', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--mpiexec', default='mpirun')
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    cases = [('dppc/VVER_NVT_BUSSI1', False),
             ('dppc/VRES_NPT_BUSSI_group1', True),
             ('nano_material/VVER_NPT_BUSSI1', False)]
    with tempfile.TemporaryDirectory(prefix='genesis-ele-existing-') as tmp:
        root = Path(tmp)
        for name, normalized in cases:
            source = repository / 'tests/regression_test/test_spdyn' / name
            lines, section = [], ''
            for line in (source / 'inp').read_text().splitlines():
                if line.lstrip().startswith('['):
                    section = line.strip().upper()
                if section == '[INPUT]' and '=' in line and not line.lstrip().startswith('#'):
                    key, value = line.split('=', 1)
                    value = value.split('#')[0].strip()
                    if value:
                        line = key + '= ' + str((source/value).resolve())
                lines.append(line)
            text = '\n'.join(lines) + '\n'
            if '[OUTPUT]' not in text:
                text = '[OUTPUT]\n' + text
            text = text.replace('[OUTPUT]', '[OUTPUT]\ndcdfile=coords.dcd')
            text = text.replace('[DYNAMICS]', '[DYNAMICS]\ncrdout_period=20')
            text = text.replace('[ENERGY]', '[ENERGY]\nefield_z=-0.00667\n'
                                'efield_normal=' + ('YES' if normalized else 'NO'))
            variants = [('enabled', args.spdyn.resolve()),
                        ('disabled', args.spdyn.resolve())]
            if args.baseline:
                variants.append(('baseline', args.baseline.resolve()))
            energies, trajectories, enabled = [], [], None
            for variant, executable in variants:
                work = root / (name.replace('/', '-') + '-' + variant)
                work.mkdir()
                content = text
                if variant == 'enabled':
                    enabled = work
                    content = content.replace('[OUTPUT]', '[OUTPUT]\nelefile=potential.csv')
                    content = content.replace('[DYNAMICS]', '[DYNAMICS]\neleout_period=20')
                (work/'inp').write_text(content)
                env = os.environ.copy()
                env['OMP_NUM_THREADS'] = '1'
                run = subprocess.run([*shlex.split(args.mpiexec), '-np', '1',
                                      str(executable), 'inp'], cwd=work, env=env,
                                     text=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, timeout=180)
                (work/'run.log').write_text(run.stdout)
                assert run.returncode == 0 and 'Output_Time>' in run.stdout, run.stdout[-6000:]
                energies.append([line for line in run.stdout.splitlines()
                                 if line.startswith('INFO:')])
                trajectories.append(dcd_numeric_records(work/'coords.dcd'))
            assert all(lines == energies[0] for lines in energies), name
            assert all(records == trajectories[0] for records in trajectories), name
            fixed, frames, data = read_csv(enabled/'potential.csv')
            assert len(frames) == 1 and int(frames[0]['frame_step']) == 20
            frame = frames[0]
            assert np.isfinite(data).all()
            assert np.max(np.abs(data[:, 7]-data[:, 5]-data[:, 6])) < 2.e-14
            assert abs(data[:, 5].mean()) < 1.e-12
            cell = np.array([np.fromstring(frame['cell_vector_'+a], sep=',') for a in 'abc'])
            volume = np.prod(cell.diagonal())
            header = energies[0][0].split()[1:]
            logged_volume = float(energies[0][-1].split()[1+header.index('VOLUME')])
            assert abs(volume-logged_volume) < .001
            field = np.fromstring(frame['effective_field'], sep=',')
            center = np.fromstring(frame['external_reference_position'], sep=',')
            assert np.max(np.abs(data[:, 6]+(data[:, 2:5]-center)@field)) < 2.e-14
            print('PASS', name, 'normalized field:', normalized, 'rows:', len(data),
                  'identical energies and DCD numerical records', flush=True)


if __name__ == '__main__':
    main()
