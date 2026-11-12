#!/usr/bin/env python3
"""End-to-end SPDYN CSV checks using a deterministic neutral charge lattice.

Requires NumPy and a double-precision MPI build. Example:
  python3 tests/regression_test/test_ele_output.py --spdyn src/spdyn/spdyn
Add --baseline /path/to/unchanged/spdyn to compare control trajectories too.
All inputs, trajectories, logs and CSVs are created in a temporary directory.
"""

import argparse
import itertools
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import tempfile

import numpy as np

from test_ele_solver import exact_potential, KCAL_PER_MOL_PER_EV


COORDS = np.array(list(itertools.product((-12., 12.), repeat=3)))
CHARGES = np.array([(-1.) ** i for i in range(8)])
SCHEMES = ('opt_1dalltoall', 'noopt_1dalltoall', 'opt_2dalltoall',
           'noopt_2dalltoall', 'opt_slab', 'noopt_slab', 'autoselect')


def fixture(path):
    top = ('[ defaults ]\n1 2 yes 0.5 0.833333333333\n'
           '[ atomtypes ]\nION 18 39.948 0.0 A 0.3 0.1\n'
           '[ moleculetype ]\nIONS 0\n[ atoms ]\n')
    for i, q in enumerate(CHARGES, 1):
        top += f'{i} ION 1 ION A{i} {i} {q} 39.948\n'
    (path / 'ions.top').write_text(
        top + '[ system ]\nNeutral charge lattice\n[ molecules ]\nIONS 1\n')
    gro = 'Neutral charge lattice\n    8\n'
    for i, (x, y, z) in enumerate(COORDS / 10, 1):
        gro += (f'{1:5d}{"ION":<5}{("A"+str(i)):>5}{i:5d}'
                f'{x:8.3f}{y:8.3f}{z:8.3f}\n')
    (path / 'ions.gro').write_text(gro + '   4.80000   4.80000   4.80000\n')


def control(scheme='opt_1dalltoall', output=True, period=2, grid=32,
            field=(0., 0., 0.), integrator='VVER',
            nsteps=2, npt=False, normalized=False, timestep=1.e-9):
    shape = (grid, grid, grid) if isinstance(grid, int) else grid
    grid_lines = '' if grid is None else '\n'.join(
        f'pme_ngrid_{axis}={n}' for axis, n in zip('xyz', shape))
    return f'''[INPUT]
grotopfile=ions.top
grocrdfile=ions.gro
[OUTPUT]
dcdfile=coords.dcd
{'elefile=potential.csv' if output else ''}
[ENERGY]
forcefield=GROAMBER
electrostatic=PME
cutoffdist=8
switchdist=7
pairlistdist=10
pme_alpha=0.25
{grid_lines}
pme_scheme={scheme}
pme_nspline=4
efield_x={field[0]}
efield_y={field[1]}
efield_z={field[2]}
efield_normal={'YES' if normalized else 'NO'}
[DYNAMICS]
integrator={integrator}
nsteps={nsteps}
timestep={timestep}
iseed=777
eneout_period=1
crdout_period=2
nbupdate_period=2
elec_long_period=2
thermostat_period=2
barostat_period=2
eleout_period={period if output else 0}
stoptr_period=0
[ENSEMBLE]
ensemble={'NPT' if npt else 'NVE'}
tpcontrol={'BUSSI' if npt else 'NO'}
temperature={300 if npt else 0}
pressure=1.0
group_tp={'YES' if npt else 'NO'}
[CONSTRAINTS]
rigid_bond=NO
[BOUNDARY]
type=PBC
box_size_x=48
box_size_y=48
box_size_z=48
'''


def run(path, executable, command, content, ranks=1, threads=1,
        expected_error=None, existing_csv=False):
    path.mkdir()
    fixture(path)
    (path / 'inp').write_text(content)
    if existing_csv:
        (path / 'potential.csv').write_text('preserve existing scientific data\n')
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(threads)
    result = subprocess.run([*command, '-np', str(ranks), str(executable), 'inp'],
                            cwd=path, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=120)
    (path / 'run.log').write_text(result.stdout)
    if expected_error is not None:
        assert result.returncode != 0, result.stdout[-5000:]
        assert expected_error.lower() in result.stdout.lower(), result.stdout[-5000:]
        assert 'Output_Time>' not in result.stdout
        return result.stdout
    assert result.returncode == 0, result.stdout[-5000:]
    assert 'Output_Time>' in result.stdout, result.stdout[-5000:]
    return result.stdout


def read_csv(path):
    fixed, frames = {}, []
    header = None
    data = []
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                key, value = line[1:].strip().split('=', 1)
                if key == 'frame_step':
                    frames.append({key: value})
                elif frames:
                    frames[-1][key] = value
                else:
                    fixed[key] = value
            elif header is None:
                header = line.strip().split(',')
            else:
                data.append(line)
    assert header == ['step', 'time', 'x', 'y', 'z', 'charge_potential',
                      'external_potential', 'total_potential'], header
    return fixed, frames, np.loadtxt(data, delimiter=',', ndmin=2)


def dcd_numeric_records(path):
    """Compare binary records excluding the human-readable dated title."""
    data = path.read_bytes()
    endian = '<' if struct.unpack('<i', data[:4])[0] == 84 else '>'
    records = []
    offset = 0
    while offset < len(data):
        size = struct.unpack_from(endian + 'i', data, offset)[0]
        payload = data[offset+4:offset+4+size]
        assert struct.unpack_from(endian+'i', data, offset+4+size)[0] == size
        records.append(payload)
        offset += size + 8
    return records[:1] + records[2:]


def check_csv(path, field=(0., 0., 0.), expected_steps=(2,),
              origin=(0., 0., 0.), oracle=True, oracle_budget=3.e-10):
    fixed, frames, data = read_csv(path)
    assert np.isfinite(data).all()
    assert [int(f['frame_step']) for f in frames] == list(expected_steps)
    assert fixed['potential_unit'] == 'V'
    assert abs(float(fixed['coulomb_prefactor']) - 332.0636930) < 1.e-12
    previous = 0
    for frame, step in zip(frames, expected_steps):
        grid = tuple(int(frame['grid_dim_' + axis]) for axis in 'xyz')
        count = int(np.prod(grid))
        rows = data[previous:previous+count]
        previous += count
        assert len(rows) == count
        assert np.all(rows[:, 0] == step)
        assert np.max(np.abs(rows[:, 1] - step * float(fixed['timestep']))) < 1.e-16
        corner = np.fromstring(frame['grid_origin'], sep=',')
        spacing = np.fromstring(frame['grid_spacing'], sep=',')
        cell = np.array([np.fromstring(frame['cell_vector_' + a], sep=',')
                         for a in 'abc'])
        assert np.max(np.abs(np.diag(cell) - spacing * grid)) < 1.e-12
        positions = np.stack(np.meshgrid(*[corner[a] + spacing[a] *
            np.arange(grid[a]) for a in range(3)], indexing='ij'), axis=-1)
        assert np.max(np.abs(rows[:, 2:5] - positions.reshape(-1, 3))) < 1.e-12
        center = np.fromstring(frame['external_reference_position'], sep=',')
        assert np.max(np.abs(center - origin)) < 1.e-12
        assert np.max(np.abs(corner - (center - np.array(grid)//2*spacing))) < 1.e-12
        effective = np.fromstring(frame['effective_field'], sep=',')
        scale = np.fromstring(frame['field_box_scale'], sep=',')
        if int(fixed['efield_normal']):
            assert np.max(np.abs(scale - 48. / np.diag(cell))) < 1.e-12
        assert np.max(np.abs(effective - np.array(field) * scale *
                             23.06 / KCAL_PER_MOL_PER_EV)) < 1.e-15
        expected_ext = -(rows[:, 2:5] - center) @ effective
        assert np.max(np.abs(rows[:, 6] - expected_ext)) < 2.e-14
        assert np.max(np.abs(rows[:, 7] - rows[:, 5] - rows[:, 6])) < 2.e-14
        assert abs(rows[:, 5].mean()) < 1.e-13
        if oracle:
            reference = exact_potential(COORDS, CHARGES, grid, corner,
                                        np.diag(cell), .25, 1., 3)
            error = np.max(np.abs(rows[:, 5].reshape(grid) - reference))
            assert error < oracle_budget, error
    assert previous == len(data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spdyn', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--mpiexec', default='mpirun')
    args = parser.parse_args()
    executable = args.spdyn.resolve()
    command = shlex.split(args.mpiexec)
    with tempfile.TemporaryDirectory(prefix='genesis-ele-output-') as tmp:
        root = Path(tmp)
        first = None
        for scheme in SCHEMES:
            case = root / scheme
            log = run(case, executable, command, control(scheme=scheme), ranks=8)
            data = check_csv(case / 'potential.csv')
            if first is None:
                first = data
            else:
                assert np.max(np.abs(first - data)) < 3.e-12
            print('PASS native grid, Fourier oracle, MPI8:', scheme)
        for ranks, threads in ((1, 1), (2, 2), (8, 2)):
            case = root / f'ranks{ranks}-threads{threads}'
            run(case, executable, command, control(period=1), ranks, threads)
            data = check_csv(case / 'potential.csv', expected_steps=(1, 2))
            assert np.max(np.abs(first - data[data[:, 0] == 2])) < 3.e-12
            print('PASS MPI/thread ownership and cadence:', ranks, threads)
        field = (.001, -.002, -.00667)
        case = root / 'field'
        run(case, executable, command, control(field=field))
        check_csv(case / 'potential.csv', field=field)
        print('PASS physical coordinates and signed external field')
        case = root / 'odd-native-grid'
        run(case, executable, command, control(grid=(48, 48, 45)))
        check_csv(case / 'potential.csv', oracle_budget=1.e-3)
        print('PASS integer half-grid origin for odd native dimensions')
        case = root / 'inactive-configured-field'
        run(case, executable, command, control(field=(0., 0., 1.e-11)))
        fixed, frames, _ = read_csv(case / 'potential.csv')
        assert np.fromstring(fixed['configured_field'], sep=',')[2] == 1.e-11
        assert fixed['use_efield'] == '0'
        assert np.max(np.abs(np.fromstring(frames[0]['effective_field'], sep=','))) == 0
        print('PASS configured field retained below activation threshold')
        case = root / 'auto-grid'
        run(case, executable, command, control(grid=None))
        # Automatic grids need not align with this charge lattice; verify
        # geometry/field separately from finite-order mesh accuracy.
        check_csv(case / 'potential.csv', oracle_budget=1.e-3)
        print('PASS automatic native grid metadata')
        case = root / 'respa'
        run(case, executable, command, control(integrator='VRES'))
        check_csv(case / 'potential.csv')
        print('PASS RESPA output state and cadence')
        for normalized in (False, True):
            case = root / f'npt-normalized{normalized}'
            run(case, executable, command, control(npt=True, normalized=normalized,
                field=field, nsteps=4, timestep=.001))
            check_csv(case / 'potential.csv', field=field, expected_steps=(2,4),
                      oracle=False)
            print('PASS changing NPT geometry and field normalization:', normalized)
        enabled = root / 'enabled'
        disabled = root / 'disabled'
        a = run(enabled, executable, command, control(field=field))
        b = run(disabled, executable, command, control(output=False, field=field))
        assert [l for l in a.splitlines() if l.startswith('INFO:')] == [
            l for l in b.splitlines() if l.startswith('INFO:')]
        assert dcd_numeric_records(enabled/'coords.dcd') == dcd_numeric_records(disabled/'coords.dcd')
        if args.baseline:
            baseline = root / 'baseline'
            c = run(baseline, args.baseline.resolve(), command,
                    control(output=False, field=field).replace('eleout_period=0\n', ''))
            assert [l for l in b.splitlines() if l.startswith('INFO:')] == [
                l for l in c.splitlines() if l.startswith('INFO:')]
            assert dcd_numeric_records(baseline/'coords.dcd') == dcd_numeric_records(disabled/'coords.dcd')
        print('PASS enabled/disabled/control energies and byte-identical DCD numerical records')
        invalid = [
            ('missing-file', control().replace('elefile=potential.csv', ''), 'elefile must be specified'),
            ('negative-period', control(period=-1), 'eleout_period must be nonnegative'),
            ('zero-period', control(period=0), 'eleout_period must be positive'),
            ('cutoff', control().replace('electrostatic=PME', 'electrostatic=CUTOFF'), 'requires PME'),
            ('respa-cadence', control(integrator='VRES', period=1), 'eleout_period must be a multiple'),
            ('respa-truncation', control(integrator='VRES', nsteps=3).replace('crdout_period=2', 'crdout_period=0'), 'nsteps must be a multiple'),
            ('leap', control(integrator='LEAP'), 'LEAP integrator is not supported'),
            ('minimize', re.sub(r'\[DYNAMICS\].*?(?=\[ENSEMBLE\])',
                '[MINIMIZE]\nmethod=SD\nnsteps=2\neneout_period=1\nnbupdate_period=1\n',
                control(), flags=re.S), 'electrostatic output is supported only for standard MD'),
        ]
        for name, content, error in invalid:
            run(root/name, executable, command, content, expected_error=error)
            print('PASS rejected invalid input:', name)
        case = root / 'existing-file'
        run(case, executable, command, control(), ranks=8,
            expected_error='cannot create potential.csv', existing_csv=True)
        assert (case / 'potential.csv').read_text() == 'preserve existing scientific data\n'
        print('PASS MPI-safe rejection preserves an existing CSV')
        print('All electrostatic CSV integration checks passed.')


if __name__ == '__main__':
    main()
