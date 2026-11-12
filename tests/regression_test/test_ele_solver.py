#!/usr/bin/env python3
"""Independent reciprocal-space oracle for the electrostatic output solver.

Requires NumPy and a configured/built double-precision MPI SPDYN tree. Example:
  python3 tests/regression_test/test_ele_solver.py --build .

The reference uses exact particle structure factors, not mesh assignment or
the production solver's B-spline deconvolution. Temporary build artifacts are
removed automatically. Fractional particle positions have finite PME error;
grid-aligned charges provide an essentially exact normalization/phase check.
"""

import argparse
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np


KCAL_PER_MOL_PER_EV = 1.602176634e-19 * 6.02214076e23 / 4184.
ELECOEF = {1: 332.0716, 2: 332.05221729, 3: 332.0636930}


def exact_potential(coords, charges, grid, origin, box, alpha, dielectric,
                    forcefield=1, spline=4):
    """Evaluate the finite reciprocal sum by exact structure factors and FFT."""
    frequencies = [np.fft.fftfreq(n) * n for n in grid]
    wave2 = sum((2 * np.pi * frequencies[a].reshape(
        tuple(grid[a] if d == a else 1 for d in range(3))) / box[a]) ** 2
                for a in range(3))
    rho = np.zeros(grid, dtype=complex)
    for position, charge in zip(coords, charges):
        phase = []
        for axis, n in enumerate(grid):
            angle = 2 * np.pi * frequencies[axis] * (
                position[axis] - origin[axis]) / box[axis]
            values = np.exp(-1j * angle)
            # The real grid represents the even +/- Nyquist average. Odd
            # assignment orders omit Nyquist, as the SPDYN kernels do.
            if n % 2 == 0:
                values[n // 2] = (0 if spline % 2 else
                                   np.cos(angle[n // 2]))
            phase.append(values)
        rho += charge * (phase[0][:, None, None] * phase[1][None, :, None]
                         * phase[2][None, None, :])
    green = np.zeros(grid)
    nz = wave2 > 0
    green[nz] = np.exp(-wave2[nz] / (4 * alpha ** 2)) / wave2[nz]
    green *= 4 * np.pi * ELECOEF[forcefield] / (
        np.prod(box) * dielectric * KCAL_PER_MOL_PER_EV)
    return np.fft.ifftn(rho * green).real * np.prod(grid)


def solver(binary, coords, charges, grid, origin, box, alpha=.25,
           dielectric=1., forcefield=1, spline=4):
    header = [*grid, spline, alpha, dielectric, len(charges), forcefield]
    lines = [' '.join(map(str, header)),
             ' '.join(map(str, [*origin, *box]))]
    lines += [' '.join(map(str, [*r, q])) for r, q in zip(coords, charges)]
    result = subprocess.run([str(binary)], input='\n'.join(lines) + '\n',
                            text=True, capture_output=True, check=True)
    values = np.fromstring(result.stdout, sep=' ')
    assert values.size == np.prod(grid), result.stdout[-1000:]
    values = values.reshape(tuple(reversed(grid))).transpose(2, 1, 0)
    assert np.isfinite(values).all()
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--fc', default='mpifort')
    args = parser.parse_args()
    build = args.build.resolve()
    source = Path(__file__).with_name('ele_solver_driver.F90')
    with tempfile.TemporaryDirectory(prefix='genesis-ele-oracle-') as tmp:
        binary = Path(tmp) / 'oracle_driver'
        subprocess.run([*shlex.split(args.fc), '-fopenmp',
                        '-I' + str(build / 'src/lib'),
                        '-I' + str(build / 'src/spdyn'), str(source),
                        str(build / 'src/spdyn/sp_electrostatic.o'),
                        str(build / 'src/lib/lib.a'), '-llapack', '-lblas',
                        '-o', str(binary)], check=True)
        subprocess.run([str(binary), 'finite-check'], check=True)
        origin = np.array([-16., -12., -20.])
        box = np.array([32., 24., 40.])
        coords = np.array([[-8., -6., -10.], [6., 4.5, 7.5]])
        charges = np.array([1., -.7])  # also exercises a non-neutral cell
        grid = (32, 32, 32)
        for spline in (3, 4, 5, 6, 8):
            for forcefield in (1, 2, 3):
                got = solver(binary, coords, charges, grid, origin, box,
                             forcefield=forcefield, spline=spline)
                expected = exact_potential(coords, charges, grid, origin,
                                           box, .25, 1., forcefield, spline)
                error = np.max(np.abs(got - expected))
                print(f'aligned spline={spline} forcefield={forcefield} '
                      f'max_abs_error={error:.3e} V')
                assert error < 2.e-11, error
                assert abs(got.mean()) < 2.e-14
        coords += np.array([.237, -.193, .311])
        errors = []
        # Budgets are 0.3 mV on the coarse mesh and 0.03 mV on the fine
        # mesh. Unlike the aligned cases, fractional sources have genuine
        # finite-order assignment error; refinement must reduce it strongly.
        for n, budget in ((32, 3.e-4), (64, 3.e-5)):
            shape = (n, n, n)
            got = solver(binary, coords, charges, shape, origin, box)
            expected = exact_potential(coords, charges, shape, origin,
                                       box, .25, 1.)
            error = np.max(np.abs(got - expected))
            errors.append(error)
            print(f'fractional grid={shape} max_abs_error={error:.3e} V')
            assert error < budget
        assert errors[1] < errors[0] / 5, errors
        for shape in ((30, 27, 45), (48, 48, 45)):
            lengths = np.array(shape, dtype=float)
            corner = -np.floor_divide(shape, 2).astype(float)
            particles = corner + np.array([[3., 4., 5.], [17., 19., 23.]])
            got = solver(binary, particles, charges, shape, corner, lengths)
            expected = exact_potential(particles, charges, shape, corner,
                                       lengths, .25, 1.)
            error = np.max(np.abs(got - expected))
            assert error < 2.e-11, error
            wrong = exact_potential(particles, charges, shape,
                                    corner + np.array([0., 0., .5]),
                                    lengths, .25, 1.)
            assert np.max(np.abs(got - wrong)) > 1.e-3
            print(f'odd/anisotropic grid={shape} max_abs_error={error:.3e} V; '
                  'half-spacing negative control detected')
        a = solver(binary, coords, charges, grid, origin, box)
        b = solver(binary, coords + box, charges, grid, origin, box)
        assert np.max(np.abs(a-b)) < 2.e-13
        b = solver(binary, coords, charges, grid, origin, box, dielectric=2.)
        assert np.max(np.abs(a/2-b)) < 2.e-13
        zero = solver(binary, coords, charges * 0, grid, origin, box)
        assert np.max(np.abs(zero)) == 0
        for invalid in (np.nan, np.inf, -np.inf):
            for field in ('coords', 'charges', 'origin', 'box', 'alpha', 'dielectric'):
                arguments = dict(coords=coords.copy(), charges=charges.copy(),
                                 grid=grid, origin=origin.copy(), box=box.copy(),
                                 alpha=.25, dielectric=1.)
                if field in ('alpha', 'dielectric'):
                    arguments[field] = invalid
                else:
                    arguments[field].flat[0] = invalid
                try:
                    solver(binary, **arguments)
                except subprocess.CalledProcessError as error:
                    assert 'Nonfinite' in error.stdout + error.stderr
                else:
                    raise AssertionError(f'Accepted nonfinite {field}: {invalid}')
        print('PASS nonfinite inputs rejected under production optimization flags.')
        print('PASS: exact Fourier normalization, phase, spline orders, '
              'forcefields, zero mode, refinement, periodicity and dielectric.')


if __name__ == '__main__':
    main()
