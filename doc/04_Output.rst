.. highlight:: bash
.. _output:

=======================================================================
Output section
=======================================================================

**GENESIS** yields trajectory data (coordinates and velocities)
in the *DCD* file format regardless of the force field or MD algorithm. 
**GENESIS** can also generate a restart file (*rstfile*) 
during or at the end of the simulation,
which can be used to restart and extend the simulation continuously.
Output frequency of each file (e.g., crdout_period and velout_period)
is specified in the **[DYNAMICS]** section in the case of the MD, 
REMD, and RPATH simulations, or **[MINIMIZE]** section in the case of
the energy minimization.

-----------------------------------------------------------------------

General output files
====================

**dcdfile**

  Filename for the coordinates trajectory data.
  Coordinates are written in the DCD format,
  which is commonly used in various MD software such as CHARMM and NAMD.
  The filename must be given in the case of ``crdout_period > 0``.
  However, if ``crdout_period = 0`` is specified in the control file, 
  no **dcdfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**dcdvelfile**

  Filename for the velocity trajectory data.
  Velocities are written in the DCD format.
  The filename must be given in the case of ``velout_period > 0``.
  However, if ``velout_period = 0`` is specified in the control file, 
  no **dcdvelfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**rstfile**

  Filename for the restart data.
  The rstfile contains coordinates, velocities, simulation box size, and so on. 
  This file can be used to extend the simulation continuously.
  In addition, it can be used to switch the simulation algorithms
  (e.g., from minimization to MD, from MD to REMD, from REMD to minimization, etc).
  The filename must be given in the case of ``rstout_period > 0``.
  However, if ``rstout_period = 0`` is specified in the control file,
  no **rstfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**pdbfile** (for ATDYN only)

  Filename for the restart PDB file. This file is updated every ``rstout_period`` steps.

**elefile** (SPDYN standard MD only)

  CSV filename for instantaneous electrostatic potential maps. Enable it with
  a positive ``eleout_period`` in **[DYNAMICS]**. The filename must not already
  exist; use a new filename for each restarted run. No file is opened when
  both options are omitted. Specifying only one option is an error.

  This output supports VVER and VRES with standard point-charge PME force
  fields, including NVT and NPT. Minimization, REMD, bias/RPATH, alchemical
  simulations, and SPICA/pSPICA electrostatics are rejected. It does not
  enable electrostatic output in ATDYN or CGDYN.

  Parallel-I/O restart input supports zero-field maps. A requested nonzero
  field with this output is rejected on that input path because the existing
  parallel-I/O setup does not apply external fields. Use ordinary input or
  restart files for simulations with an applied field.

  Each selected positive step contains every point of the actual PME grid
  selected by the FFT backend, with no duplicate periodic endpoints::

    step,time,x,y,z,charge_potential,external_potential,total_potential

  Time is in ps, coordinates in Angstrom, and potentials in volts. Step and
  time match the engine's run-local values; restarted runs reset this
  chronology. The ``initial_time`` setting is recorded separately.

  ``charge_potential`` is the Gaussian-smoothed reciprocal Ewald potential
  of all atomic charges, as in VMD PMEPot. It includes the instantaneous
  ionic and solvent response. It omits the direct-space term, intramolecular
  exclusions, and atom self-energy subtraction. The zero wavevector is
  omitted, corresponding to a uniform neutralizing background for a charged
  cell, and the spatial mean is zero. The Ewald factor controls smoothing;
  finite B-spline/grid errors should be checked by grid refinement.

  ``external_potential`` is ``-dot(effective_field,position-cell_center)``.
  The effective field includes ``efield_normal`` box scaling and the legacy
  field-force conversion factor 23.06 divided by the SI conversion from eV
  to kcal/mol. The latter is computed from the elementary charge, Avogadro
  constant, and thermochemical calorie. ``total_potential`` is the sum.
  The external reference at the cell center differs from the legacy
  external-energy routine's reference by a spatially constant offset.
  Potential differences and fields are unaffected by that gauge choice.

  Fixed metadata is written once in ``#key=value`` comments: units,
  potential/reference definitions, numerical precision, Ewald factor,
  dielectric constant, Coulomb prefactor, spline order, configured field,
  and conversion factors. Each ``#frame_step`` block records the current
  grid dimensions, grid origin/spacing, cell vectors, reference position,
  effective field, and net charge. Rows are ordered with x outermost and z
  innermost. Skip comment lines when reading the CSV.

  Arithmetic follows the build precision: mixed and single builds compute
  this map in single precision, while double builds use double precision.
  The printed digits preserve the stored values; they do not imply greater
  accuracy. Use double precision and check spatial/sampling convergence for
  demanding quantitative comparisons.

  Coordinates use the input-coordinate frame. For axis length L and native
  grid count N, the first point is ``cell_center-floor(N/2)*L/N``. Thus an
  odd grid starts half a spacing above the geometric lower cell boundary.
  DCD coordinates use the internally centered frame; add the recorded
  external reference position when comparing a DCD with the CSV.

  The calculation uses separate analysis arrays and does not change forces
  or advance the random-number generator. It gathers atoms and computes a
  serial FFT on the main rank, with memory and time costs proportional to
  the grid size and atom count. CSV files can be large. Average equilibrated
  frames and extract slices/profiles during analysis. Interpolate the
  periodic charge potential before adding the analytic external ramp;
  the total potential is not periodic across the applied voltage jump.

  Example::

    [OUTPUT]
    elefile = potential.csv

    [DYNAMICS]
    eleout_period = 1000


Output files in REMD and RPATH simulations
=====================================================================
When the user performs REMD or RPATH simulations,
the user must include '{}' in the output filename.
This {} is automatically replaced with the replica index.

**remfile** (only for REMD simulations)

  This file contains parameter index data from the REMD simulation, which is
  written for each replica every ``exchange_period`` steps.
  This is used as an input file for the *remd_convert* tool
  to sort the coordinates trajectory data by parameters. 
  The filename must contain '{}', which is automatically replaced with the replica index.
  Note that the information about the parameter index as well as replica index
  in the entire REMD simulation is written in the standard (single) output file
  (see online Tutorials).

**logfile** (only for REMD and RPATH simulations)

  This file contains the energy trajectory data from the REMD or RPATH simulations, which is
  written for each replica every ``exchange_period`` steps.
  This is used as an input file for the *remd_convert* tool
  to sort the coordinates trajectory data by parameters.
  The filename must contain '{}', which is automatically replaced with the replica index.

**rpathfile** (only for RPATH simulations)

  This file contains the trajectory of image coordinates in the string method,
  which are reference values used in the restraint functions.
  Columns correspond to the collective variables, and rows are time steps.
  This data is written with the same timing as the *dcdfile*.
  For details, see :ref:`RPATH`.


Output file in GaMD simulations
===================================

**gamdfile**

  This file provides GaMD parameters determined during the GaMD simulation.
  The filename must be given in the case of ``update_period > 0`` in **[GAMD]** section.
  The GaMD simulation updates its parameters every ``update_period`` steps, and then the updated parameters
  are output to **gamdfile**.
  This file includes the maximum, minimum, average, and deviation of the total potential or dihedral potential,
  which are calculated within the interval ``update_period``.


Output file in Vibrational analysis
===================================

**minfofile**

  This file provides the coordinates and normal mode vectors of the molecules
  specified for vibrational analysis. It is used in **SINDO** for visualizing
  the vibrational motion. It is also an input file to start anharmonic
  vibrational calculations. See :ref:`Vibration` for the vibrational 
  analysis.


Output file in FEP simulations
===================================

**fepfile**

  This file provides energy differences between adjacent windows in FEP simulations.
  The filename must be given in the case of ``fepout_period > 0`` in **[ALCHEMY]** section.
  However, if ``fepout_period = 0`` is specified in the control file, 
  no **fepfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.
  If FEP/:math:`\lambda`-REMD simulations are performed
  (i.e., both **[REMD]** and **[ALCHEMY]** sections are specified in the control file),
  the filename must contain '{}', which is automatically replaced with the replica index.


Output files in ABMD simulations
================================
When the user performs RE-ABMD, gREST-ABMD or BE-ABMD simulations,
the user must include '{}' in the output filename.
This {} is automatically replaced with the replica index.

**abmdfile**

  This file provides the time-dependent biassing potential along predefined collective variables (CVs).
  The filename must be given in the case of ``monitor_period > 0`` in **[ABMD]** section.
  ABMD and MW-ABMD provide a single file, whereas RE-ABMD, gREST-ABMD, and BE-ABMD provide multiple files for each replica.

**splinefile**

  This file provides the spline coefficients used to reconstruct the bias potential when restarting ABMD simulations.
  The filename must be given in the case of ``monitor_period_spline > 0`` in **[ABMD]** section.
  ABMD and MW-ABMD provide a single file, whereas RE-ABMD, gREST-ABMD, and BE-ABMD provide multiple files for each replica.

**reacfile**

  This file provides the computed values of the collective variables (CVs) during the simulation.
  The filename must be given in the case of ``monitor_period > 0`` in **[ReacCoord]** section.
  ABMD provides a single file, whereas MW-ABMD, RE-ABMD, gREST-ABMD, and BE-ABMD provide multiple files for each replica.


Examples
========

For normal MD simulations:
::
  
  [OUTPUT]
  dcdfile = run.dcd
  rstfile = run.rst

For REMD simulations:
::

  [OUTPUT]
  logfile = run_rep{}.log
  dcdfile = run_rep{}.dcd
  remfile = run_rep{}.rem
  rstfile = run_rep{}.rst
