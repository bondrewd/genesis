!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_electrostatic_mod
!> @brief   Smoothed reciprocal electrostatic potential on the PME grid.
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_electrostatic_mod

  use constants_mod
  use math_libs_mod, only: b_spline_coef
  use messages_mod,  only: ErrOut
  use, intrinsic :: iso_fortran_env, only: int32, int64
#ifdef HAVE_MPI_GENESIS
  use mpi_parallel_mod, only: mpi_comm_country
  use mpi
#endif

  implicit none
  private

  ! One electron volt in kcal/mol, using the SI-defined elementary charge,
  ! Avogadro constant, and thermochemical calorie.  In particular this is not
  ! the rounded 23.06 used by the legacy external-field force calculation.
  real(wp), public, parameter :: kcal_per_mol_per_ev = &
       1.602176634e-19_wp*6.02214076e23_wp/4184.0_wp

  public :: compute_electrostatic_potential
  public :: electrostatic_is_finite

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_electrostatic_potential
  !> @brief        Calculate a Gaussian-smoothed charge potential in volts.
  !! @param[in]    coord      : physical atomic coordinates, Angstrom
  !! @param[in]    charge     : atomic charges, elementary charge units
  !! @param[in]    ngrid      : actual simulation PME grid dimensions
  !! @param[in]    origin     : physical lower corner of the grid, Angstrom
  !! @param[in]    box        : orthorhombic cell lengths, Angstrom
  !! @param[in]    alpha      : actual Ewald factor, inverse Angstrom
  !! @param[in]    dielectric : force-field dielectric constant
  !! @param[in]    nspline    : actual simulation PME B-spline order
  !! @param[out]   potential  : charge potential with zero cell mean, volts
  !
  ! This analysis calculation never touches a force-kernel workspace.  It
  ! approximates the point-charge structure factor S(k)=sum q exp(-i k.r)
  ! using one charge-assignment deconvolution and computes
  !
  !   phi(r) = 4*pi*ELECOEF/(V*dielectric*kcal_per_mol_per_ev)
  !            * sum[k/=0] exp(-k^2/(4*alpha^2))/k^2 S(k) exp(i k.r).
  !
  ! Coordinates in this formula are relative to origin.  Let w_l(0) be the
  ! normalized B-spline weights for an atom exactly on grid index u=0.
  ! The assignment used here puts w_l(0) at index 1-l (zero-based), so its
  ! Fourier transfer function is H(m)=sum_l w_l(0) exp(+2*pi*i*m*(l-1)/N).
  ! Therefore S(m) is approximated by FFT(Q)/H(m), INCLUDING H's phase.
  ! The force kernel's |H|^-2 factor instead compensates both assignment
  ! and interpolation back to atoms and must not be used for this map.
  !
  ! zfft1d uses unnormalized transforms in both directions in GENESIS;
  ! consequently no division by the number of grid points is applied.
  ! The zero Fourier mode is omitted.  For a charged cell this is the
  ! conventional uniform neutralizing background.  No direct-space term,
  ! exclusions, or atom self-energy subtraction belongs to this map.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_electrostatic_potential(coord, charge, ngrid, origin, &
                                             box, alpha, dielectric, &
                                             nspline, potential)

    real(wp),                intent(in)    :: coord(:,:), charge(:)
    integer,                 intent(in)    :: ngrid(3), nspline
    real(wp),                intent(in)    :: origin(3), box(3)
    real(wp),                intent(in)    :: alpha, dielectric
    real(wp), allocatable,   intent(out)   :: potential(:,:,:)

    complex(wp), allocatable :: grid(:,:,:), correction(:,:)
    real(wp),    allocatable :: wave2(:,:)
    real(wp)                 :: weights(nspline,3), integer_weights(nspline)
    real(wp)                 :: u, fraction, angle, bnormal, k2, factor
    real(wp)                 :: prefactor, grid_mean, r_scale(3)
    complex(wp)              :: transfer
    integer                  :: i, j, k, l, axis, atom, m, nmax
    integer                  :: base(3), ix, iy, iz, gx, gy, gz, nfactor, alloc_stat


    if (size(coord,1) /= 3 .or. size(coord,2) /= size(charge)) &
      call solver_error('Compute_Electrostatic_Potential> Invalid atom arrays.')
    if (any(ngrid < 1) .or. nspline < 3) &
      call solver_error('Compute_Electrostatic_Potential> Invalid PME grid or order.')
    if (.not. all(electrostatic_is_finite(box)) .or. &
        .not. electrostatic_is_finite(alpha) .or. &
        .not. electrostatic_is_finite(dielectric)) &
      call solver_error('Compute_Electrostatic_Potential> Nonfinite cell or Ewald parameters.')
    if (any(box <= 0.0_wp) .or. alpha <= 0.0_wp .or. dielectric <= 0.0_wp) &
      call solver_error('Compute_Electrostatic_Potential> Invalid cell or Ewald parameters.')
    if (.not. all(electrostatic_is_finite(origin)) .or. &
        .not. all(electrostatic_is_finite(coord)) .or. &
        .not. all(electrostatic_is_finite(charge))) &
      call solver_error('Compute_Electrostatic_Potential> Nonfinite coordinates or charges.')

    ! FFTE accepts dimensions factorizable into 2, 3, and 5 only.
    do axis = 1, 3
      nfactor = ngrid(axis)
      do i = 2, 5
        if (i == 4) cycle
        do while (mod(nfactor,i) == 0)
          nfactor = nfactor/i
        end do
      end do
      if (nfactor /= 1) &
        call solver_error('Compute_Electrostatic_Potential> Grid is not a 2/3/5 FFT size.')
    end do

    nmax = maxval(ngrid)
    allocate(grid(ngrid(1),ngrid(2),ngrid(3)), &
             correction(nmax,3), wave2(nmax,3), stat=alloc_stat)
    if (alloc_stat /= 0) call solver_error('Compute_Electrostatic_Potential> Allocation failed.')
    grid = cmplx(0.0_wp,0.0_wp,kind=wp)

    r_scale = real(ngrid,wp)/box
    call b_spline_coef(nspline, 0.0_wp, integer_weights)
    integer_weights = integer_weights/sum(integer_weights)

    ! Deposit charges in deterministic atom order, wrapping indices rather
    ! than moving or modifying the caller's coordinates.
    do atom = 1, size(charge)
      do axis = 1, 3
        u = modulo((coord(axis,atom)-origin(axis))*r_scale(axis), &
                   real(ngrid(axis),wp))
        if (u >= real(ngrid(axis),wp)) u = 0.0_wp
        base(axis) = int(u)
        fraction = u - real(base(axis),wp)
        call b_spline_coef(nspline, fraction, weights(:,axis))
        bnormal = sum(weights(:,axis))
        weights(:,axis) = weights(:,axis)/bnormal
      end do
      do iz = 1, nspline
        gz = modulo(base(3)-iz+1,ngrid(3)) + 1
        do iy = 1, nspline
          gy = modulo(base(2)-iy+1,ngrid(2)) + 1
          do ix = 1, nspline
            gx = modulo(base(1)-ix+1,ngrid(1)) + 1
            grid(gx,gy,gz) = grid(gx,gy,gz) + charge(atom) &
                              *weights(ix,1)*weights(iy,2)*weights(iz,3)
          end do
        end do
      end do
    end do

    call transform_grid(grid, ngrid, -1)

    do axis = 1, 3
      do i = 1, ngrid(axis)
        m = i - 1
        if (m > ngrid(axis)/2) m = m - ngrid(axis)
        wave2(i,axis) = (2.0_wp*PI*real(m,wp)/box(axis))**2
        transfer = cmplx(0.0_wp,0.0_wp,kind=wp)
        do l = 1, nspline
          angle = 2.0_wp*PI*real(m,wp)*real(l-1,wp)/real(ngrid(axis),wp)
          transfer = transfer + integer_weights(l) &
                     *cmplx(cos(angle),sin(angle),kind=wp)
        end do

        ! Odd-order assignment has no Nyquist component.  The simulation's
        ! PME kernels omit this component as well; never divide by its zero.
        if (mod(nspline,2) == 1 .and. mod(ngrid(axis),2) == 0 .and. &
            m == ngrid(axis)/2) then
          correction(i,axis) = cmplx(0.0_wp,0.0_wp,kind=wp)
        else
          correction(i,axis) = 1.0_wp/transfer
          if (m == 0 .or. &
              (mod(ngrid(axis),2) == 0 .and. m == ngrid(axis)/2)) &
            correction(i,axis) = cmplx(real(correction(i,axis),wp), &
                                      0.0_wp,kind=wp)
        end if
      end do
    end do

    prefactor = 4.0_wp*PI*ELECOEF &
                 /(product(box)*dielectric*kcal_per_mol_per_ev)
    do k = 1, ngrid(3)
      do j = 1, ngrid(2)
        do i = 1, ngrid(1)
          k2 = wave2(i,1) + wave2(j,2) + wave2(k,3)
          if (i == 1 .and. j == 1 .and. k == 1) then
            grid(i,j,k) = cmplx(0.0_wp,0.0_wp,kind=wp)
          else
            factor = prefactor*exp(-k2/(4.0_wp*alpha*alpha))/k2
            grid(i,j,k) = grid(i,j,k)*factor &
                          *correction(i,1)*correction(j,2)*correction(k,3)
          end if
        end do
      end do
    end do

    call transform_grid(grid, ngrid, 1)

    if (.not. all(electrostatic_is_finite(real(grid,wp))) .or. &
        .not. all(electrostatic_is_finite(aimag(grid)))) &
      call solver_error('Compute_Electrostatic_Potential> Nonfinite Fourier grid.')

    allocate(potential(ngrid(1),ngrid(2),ngrid(3)), stat=alloc_stat)
    if (alloc_stat /= 0) call solver_error('Compute_Electrostatic_Potential> Allocation failed.')
    potential = real(grid,wp)
    grid_mean = sum(potential)/product(real(ngrid,wp))
    potential = potential - grid_mean
    if (.not. all(electrostatic_is_finite(potential))) &
      call solver_error('Compute_Electrostatic_Potential> Nonfinite grid potential.')

  end subroutine compute_electrostatic_potential

  ! Integer exponent tests remain effective under the -ffast-math flags
  ! used by GENESIS, which otherwise optimize IEEE_IS_FINITE to true.
  ! GENESIS wp is IEEE binary32 or binary64 in the supported builds.
  elemental logical function electrostatic_is_finite(value)
    real(wp), intent(in) :: value
    integer(int32) :: bits32
    integer(int64) :: bits64
    if (storage_size(value) == 32) then
      bits32 = transfer(value, 0_int32)
      electrostatic_is_finite = &
        iand(bits32,int(z'7F800000',int32)) /= int(z'7F800000',int32)
    else if (storage_size(value) == 64) then
      bits64 = transfer(value, 0_int64)
      electrostatic_is_finite = &
        iand(bits64,int(z'7FF0000000000000',int64)) /= &
             int(z'7FF0000000000000',int64)
    else
      electrostatic_is_finite = .false.
    end if
  end function electrostatic_is_finite

  subroutine solver_error(message)
    character(*), intent(in) :: message
#ifdef HAVE_MPI_GENESIS
    logical :: initialized
    integer :: status
#endif
    write(ErrOut,'(A)') trim(message)
    flush(ErrOut)
#ifdef HAVE_MPI_GENESIS
    call mpi_initialized(initialized, status)
    if (initialized) then
      call mpi_abort(mpi_comm_country, 1, status)
    end if
#endif
    error stop 1
  end subroutine solver_error

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    transform_grid
  !> @brief        Deterministic separable, unnormalized 3D complex FFT.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transform_grid(grid, ngrid, direction)

    complex(wp),             intent(inout) :: grid(:,:,:)
    integer,                 intent(in)    :: ngrid(3), direction

    complex(wp), allocatable :: line(:), work(:)
    integer                  :: i, j, k, nmax, alloc_stat


    nmax = maxval(ngrid)
    allocate(line(nmax), work(2*nmax), stat=alloc_stat)
    if (alloc_stat /= 0) call solver_error('Transform_Grid> Allocation failed.')

    call zfft1d(line, ngrid(1), 0, work)
    do k = 1, ngrid(3)
      do j = 1, ngrid(2)
        line(1:ngrid(1)) = grid(:,j,k)
        call zfft1d(line, ngrid(1), direction, work)
        grid(:,j,k) = line(1:ngrid(1))
      end do
    end do

    call zfft1d(line, ngrid(2), 0, work)
    do k = 1, ngrid(3)
      do i = 1, ngrid(1)
        line(1:ngrid(2)) = grid(i,:,k)
        call zfft1d(line, ngrid(2), direction, work)
        grid(i,:,k) = line(1:ngrid(2))
      end do
    end do

    call zfft1d(line, ngrid(3), 0, work)
    do j = 1, ngrid(2)
      do i = 1, ngrid(1)
        line(1:ngrid(3)) = grid(i,j,:)
        call zfft1d(line, ngrid(3), direction, work)
        grid(i,j,:) = line(1:ngrid(3))
      end do
    end do

  end subroutine transform_grid

end module sp_electrostatic_mod
