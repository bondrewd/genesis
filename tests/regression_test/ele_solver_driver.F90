! Numerical test driver; link against the built SPDYN analysis solver and lib.a.
program ele_solver_driver
  use constants_mod
  use sp_electrostatic_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64
  implicit none
  integer :: ngrid(3), nspline, natom, i, ix, iy, iz, forcefield
  real(wp) :: origin(3), box(3), alpha, dielectric
  real(wp), allocatable :: coord(:,:), charge(:), potential(:,:,:)
  real(wp) :: nonfinite
  character(32) :: mode
  integer(int32), parameter :: invalid32(4) = &
    [int(z'7FC00001',int32), int(z'7F800001',int32), &
     int(z'7F800000',int32), int(z'FF800000',int32)]
  integer(int64), parameter :: invalid64(4) = &
    [int(z'7FF8000000000001',int64), int(z'7FF0000000000001',int64), &
     int(z'7FF0000000000000',int64), int(z'FFF0000000000000',int64)]

  call get_command_argument(1,mode)
  if (trim(mode) == 'finite-check') then
    do i = 1, 4
      if (storage_size(nonfinite) == 32) then
        nonfinite = transfer(invalid32(i),0.0_wp)
      else
        nonfinite = transfer(invalid64(i),0.0_wp)
      end if
      if (electrostatic_is_finite(nonfinite)) error stop 1
    end do
    if (.not. electrostatic_is_finite(0.0_wp)) error stop 1
    if (.not. electrostatic_is_finite(1.0_wp)) error stop 1
    write(*,*) 'PASS quiet/signaling NaN and both infinities rejected.'
    stop
  end if

  read(*,*) ngrid, nspline, alpha, dielectric, natom, forcefield
  read(*,*) origin, box
  select case (forcefield)
  case (1)
    ELECOEF = ELECOEF_CHARMM
  case (2)
    ELECOEF = ELECOEF_AMBER
  case (3)
    ELECOEF = ELECOEF_GROMACS
  end select
  allocate(coord(3,natom), charge(natom))
  do i = 1, natom
    read(*,*) coord(:,i), charge(i)
  end do
  call compute_electrostatic_potential(coord, charge, ngrid, origin, box, &
                                      alpha, dielectric, nspline, potential)
  do iz = 1, ngrid(3)
    do iy = 1, ngrid(2)
      do ix = 1, ngrid(1)
        write(*,'(ES25.16E3)') potential(ix,iy,iz)
      end do
    end do
  end do
end program ele_solver_driver
