!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_electrostatic_output_mod
!> @brief   Reproducible electrostatic potential snapshots for SPDYN MD.
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_electrostatic_output_mod

  use sp_electrostatic_mod
  use sp_energy_pme_mod, only: get_pme_grid, pme_nspline
  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_boundary_str_mod
  use sp_dynvars_str_mod
  use sp_domain_str_mod
  use sp_enefunc_str_mod
  use constants_mod
  use messages_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  public :: open_electrostatic_output
  public :: close_electrostatic_output
  public :: output_electrostatic

contains

  subroutine open_electrostatic_output(output)

    type(s_output), intent(inout) :: output
    integer :: status

    status = 0
    output%out_ele_header = .true.
    if (main_rank) &
      open(newunit=output%eleunit, file=trim(output%elefile), &
           status='new', action='write', form='formatted', iostat=status)
#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(status, 1, mpi_integer, 0, mpi_comm_country, ierror)
#endif
    if (status /= 0) &
      call electrostatic_error('Open_Electrostatic_Output> cannot create '//&
                               trim(output%elefile)//'; file may exist.')

  end subroutine open_electrostatic_output

  subroutine close_electrostatic_output(output)

    type(s_output), intent(in) :: output
    integer :: status

    status = 0
    if (main_rank) close(output%eleunit, iostat=status)
#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(status, 1, mpi_integer, 0, mpi_comm_country, ierror)
#endif
    if (status /= 0) &
      call electrostatic_error('Close_Electrostatic_Output> file I/O error.')

  end subroutine close_electrostatic_output

  ! Gather only owned atoms into global atom-id order.  A single contributing
  ! owner per atom makes the reduction independent of MPI decomposition and
  ! avoids halo duplicates.  The analysis never modifies simulation arrays.
  subroutine output_electrostatic(output, dynamics, boundary, dynvars, &
                                   domain, enefunc)

    type(s_output),   intent(inout) :: output
    type(s_dynamics), intent(in)    :: dynamics
    type(s_boundary), intent(in)    :: boundary
    type(s_dynvars),  intent(in)    :: dynvars
    type(s_domain),   intent(in)    :: domain
    type(s_enefunc),  intent(in)    :: enefunc

    real(wp), allocatable :: local_atoms(:,:), global_atoms(:,:)
    real(wp), allocatable :: potential(:,:,:)
    integer,  allocatable :: local_owners(:), global_owners(:)
    real(wp) :: box(3), center(3), origin(3), spacing(3), position(3)
    real(wp) :: configured_field(3), effective_field(3), scale(3)
    real(wp) :: applied, total
    integer :: natom, icell, iatom, id, ngrid(3), ix, iy, iz, status
    integer :: file, alloc_stat

    if (domain%num_atom_all > int(huge(natom)/4,iintegers)) &
      call electrostatic_error('Output_Electrostatic> atom count exceeds '//&
                               'MPI count range.')
    natom = int(domain%num_atom_all)
    allocate(local_atoms(4,natom), global_atoms(4,natom), &
             local_owners(natom), global_owners(natom), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call electrostatic_error('Output_Electrostatic> allocation failed.')
    local_atoms = 0.0_wp
    local_owners = 0

    do icell = 1, domain%num_cell_local
      do iatom = 1, domain%num_atom(icell)
        id = domain%id_l2g(iatom,icell)
        if (id < 1 .or. id > natom) &
          call electrostatic_error('Output_Electrostatic> invalid atom ID.')
        local_atoms(1:3,id) = real(domain%coord(1:3,iatom,icell),wp)
        local_atoms(4,id) = domain%charge(iatom,icell)
        local_owners(id) = local_owners(id) + 1
      end do
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(local_atoms, global_atoms, 4*natom, mpi_wp_real, &
                    mpi_sum, 0, mpi_comm_country, ierror)
    if (ierror /= mpi_success) &
      call electrostatic_error('Output_Electrostatic> atom reduction failed.')
    call mpi_reduce(local_owners, global_owners, natom, mpi_integer, &
                    mpi_sum, 0, mpi_comm_country, ierror)
    if (ierror /= mpi_success) &
      call electrostatic_error('Output_Electrostatic> owner reduction failed.')
#else
    global_atoms = local_atoms
    global_owners = local_owners
#endif
    deallocate(local_atoms, local_owners)
    if (.not. main_rank) return
    if (any(global_owners /= 1)) &
      call electrostatic_error('Output_Electrostatic> atoms must have '//&
                               'exactly one owner.')
    deallocate(global_owners)

    box = [real(boundary%box_size_x,wp), &
           real(boundary%box_size_y,wp), real(boundary%box_size_z,wp)]
    center = [boundary%origin_x, boundary%origin_y, boundary%origin_z]
    call get_pme_grid(domain, ngrid)
    spacing = box / real(ngrid,wp)
    ! Match the force kernels' integer half-grid shift. An odd grid axis
    ! starts half a grid spacing above the geometric cell lower corner.
    origin = center - real(ngrid/2,wp)*spacing
    do id = 1, natom
      global_atoms(1:3,id) = global_atoms(1:3,id) + center
    end do
    configured_field = enefunc%efield
    effective_field = 0.0_wp
    scale = 1.0_wp
    if (enefunc%use_efield) then
      if (enefunc%efield_normal) &
        scale = domain%system_size_ini / domain%system_size
      ! Match the existing field force exactly in physical volts: the force
      ! routine uses rounded 23.06 kcal/mol per configured V/Angstrom.
      effective_field = configured_field*scale*23.06_wp/kcal_per_mol_per_ev
    end if
    if (.not. all(electrostatic_is_finite(configured_field)) .or. &
        .not. all(electrostatic_is_finite(effective_field))) &
      call electrostatic_error('Output_Electrostatic> Nonfinite electric field.')

    call compute_electrostatic_potential(global_atoms(1:3,:), &
                                         global_atoms(4,:), ngrid, origin, &
                                         box, enefunc%pme_alpha, &
                                         enefunc%dielec_const, pme_nspline, &
                                         potential)
    file = output%eleunit
    if (output%out_ele_header) then
      write(file,'(A)',iostat=status) &
        '#length_unit=angstrom', &
        '#time_unit=ps', &
        '#potential_unit=V'
      call check_electrostatic_io(status)
      call comment_real(file, 'ewald_alpha', real(enefunc%pme_alpha,dp))
      call comment_real(file, 'dielectric', real(enefunc%dielec_const,dp))
      write(file,'(A)',iostat=status) &
        'step,time,x,y,z,charge_potential,external_potential,total_potential'
      call check_electrostatic_io(status)
      output%out_ele_header = .false.
    end if

    do ix = 1, ngrid(1)
      position(1) = origin(1) + real(ix-1,wp)*spacing(1)
      do iy = 1, ngrid(2)
        position(2) = origin(2) + real(iy-1,wp)*spacing(2)
        do iz = 1, ngrid(3)
          position(3) = origin(3) + real(iz-1,wp)*spacing(3)
          applied = -dot_product(effective_field,position-center)
          total = potential(ix,iy,iz) + applied
          write(file,'(I0,7(",",ES24.16E3))',iostat=status) &
            dynvars%step, dynvars%time, position, potential(ix,iy,iz), &
            applied, total
          call check_electrostatic_io(status)
        end do
      end do
    end do
    flush(file,iostat=status)
    call check_electrostatic_io(status)

  end subroutine output_electrostatic

  subroutine comment_real(file, key, value)
    integer, intent(in) :: file
    character(*), intent(in) :: key
    real(dp), intent(in) :: value
    integer :: status
    write(file,'(A,A,A,ES24.16E3)',iostat=status) '#',key,'=',value
    call check_electrostatic_io(status)
  end subroutine comment_real

  subroutine check_electrostatic_io(status)
    integer, intent(in) :: status
    if (status /= 0) &
      call electrostatic_error('Output_Electrostatic> CSV file I/O error.')
  end subroutine check_electrostatic_io

  ! The generic error handler exits only the calling rank.  Root-only I/O
  ! failures must instead terminate the communicator so peers cannot hang.
  subroutine electrostatic_error(message)
    character(*), intent(in) :: message
#ifdef HAVE_MPI_GENESIS
    write(ErrOut,'(A)') trim(message)
    flush(ErrOut)
    call mpi_abort(mpi_comm_country, 1, ierror)
#endif
    call error_msg(message)
  end subroutine electrostatic_error

end module sp_electrostatic_output_mod
