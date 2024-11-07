program reorder_xyz
  implicit none
  integer, parameter :: natoms = 3192
  integer, parameter :: nSiOH = 96, nSi_fw = 128, nO_fw = 304, nCO2 = 90
  integer :: i, j, t, set, ios
  character(len=10) :: atom_type
  real(8) :: x, y, z
  character(len=100) :: line
  integer :: timestep
  character(len=20) :: filename_in, filename_out
  character(len=4), dimension(natoms) :: atom_types
  real(8), dimension(natoms, 3) :: coordinates
  integer :: nSiOH_count, nSi_fw_count, nO_fw_count, nCO2_count

  ! File names
  filename_in = "nvt.xyz"
  filename_out = "reordered_nvt.xyz"

  open(unit=10, file=filename_in, status="old")
  open(unit=20, file=filename_out, status="unknown")

  ! Loop over timesteps
  do t = 1, 80001
     read(10, *)
     read(10, *)

     ! Write the timestep data in reordered format for this set
     write(20, *) natoms
     write(20, '(A)', advance="no") "Atoms. Timestep: "
     write(20, '(I0)') 

     ! Read all atoms in the current set of the timestep
     do i = 1, natoms
        read(10, *) atom_type, x, y, z
        atom_types(i) = trim(adjustl(atom_type))
        coordinates(i, :) = (/ x, y, z /)
     end do

     print*, atom_types(1), coordinates(1,:)
     print*, atom_types(799), coordinates(799,:)

     ! Write SiOH atoms in 4 sets with a shift in `nSiOH_count` after each set
     nSiOH_count = 0

     do j = 1, 4  ! Loop over 4 replicates
        do i = 1, nSiOH
           nSiOH_count = nSiOH_count + 1
           write(20, '(A,3F10.5)') atom_types(nSiOH_count), coordinates(nSiOH_count, 1), coordinates(nSiOH_count, 2), coordinates(nSiOH_count, 3)
        end do

        print*,nSiOH_count
        nSiOH_count = nSiOH_count + (798 - nSiOH)  ! Apply the offset for the next set
     end do


     ! Write framework Si atoms in 4 sets with a shift in `Si` after each set
     nSi_fw_count = nSiOH

     do j = 1, 4 ! loop over 4 replicates
        do i = 1, nSi_fw
           nSi_fw_count = nSi_fw_count + 1
           write(20, '(A,3F10.5)') atom_types(nSi_fw_count), coordinates(nSi_fw_count, 1), coordinates(nSi_fw_count, 2), coordinates(nSi_fw_count, 3)
        end do
        print*,nSi_fw_count
        nSi_fw_count = nSi_fw_count + (798 - nSi_fw)
     end do


     ! Write framework O atoms in 4 sets with a shift in `O` after each set
     nO_fw_count = nSiOH + nSi_fw

     do j = 1, 4 ! loop over the 4 replicates
        do i = 1, nO_fw
           nO_fw_count = nO_fw_count + 1
           write(20, '(A,3F10.5)') atom_types(nO_fw_count), coordinates(nO_fw_count, 1), coordinates(nO_fw_count, 2), coordinates(nO_fw_count, 3)
        end do
        nO_fw_count = nO_fw_count + (798 - nO_fw)
     end do

     ! Write CO2 molecules 
     nCO2_count = nSiOH + nSi_fw + nO_fw

     do j = 1, 4 ! loop over the 4 replicates
        do i = 1, 3 * nCO2  ! Each CO2 has 3 atoms (C and 2 O)
           nCO2_count = nCO2_count + 1
           write(20, '(A,3F10.5)') atom_types(nCO2_count), coordinates(nCO2_count, 1), coordinates(nCO2_count, 2), coordinates(nCO2_count, 3)
        end do
        nCO2_count = nCO2_count + (798 - (3*nCO2))
     end do

  end do

  close(10)
  close(20)
end program
