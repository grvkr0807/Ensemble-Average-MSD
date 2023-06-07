

Program EnsembleAverageMSD
 implicit none

 integer, parameter :: fh = 15
 character(len=300) :: buffer, label
 integer :: pos
 integer :: ios=0
 integer :: line=0

 ! To be read from CONTROL file-------------------------------------------------------
 integer(kind=4) :: concentration
 integer(kind=4) :: AtomsPerMolecule, MoleculePerBox
 integer(kind=4) :: FirstStep, LastStep, TimeStep, Step_Equilibration, LastStep_MSD
 character(len=300) :: InFile, OutFile, InTemp, OutTemp
 ! ------------------------------------------------------------------------------------

 integer(kind=4) :: Natoms
 integer(kind=4) :: i_LastStep, i_EquilibrationStep, Data_Steps, MSD_Steps
 integer(kind=4) :: iat, istep, ref_step, index2, index1   ! Loop variables
 integer(kind=4) :: AtomID, mol, type
 real(kind=8) :: q, pos_junk
 character(len=50) :: junk
 real(kind=8), allocatable :: coord(:,:,:), msd(:)

 !----------- READ CONTROL FILE -------------------

 open(fh, file='CONTROL.txt')
 do while (ios == 0)
    read(fh, '(A)', iostat=ios) buffer
    if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('concentration')
            read(buffer, *, iostat=ios) concentration
        case ('AtomsPerMolecule')
            read(buffer, *, iostat=ios) AtomsPerMolecule
        case ('MoleculePerBox')
            read(buffer, *, iostat=ios) MoleculePerBox
        case ('FirstStep')
            read(buffer, *, iostat=ios) FirstStep
        case ('TimeStep')
            read(buffer, *, iostat=ios) TimeStep
        case ('LastStep')
            read(buffer, *, iostat=ios) LastStep
        case ('InFile')
            read(buffer, *, iostat=ios) InTemp
        case ('OutFile')
            read(buffer, *, iostat=ios) OutTemp
        case ('Step_Equilibration')
            read(buffer, *, iostat=ios) Step_Equilibration
        case ('LastStep_MSD')
            read(buffer, *, iostat=ios) LastStep_MSD
        case default
            print *, 'Skipping invalid label at line', line
        end select
    end if
 end do
 close(fh)

 !----------FINISHED READING CONTROL FILE ------------------

 InFile= trim(adjustl(InTemp))
 OutFile= trim(adjustl(OutTemp))


 Natoms= AtomsPerMolecule*MoleculePerBox

 i_LastStep= ((LastStep - FirstStep)/TimeStep) + 1
 i_EquilibrationStep= ((Step_Equilibration - FirstStep)/TimeStep)
 Data_Steps= i_LastStep - i_EquilibrationStep
 MSD_Steps= LastStep_MSD/TimeStep - Step_Equilibration/TimeStep

 allocate(coord(Natoms,3,Data_Steps))
 allocate(msd(MSD_Steps))

 !---------Ignore coordinates of equilibration part--------

 open(22, file=InFile, status='old')
 do istep= 1, i_EquilibrationStep

    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    do iat= 1,Natoms

        read(22, *) AtomID, mol, type, q, pos_junk, pos_junk, pos_junk
        
    enddo

 enddo

! Read the coordinates post-equilibration -------------
 do istep= 1, Data_Steps

    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    read(22, *) junk
    do iat= 1,Natoms

        read(22, *) AtomID, mol, type, q, coord(iat,1,istep), coord(iat,2,istep), coord(iat,3,istep)
        
    enddo

 enddo
 close(22)
 ! ------Finished reading LAMMPS data -----------

 open(23, file=OutFile, status='new')

 do istep= 1, MSD_Steps
    msd(istep)= 0.0d0
    index1= 1
    index2= istep+index1
    do while (index2.le.Data_Steps)
        do iat= 1, Natoms
            msd(istep)= msd(istep) + (coord(iat,1,index2) - coord(iat,1,index1))**2 + &
            (coord(iat,2,index2) - coord(iat,2,index1))**2 + (coord(iat,3,index2) - coord(iat,3,index1))**2
        enddo
        index1= index1 + 1
        index2= istep + index1
    enddo

    msd(istep)= msd(istep)/Natoms
    msd(istep)= msd(istep)/(index1 - 1)
    write(23, '(I8, F18.4)') (istep)*TimeStep, msd(istep)

 enddo
 close(23)

end program EnsembleAverageMSD
