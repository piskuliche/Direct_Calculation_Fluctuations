Program msd_rot_calc

    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, &
        error_unit, iostat_end, iostat_eor
    implicit none
    integer :: i, j, k, a, it, orig,bl
    integer :: ntimes, or_int
    integer :: natms, nmol, atms_per_mol
    integer :: startskip, endskip, startconfig, endconfig, nconfigs, ntos
    integer :: nblocks, shamblock, nperblock
    integer :: msdindex, c2index1, c2index2
    integer :: com, counti
    integer :: tospblock
    real    :: dt, eavg, volume
    real    :: dp, dsq
    real    :: MASS, cm_tmpmsd, norm
   
    integer    :: ioerr
    logical    :: needblock

    real, dimension(3) :: L
    real, dimension(3) :: etmp
    real, dimension(10) :: M, tmpmsd
    real, dimension(5000,1000,10,3) :: r
    real, dimension(5000,1000,3) :: r_cm
    real, dimension(1000,10,3) :: r_old, shift
    real, dimension(1000,3) :: r_cm_old, shift_cm
    real, dimension(5000,1000,3) :: e
    real, dimension(0:5000,10) :: msd
    real, dimension(0:5000) :: c1, c2
    real, dimension(0:5000) :: msd_cm
    

    character(len=10) nfile, mol_name
    character(len=2) ctmp
    character(len=2) blstr, shamstr

    ! Read in the input parameters for the calculation
    open(10,file='corr_calc.in',status='old')

    read(10,*)
    read(10,*) nfile
    read(10,*)
    read(10,*) ntimes, or_int
    read(10,*)
    read(10,*) volume
    read(10,*)
    read(10,*) mol_name
    ! mol_type = 1 for water, acn
    ! mol_type = 2 for co2
    close(10)

    M = 0.0
    MASS = 0.0
    open(11,file='../../'//trim(mol_name)//'.txt',status='old')
    read(11,*) atms_per_mol
    read(11,*) nmol
    read(11,*) startskip, endskip
    read(11,*) c2index1, c2index2
    do i = 1, atms_per_mol
        read(11,*) M(i)
        MASS = MASS + M(i)
    enddo
    close(11)
    write( unit=output_unit,*) atms_per_mol 

    L(1) = volume ** (1.0/3.0)
    L(2) = L(1)
    L(3) = L(1)


    ! Read the trajectory file
    open(12,file='traj_'//trim(nfile)//'_'//trim(mol_name)//'.xyz',status='old') ! Open trajecotry file
    write( unit=output_unit,*) 'Reading file ', trim(nfile)
    
    r_cm = 0.0
    r = 0.0
    c1 = 0.0; c2 = 0.0
    ! This reads nconfigs frames
    do i=1, ntimes
        if (MOD(i,5000) == 0) then
            write(*,*) 'Configuration ',i,'read in.'
        end if
        do j=1, startskip
            read(12,*)
        enddo
        do j=1, nmol
            do k=1, atms_per_mol
                read(12,*) ctmp, (r(i,j,k,a), a=1,3)
                if ( k /= 1 ) then
                    ! Read in w/ periodic boundary conditions
                    r(i,j,k,:) = r(i,j,k,:) - L(:)*anint(( r(i,j,k,:) - &
                        r(i,j,1,:) )/L(:))
                end if
                ! Center of MAss Calc
                r_cm(i,j,:) = r_cm(i,j,:)+r(i,j,k,:)*M(k)/MASS
            enddo
            r_cm(i,j,:) = r_cm(i,j,:)+r(i,j,k,:)*M(k)/MASS
            etmp(:) = r(i,j,c2index2,:) - r(i,j,c2index1,:)
            norm = dot_product(etmp, etmp)
            e(i,j,:) = etmp(:)/sqrt(norm)
        enddo
        do j=1, endskip
            read(12,*)
        enddo
    enddo
    close(12)

    c1(0) = real(nmol)
    c2(0) = real(nmol)
    
    ! Loop over time origins
    msd=0.0
    msd_cm=0.0
    shift = 0.0
    shift_cm = 0.0
    ! Loop over the timesteps in each trajectory
    do it=2,ntimes
        i=1
        ! Loop over molecules
        do j = 1, nmol
            do k=1, atms_per_mol
                ! Calculate the shift if it occurs.
                shift(j,k,:) = shift(j,k,:) - L(:)*anint((r(it,j,k,:) - &
                    r(it-1,j,k,:) )/L(:))
                ! Calculate the square displacements
                dsq = ( r(it,j,k,1) + shift(j,k,1) - r(i,j,k,1) ) ** 2. &
                     +( r(it,j,k,2) + shift(j,k,2) - r(i,j,k,2) ) ** 2. &
                     +( r(it,j,k,3) + shift(j,k,3) - r(i,j,k,3) ) ** 2.
                msd(it-1, k) = msd(it-1, k) + dsq
                ! Calculate the contribution to the c1,c2
            enddo ! End Atoms Loop (k)
            ! Calculate the shift if it occurs.
            shift_cm(j,:) = shift_cm(j,:) - L(:)*anint((r_cm(it,j,:) - &
                            r_cm(it-1,j,:) )/L(:))
            ! Calculate the square displacements
            dsq = ( r_cm(it,j,1) + shift_cm(j,1) - r_cm(i,j,1) ) ** 2. &
                +( r_cm(it,j,2) + shift_cm(j,2) - r_cm(i,j,2) ) ** 2. &
                +( r_cm(it,j,3) + shift_cm(j,3) - r_cm(i,j,3) ) ** 2.
            msd_cm(it-1) = msd_cm(it-1) + dsq
            dp  = e(it,j,1)*e(i,j,1) + e(it,j,2)*e(i,j,2) + e(it,j,3)*e(i,j,3)
            c1(it-1) = c1(it-1) + dp
            c2(it-1) = c2(it-1) + 0.5*(3.0*dp**2-1.0)
        enddo ! End Molecules Loop (j)
        r_old(:,:,:) = r(it,:,:,:)
        r_cm_old(:,:) = r_cm(it,:,:)
   enddo ! End t's loop (it)

    ! Write out the MSD, C1, C2
    open(21,file='c1_'//trim(nfile)//'_'//trim(mol_name)//'.dat')  !open C1(t) file for this trajectory
    open(22,file='c2_'//trim(nfile)//'_'//trim(mol_name)//'.dat')  !open C2(t) file for this trajectory
    open(23,file='msd_'//trim(nfile)//'_'//trim(mol_name)//'.dat') !open MSD(t) file for this trajectory

    do it = 0, ntimes - 1
        write(21,'(2F12.5)') real(it), c1(it)/real(nmol)
        write(22,'(2F12.5)') real(it), c2(it)/real(nmol)
        write(23,'(5F12.5)') real(it), msd_cm(it)/real(nmol),(msd(it,k)/real(nmol), k=1,3)
    enddo
   ! close(21)
    !close(22)
    close(23) 
END Program msd_rot_calc 
