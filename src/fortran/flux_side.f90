!   Program that calculates the flux-side correlation function
!       and averages it to get the total.
!


Program flux_side
    
    implicit none
    integer :: i, j, k, it
    integer :: heaviside
    integer :: ntimes
    integer :: counter
    integer :: natms, nmol, atms_per_mol
    integer :: progkey
    real :: dt, volume
    real :: dx, dy, dz, dr
    real :: dvx, dvy, dvz, dvs
    real :: vr,vs, eval
    real :: constraint

    real, dimension(3) :: L
    real, dimension(0:5000) :: fsc
    real, dimension(0:5000) :: r
    real, dimension(2,3) :: iro,viro
    real, dimension(2,3) :: ir

    character(len=10) nfile, mol_name
    character(len=2) ctmp
    
    atms_per_mol = 3 ! Set # of atoms per mol

    ! Read in the input parameters for the calculation
    open(10,file='corr_calc.in',status='old')

    read(10,*)
    read(10,*) nfile
    read(10,*)
    read(10,*) ntimes, dt
    read(10,*)
    read(10,*) volume
    read(10,*)
    read(10,*) mol_name
    read(10,*)
    read(10,*) constraint
    read(10,*)
    read(10,*) progkey
    close(10)

    L(1) = volume ** (1.0/3.0)
    L(2) = L(1)
    L(3) = L(1)
    
    ! Open Trajectory File
    open(11, file='traj_'//trim(nfile)//'_'//trim(mol_name)//'.xyz',status='old')
    open(12, file='vel_'//trim(nfile)//'_'//trim(mol_name)//'.vxyz',status='old')

    ! Zero the results for the NVE trajectory
    fsc=0.0

    ! Read in the first configuration and set it as the zero.
    read(11,*) natms
    read(11,*)
    nmol = natms/atms_per_mol
    write(*,*) nmol
    do j = 1, nmol
        do k = 1, atms_per_mol
            read(11,*)
        enddo
    enddo
    read(11,*) ctmp, (iro(1,k), k=1,3)
    read(11,*) ctmp, (iro(2,k), k=1,3)
    write(*,*) iro(1,1)
    ! Calculate separation
    dx = iro(1,1) - iro(2,1)
    dy = iro(1,2) - iro(2,2)
    dz = iro(1,3) - iro(2,3)
    dx = dx - L(1)*anint(dx/L(1))
    dy = dy - L(2)*anint(dy/L(2))
    dz = dz - L(3)*anint(dz/L(3))

    dr = sqrt(dx**2+dy**2+dz**2)
    r(0) = dr

    ! Read in first set of velocities
    if (progkey == 1) then
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        write(*,*) nmol
        write(*,*) atms_per_mol
        do j = 1, nmol
            do k = 1, atms_per_mol
                read(12,*)
            enddo
        enddo
        read(12,*) (viro(1,k), k=1,3)
        read(12,*) (viro(2,k), k=1,3)
    else if (progkey == 2) then
        read(12,*)
        read(12,*)
        write(*,*) nmol
        write(*,*) atms_per_mol
        do j = 1, nmol
            do k = 1, atms_per_mol
                read(12,*)
            enddo
        enddo
        read(12,*) ctmp, (viro(1,k), k=1,3)
        read(12,*) ctmp, (viro(2,k), k=1,3)
    end if
    close(12)
    ! Calculates the distances
    dvx = viro(1,1) - viro(2,1)
    dvy = viro(1,2) - viro(2,2)
    dvz = viro(1,3) - viro(2,3)
    vr = dvx*dx/dr + dvy*dy/dr + dvz*dz/dr
    vs = -vr
    eval = dr - constraint
    write(*,*) dr, constraint
    write(*,*) eval
    write(*,*) vr, vs
    if (eval .ge. 0.0) then
        fsc(0) = 0.0
    else
        fsc(0) = 1.0
    end if
    open(22, file='kappa.dat')
    if (vs .ge. 0.0) then
        write(22,'(2E15.5)') fsc(0)*vs, 0.0
    else
        write(22,'(2E15.5)') 0.0, fsc(0)*vs
    end if
    close(22)

    
    ! Loop over the timesteps in each trajectory
    do it = 1, ntimes-1

        ! Read in configuration at timestep it
        read(11,*)
        read(11,*)
        counter = 0
        do j = 1, nmol
            do k = 1, atms_per_mol
                read(11,*)
            enddo
        enddo
        read(11,*) ctmp, (ir(1,k), k=1,3)
        read(11,*) ctmp, (ir(2,k), k=1,3)
        ! Calculates the distances
        dx = ir(1,1) - ir(2,1)
        dy = ir(1,2) - ir(2,2)
        dz = ir(1,3) - ir(2,3)
        dx = dx - L(1)*anint(dx/L(1))
        dy = dy - L(2)*anint(dy/L(2))
        dz = dz - L(3)*anint(dz/L(3))
        dr = sqrt(dx**2+dy**2+dz**2)
        r(it) = dr
        eval = dr - constraint
        if (eval .ge. 0) then
            fsc(it) = 0.0
        else
            fsc(it) = 1.0 
        end if
    enddo ! end loop over timesteps
    close(11)
    fsc = fsc*vs
    ! Write out the FS Correlation Function
    open(21,file='fsc_'//trim(nfile)//'_'//trim(mol_name)//'.dat')

    do it = 0, ntimes - 1
        write(21,'(3E15.5)') real(it), fsc(it), r(it)
    enddo
    close(21)
    open(22, file='kappa.dat')
    do it = 0, ntimes - 1 
        if (vs .ge. 0.0) then
            write(22,'(2E15.5)') vs, 0.0
        else
            write(22,'(2E15.5)') 0.0, vs
        end if
    enddo
    close(22)

End Program flux_side

