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
    real :: dt, volume
    real :: dx, dy, dz, dr
    real :: dvx, dvy, dvz, dvs
    real :: vr,eval

    real, dimension(3) :: L
    real, dimension(0:5000) :: fsc
    real, dimension(0:5000) :: r
    real, dimension(2,3) :: iro,viro
    real, dimension(2,3) :: ir

    character(len=10) nfile, mol_name
    character(len=2) ctmp
    
    atms_per_mol = 3 ! Set # of atoms per mol

    ! Read in the input parameters for the calculation
    open(10,file='msd_rot_calc.in',status='old')

    read(10,*)
    read(10,*) nfile
    read(10,*)
    read(10,*) ntimes, dt
    read(10,*)
    read(10,*) volume
    read(10,*)
    read(10,*) mol_name
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
    close(12)
    ! Calculates the distances
    dvx = viro(1,1) - viro(2,1)
    dvy = viro(1,2) - viro(2,2)
    dvz = viro(1,3) - viro(2,3)
    vr = dvx*dx/dr + dvy*dy/dr + dvz*dz/dr
    eval = dr - 2.6
    if (eval .ge. 0) then
        fsc(0) = 0
    else
        fsc(0) = vr 
    end if
    

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
        eval = dr - 2.6
        if (eval .ge. 0) then
            fsc(it) = 0
        else
            fsc(it) = vr 
        end if
    enddo ! end loop over timesteps
    close(11)

    ! Write out the FS Correlation Function
    open(21,file='fsc_'//trim(nfile)//'_'//trim(mol_name)//'.dat')

    do it = 0, ntimes - 1
        write(21,'(3F12.5)') real(it), fsc(it), r(it)
    enddo
    close(21)

End Program flux_side

