Program checkhbond

    implicit none
    integer :: i, j, k, ntimes, t
    integer :: natms, nmol, atms_per_mol,cntj
    real :: dt, volume
    real, dimension(3) :: L
    real, dimension(500,3) :: rO, r1, r2
    integer, dimension(1000) :: OH,OH1
    integer, dimension(1000,6) :: Hbonds
    integer, dimension(500) :: jumps,jumps_old

    character(len=10) nfile, mol_name
    character(len=2) ctmp

    atms_per_mol = 3  ! Set # of atoms per mol !!!!!!!!!!!!

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

    close(10)
    L(1)=volume ** (1.0/3.0)
    L(2)=L(1)
    L(3)=L(1)


    open(11,file='traj_'//trim(nfile)//'_'//trim(mol_name)//'.xyz',status='old')
    open(12,file="hbonds.out")
    open(13,file="Ocoords.dat")
    open(14,file='grpjump_'//trim(nfile)//'_'//trim(mol_name)//'.dat')
    rO=0.0;r1=0.0;r2 = 0.0
    OH = 0; Hbonds = 0; OH1 = 0
    ! Read in the first configuration and set it as the zero.
    jumps_old=1
    do t = 1, ntimes
        if (mod(t-1,500) .eq. 0) then
            write(*,*) "STEP: ", t-1
        endif
        read(11,*) natms
        read(11,*)
        nmol = natms/atms_per_mol
        do j = 1, nmol
            read(11,*) ctmp, (rO(j,k),k=1,3)
            read(11,*) ctmp, (r1(j,k),k=1,3)
            read(11,*) ctmp, (r2(j,k),k=1,3)
            r1(j,:) = r1(j,:) - L(:)*anint(( r1(j,:) - rO(j,:) )/L(:))
            r2(j,:) = r2(j,:) - L(:)*anint(( r2(j,:) - rO(j,:) )/L(:))
        enddo
        call findhbonds(nmol, L, rO, r1, r2, OH)
        if (t .eq. 1) then
            OH1 = OH(:)
        endif
        jumps=jumps_old(:)
        call findjump(nmol, OH1, OH, jumps)
        jumps_old=jumps(:)
        call sorthbonds(nmol, OH, Hbonds)
        cntj=0
        do j=1,nmol
            if (jumps(j) .eq. 0) then
                cntj = cntj + 1
            endif
        enddo
        write(14,*) t, 1.0-cntj/real(nmol)
        do j=1,nmol
            write(12,*) (Hbonds(j,k), k=1,6)
            write(13,*) (rO(j,k), k=1,3)
        enddo
    enddo
    close(11)
    close(12)
    close(13)
    close(14)

End Program checkhbond

Subroutine findjump(nmol, OH1, OH,jumps)
    integer :: nmol, i, j, k, o1, o2, OHCOUNT, chk1, chk2
    integer, dimension(1000) :: OH1, OH
    integer, dimension(500) :: jumps

    ! Note: this subroutine takes the following as arguments
    ! nmol: number of water molecules
    ! OH1: array of nmol*2 that has t=0 h-bond acceptors for each OH (or zero)
    ! OH: same as above for time t
    ! jumps: nmol array that starts as 1, but then elements are flipped to 0 by
    !        this subroutine. Nothing ever increments it back to 1 - so even if 
    !        it exchanges back, it won't cause it to go back to 1.


    OHCOUNT=nmol*2
    ! Loop over water molecules
    do j=1,nmol
        chk1=0
        chk2=0
        ! sets the indices for the OHs in the OH vector
        o1=2*j-1; o2 = 2*j
        ! This checks that if at time zero an OH is not H-bonded
        ! if it isn't, it waits until it finds an H-bond at time t, and then 
        ! it resets the H-bond at time zero to be the one at time t
        if ( OH1(o1) .eq. 0 .and. OH(o1) .ne. 0) then
            OH1(o1) = OH(o1)
        endif
        ! Repeats for the second OH in the H-bond
        if ( OH1(o2) .eq. 0 .and. OH(o2) .ne. 0) then
            OH1(o2) = OH(o2)
        endif
        ! Checks if the OH was H-bonded at time zero and is also H-bonded at
        ! time t
        if ( OH1(o1) .ne. 0 .and. OH(o1) .ne. 0) then
            ! If the OH is hydrogen bonded to the same acceptor at time 0 and t, 
            ! adds 1 to chk, otherwise doesn't add anything.
            if (OH(o1) .eq. OH1(o1)) then
                chk1 = chk1 + 1
            endif
        else ! transient break condition adds 1 to track
            chk1 = chk1 + 1
        endif
        
        ! Repeated for OH2 
        ! Checks of o2 is h-bonded at time 0 and also at time t
        if ( OH1(o2) .ne. 0 .and. OH(o2) .ne. 0) then
            ! if so, checks if they are hbonded to the same molec
            ! adds 1 to chk, otherwise don't add anything
            if (OH(o2) .eq. OH1(o2)) then
                chk2 = chk2 + 1
            endif
        else ! transient break condition adds 1 to track
            chk2 = chk2 + 1
        endif

        ! chk=0 means exchange ocurred for the first oh in mol j
        ! increments jumps by setting the donor molecule j, 
        ! the original acceptor OH1(o1),
        ! and the new acceptor OH(o1) 
        ! all to zero.
        if (chk1 .ne. 1) then
            jumps(j)=0
            jumps(OH(o1))=0
            jumps(OH1(o1))=0
        endif
        ! Repeated for the second OH of molecule j. 
        if (chk2 .ne. 1) then
            jumps(j)=0
            jumps(OH(o2))=0
            jumps(OH1(o2))=0
        endif
    enddo

End Subroutine
    


Subroutine sorthbonds(nmol, OH, Hbonds)
    integer :: OHCOUNT, nmol, oh1, oh2, dcnt
    integer, dimension(1000) :: OH
    integer, dimension(1000,6) :: Hbonds
    ! mol donates donates accepts accepts accepts
    OHCOUNT=nmol*2
    Hbonds=0
    do j = 1, nmol
        oh1 = 2*j-1; oh2 =2*j
        Hbonds(j,1) = j
        Hbonds(j,2) = OH(oh1)
        Hbonds(j,3) = OH(oh2)
        dcnt = 0
        do k=1, OHCOUNT
            if ( OH(k) == j .and. dcnt+4 .le. 6 ) then
                Hbonds(j,4+dcnt) = k
                dcnt = dcnt + 1
                if ( dcnt > 3 ) then
                    write(*,*) "Error: dcnt too large"
                endif
            endif
        enddo
    enddo
End Subroutine

Subroutine findhbonds(nmol,L, RO, R1, R2, OH)

    implicit none
    integer :: j, k, nmol, OHCOUNT, ohtmp,hbnds
    real :: dOO, dOH1, dOH2, rOOmax, rOHmax, angmax, norm, ang
    real,dimension(3) :: L, eootmp, eohtmp, eoo, eoh,dOOtmp, dOH1tmp, dOH2tmp
    real, dimension(500,3) :: RO, R1, R2
    integer, dimension(1000) :: OH
   
    rOOmax = 3.1
    rOHmax = 2.05
    angmax = 15
    eoo = 0.0; eoh = 0.0; dOOtmp = 0.0; dOH1tmp = 0.0; dOH2tmp =0.0
    dOO = 0.0; dOH1 = 0.0; dOH2 =0.0; ang = 0.0
    OHCOUNT=nmol*2
    OH=0
    hbnds = 0
    do j = 1, nmol
        ohtmp = 2*j-1
        do k = 1, nmol
            if ( k .NE. j ) then 
                dOOtmp(:) = RO(k,:)-RO(j,:) - L(:)*anint(( RO(k,:) - RO(j,:) )/L(:))
                dOO = sqrt(dot_product(dOOtmp,dOOtmp))
                ! Checks the OO distance
                if (dOO < rOOmax) then
                    dOH1tmp(:) = RO(k,:)-R1(j,:) - L(:)*anint(( RO(k,:) - R1(j,:) )/L(:))
                    dOH2tmp(:) = RO(k,:)-R2(j,:) - L(:)*anint(( RO(k,:) - R2(j,:) )/L(:))
                    dOH1 = sqrt(dot_product(dOH1tmp,dOH1tmp))
                    dOH2 = sqrt(dot_product(dOH2tmp,dOH2tmp))
                    ! Checks the OH1 distance
                    if ( dOH1 < rOHmax ) then
                        eoh(:) = dOH1tmp(:)/dOH1
                        eoo(:) = dOOtmp(:)/dOO
                        ang = acosd(dot_product(eoo,eoh))
                        ! Checks the angle
                        if ( ang < angmax) then
                            OH(ohtmp) = k
                            !write(*,*) hbnds, j, k
                            hbnds = hbnds + 1
                        endif
                    ! Checks the OH2 distance
                    else if ( dOH2 < rOHmax ) then
                        eoh(:) = dOH2tmp(:)/dOH2
                        eoo(:) = dOOtmp(:)/dOO
                        ang = acosd(dot_product(eoo,eoh))
                        ! Checks the angle
                        if ( ang < angmax) then
                            OH(ohtmp+1) = k
                            !write(*,*) hbnds, j, k
                            hbnds = hbnds +1 
                        endif
                    endif
                endif 
            endif
        enddo
    enddo

    !hbnds = 0
    !do j=1, 686
    !    if ( OH(j) .ne. 0 ) then
    !        hbnds = hbnds + 1
    !    endif
    !enddo
        
    !write(*,*) hbnds


End Subroutine
