Program jumpcalc
    implicit none
    integer :: i, j, k, ntimes, t,flag
    integer :: natms, nmol, atms_per_mol
    real :: dt, volume,cntzero
    real, dimension(3) :: L
    real, dimension(2000,500,3) :: rO, r1, r2
    real, dimension(2000) :: CRP
    integer, dimension(2000,1000,6) :: Hbonds
    real,dimension(2000) :: DMSD, cnt, cMSD
    character(len=10) nfile, mol_name
    character(len=2) ctmp
    
    atms_per_mol = 3  ! Set # of atoms per mol !!!!!!!!!!!!
    nmol=343
    Hbonds = 0; rO = 0.0
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
    Hbonds=0
    open(12,file="hbonds.out",status="old")
    open(13,file="Ocoords.dat",status="old")
    do t=1,ntimes
        do j=1,nmol
            read(12,*) (Hbonds(t,j,k), k=1,6)
            read(13,*) (rO(t,j,k),k=1,3)
        enddo
    enddo
    close(12)
    close(13)
    call JumpCalculation(nmol, ntimes, Hbonds, CRP)
    call JumpD(nmol,ntimes,L,Hbonds,rO,cMSD,DMSD,cntzero,cnt)
    open(14,file='jumpD_'//trim(nfile)//'_'//trim(mol_name)//'.dat')
    open(15,file='msdj_'//trim(nfile)//'_'//trim(mol_name)//'.dat')
    open(16,file='cntj_'//trim(nfile)//'_'//trim(mol_name)//'.dat')
    write(14,*) 0, 0.0
    write(15,*) 0, 0.0
    write(16,*) 0, cntzero/cntzero, cntzero
    do t=1,ntimes
        write(14,*) t, DMSD(t)
        write(15,*) t, cMSD(t)
        write(16,*) t, cnt(t)/cntzero, cnt(t)
    enddo
    close(14) 
    close(15)
    close(16)

End Program

Subroutine JumpD(nmol,ntimes,L,Hbonds,rO,cMSD,DMSD,cntzero,cnt)
    implicit none
    integer :: nmol, ntimes, j, k, t, flag
    real,dimension(3) :: L
    integer, dimension(2000,500) :: intact
    integer, dimension(2000,1000,6) :: Hbonds
    real :: cntzero
    real,dimension(2000,500) :: MSD
    real,dimension(2000) :: DMSD, cnt, cMSD
    real, dimension(2000,500,3) :: rO
    
    DMSD = 0.0; cMSD(t) = 0.0; MSD = 0.0;cnt=0.0; cntzero=0.0

    call CheckHbonds(nmol,ntimes,Hbonds,intact)
    call CalcMSD(nmol,ntimes,L,rO,intact,MSD,cnt,cntzero)
    do j=1,nmol
        do t=1,ntimes
            DMSD(t) = DMSD(t) + MSD(t,j)*intact(t+1,j)
            cMSD(t) = cMSD(t)  + MSD(t,j)
        enddo
    enddo
    do t=1,ntimes
        if ( cnt(t) .ne. 0.0 ) then
            DMSD(t) = DMSD(t)/cnt(t)
        endif
        cMSD(t) = cMSD(t)/real(nmol)
    enddo

End Subroutine

Subroutine CalcMSD(nmol,ntimes, L, rO,intact,MSD,cnt,cntzero)
    implicit none
    integer :: nmol, ntimes, j, k, t, flag
    real,dimension(3) :: L
    real :: dOsq,cntzero
    integer, dimension(2000,500) :: intact
    real,dimension(2000) :: cnt,msdk
    real,dimension(2000,500) :: MSD
    real,dimension(500,3) :: shift
    real, dimension(2000,500,3) :: rO
    character(len=2) ctmp
    shift = 0.0; MSD = 0.0; cnt=0.0; msdk = 0.0
    do j=1,nmol
        if ( intact(1,j) .eq. 1 ) then
            cntzero = cntzero + 1.0
        endif
    enddo
    do t=2,ntimes
        do j=1,nmol
            shift(j,:) = shift(j,:) - L(:)*anint(( rO(t,j,:) - rO(t-1,j,:) )/L(:))
            dOsq = ( rO(t,j,1) + shift(j,1) - rO(1,j,1) )**2 &
                 + ( rO(t,j,2) + shift(j,2) - rO(1,j,2) )**2 &
                 + ( rO(t,j,3) + shift(j,3) - rO(1,j,3) )**2
            MSD(t-1,j) = dOsq
            if ( intact(t,j) .eq. 1 ) then
                cnt(t-1) = cnt(t-1) + 1.0
            endif
        enddo
    enddo

End Subroutine

Subroutine CheckHbonds(nmol, ntimes, Hbonds, intact)
    implicit none
    integer :: nmol, ntimes, j, k, t, flag
    integer, dimension(2000,500) :: intact
    integer, dimension(2000,1000,6) :: Hbonds
    real, dimension(2000) :: valcnt
    valcnt=0
    intact=0
    do j=1,nmol
        flag=1
        ! Check that there is no change in partners
        do t=1,ntimes
            do k=2,6
                if ( Hbonds(t,j,k) .ne. Hbonds(1,j,k) ) then
                    if ( Hbonds(t,j,k) .ne. 0 ) then
                        flag=0
                    endif
                endif
            enddo
            intact(t,j) = flag
            valcnt(t) = valcnt(t) + flag
        enddo
    enddo
    open(13,file="jumpD_cnt.dat")
    do t=1,200
        write(13,*) (t-1)*10.0/1000.0, valcnt(t)/valcnt(1)
    enddo
    close(13)

End Subroutine


Subroutine JumpCalculation(nmol, ntimes, Hbonds, CRP)

    integer :: nmol, ntimes, flag, t, j
    real, dimension(2000) :: CRP
    integer, dimension(2000,1000,6) :: Hbonds

    CRP=0
    do j=1,nmol
        flag=1
        if (Hbonds(1,j,2) .ne. 0) then
            do t=1,ntimes
                if ( Hbonds(t,j,2) .ne. Hbonds(1,j,2) ) then
                    if ( Hbonds(t,j,2) .ne. 0 ) then
                        flag=0
                    endif
                endif
                CRP(t) = CRP(t) + flag
            enddo
        endif
        flag=1
        if (Hbonds(1,j,3) .ne. 0) then
            do t=1,ntimes
                if ( Hbonds(t,j,3) .ne. Hbonds(1,j,3) ) then
                    if ( Hbonds(t,j,3) .ne. 0 ) then
                        flag=0
                    endif
                endif
                CRP(t) = CRP(t) + flag
            enddo
        endif
    enddo
    CRP(:) = CRP(:)/CRP(1)
    open(13,file="jumpfort.dat")
    do t=1,ntimes
        write(13,*) t, CRP(t)
    enddo
    close(13)

End Subroutine




