!   Program to calculate the individual MSDs and rotational correlation functions
!    and average them to get the total.
!

  Program collect_msd_rot

    implicit none
    integer :: i, j, k, it, ibl
    integer :: ntimes, nfiles, nbl, bl_size
    integer :: nmol
    real :: dt, eavg, atmp, keavg, potavg, ljavg, coulavg, tval

    real, dimension(0:2000,3) :: msd, msd_tot, emsd_tot, e2msd_tot
    real, dimension(0:2000,3) :: err_msd, err_emsd, err_enermsd
    real, dimension(0:2000,3) :: kemsd_tot, potmsd_tot, ljmsd_tot, coulmsd_tot,
    volmsd_tot
    real, dimension(0:2000) :: c1, c1_tot, c2, c2_tot
    real, dimension(0:2000) :: ec1_tot, ec2_tot, err_ec2
    real, dimension(0:2000) :: kec2_tot, potc2_tot, ljc2_tot, coulc2_tot,
    volc2_tot
    real, dimension(0:2000,20) :: ec1_bl,ec2_bl,c2_bl
    real, dimension(0:2000,3,20) :: msd_bl, emsd_bl

    real, dimension(3000) :: ener, dener
    real, dimension(3000) :: ke, dke, pot, dpot, lj, dlj, coul, dcoul, ewald,
    vol, dvol
    real, dimension(3000,20) :: dener_bl
    real, dimension(20) :: eavg_bl, cnt_bl

    character(len=10), dimension(3000) :: ext

    ! Read in the input parameters for the calculation
    open(10,file='collect_msd_rot.in',status='old')

    read(10,*)
    read(10,*) nfiles
    read(10,*) 
    read(10,*) ntimes, dt
    read(10,*)
    read(10,*) nbl

    bl_size = nfiles/nbl
    call tvalue(nbl, tval)

    ! First read in the extensions for the directories and trajectory files
    open(9,file='file_names',status='old')

    do i = 1, nfiles
       read(9,*) ext(i)
    enddo
    close(9)

    ! Read in the energies for each trajectory
    open(8,file='e_init.out',status='old')
    open(81,file='ke_init.out',status='old')
    open(82,file='v_init.out',status='old')
    open(83,file='lj_init.out',status='old')
    open(84,file='coul_init.out',status='old')
    open(85,file='ewald_init.out',status='old')
    open(86,file='vol_init.out',status='old')

    eavg = 0.0; eavg_bl = 0.0; keavg = 0.0; potavg = 0.0; ljavg = 0.0; coulavg = 0.0, volavg = 0.0
    do i = 1, nfiles
       read(8,*) ener(i)
       eavg = eavg + ener(i)

       read(81,*) ke(i)
       keavg = keavg + ke(i)
       read(82,*) pot(i)
       potavg = potavg + pot(i)
       read(83,*) lj(i)
       ljavg = ljavg + lj(i)
       read(84,*) coul(i)
       read(85,*) ewald(i)
       coul(i) = coul(i) + ewald(i)
       coulavg = coulavg + coul(i) 
       read(86,*) vol(i)
       volavg = volavg + vol(i)

       ibl = int(real(i-1)/real(bl_size)) + 1
       eavg_bl(ibl) = eavg_bl(ibl) + ener(i)
    enddo
    close(8)
    close(81)
    close(82)
    close(83)
    close(84)
    close(86)
    eavg = eavg/real(nfiles)
    eavg_bl = eavg_bl/real(bl_size)
    keavg = keavg/real(nfiles)
    potavg = potavg/real(nfiles)
    ljavg = ljavg/real(nfiles)
    coulavg = coulavg/real(nfiles)
    volavg = volavg/real(nfiles)

    write(6,*) ' nfiles = ',nfiles
    write(6,*) ' eavg   = ',eavg
    write(6,*) ' keavg   = ',keavg
    write(6,*) ' potavg   = ',potavg
    write(6,*) ' ljavg   = ',ljavg
    write(6,*) ' coulavg   = ',coulavg
    write(6,*) ' volavg   = ', volavg
    do ibl = 1, nbl
       write(6,*) ibl, eavg_bl(ibl)
    enddo

    ! Calculate the energy fluctuation for each trajectory
    do i = 1, nfiles
       dener(i) = ener(i) - eavg
       dke(i) = ke(i) - keavg
       dpot(i) = pot(i) - potavg
       dlj(i) = lj(i) - ljavg
       dcoul(i) = coul(i) - coulavg
       dvol(i) = vol(i) - volavg

       ibl = int(real(i-1)/real(bl_size)) + 1
       dener_bl(i,ibl) = ener(i) - eavg_bl(ibl)
    enddo
    open(7,file='deners.dat')

    do i = 1, nfiles
       ibl = int(real(i-1)/real(bl_size)) + 1
       write(7,*) i, dener(i), dener_bl(i,ibl)
    enddo
    close(7)

    ! Zero the total results
    msd_tot = 0.0; c1_tot = 0.0; c2_tot = 0.0
    emsd_tot = 0.0; e2msd_tot = 0.0; ec1_tot = 0.0; ec2_tot = 0.0
    kemsd_tot = 0.0; potmsd_tot = 0.0; ljmsd_tot = 0.0; coulmsd_tot = 0.0
    kec2_tot = 0.0; potc2_tot = 0.0; ljc2_tot = 0.0; coulc2_tot = 0.0
    msd_bl = 0.0; emsd_bl = 0.0; ec1_bl = 0.0; ec2_bl = 0.0; c2_bl = 0.0
    cnt_bl = 0.0

    ! Then read in each file separately
    
    
    do i = 1, nfiles
       open(11,file=trim(ext(i))//'/c1_'//trim(ext(i))//'.dat',status='old')  !open c1 file
       open(12,file=trim(ext(i))//'/c2_'//trim(ext(i))//'.dat',status='old')  !open c2 file
       open(13,file=trim(ext(i))//'/msd_'//trim(ext(i))//'.dat',status='old')  !open msd file
       ibl = int(real(i-1)/real(bl_size)) + 1
       cnt_bl(ibl) = cnt_bl(ibl) + 1.0
       

       ! Read in the correlation functions
       do j = 0, ntimes - 1
          read(11,*) atmp, c1(j)
          read(12,*) atmp, c2(j)
          read(13,*) atmp, (msd(j,k), k=1,3)
       enddo
       close(11)
       close(12)
       close(13)

       c1_tot = c1_tot + c1
       c2_tot = c2_tot + c2
       msd_tot = msd_tot + msd

       ec1_tot = ec1_tot + dener(i)*c1
       ec2_tot = ec2_tot + dener(i)*c2
       emsd_tot = emsd_tot + dener(i)*msd
       e2msd_tot = e2msd_tot + dener(i)**2*msd

       ! Decomposition parts of the MSD
       kemsd_tot = kemsd_tot + dke(i)*msd
       potmsd_tot = potmsd_tot + dpot(i)*msd
       ljmsd_tot = ljmsd_tot + dlj(i)*msd
       coulmsd_tot = coulmsd_tot + dcoul(i)*msd

       ! Decomposition parts of C2
       kec2_tot = kec2_tot + dke(i)*c2
       potc2_tot = potc2_tot + dpot(i)*c2
       ljc2_tot = ljc2_tot + dlj(i)*c2
       coulc2_tot = coulc2_tot + dcoul(i)*c2

       msd_bl(:,:,ibl) = msd_bl(:,:,ibl) + msd(:,:)

       c2_bl(:,ibl) = c2_bl(:,ibl) + c2(:)

       ec1_bl(:,ibl) = ec1_bl(:,ibl) + dener_bl(i,ibl)*c1(:)
       ec2_bl(:,ibl) = ec2_bl(:,ibl) + dener_bl(i,ibl)*c2(:)
       emsd_bl(:,:,ibl) = emsd_bl(:,:,ibl) + dener_bl(i,ibl)*msd(:,:)

    enddo  ! end loop over files

    ! Normalize
    e2msd_tot = e2msd_tot/real(nfiles)
    emsd_tot = emsd_tot/real(nfiles)
    kemsd_tot = kemsd_tot/real(nfiles)
    potmsd_tot = potmsd_tot/real(nfiles)
    ljmsd_tot = ljmsd_tot/real(nfiles)
    coulmsd_tot = coulmsd_tot/real(nfiles)
    msd_tot = msd_tot/real(nfiles)

    c1_tot = c1_tot/real(nfiles)
    ec1_tot = ec1_tot/real(nfiles)
    c2_tot = c2_tot/real(nfiles)
    ec2_tot = ec2_tot/real(nfiles)
    kec2_tot = kec2_tot/real(nfiles)
    potc2_tot = potc2_tot/real(nfiles)
    ljc2_tot = ljc2_tot/real(nfiles)
    coulc2_tot = coulc2_tot/real(nfiles)

    do ibl = 1, nbl
       ec1_bl(:,ibl) = ec1_bl(:,ibl)/cnt_bl(ibl)
       ec2_bl(:,ibl) = ec2_bl(:,ibl)/cnt_bl(ibl)
       c2_bl(:,ibl) = c2_bl(:,ibl)/cnt_bl(ibl)
       emsd_bl(:,:,ibl) = emsd_bl(:,:,ibl)/cnt_bl(ibl)
       msd_bl(:,:,ibl) = msd_bl(:,:,ibl)/cnt_bl(ibl)
    enddo


    ! Calculate error bars
    err_ec2 = 0.0; err_msd = 0.0; err_emsd = 0.0
    do it = 0, ntimes - 1
       do ibl = 1, nbl
          err_ec2(it) = err_ec2(it) + (ec2_bl(it,ibl) - ec2_tot(it))**2
          err_msd(it,:) = err_msd(it,:) + (msd_bl(it,:,ibl) - msd_tot(it,:))**2
          err_emsd(it,:) = err_emsd(it,:) + (emsd_bl(it,:,ibl) - emsd_tot(it,:))**2
       enddo
       err_ec2(it) = tval*sqrt(err_ec2(it))/real(nbl)
       err_msd(it,:) = tval*sqrt(err_msd(it,:))/real(nbl)
       err_emsd(it,:) = tval*sqrt(err_emsd(it,:))/real(nbl)
    enddo

    err_enermsd = 0.0
    do it = 1, ntimes - 1
       do ibl = 1, nbl
          err_enermsd(it,1) = err_enermsd(it,1) + ( (emsd_bl(it,1,ibl)/msd_bl(it,1,ibl)) & 
               - (emsd_tot(it,1)/msd_tot(it,1)) )**2
          err_enermsd(it,2) = err_enermsd(it,2) + ( (emsd_bl(it,2,ibl)/msd_bl(it,2,ibl)) & 
               - (emsd_tot(it,2)/msd_tot(it,2)) )**2
          err_enermsd(it,3) = err_enermsd(it,3) + ( (emsd_bl(it,3,ibl)/msd_bl(it,3,ibl)) & 
               - (emsd_tot(it,3)/msd_tot(it,3)) )**2
       enddo
       err_enermsd(it,:) = tval*sqrt(err_enermsd(it,:))/real(nbl)
    enddo



    !write out the results for the total of all trajectories
    open(21,file='coll_c1_tot.dat')  !open C1(t) file 
    open(22,file='coll_c2_tot.dat')  !open C2(t) file 
    open(23,file='coll_msdo_tot.dat') !open MSD(t) file 
    open(24,file='coll_msd1_tot.dat') !open MSD(t) file 
    open(25,file='coll_msd2_tot.dat') !open MSD(t) file 
    open(31,file='coll_ener_c1_tot.dat')  !open dE-C1(t) file 
    open(32,file='coll_ener_c2_tot.dat')  !open dE-C2(t) file 
    open(33,file='coll_ener_msdo_tot.dat') !open dE-MSD(t) file 
    open(38,file='coll_ener_msd1_tot.dat') !open dE-MSD(t) file 
    open(39,file='coll_ener_msd2_tot.dat') !open dE-MSD(t) file 
    open(34,file='coll_ke_msd_tot.dat') !open dKE-MSD(t) file 
    open(35,file='coll_pot_msd_tot.dat') !open dV-MSD(t) file 
    open(36,file='coll_lj_msd_tot.dat') !open dLJ-MSD(t) file 
    open(37,file='coll_coul_msd_tot.dat') !open dCoul-MSD(t) file 
    open(41,file='coll_de_c1_tot.dat')  !open dE-C1(t) file 
    open(42,file='coll_de_c2_tot.dat')  !open dE-C2(t) file 
    open(43,file='coll_de_msdo_tot.dat') !open dE-MSD(t) file 
    open(44,file='coll_de_msd1_tot.dat') !open dE-MSD(t) file 
    open(45,file='coll_de_msd2_tot.dat') !open dE-MSD(t) file 
    open(46,file='coll_de2_msd_tot.dat') !open dE^2-MSD(t) file 
    open(50,file='bl_c2_tot.dat')  !open block C2(t) file 
    open(51,file='bl_de_c1_tot.dat')  !open block dE-C1(t) file 
    open(52,file='bl_de_c2_tot.dat')  !open block dE-C2(t) file 
    open(53,file='bl_de_msd_tot.dat') !open block dE-MSD(t) file 
    open(54,file='bl_ener_msdo_tot.dat') !open ener-MSDo(t) file 
    open(55,file='coll_ener2_msd_tot.dat') !open ener^2-MSD(t) file 

    open(61,file='coll_ke_c2_tot.dat')  !open dKE-C2(t) file 
    open(62,file='coll_pot_c2_tot.dat')  !open dV-C2(t) file 
    open(63,file='coll_lj_c2_tot.dat')  !open dLJ-C2(t) file 
    open(64,file='coll_coul_c2_tot.dat')  !open dCoul-C2(t) file 

    do it = 0, ntimes - 1
       write(21,'(2F12.5)') real(it)*dt, c1_tot(it)
       write(22,'(2F12.5)') real(it)*dt, c2_tot(it)
       write(23,'(5F12.5)') real(it)*dt, msd_tot(it,1), err_msd(it,1)
       write(24,'(5F12.5)') real(it)*dt, msd_tot(it,2), err_msd(it,2)
       write(25,'(5F12.5)') real(it)*dt, msd_tot(it,3), err_msd(it,3)

       write(31,'(2F12.5)') real(it)*dt, ec1_tot(it)/c1_tot(it)
       write(32,'(2F12.5)') real(it)*dt, ec2_tot(it)/c2_tot(it)

       write(41,'(2F12.5)') real(it)*dt, ec1_tot(it)
       write(42,'(5F12.5)') real(it)*dt, ec2_tot(it), err_ec2(it)
       write(43,'(5F12.5)') real(it)*dt, emsd_tot(it,1), err_emsd(it,1)
       write(44,'(5F12.5)') real(it)*dt, emsd_tot(it,2), err_emsd(it,2)
       write(45,'(5F12.5)') real(it)*dt, emsd_tot(it,3), err_emsd(it,3)
       write(46,'(5F12.5)') real(it)*dt, (e2msd_tot(it,k), k=1,3)

       write(61,'(5F12.5)') real(it)*dt, kec2_tot(it)
       write(62,'(5F12.5)') real(it)*dt, potc2_tot(it)
       write(63,'(5F12.5)') real(it)*dt, ljc2_tot(it)
       write(64,'(5F12.5)') real(it)*dt, coulc2_tot(it)

       write(50,'(21F12.5)') real(it)*dt, (c2_bl(it,ibl), ibl=1,nbl)
       write(51,'(21F12.5)') real(it)*dt, (ec1_bl(it,ibl), ibl=1,nbl)
       write(52,'(21F12.5)') real(it)*dt, (ec2_bl(it,ibl), ibl=1,nbl)
       write(53,'(21F12.5)') real(it)*dt, (emsd_bl(it,1,ibl), ibl=1,nbl)
    enddo

    do it = 1, ntimes - 1
       write(33,'(5F12.5)') real(it)*dt, emsd_tot(it,1)/msd_tot(it,1), err_enermsd(it,1)
       write(38,'(5F12.5)') real(it)*dt, emsd_tot(it,2)/msd_tot(it,2), err_enermsd(it,2)
       write(39,'(5F12.5)') real(it)*dt, emsd_tot(it,3)/msd_tot(it,3), err_enermsd(it,3)
       write(34,'(5F12.5)') real(it)*dt, (kemsd_tot(it,k)/msd_tot(it,k), k=1,3)
       write(35,'(5F12.5)') real(it)*dt, (potmsd_tot(it,k)/msd_tot(it,k), k=1,3)
       write(36,'(5F12.5)') real(it)*dt, (ljmsd_tot(it,k)/msd_tot(it,k), k=1,3)
       write(37,'(5F12.5)') real(it)*dt, (coulmsd_tot(it,k)/msd_tot(it,k), k=1,3)
       write(54,'(21F12.5)') real(it)*dt, (emsd_bl(it,1,ibl)/msd_bl(it,1,ibl), ibl=1,nbl)
!       write(55,'(5F12.5)') real(it)*dt, (e2msd_tot(it,k)/msd_tot(it,k) + & 
!            (emsd_tot(it,k)/msd_tot(it,k))**2, k=1,3)
       write(55,'(5F12.5)') real(it)*dt, -e2msd_tot(it,1)/msd_tot(it,1),(emsd_tot(it,1)/msd_tot(it,1))**2
    enddo


    close(21)
    close(22)
    close(23)
    close(24)
    close(25)

    close(31)
    close(32)
    close(33)
    close(34)
    close(35)
    close(36)
    close(37)
    close(38)
    close(39)
    close(41)
    close(42)
    close(43)
    close(51)
    close(52)
    close(53)
    close(54)
    close(55)
    close(61)
    close(62)
    close(63)
    close(64)


  End Program collect_msd_rot


!     Subroutine to Determine the value of t for 95% Confidence Limits                                                                                  

      Subroutine tvalue(nblocks,t)
      implicit none

      integer :: nblocks,nu
      real :: t

!     Calculate the Degrees-of-Freedom (average)

      nu = nblocks - 1

      if(nu.eq.1) then
         t = 12.7d0
      elseif(nu.eq.2) then
         t = 4.30d0
      elseif(nu.eq.3) then
         t = 3.18d0
      elseif(nu.eq.4) then
         t = 2.78d0
      elseif(nu.eq.5) then
         t = 2.57d0
      elseif(nu.eq.6) then
         t = 2.45d0
      elseif(nu.eq.7) then
         t = 2.36d0
      elseif(nu.eq.8) then
         t = 2.31d0
      elseif(nu.eq.9) then
         t = 2.26d0
      elseif(nu.eq.10) then
         t = 2.23d0
      else
         t = 2d0   ! this is approximate, but only off by 10% at most.                                                                                     
      endif

    End Subroutine tvalue
