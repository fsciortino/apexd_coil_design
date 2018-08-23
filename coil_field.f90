!    =================================================================
!    Subroutine to construct a Python module via f2py. 
!      Computation of magnetic field by ring and line current in APEX-D.
!      Consideres line loop current, coil volume integral, line (element)
!      current options. 
!
!    History:
!      - Computation originally in mag1.f90, written by H. Saitoh, 28.11.2017
!        Compiled with "gfortran -o mag1 mag1.f90"
!      - Modified to have subroutine apexd_field with f2py setup.
!        Compile into a Python module with
!            f2py -c -m apexd_field apexd_field.f90
!        (additional options "--fcompiler=gfortran --compiler=gcc" may be
!        needed on some systems)
!      
!    How it works:
!      According to set values (coil position, size, etc.), magnetic field
!      and flux function are calculated using complete elliptic function.
!      Can be used for calculations of B field inside SC coil.
!
!    How to use:
!      In "set parameters", decide coil parameters.
!      Only for the first run, set itf=1 so that elliptic function data
!      el2.dat are generated. As far as this file exits, you may set itf=0.
!      Then elliptic function data are loaded from this file, which reduces
!      calculation time.
!      Data files are saved as b.txt and bl.txt. See comments.

!    =================================================================
     subroutine apexd_field(ci1,rs1,rs2,zs1,zs2, xm,zm, Bxm,Bzm,Bs,rAm)
!    =================================================================
       implicit none
       double precision, intent(inout) :: ci1,rs1,rs2,zs1,zs2
!f2py intent(inout) :: ci1, rs1, rs2, zs1, zs2
       double precision, dimension(0:205), intent(out) :: xm, zm
!f2py intent(out) :: xm, zm
       double precision, dimension(0:205,0:205), intent(out) :: Bxm, Bzm, Bs, rAm
!f2py intent(out) :: Bxm,Bzm,Bs,rAm
       
       double precision, parameter :: pi = 3.14159265358979323846d0;
       integer :: nxmax,nzmax,nzmid,i,j,k,is,nx,nz,itf,itfo,ick
       double precision,dimension(0:100005) :: elik2,elie2
  
       double precision :: Bx,By,Bz,Bxt,Byt,Bzt,rAmt
       double precision :: x1,y1,z1,x2,y2,z2,ci,rc,zc,xc,dl,x,y,z,r
       double precision :: xmin,xmax,zmin,zmax,dx,dz,b,xi,at
       double precision :: ci11,rc11,zc11,xc11
       double precision :: rsd1,zsd1,zds1,rds1,csq1,xc1
       integer :: nrc1,nzc1,ic,jc

      
!      ----------------------------------
!      Obtain elliptic functions
       LOGICAL :: file_exists 
       INQUIRE( FILE="el2.dat", EXIST=file_exists ) 
       IF ( file_exists ) then
          OPEN (50,File="el2.dat")
          DO I=0,99999
               READ(50,*) xi,elik2(i),elie22(i)
          end do
          CLOSE(50)
       else
          call efile2(elik2,elie2) ! calculated K and E are saved as arrays
          OPEN (2,File="el2.dat")
          DO i=0,99999
             write(2,703) I/100000.0d0,elik2(i),elie2(i)
          end do
          CLOSE(2)
       endif

       if(elik(1)-1.57079d0.gt.0.01d0) then ! Warning when E() is empty.
          write(6,*) ""
          write(6,*) "##### Check eliptic functions! Data is null! #####"
          write(6,*) ""
       endif

!      -----------------------------------
!      Calculation region grid number and area
!      nx and nz must be even.
       nx=50; nz=50;
       xmin=-200.0d-3; xmax=200.0d-3;
       zmin=-200.0d-3; zmax=200.0d-3;

       nrc1 = 100; nzc1 = 100; ! calculation line numbers
       rsd1 = (rs2-rs1)/dble(nrc1); zsd1 = (zs2-zs1)/dble(nzc1); 
       csq1 = 1.0d0 / dble(nrc1+1) / dble(nzc1+1);
       xc1 = 0.0d0; ! offset in x

!      -----------------------------------
!      Calculation region, start
!      Set mesh size in r(x) and z, and make arrays for r and z
       nxmax=nx+1; nzmax=nz+1; nzmid=nz/2+1;

       dx=(xmax-xmin)/dble(nx)
       dz=(zmax-zmin)/dble(nz)

       do i=1,nxmax
          xm(i)=xmin+dble(i-1)*dx
       end do
       
       do j=1,nzmax
          zm(j)=zmin+dble(j-1)*dz
       end do

      
!      ##### B calculation, start ####
!      calculate B and 2 pi r Atheta at each point
       Bxt=0.0d0; Byt=0.0d0; Bzt=0.0d0; rAmt=0.0d0
      
!      calculate B and rAtheta at each (x,z)
       do i = 1,nxmax
          x = xm(i)
          do j = 1,nzmax
             z = zm(j)

!            ## superpose coil currents for B, start ##
             Bxt=0.0d0; Byt=0.0d0; Bzt=0.0d0;
               
             do ic=0,nzc1
                zds1 = zs1 + zsd1*dble(ic)
                do jc=0,nrc1
                   rds1 = rs1 + rsd1*dble(jc)
                   call bloop(elik2,elie2,rds1,dabs(x-xc1),z-zds1,ci1*csq1,Bx,Bz)
                   Bxt=Bxt+Bx; Bzt=Bzt+Bz;
                end do
             end do

             Bxm(i,j) = Bxt; Bzm(i,j) = Bzt; Bs(i,j) = dsqrt( Bxt*Bxt + Bzt*Bzt )

!            ## superpose coil currents for rAtheta, start ##
             rAmt=0.0d0

!            # coil volume #
             do ic=0,nzc1
                zds1 = zs1 + zsd1*dble(ic)
                do jc=0,nrc1
                   rds1 = rs1 + rsd1*dble(jc)
                   call rathe(elik2,elie2,rds1,dabs(x-xc1),z-zds1,ci1*csq1,at)
                   rAmt=rAmt+2.0d0*pi*at
                end do
             end do
      
             rAm(i,j) = rAmt

          end do
       end do

     703 format(1p3d23.15)

     end subroutine apexd_field




!     =================================================================
      subroutine bline(x1,y1,z1,x2,y2,z2,ci,x,y,z,dl,Bx,By,Bz)
!     current element at r1, directed in r2, makes B at r
!     =================================================================
        implicit none
        double precision, intent(in) :: x1,y1,z1,x2,y2,z2,ci,x,y,z,dl
        double precision, intent(out) :: Bx,By,Bz
        double precision :: Bx1,By1,Bz1,Bx2,By2,Bz2,r2,r1r2,r1r,rcos,rsin
        double precision :: xcc = 1.0d-7
      
        Bx1 = y2*(z-z1) - z2*(y-y1)
        By1 = z2*(x-x1) - x2*(z-z1)
        Bz1 = x2*(y-y1) - y2*(x-x1)
        
        Bx2 = Bx1 / dsqrt( Bx1*Bx1 + By1*By1 + Bz1*Bz1 )
        By2 = By1 / dsqrt( Bx1*Bx1 + By1*By1 + Bz1*Bz1 )
        Bz2 = Bz1 / dsqrt( Bx1*Bx1 + By1*By1 + Bz1*Bz1 )
      
        r2 = dsqrt( x2*x2 + y2*y2 + z2*z2 )
        r1r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)
        r1r = dsqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) )
        rcos = ( x2*(x-x1) + y2*(y-y1) + z2*(z-z1) ) / r2 / r1r
        rsin = dsqrt(1-rcos*rcos)
      
!       write(6,*) rcos,rsin

        Bx = Bx2 * xcc * dl * ci * rsin / r1r2
        By = By2 * xcc * dl * ci * rsin / r1r2
        Bz = Bz2 * xcc * dl * ci * rsin / r1r2

      end subroutine bline


!     =================================================================
      SUBROUTINE BLOOP(elik2,elie2,RC,R,Z,CI,BR,BZ)
!     calculate the Br and Bz produced by a loop current
!     RC:loop radius R,Z:position CI:coil current
!     =================================================================
!                                                   Biot-Savalt formula
!                                                   ===================
        implicit none
        double precision, dimension(0:100005) :: elik2,elie2
        DOUBLE precision :: ZK,ZZK,RC,R,Z,CI,BR,BZ,G1,FK,FE,A,G,FF,E,H,xcc,fk1,fe1,fk2,fe2
        integer :: i,i2

        xcc=2.0d-7

        ZK=4.0D0*RC*R/((RC+R)*(RC+R)+Z*Z)
        ZZK=dsqrt(ZK)

!       write(6,*) "k",ZZK
!       write(6,*) "K",elik2(100)
!       write(6,*) "E",elie2(ZK)

        IF(ZK.GE.0.999999D0) GO TO 10

        G1=dsqrt(1.0D0-ZK)
        IF(ZK.GT.0.9998D0) FK=DLOG(4.0D0/G1)+0.25D0*(DLOG(4.0D0/G1)-1.0D0)*G1*G1
        IF(ZK.GT.0.9998D0) FE=1.0D0+0.5D0*(DLOG(4.0D0/G1)-0.5D0)*G1*G1
        IF(ZK.GT.0.9998D0) GO TO 20

        I=IDINT(ZZK*100000.0D0)

!       Linear interpolation of FK and FE start
        i2=i+1;
        fk1=elik2(i); fe1=elie2(i);
        fk2=elik2(i2); fe2=elie2(i2);

        fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
        fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);

!       no interpolation of FK and FE
!       FK=elik2(I)
!       FE=elie2(I)

20      A=XCC*CI/dsqrt((RC+R)*(RC+R)+Z*Z)
        G=RC*RC+R*R+Z*Z
        FF=(RC-R)*(RC-R)+Z*Z
        E=A*Z/R
        H=RC*RC-R*R-Z*Z
        BZ=A*(H*FE/FF+FK)
        BR=E*(G*FE/FF-FK)

        IF(I.EQ.0) THEN
           BR=0.0d0
        ENDIF

        RETURN

10      BZ=0.0D0
        BR=0.0D0
        RETURN
      END SUBROUTINE BLOOP


!     =================================================================
      SUBROUTINE rathe(elik2,elie2,rc,x,z,ci,at)
!     calculate r A_theta produced by a loop current
!     RC:loop radius R,Z:position CI:coil current
!     Contour plot of this (r A_theta) gives magnetic surfaces
!     =================================================================
!                                                   Biot-Savalt formula
!                                                   ===================
        implicit none
        double precision, intent(in), dimension(0:100005) :: elik2,elie2
        double precision ZK,ZZK,A0,A1,A2,FK,FE,xcc,fk1,fe1,fk2,fe2
        double precision, intent(in) :: rc,x,z,ci
        double precision, intent(out) :: AT
        integer i,i2

        DATA XCC/2.d-7/

!       ZK is k^2
        ZK=4.0d0*rc*x/((rc+x)*(rc+x)+z*z)
        ZZK=dsqrt(ZK)
        
        IF(ZK.GE.0.999999d0) GO TO 20
        IF(ZK.GT.0.9998d0) GO TO 10
        
        A0=2.0d0*XCC*CI*x/dsqrt(ZK)
        A1=dsqrt(RC/x)
        A2=1.0d0-0.5d0*ZK
        I=IDINT(ZZK*100000.0d0)
        
!       Linear interpolation of FK and FE start
        i2=i+1;
        fk1=elik2(i); fe1=elie2(i);
        fk2=elik2(i2); fe2=elie2(i2);

        fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
        fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);

!       no interpolation of FK and F
!       FK=elik2(I)
!       FE=elie2(I)

        AT=A0*A1*(A2*FK-FE)

        IF(I.EQ.0) THEN
           AT=0.0d0
        ENDIF
        
        RETURN

10      A1=XCC*CI*RC
        A2=dlog(8.0d0*x/dsqrt((x-RC)*(x-RC)+Z*Z))
        AT=A1*(A2-2.0d0)
        RETURN

20      AT=0.0d0
        RETURN
      END SUBROUTINE rathe



!     =================================================================
      SUBROUTINE EFILE2(elik2,elie2)
!     complete elliptic fuctions
!     calculation speed may be adjusted by setting nmax x2
!     =================================================================
!                                                    elliptic functions
!                                                    ==================
        implicit none
        double precision,dimension(0:100005) :: elik2,elie2
        double precision ZZZAAA,FEF1,FEF2
        integer ji
        
        write(6,*) "start: elliptic fuction calculation"
        write(6,*) "... may takes several seconds ..."

        DO JI=0,99999
           ZZZAAA=DBLE(JI)/100000.0d0
           ELIK2(JI)=FEF1(ZZZAAA)
           ELIE2(JI)=FEF2(ZZZAAA)

!          just to show how calculation is going
           IF ( MOD(JI,20000).EQ.0 ) THEN
              WRITE(6,*) DBLE(JI)/100000.0*100.0,'% '
           ENDIF

        end do

        WRITE(6,*) DBLE(1)/1.0*100.0,'% '
        write(6,*) "end: elliptic fuction calculation"

        RETURN
        
      END SUBROUTINE EFILE2

!     ---------------------------------------
      double precision function FEF1(k)  ! Complete elliptic integral of the first kind 
        implicit none
        double precision, intent(in) :: k  !
        double precision :: pi, m, dt, t, tmin, tmax
        integer :: i
        integer, parameter :: nmax=2000
        double precision :: f, x
        
        f(m,x) = 1.0d0/dsqrt(1.0d0-(m*dsin(x))**2)
        
        if(k.ge.1.0d0)then
           write(*,*) "(error ! : k must 0=<k<1.)"
           return
        end if

        pi = 3.14159265358979323846d0
        
        tmin = 0.0d0
        tmax = pi/2.0d0
        dt = (tmax-tmin)/dble(nmax-1)
        
        FEF1 = 0.5d0*dt*(f(k,tmin)+f(k,tmax))
        do i=1,nmax-2
           t = tmin+dt*dble(i)
           FEF1 = FEF1+dt*f(k,t)
        end do
        
        return
      end function FEF1

!     -------------------------------------
      double precision function FEF2(k)  ! Complete elliptic integral of the second kind 
        implicit none
        double precision, intent(in) :: k  !
        double precision :: pi, m, dt, t, tmin, tmax
        integer :: i
        integer, parameter :: nmax=2000
        double precision :: f, x
        
        f(m,x) = dsqrt(1.0d0-(m*dsin(x))**2)
        
        pi = 3.14159265358979323846d0
        
        if(k.gt.1.0d0)then
           write(*,*) "(error) ! : k must 0=<k=<1."
           return
        end if
        
        tmin = 0.0d0
        tmax = pi/2.0d0
        dt = (tmax-tmin)/dble(nmax-1)
        
        FEF2 = 0.5d0*dt*(f(k,tmin)+f(k,tmax))
        do i=1,nmax-2
           t = tmin+dt*dble(i)
           FEF2 = FEF2+dt*f(k,t)
        end do
        
        return
      end function FEF2
      
      



