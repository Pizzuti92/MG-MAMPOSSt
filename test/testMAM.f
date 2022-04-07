       program testMAM
!>  test some routines in libMAM.a
!>  Usage (assuming testMAM.f is in the main folder): 
!>  f95 -o testMAM.e testMAM.f -L build/ -lMAM -L build/GamI/ -lGamI -L build/Newuoa/ -lNewuoa -L build/JJin/ -lJJin -L build/Utili/ -lUtili -L build/Powell/ -lPowell
       
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       

       include 'paramsoptS.i' 
       
       h0=70
       Omegam=0.3
       Omegal=0.7
       kmp=9
       r200=1.2
       rs=0.5
       tmass=100.23
       nhone=10
       screen=0.3
       write(*,*) 'test dphidr= ', dphidr(2.12d0)
       
       
c************ test field ***********************************************
       h0=70
       Omegam=0.3
       Omegal=0.7
       za=0.3
       r200=1.2
       rs=0.5
       tmass=0.8
       write(*,*) 'test field= ', field(2.12d0)       
       
C********* test sr2int *************************************************

       nbs=1
 !  Values of the parameters of Hansen & Moore's (2006)
 !     beta = a + b dln(rho)/dln(r) relation
       ahm=-0.15
       bhm=-0.192
 
       r200=1.41
       rc=0.3
       rs=0.3
       cbe=1.2
       tmass=0.10
       screen=0.4
       kmp=7
       v200=r200*10*85
       knfit=1
       kani=4
       kscr=1
       rcut=1.0
       zb=0.3
       nhs=2

       alr=dlog(2.5d0) 
       sint=sr2int(alr)  
       write(*,*) 'test sr2int= ', sint
       
 
C********* test sr2out ************************************************* 
       nbs=1
 !  Values of the parameters of Hansen & Moore's (2006)
 !     beta = a + b dln(rho)/dln(r) relation
       ahm=-0.15
       bhm=-0.192
 
       r200=1.41
       rc=0.3
       rs=0.3
       cbe=1.2
       tmass=1.0
       screen=0.4
       kmp=7
       knfit=1
       kani=4
       kscr=1
       rcut=1.0
       zb=0.3
       nhs=2


       tt=0.5d0 
       sint=sr2out(tt)  
       write(*,*) 'test sr2out= ', sint
C******** test gt_mod ************************************************** 
 
        r2=1.41
        rss=0.3
        ya=0.1
        yb=0.1
        write(*,*) 'test gt_mod= ', gt_mod(1.2d0,r2,rss,ya,yb) 
 
C******** test gt_true ************************************************* 
 
        r200t=1.41
        rst=0.3

        write(*,*) 'test gt_true= ', gt_true(1.2d0)
       
C******* test Likelens_bh **********************************************      
      r200=1.41
      rs=0.3 
      cbe=1.2
      tmass=1.2
      screen=0.4
!      Olam=0.7
!      Omegam=0.3
!      h0=70
      rupin=r200
      rst=0.3    !lensing "true" value of the scale radius
      r200t=1.41 !lensing "true value of the virial raidius
      wr2=0.25  !intrinsic ellipticity
      wr=0.005  !lss noise
      wx=30     !number of galaxies per arcmin^2
      call  Likelens_bh(10,plens)
      write(*,*) 'test Likelens_bh= ', plens
       end
