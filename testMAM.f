       program testMAM
!>  test some routines in libMAM.a
!>  Usage (assuming testMAM.f is in the main folder): 
!>  f95 -o testMAM.e testMAM.f -L build/ -lMAM -L build/GamI/ -lGamI -L build/Newuoa/ -lNewuoa -L build/JJin/ -lJJin -L build/Utili/ -lUtili -L build/Powell/ -lPowell
       
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       
       include 'datavdp.i'
       include 'paramsoptS.i' 
       include 'barions.i'
       
       parameter (iii=5000, grav=4.302e-9, pig=3.141592653589793d0)
       dimension rvdp(iii),vdp(iii),sig2i(iii)
       dimension bcgr(iii),bcgvdp(iii)
       character(len = 16) :: filebcg
       character(len = 41) :: fname
       parameter  (nout=200)
       

       open(12,file='data/massbar.dat',status='old')
       i=0
 723   continue
       read(12,*,end=929) rgasr,xmgasr
       i=i+1
       rgas(i)=rgasr
       xmgas(i)=xmgasr*1.e+13
       
       goto 723
 929   continue
       close(12)
       ngas=i
       
       rjaf=0.03d0
       xlumbcg=492000000000.0d0
       xmstar=4.3

       !Cosmology
       ncutof = 1
     
       
c       h0=70
c       Omegam=0.3
c       Olam=0.7
c       za=0.44
c       !RG
       
c       r200 = 1.5
c       eim = 1.4
c       ncutof = 0
c       kmp = 13
c       nhone = 2
       
c       rhos = 1e15
c       rss = 0.5
c       rs = 0.5
c       eim = 1.0
c       Q = 1
c       screen = Q
c       tmass=10.0d0
c       phinf = 4e-5
c       b = 2.3
      

! ======================================================================      
       write(*,*) "Test of general screening in Chameleon gravity"
       
       !all parameters needed
       rs = 0.87
       r200 = 2.0
       tmass = 10.0d0
       screen = 1.0d0
       dA3 = 0.73
       xlumbcg = 4.9e11
       xmstar = 4.46
       rjaf = 0.039
       rtbcg=20.
       rabcg=20.2
       
       rc = 0.46
       h0=70
       Omegam=0.3
       Olam=0.7
       za=0.44
       rhog1 = 25.7e13
       rsga1  = 0.369
       rhost1 = 3.547e13
       rst1 = 0.369
       cbe = 3.16
       cbe0 = 1.41
       xlim = 1.0e-6
       
       
       xsc = ScreenGen(xlim)
       write(*,429)  xsc
 429  format("Value of the screeining:",1(f6.3))       
       xscCH = xsc
       write(*,*) "Value of the effective mass:"
       write(*,430) dphiGen(2.0d0)
 430  format(1(e14.4))
 
       write(*,*) "Value of the total mass:"
       write(*,430)  tmass_Gen(1.2d0)
       
       kmp = 17
       
       call vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)
       ifin = 6
       write(*,*) iout
       rbcg(1) = 1.000e-03
       rbcg(2) = 2.154e-03
       rbcg(3) = 4.642e-03
       rbcg(4) = 1.000e-02
       rbcg(5) = 2.154e-02
       rbcg(6) = 4.642e-02
      
      rlim = 2.*rbcg(6) 
      do k=1,iout
        if (bcgr(k).le.rlim) imax=k
      enddo
      write(*,*) imax
       do i=1,ifin
        rbcgj=rbcg(i)             
        call hunt(rvdp,imax,rbcgj,jlo)
        write(*,*) rbcgj, vdp(jlo)
       enddo
       
c       stop('DONE')
      
c ======================================================================
        
C       nstep = 500
c       !fac200 = 1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
C       fname = '/Users/lorenzopizzuti/Desktop/testbar.dat'
c       open(20,file=fname,status='unknown')
c       do i=1,nstep
c        x=0.01*i
        
        !xm = (dlog(1.d0+tt)-tt/(1.d0+tt))*fac200
c        write(20,*) x, tmass_Gen(x) 
c        write(*,*) x, tmass_Gen(x)  !, r200, rs, tmass, screen
c       enddo
c       close(20)
       

       
c       stop('fa che funzioni')

       
       
c       kbcg=1
       
c       if (kbcg.eq.1) then
c            !read data file of BCG
c       filebcg='data/vel_bcg.dat'
c       open (12, FILE=filebcg, STATUS='OLD')
c        jj=0
c        j0=-1
c        dimin=1.e12
c 406    continue
                
c        read(12,*,end=176) dmpc,dvr,dev 
c        jj=jj+1
c        rbcg(jj)=dmpc
c        sbcg(jj)=dvr
c        esbcg(jj)=dev
cc       write(*,*) rbcg(jj), sbcg(jj)
c        goto 406
c 176    continue
c        close(12)
c        nbcg=jj
c       endif
       
c      open(12,file='data/massbar.dat',status='old')
c      i=0
c 777  continue
c      read(12,*,end=999) rgasr,xmgasr
c      i=i+1
c      rgas(i)=rgasr
c      xmgas(i)=xmgasr*1.e+13

c      goto 777
c 999  continue
c      close(12)
c      ngas=i
cc      do i=1,1000
cc       x=0.02*i
cc       write(*,*) x, fbars(x)
cc      enddo
      
c       h0=70
c       Omegam=0.3
c       Omegal=0.7
c       kmp=9
c       r200=1.2
c       rs=0.5
c       tmass=100.23
c       nhone=10
c       screen=0.3
c       write(*,*) 'test dphidr= ', dphidr(2.12d0)
       
       
cc************ test field ***********************************************
c       h0=70
c       Omegam=0.3
c       Omegal=0.7
c       za=0.3
c       r200=1.2
c       rs=0.5
c       tmass=0.8
c       write(*,*) 'test field= ', field(2.12d0)       


cC****** test frlin *****************************************************
c       kscr=-1  !screening parameter. kscr=-1 corresponds to linear Horndeski gravity
c       screen=1.2
c       r200=1.2
c       rs=0.5
c       tmass=0.8
c       write(*,*) 'test frlin= ', frlin(2.12d0)      
cC********* test sr2int *************************************************

       nbs=1
 !     Values of the parameters of Hansen & Moore's (2006)
 !     beta = a + b dln(rho)/dln(r) relation
       ahm=-0.15
       bhm=-0.192
       h0=70
       Omegam=0.3
       Olam=0.7
       r200=2.0
       rc=0.75
       rs=0.3
       cbe=2.7
       cbe0=1.0 !not used unless kani=21, 41
       tmass=0.10 
       screen=2.
       kmp=1    !mass profile model. kmp=7 corresponds to mNFW_LH
       v200=r200*10*70
       knfit=1  !Number density profile model. knfit=1 Corresponds to pNFW
       kani=8   !Anisotropy model. kani=4 corresponds to BP profile
       kscr=1   !screening option
       rcut=1.0
       zb=0.3
       za=0.3
       nhs=2
       eim=5.   
       kbcg=1   
       rcut = 0.7
       xmstar=4.5d0
       xlumbcg=492000000000.0d0
       rjaf=0.039
       rtbcg=20.
       rabcg=20.2
       alr=dlog(2.5d0) 
       sint=sr2int(alr)  
       
       rc=0.039d0
       kdbcg=0
       ncutof=0
!************************** TEST BURKERT MG ****************************
       alimi=0.001
       xsccH=alimi
       alimit=alimi
       tmass=10.
       r200=2.2d0
       rs=0.3d0
       call find_xsc(xscCH,alimit)       
       write(*,*) xsccH
    
       
       Hz2=h0**2*(Omegam*(1+za)**3+Olam)
       phinf=tmass*1e-5
       
       dm200=100*Hz2/grav*r200**3
       
       t=2.0d0
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm1=fac200*gip
         
         write(*,12) (xm1+dphiE(t/rs))*dm200 
   12  format ('test dphib = ', E13.4 )
       
       write(*,12) dm200
       

! *********************************************************************       
       
       
      write(*,*) 'inclusion of gas and galaxies from a separate file'
      open(12,file='data/massbar.dat',status='old')
      i=0
 717  continue
      read(12,*,end=999) rgasr,xmgasr
      i=i+1
      rgas(i)=rgasr
      xmgas(i)=xmgasr*1.e+13

      goto 717
 999  continue
      close(12)
      ngas=i
       
       write(*,*) 'test RF', fmodRG(0.5d0)
       
       
       write(*,*) 'test frhoproj', frhoproj(0.4d0)
       write(*,*) 'test frhoprojH', frpHer(0.4d0)
       write(*,*) 'test fjH',sigmar2(0.4d0)
       
       write(*,*) 'test sr2int= ', sint

       write(*,*) 'test fmass= ', fmass(2.5d0) 

       alimit=alimi
       kmp=17
       xlim = 1.0e-3
       xscCH = ScreenGen(xlim)       
       
       
       
       !write(*,*) 'test dphiE=', dphiE(3.0d0)

c       write(*,*) 'test FINTEGRAND= ', fintegrand(alr)
   
c       call vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)
c       open(20,file='data/vdp_theo.txt',status='unknown')
c       do i=1,nout
c        write(20,*) rvdp(i), vdp(i)
c       enddo
c       close(20)
c       call vdplik(chibcg)
       
c       write(*,*) 'chisqure is ', chibcg
       
cC********* test sr2out ************************************************* 
c       nbs=1
c !     Values of the parameters of Hansen & Moore's (2006)
c !     beta = a + b dln(rho)/dln(r) relation
c       ahm=-0.15
c       bhm=-0.192
 
c       r200=1.41
c       rc=0.3
c       rs=0.3
c       cbe=1.2
c       cbe0=1.0 !not used unless kani=21, 41   
c       knfit=1  !number density profile model. knfit=2 corresponds to pHer
c       kani=4   !velocity anisotropy model. kani=4 corresponds to Tiret profile

c       tt=0.5d0 
c       sint=sr2out(tt)  
c       write(*,*) 'test sr2out= ', sint
cC******* test fmg,gmg, f_nfw, g_nfw ************************************
      
c      y1=0.5
c      y2=0.2
c      tt=1.2d0
c      write(*,*) 'test f_nfw= ', f_nfw(tt)
c      write(*,*) 'test g_nfw= ', g_nfw(tt)
c      write(*,*) 'test fmg= ', fmg(tt,y1,y2)  
c      write(*,*) 'test gmg= ', gmg(tt,y1,y2) 


cC******** test gt_mod ************************************************** 
 
c        r2=1.41
c        rss=0.3
c        ya=0.1
c        yb=0.1
c        write(*,*) 'test gt_mod= ', gt_mod(1.2d0,r2,rss,ya,yb) 
 
cC******** test gt_true ************************************************* 
 
c        r200t=1.41
c        rst=0.3

c        write(*,*) 'test gt_true= ', gt_true(1.2d0)
       
cC******* test Likelens_bh **********************************************      
c      r200=1.41
c      rs=0.3 
c      cbe=1.2
c      tmass=1.2
c      screen=0.4
c!      Olam=0.7
c!      Omegam=0.3
c!      h0=70
c      rupin=r200
c      rst=0.3    !lensing "true" value of the scale radius
c      r200t=1.41 !lensing "true value of the virial raidius
c      wr2=0.25  !intrinsic ellipticity
c      wr=0.005  !lss noise
c      wx=30     !number of galaxies per arcmin^2
c      call  Likelens_bh(10,plens)
c      write(*,*) 'test Likelens_bh= ', plens
       end
