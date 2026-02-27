      program gomamposstoptombcgrvsph

c     derived from gomamposst.f ...'opt' for optimization
c     derived from gomamposstopt.f ...'om' for omegas

c     The pgm assumes interlopers have been removed
c     unless one requires a joint halo+interlopers MAMPOSSt
c     analysis

c     Distances must be given in kpc from cluster center

c     Velocities must be given in km/s, rest-frame
c     (hence after 1/(1+<z>) correction, with <v>=0)

c     History of modifications:

c     Gap & Interloper removal part removed

c     H0 is an input parameter (in place of the interloper rejection parameter)

c     PIEMD M(r) inserted

c     Modified Tiret anisotropy model inserted (it has non-zero anisotropy 
c     at r=0)

c     Minimization program inserted

c     beta-model with fixed (free choice) alpha for N(R) inserted

c     Burkert and Softened IS M(r) inserted

c     Hansen+Moore beta(r)=a+b*dln(rho)/dln(r) inserted

c     Handling optimization with only 1 free parameter

c     Inclusion of velocity errors in the analysis

c     Corrected error in normalisation of Hernquist M(r)

c     Removed file of 'true' values from the input files

c     Added a Ngals-dependent scaling for the grid width
c     so that samples with fewer galaxies use wider grids

c     Added a fixed-m Einasto model

c     Forcing beta<1 

c     Pushing the determination of los VDP to lower distances
c     from the center (5 kpc)

c     weights added in input file

c     omega0, omegal now free params read in the params file
c     and passed through an include among different subroutines

c     modify mass profiles by including the baryonic mass profile
c     of the BCG, with a single free parameter, M*/L (with B. Sartoris)

c     modify the cluster mass profile by adding the inner slope free 
c     to the NFW model (i.e. implementing the gNFW model) (with B. Sartoris)

c     add an anisotropy model, Tiret_TAL, where r_beta=r_tracer

c     consider case of only real members by allowing integration to r200

c     add two anisotropy models, Tiret with beta0 <> 0      

      
c     Last modif: Paris, June 2012
c                 Trieste, September 2012
c                 Trieste, October 2012
c                 Milano, November 2012
c                 Trieste, December 2012
c                 Trieste, January 2013
c                 Paris, February 2013
c                 Trieste, July 2013
c                 Munchen, December 2014
c                 Trieste, January 2015
c                 Trieste, December 2015
c                 Trieste, May-June-July 2016 (with B. Sartoris)
c                 Paris, September 2018
c                 Trieste, Feb/Mar 2021      

c     if needed, re-create the archives:

c     cd ~/FORT/Utili
c     ifort -c *.f                       [f95 instead of ifort on laptop]
c     ar -r libUtil.a *.o

c     cd ~/FORT/GamI
c     ifort -c *.f
c     ar -r libGAM.a *.o

c     cd ~/FORT/Newuoa
c     ifort -c *.f
c     ar -r libNewuoa.a *.o

c     cd ~/FORT/Powell
c     ifort -c *.f
c     ar -r libPowell.a *.o

c     cd ~/FORT/JJin
c     ifort -c *.for
c     ifort -c *.f
c     ar -r libJin.a *.o

c     cd ~/FORT/Hyper
c     ifort -c *.f
c     ar -r libHyper.a *.o

c     compile with (one-line command; use gfortran if ifort not available):
c
c     cd ~/MAMPOSST/Code_AB
c     ifort -o gomamposstoptombcgrvsph.e gomamposstoptombcgrvsph.f 
c     -L/home/biviano/FORT/GamI/ -lGAM 
c     -L/home/biviano/FORT/Utili/ -lUtil 
c     -L/home/biviano/FORT/Newuoa -lNewuoa 
c     -L/home/biviano/FORT/Powell -lPowell
c     -L/home/biviano/FORT/Hyper -lHyper
c     -L/home/biviano/FORT/JJin -lJin

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,clight=2.99792458d5)
      dimension di(iga),ve(iga),eve(iga),dns(iga),ra(iga),de(iga),
     &     rso(iga),vso(iga),bins(400),pars(28),
     &     syb(1000),sybs(1000),vkb(1000),vkbs(1000),
     &     xx(iga),yy(iga),yboo(iga),wboo(iga),ww(iga),
     &     iw(iga)
      character*75 fgroup,frnvn,fbfn,fbn,ffn,fmlv,fprop,line,
     &     fsb,fsf,fsft,fkb,fkf,fkft,fparmam
      include 'datarv.i'
      include 'paramsoptbcg.i'
      include 'units.i'
      include 'sr.i'
      include 'vkr.i'
      include 'vlos.i'
      include 'tsallis.i'
      include 'omegas.i'
      external fcn1,fcn2,fcn3,fa,fk
      external sr2int,sr2out,sigmar1,sigmar2,sigmar3,gwenu
      external vr4nuint
      external sigmarnorm
      external gammarec

      call uerset(1,levold)
      icsr = 100
      icsrk = 100

c     number of bootstrap resamplings for the errors on sigma_los

      nboo=500

c     User's inputs

      qtsa=1.     ! forcing Gaussian MAMPOSSt

c     Tsallis constants and exponent

      if (abs(qtsa-1.).gt.0.01) then
         ktsa=1
         atsa=qtsa/(1.-qtsa)
         batsa=1.d0/(2.5d0+atsa)
         if (atsa.gt.0) then
            atatsa=gammarec(2.5d0+atsa)/gammarec(1.d0+atsa)*batsa**1.5
         elseif (atsa.lt.-2.5) then
            atatsa=gammarec(-atsa)/gammarec(-1.5-atsa)*(-batsa)**1.5
         endif
      else
         ktsa=0
      endif
c
c     Values of the parameters of Hansen & Moore's (2006)
c     beta = a + b dln(rho)/dln(r) relation
c
      ahm=-0.15
      bhm=-0.192

c
      
      write(*,*) ' '
      write(*,*) ' '
      write(*,*)    ' ***************************************'
      write(*,*)    ' ***                                 ***'
      if (ktsa.gt.0.5) then
         write(*,*) ' *** Welcome into Tsallis MAMPOSSt!  ***'
      else
         write(*,*) ' *** Welcome into Gaussian MAMPOSSt! ***'
      endif
      write(*,*)    ' ***                                 ***'
      write(*,*)    ' ***************************************'
      write(*,*) ' '
      write(*,*) ' '

      write(*,100)
 100  format(' Input files: '/ 
     &     ' (1) data (R [kpc], Vrf [km/s], errVrf [km/s], weights',/
     &     ' (2) input parameters for MAMPOSSt ')
      read(*,9) fgroup
      read(*,9) fparmam
 9    format(a)
      open(10,file=fgroup,status='old')
      open(29,file=fparmam,status='old')

      write(*,196) fgroup
 196  format(//' Data-set is ',a50,//)

c     read file with parameters for MAMPOSSt

      do i=1,28
         read(29,*) pars(i)
      enddo
      close(29)

      nr200=nint(pars(1))   ! number of steps for r200 fit
      nrc=nint(pars(2))     ! number of steps for rc fit, scale radius of N(R)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces LfM, c_nu=c]
                            !   [if = -2 fit N(R) outside MAMPOSSt]
      nrs=nint(pars(3))     ! number of steps for rs fit, scale radius of M(r)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces MfL, c=c_nu]
                            !   [if = -2 forces LCDM, c=c(M)]
      nbs=nint(pars(4))     ! number of steps for anisotropy parameter
                            !   [if = -1 forces a_ML=r_s]
                            !   [if = -2 forces Hansen+Moore]
      ngam=pars(5)          ! number of steps for inner slope of gNFW, gamma
      nml=pars(6)          ! number of steps for M*/L

      h0=pars(7)          ! Hubble constant at z=0
      za=pars(8)     ! average redshift of the cluster (needed to evaluate Hz)
                                        ! since velocities are given as rest-frame in input file
      rlowin=pars(9)        ! Inner radius for sample selection (Mpc) 
      rupin=pars(10)         ! Outer radius for sample selection (Mpc) 
      kintd=nint(pars(11))   ! 0: Interlopers removed; 1: remove interlopers statistically, 2: only real true members are present 
      r200g=pars(12)        ! r200 initial guess (Mpc)
      r200=r200g             ! set r200 initial value to guess value
      knfit=nint(pars(13))   ! N(R) model, projected NFW / projected Hernquist / beta-model (1/2/3)
                            !             
      rcg=pars(14)            ! N(R) scale radius initial guess (Mpc)
      al=pars(15)            ! N(R) negative exponent (only if knfit=3)
      kmp=nint(pars(16))     ! rho(r) model: gNFW/Hernquist/PIEMD/Burkert/SoftIS/Einasto_m=5 (1/2/3/4/5/6)
      rsg=pars(17)           ! rho(r) scale radius initial guess (Mpc)
      gam=pars(18)           ! inner slope of gNFW M(r) >0
      kani=nint(pars(19))    ! Anisotropy model, beta'=constant, MamLok, OsiMer, 
                             !   simplified Wojtak, simplified Tiret, modified Tiret, Tiret_TAL, Tiret_tan, Tiret_rad (0,1,2,3,4,5,6,7,8)
      if (nbs.eq.-1.) kani=1 !   forced to MamLok if requested
      if (nbs.eq.-2.) kani=-1 !   if Hansen&Moore, beta(r) depends on rho(r)

      cbeg=pars(20)         ! Anisotropy initial guess, beta', a_ML, a_OM, a_W, beta'_inf
      rcut=pars(21)  ! PIEMD model rcut in Mpc (only used if kmp=3)
      
      kbsp=nint(pars(22))   ! run MAMPOSSt in fast mode? N/Y=0/1

      kopt=nint(pars(23))   ! optimization algorithm: 0/1/2=bobyqa/newuao/powell
                                            ! -1 skip optimization
      omega0=pars(24)       ! omega matter
      omegal=pars(25)       ! omega lambda

      rjaf=pars(26)         ! scale radius of BCG Jaffe lum dens. prof (Mpc)
      xlumbcg=pars(27)      ! BCG total luminosity
      xmstarlg=pars(28)     ! free parameter: mass-to-light ratio of the
                            ! stellar population of te BCG (1st guess)


cccc
cccc  force kintd=2, no interlopers, integration on the r200 sphere
cccc
cccc      kintd=2

cccc
cccc  force anisotropy Tiret with beta0=-1
cccc
cccc      kani=7

cccc
cccc  force anisotropy Tiret with beta0=0.4
cccc
cccc      kani=8

      
      write(*,346) nr200,nrc,nrs,nbs,ngam,nml
 346  format(' Grid steps in r200, rtr, rs, anis, gamma, M*/L: ',6(i4))

c
c     Stop if H&M beta(r) required and M(r) difft from NFW
c
cc      if (kani.lt.0.and.kmp.ne.1) then
      if (kani.lt.0.and.knfit.ne.1) then
         write(*,*) ' Hansen & Moore beta(r) currently '
cc         write(*,*) ' implemented only for NFW M(r) model '
         write(*,*) ' implemented only for pNFW N(R) model '
         stop
      endif
         
      write(*,800)
 800  format(/' Output files: '/
     &     ' (1) normalized distances and velocities '/
     &     ' (2) best-fit result parameters for N(R) '/
     &     ' (3) binned N(R) '/
     &     ' (4) fitted N(R) '/
     &     ' (5) fit Max Lik values for velocity distribution '/
     &     ' (6) binned sigma_los(R)'/
     &     ' (7) fitted sigma_los(R) for MAMPOSSt solution')
      read(*,9) frnvn
      read(*,9) fbfn
      read(*,9) fbn
      read(*,9) ffn
      read(*,9) fmlv
      read(*,9) fsb
      read(*,9) fsf
      open(20,file=frnvn,status='unknown')
      open(30,file=fbfn,status='unknown')
      open(40,file=fbn,status='unknown')
      open(50,file=ffn,status='unknown')
      open(60,file=fmlv,status='unknown')
      open(70,file=fsb,status='unknown')
      open(80,file=fsf,status='unknown')

c     read system properties



c     read radial positions, velocities and velocity errors;
c     positions are in kpc, vels are in km/s, assumed rest-frame
c     errors are in km/s

      read(10,9) line
      read(10,9) line
      j=0
      j0=-1
      dimin=1.e12
 222  continue
      read(10,*,end=111) dkpc,vkms,evkms,wei
      j=j+1
      di(j)=dkpc/1.e3
      ve(j)=vkms
      eve(j)=evkms
cc      w(j)=1.0
      w(j)=wei
      if (di(j).lt.dimin) then
         dimin=di(j)
         j0=j
      endif
      iw(j)=j
      goto 222
 111  continue
      npg=j
      close(10)
      write(*,*) npg,' galaxies in the sample'

c     get an estimate of the velocity dispersion

      ibwt=1
      call robusti(ve,npg,ibwt,va,sv)

c     First guess at m200, v200; sigma_v is determined from r200 
c     using Mauduit & Mamon's scaling (a poly approx to it)
c     and the concentration appropriate for the m200 found
c     adopting a scaling intermediate between Duffy et al.'s
c     and (Dolag et al.'s) Gao et al.'s [at z=0]

c     These cosmo params are now in the params file
cc      omegal=0.7
cc      omega0=0.3
      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
      rv=r200
      rm200=100.*hz*hz/grav*r200**3
      cduffy=5.78*(rm200/2.e12)**(-0.089)*1.1**(-0.52)
      cgao=10.**(-0.138*dlog10(rm200*0.7)+2.646)
      cmean=(cduffy+cgao)/2.
      rsmean=rv/cmean
      v200=10.*hz*r200
      sv=(429.-0.6*cmean+0.4*cmean*cmean)*r200

      write(*,182) za,rv,v200,cmean 
 182  format(/,' Initial estimates of <z>, r200, v200 and c are: ',
     ,     f7.4,2x,f5.2,2x,f5.0,2x,f5.1,//)
     
c     call the mamposst procedure subroutine

      iu20=20
      iu30=30
      iu60=60

c     MAMPOSSt subroutine

      call mamposst(di,ve,eve,rso,vso,npg)

      close(30)
      close(60)

      if (r200.ne.r200) r200=r200g

      write(*,*) ' '
      write(*,*) ' After MAMPOSSt: '
      write(*,*) ' '
      write(*,*) '    build output files of binned N(R), VDP'
      write(*,*) '    of best-fit N(R), VDP and of'
      write(*,*) '    input N(R), VDP solutions (''true'' values)'
      write(*,*) ' '

c     output a binned number density profile N(R)

      ibwt=0
      nbins=int(dsqrt(dfloat(nga)))
      do j=1,nga
         dns(j)=rso(j)
      enddo
      call sort(nga,dns)
      do j=1,nbins+1
         if (j.eq.nbins+1) then
            do l=1,nga-nbins*nbins
               xx(l)=dns((j-1)*nbins+l)
            enddo
         else
            do l=1,nbins
               xx(l)=dns((j-1)*nbins+l)
            enddo
         endif
         lx=l-1
         call robusti(xx,lx,ibwt,cx,zero)
         area=pig*(xx(lx)*xx(lx)-xx(1)*xx(1))
         bj=cx
         dj=dfloat(lx)/area
         ej=dsqrt(dfloat(lx))/area
         write(40,*) bj,dj,ej
      enddo
      close(40)

      write(*,*) ' Binned N(R) computed'

c     output the result of the best-fit for N(R):
c     knfit=1,2,3 selects NFW, Hernquist, beta-model, resp.

      call sigmarnorm(rc,fnorm)

      if (knfit.eq.3) then

c     beta model

         do k=1,3000
            xk=(k-1)*0.002+0.001
            fr=(1.+(xk/rc)**2)**al
            write(50,*) xk,fr/fnorm*nga
         enddo

      elseif (knfit.eq.1.or.(nrc.eq.-1.and.kmp.eq.1)) then

c     NFW

         c=r200/rc
         gc=1./(dlog(c+1.)-c/(c+1.))
cc         facnorm=gc/(2.*pig)       missing the term r200^2 ? Trieste, 21/12/12
         facnorm=gc/(2.*pig)/(r200*r200)
         
         do k=1,3000
            xk=(k-1)*0.002+0.001
            cx=c*xk
            uu=1./cx
            if (uu.lt.1) then
               cm1=dacos(uu)
               fr=c*c*
     &              (1.-1./dsqrt(dabs(cx*cx-1.))*cm1)/(cx*cx-1.)
            elseif (uu.gt.1) then
               cm1=dacosh(uu)
               fr=c*c*
     &              (1.-1./dsqrt(dabs(cx*cx-1.))*cm1)/(cx*cx-1.)
            else
               fr=c*c/3.
            endif
            fr=fr*facnorm
            write(50,*) xk,fr/fnorm*nga
         enddo
         
      else

c     Hernquist

         do k=1,3000
            xk=(k-1)*0.002+0.001
            s=xk/rc
            if (s.gt.1.d0) then
               xs=dacos(1./s)/dsqrt(s*s-1.)
               fr=((2.+s*s)*xs-3.)/(2.*pig*rc*rc*(1.-s*s)**2)
            elseif (s.lt.1.d0) then
               xs=dlog((1.+dsqrt(1.-s*s))/s)/dsqrt(1.-s*s)
               fr=((2.+s*s)*xs-3.)/(2.*pig*rc*rc*(1.-s*s)**2)
            endif
            if (abs(s-1.d0).lt.0.001d0) then
               fr=2./(15.*pig*rc*rc)
            endif

            write(50,*) xk,fr/fnorm*nga
         enddo

      endif
      close(50)

ccc      write(*,*) ' Best-fit N(R) computed'

c     Here computes sigma_los 

      call sort2(nga,rso,vso)
      nbins=int(sqrt(float(nga)))/2.
      npbin=nga/nbins
      do j=1,nbins
         bins(j)=rso((j-1)*npbin)
         bins(j+1)=rso(j*npbin)
      enddo

ccc      write(*,*) ' Using ',nbins,' bins for the VDP(R) '

      do j=1,nbins
         npbin=0
         do ii=1,nga
            if (rso(ii).ge.bins(j).and.rso(ii).lt.bins(j+1)) then
               npbin=npbin+1
               xx(npbin)=rso(ii)
               yy(npbin)=vso(ii)
               ww(npbin)=1.d0
            endif
         enddo
         if (npbin.gt.1) then
            ibwt=0
            call robusti(xx,npbin,ibwt,cxd,sxd)
            ibwt=1
            call robusti(yy,npbin,ibwt,cyd,syd)

c     Bootstrap errors

            do iboo=1,nboo
               call boot(yy,ww,yboo,wboo,npbin)
               call robusti(yboo,npbin,1,cy,sy)
               syb(iboo)=sy
            enddo
            call shell(nboo,syb)
            j16=nint((nboo-1)*0.16+1)
            j84=nint((nboo-1)*0.84+1)
            esy=0.5d0*(syb(j84)-syb(j16))
            write(70,713) npbin,cxd,syd,-esy,esy
 713        format(i4,2x,f8.3,3(2x,f9.3))
         endif
      enddo
      close(70)

c     Jeans solution M(<r)+nu(r)+beta -> Sigma(R)+VDP(R),
c     following Mamon & Lokas (2005, MNRAS, 363, 705)
c     using best-fit parameters 

      errabs=0.
      rismin=1.d-190
      rismax=1.d190
      ninterp=21*2
      rinfinity=20.

      ius=80

      write(*,*) ' Evaluating expected VDP for '
      write(*,*) ' Max Lik solution values: '

 797  continue

      v200=10.*hz*r200

      write(*,516) r200,rc,r200/rc,rs,r200/rs,gam,cbe,xmstarl
 516  format('   r_200    = ',f6.3,
     &       /'   r_tracer = ',f6.3,' (c_tracer= ',f7.2,')'
     &       /'   r_mass   = ',f6.3,' (c_mass=   ',f7.2,')'
     &       /'   gamma    = ',f6.3,
     &       /'   Anisotropy parameter = ',f10.3,
     &       /'   M*/L = ',f6.3,/)

      errrel=0.0001d0

c     evaluate sigma_r at some points (it will then interpolate)

      rlow=0.005d0    ! we want the profile in the inner region (12 Dec. 2014)

      do i=1,ninterp
         xx2 = dlog(2.d0*rinfinity)
         xx1 = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         risl = dcadre(sr2int,xx1,xx2,errabs,errrel,errest,ier)
         if (ier .gt. 128) then
            print *,'MAIN: rmin=',xx1,' rmax=',xx2
            stop
         endif
         xmin=dexp(xx1)
         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         xris(i)=xx1
         yris(i)=dlog(risok)
      enddo
cc      write(*,*) ' sigma_r evaluated'

c
c compute spline coeffs for later interpolation of sigma_r
c
      call icsccu(xris,yris,ninterp,csr,icsr,ier)
      if (ier .gt. 128) then
         print *,' after sigma_r evaluated...'
         do i = 1, ninterp
            print *,' i xris yris = ', i, xris(i), yris(i)
         enddo
      endif

      errrel=0.0001d0
      errabs=0.

      do i=1,50
         xx1 = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(50-1)
         xx2=dlog(2.*rinfinity)
         xmin=dexp(xx1)
         ris2n = dcadre(fa,xx1,xx2,errabs,errrel,errest,ier)
         if (ier .gt. 128) then
            print *,'MAIN-FA: xmin=',xmin
            ris2n=0.d0
         endif
         if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then 
            pnr=sigmar1(xmin)
         elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
            pnr=sigmar2(xmin)
         else
            pnr=sigmar3(xmin)
         endif
         ris=dsqrt(ris2n/pnr)
         if (ris.eq.ris) write(ius,*) xmin,ris
      enddo
      close(ius)


      end


c
c
c
c------------------------------------------------------------------------------
        SUBROUTINE SORT(N,RA)
c------------------------------------------------------------------------------
 
c--- Routine to do a heapsort of a data array RA
c    Stolen (unabashedly) from NUMERICAL RECIPES
c
c**************************************************************************** 

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      dimension ra(n)
 
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end

c------------------------------------------------------------------------------
      SUBROUTINE SORT2(N,RA,RB)
c------------------------------------------------------------------------------
 
c--- Routine to do a heapsort of a data array RA and RB
c    Stolen (unabashedly) from NUMERICAL RECIPES
c
c**************************************************************************** 
 
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension ra(n),rb(n)

      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb=rb(l)
        else
          rra=ra(ir)
          rrb=rb(ir)
          ra(ir)=ra(1)
          rb(ir)=rb(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1)=rrb
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(i)=rb(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
        rb(i)=rrb
      go to 10
      end

c----------------------------------------------------------------------------- 
        SUBROUTINE XMIDMEAN(XDATA,N,XMID)
c----------------------------------------------------------------------------- 
 
c--- This subroutine calculates the MIDMEAN for a set of N ORDERED
c    statistics stored in XDATA. The MIDMEAN is defined to be the
c    mean of the central 50% of the data set. This corresponds to 
c    a 25% Trimmed Mean. The value of the MIDMEAN is returned as 
c    XMID and is defined:
c
c                   XMID = TRIM(.25) 
c
c    where TRIM(.25) is the 25% trimmed mean as defined above. For
c    more information on the MIDMEAN or TRIMMED MEAN see pages 312,
c    313 and pages 307, 308 in UREDA.
 
c**************************************************************************** 
 
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
        dimension xdata(n)
        data n1,n2,d1,d2,df,zero/1,2,1.0,2.0,0.25,0.0/

        ig = int(df*n)
        r = (df*n) - dfloat(ig)
        sum1 = zero
        do 11 i=ig + n2,n-ig-n1
        sum1 = sum1 + xdata(i)
11      continue
        sum3 = (d1-r)*(xdata(ig+n1) + xdata(n-ig))
        xmid = (d1/(dfloat(n)*(d1-(d2*df))))*(sum3+sum1)

        return
        end


      subroutine robusti(xin,n,ibwt,c,s)

C     Routine that uses robust techniques to estimate the central location C
C     and the spread S, for a distribution of N sorted values X.
C     Based on the work of Beers, Flynn and Gebhardt, AJ 100, 32, 1990.
C     When necessary, use is made of Numerical Recipes routines.

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter (pi=3.1415926535897932,nmax=20000)

      dimension xin(n),x(nmax),xred(nmax)

      data      t1,t2 /6.0,9.0/

C     Reset output variables

      c=0.0
      s=0.0

C     Copy input data into work array

      do i=1,n
         x(i)=xin(i)
         c=c+xin(i)
      end do
      c=c/float(n)

C     If N<3 calculation of spread is impossible

      if (n.lt.3) return

C     Calculate the median M, given the sorted distribution of X

      call sort(n,x)
      n2=n/2
      if(2*n2.eq.n)then
        xm=0.5*(x(n2)+x(n2+1))
      else
        xm=x(n2+1)
      endif

C     Calculate the Median Absolute Deviation

      do i=1,n
         xred(i)=abs(x(i)-xm)
      end do                 
      call sort(n,xred)
      n2=n/2
      if(2*n2.eq.n)then
         xmad=0.5*(xred(n2)+xred(n2+1))
      else
         xmad=xred(n2+1)
      endif
      if (xmad.eq.0.0) then
         write(*,*) 'robust: no differentiation in data!'
         return
      end if

C     Calculate the biweight location estimator

      fact=1.0/(t1*xmad)
      sum1=0.0
      sum2=0.0
      do i=1,n
         u=(x(i)-xm)*fact
         if (abs(u).lt.1.0) then
            uhlp=(1.0-u**2)**2
            sum2=sum2+uhlp
            sum1=sum1+(x(i)-xm)*uhlp
         end if
      end do
      cbi=xm+sum1/sum2
      c=cbi

C     Calculate the biweight scale if requested

      if (ibwt.gt.0.5) then
         fact=1.0/(t2*xmad)
         sum1=0.0
         sum2=0.0
         do i=1,n
            u=(x(i)-xm)*fact
            if (abs(u).lt.1.0) then
               uhlp1=1.0-u**2
               uhlp2=uhlp1**4
               uhlp3=1.0-5.0*u**2
               sum1=sum1+(x(i)-xm)**2*uhlp2
               sum2=sum2+uhlp1*uhlp3
            end if
         end do
         s=sqrt(float(n))*sqrt(sum1)/abs(sum2)
      else

C     ... and use a gapper algorithm otherwise

         sum=0.0
         do i=1,n-1
            g=x(i+1)-x(i)
            w=float(i*(n-i))
            sum=sum+w*g
         end do
         s=sqrt(pi)/(n*(n-1))*sum
      end if

      return
      end

c
c     N(R) - beta-function
c
      function sigmar3(tt)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      t=tt/rc
ccc      sigmar3=(1.d0+t*t)**al/(r200*r200)
      sigmar3=(1.d0+t*t)**al
      return
      end
c
c     N(R) - projected NFW; use eq. (41), (42) in Lokas & Mamon (2001)
c
      function sigmar1(tt)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      t=tt/rc


c     no interlopers considered

         if (t.gt.1) then
            cm1x=dacos(1./t)
            ft=(1.-1./dsqrt(dabs(t*t-1))*cm1x)/(t*t-1.)
         elseif (t.lt.1) then
            cm1x=dacosh(1./t)
            ft=(1.-1./dsqrt(dabs(t*t-1))*cm1x)/(t*t-1.)
         else
            ft=1./3.
         endif

         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         sigmar1=ft*c*c/(2.*pig)*gc/(r200*r200)

      return
      end
c
c     N(R) - projected Hernquist - see Hernquist (1990)
c
      function sigmar2(tt)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      s=tt/rc

      if (s.gt.1) then
         xs=acos(1./s)/sqrt(s*s-1.)
         fr=((2.+s*s)*xs-3.)/(2.*pig*rc*rc*(1.-s*s)**2)
      elseif (s.lt.1) then
         xs=dlog((1.+sqrt(1.-s*s))/s)/sqrt(1.-s*s)
         fr=((2.+s*s)*xs-3.)/(2.*pig*rc*rc*(1.-s*s)**2)
      else
         xs=1.
         fr=2./(15.*pig*rc*rc)
      endif

      sigmar2=fr

      return
      end
c
c     Integrand function for the determination of
c     nu(r)*sigma_r^2(r), given beta(r) and M(r)
c     from eqs. A3, A4, A6 in Mamon & Lokas (2005)
c     or from eq.(A3) in Mamon, Biviano & Boue' (2013)
c     
c       modified by GAM to handle integration of ln r
c
      function sr2int(alr)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      dimension rsvalues(28),r100values(28)
      include 'paramsoptbcg.i'

      data rsvalues/0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,
     ,     0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,
     ,     0.30,0.35,0.40,0.45,0.50,1.00/
      data r100values/1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,
     ,     1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,
     ,     1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476/
      external gammarec
c
      gm200=r200*v200*v200
c
      t = dexp(alr)
      if (kmp.eq.1) then
c
c     DM mass profile is genealized-NFW
c
         if (kgas.eq.1) then
c
c     If requested, I consider the total mass as the sum of
c     the gas and DM mass profiles; beyond r100, the fraction
c     Mgas/Mtot is assumed to be constant, hence the shape of
c     Mtot becomes that of the DM mass profile, rescaled by
c     the appropriate factor (1.-0.167). The DM profile is
c     normalised in such a way as to give 0.857 (1-0.123-0.02) of the 
c     total mass at r200, which is what is generally found in groups 
c     and clusters. 0.02 is provided by the central BGG. The gas mass fraction 
c     becomes ~0.167 (the Universal baryon ratio) at ~r100.
c     Note that the radii have to be changed accordingly, r200 is
c     not r200 of the total mass: r200mod=r200*(1-0.143)^(1/3)
c
            ficm=0.123
            fbgg=0.020
            fac200=1.*(1.d0-ficm-fbgg)/
     &           (dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
            tm=t*(1.d0-ficm-fbgg)**(-1./3.)
            xmdm=(dlog(1.d0+tm/rs)-tm/rs/(1.d0+tm/rs))*fac200
c
c     Determine r100 by interpolation
c
            j1=0
            do j=1,28
               if (rsvalues(j).le.rs) j1=j
            enddo
            if (j1.eq.0) r100=r100values(1)
            if (j1.eq.28) r100=r100values(28)
            if (j1.gt.0.and.j1.lt.28) r100=(r100values(j1+1)-
     -           r100values(j1))/
     /           (rsvalues(j1+1)-rsvalues(j1))*(rs-rsvalues(j1))+
     +           r100values(j1)
c
c     IC gas mass profile: it is computed from integration
c     of a beta profile with r_c/r_200=0.01 and beta=0.5
c     i.e. close to the average values for the Osmond &
c     Ponman sample. The normalization from Sanderson
c     et al. (2003) but it is not very accurate. It indicates
c     Mgas/Mtot=0.02 at 0.3*r_200, a value representative
c     of groups. Such a value makes Mgas/Mtot=0.123 at r200,
c     and 0.167 (the universal baryon fraction, see Zhang et al. 06)
c     at 1.23*r_200, which is close to r_100.
c     I then adopt Mgas/Mtot=0.123 at r_200, for simplicity.
c     I use Mamon+Lokas 05b expression for Mgas(r) which proves 
c     accurate to within +-2.5% over the radial range 0.01-20 r_200.
c
            ga=-1.09051d0
            z=t/0.01d0
            xmg=1.87277d-4*((z**3.d0/3.d0)**ga+
     +           (2.d0/3.d0*z**(1.5d0))**ga)**(1.d0/ga)
            if (t.le.r100) then
               xm=xmdm+xmg+fbgg
            else
               xm=xmdm/(1.d0-ficm-fbgg)
            endif
         else
c
c     If not requested, the mass profile becomes
c     identical to NFW with the normalisation
c     set to v200^2*r200/G at r/r200=1
c

            bhyp1 = 3.d0-gam
            bhyp2 = 3.d0-gam
            ahyp = 4.d0-gam

            chyp = -r200/rs
            thyp = -t/rs

            hyptmp1=HYGFX(bhyp1,bhyp2,ahyp,chyp,hypgeoc)
            hyptmp2=HYGFX(bhyp1,bhyp2,ahyp,thyp,hypgeot)
      
            fac200g=(r200/rs)**(3.d0-gam) * hypgeoc
            gnfw   =(t/rs)**(3.d0-gam) * hypgeot

            xm=gnfw/fac200g

c$$$            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
c$$$            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200
         endif
      elseif (kmp.eq.2) then
c
c     M(r) is Hernquist; no gas contribution allowed for the time being
c     Normalisation set to 1 at r/r200=1
c
         fac200=(r200+rs)*(r200+rs)/(r200*r200)
         xm=t*t/(t+rs)/(t+rs)*fac200
      elseif (kmp.eq.3) then
c
c     M(r) is PIEMD; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c     
         fac200=1./(rcut*datan(r200/rcut)-rs*datan(r200/rs))
         xm=fac200*(rcut*datan(t/rcut)-rs*datan(t/rs))
      elseif (kmp.eq.4) then
c
c     M(r) is Burkert; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         trs=t/rs
         rvrs=r200/rs
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))
         xm=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))
      elseif (kmp.eq.5) then
c
c     M(r) is Soft Isoth Sph; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         fac200=1./(r200-rs*datan(r200/rs))
         xm=fac200*(t-rs*datan(t/rs))
      else
c
c     M(r) is Einasto m=5; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm=fac200*gip
cc         write(*,999) r200,rs,t,3.*eim,2.*eim*(r200/rs)**(1./eim),
cc     &        2.*eim*(t/rs)**(1./eim),gip200,gip,xm
cc 999     format(4(1x,f5.2),2(1x,f6.2),3(1x,d11.3))
cc         eim=5.
cc         fac200=1./dgammp(3.*eim,2.*eim*(r200/rs)**(1./eim))
cc         xm=fac200*dgammp(3.*eim,2.*eim*(t/rs)**(1./eim))
      endif

c     add the mass of the BCG in units of G*M200
c     (set the free param xmstarl to zero if this is not required)

      xmbcg=xmstarl*xlumbcg*t/rjaf/(1.+t/rjaf)*grav/gm200

      xm=xm+xmbcg

c
c     nu(r) is beta-model, NFW or Hernquist
c
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
         xnu=rc/(2.*pig*t)/(t+rc)**3.
      else
         xnu=-1.d0/dsqrt(pig*1.d0)*
     &        gammarec(-al+0.5d0)/gammarec(-al+1.d0)*
     &        al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
      endif

c
c     beta(r)
c 
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
      if (kani.eq.0) then
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         app=1.d-2                                         
c
c     Radial anisotropy
c
         if(dabs(bec-1.d0).lt.app) then
            sr2int=xnu*xm
c
c     Isotropic case
c
         elseif (dabs(bec).lt.app) then
            sr2int=xnu*xm/(t*t)
c
c     General constant anisotropy
c
         else
            sr2int=xnu*xm*t**(2.*bec-2.)
         endif
c
c     Osipkov Merritt
c
      elseif (kani.eq.2) then
         sr2int=xnu*xm*(t*t+cbe*cbe)/(t*t)
c
c     Mamon Lokas anisotropy profile
c
      elseif (kani.eq.1) then
         sr2int=xnu*xm*(t+cbe)/(t*t)
c
c     Simplified Wojtak anisotropy profile
c
      elseif (kani.eq.3) then
         sr2int=xnu*xm*(t**1.5+cbe**1.5)**(4./3.)/(t*t)
c
c     Simplified Tiret (modified ML profile)
c
      elseif (kani.eq.4) then
         rm2=rs                 ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(2.*bec)/(t*t)
c
c     modified Tiret ("opposite" model)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(4.*bec)*t**(-2.-2.*bec)
c
c     Tiret with Tied Anisotropy Light (TAL)
c     
      elseif (kani.eq.6) then
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
         sr2int=xnu*xm*(t+rc)**(2.*bec)/(t*t)
c
c     Tiret with inner tangetial orbits, sigmar/sigmat<1
c     (0.7 for bec0=-1, 0.4 for bec0=-5.25)
c     
      elseif (kani.eq.7) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
c     bec0=-1.
         bec0=-5.25
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(2.*(bec-bec0))*t**(-2.+2.*bec0)
c
c     Tiret with inner radial orbits, sigmar/sigmat>1
c     (1.3 for bec0=0.4, 1.6 for bec0=0.6)
c     
      elseif (kani.eq.8) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
c     bec0=0.4
         bec0=0.6
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(2.*(bec-bec0))*t**(-2.+2.*bec0)
c
c     Hansen+Moore relation for NFW mass profile
c
      else
         rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
         sr2int=xnu*xm*t**(2.*(ahm-bhm-1.))*(rhm+t)**(-4.*bhm)
      endif

      sr2int = t*sr2int*gm200    ! integral in dlog and G*M scaling
c
      return
      end

c
c     The factor outside the integral for the
c     sigma_r^2 formulae
c
      function sr2out(t)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      external gammarec

      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
         xnu=rc/(2.*pig*t)/(t+rc)**3.
      else
         xnu=-1.d0/dsqrt(pig*1.d0)*
     &        gammarec(-al+0.5d0)/gammarec(-al+1.d0)*
     &        al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
      endif

      if (kani.eq.0) then
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         app=1.d-2                                         
c
c     Radial anisotropy
c
         if(dabs(bec-1.d0).lt.app) then
            sr2out=1./(t*t)
c
c     Isotropic case
c
         elseif (dabs(bec).lt.app) then
            sr2out=1.
c
c     General constant anisotropy
c
         else
            sr2out=t**(-2.*bec)
         endif
c
c     Osipkov Merritt
c
      elseif (kani.eq.2) then
         sr2out=1./(t*t+cbe*cbe)
c
c     Mamon Lokas anisotropy profiles
c
      elseif (kani.eq.1) then
         sr2out=1./(t+cbe)
c
c     Simplified Wojtak anisotropy profile
c
      elseif (kani.eq.3) then
         sr2out=(t**1.5+cbe**1.5)**(-4./3.)
c
c     Simplified Tiret (modified ML)
c
      elseif (kani.eq.4) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-2.*bec)
c
c     modfied Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-4.*bec)*t**(2.*bec)
c
c     Tiret with Tied Anisotropy Light (TAL)
c     
      elseif (kani.eq.6) then
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
         sr2out=(t+rc)**(-2.*bec)
c
c     Tiret with inner tangetial orbits, sigmar/sigmat=0.7
c     
      elseif (kani.eq.7) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=-5.25
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-2.*(bec-bec0))*t**(-2.*bec0)
c
c     Tiret with inner radial orbits, sigmar/sigmat=1.3
c     
      elseif (kani.eq.8) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=0.6
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-2.*(bec-bec0))*t**(-2.*bec0)
c
c     Hansen+Moore relation for NFW mass profile
c
      else
         rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
         sr2out=t**(2.*(bhm-ahm))*(rhm+t)**(4.*bhm) ! My solution
      endif
c
c     Divide by the density profile
c
      sr2out=sr2out/xnu
      return
      end

c
c     Integrand to derive the f(R,vlos)
c     distribution function for given velocity
c     vlos, given anisotropy and mass profile
c     (from Gwenael Boue')
c
      double precision function gwen()
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 rjl(1), yrjl(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'tsallis.i'
      external gammarec
c
c     nu(t)
c
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
         xnu=rc/(2.*pig*t)/(t+rc)**3.
      else
         xnu=-1.d0/dsqrt(pig*1.d0)*
     &        gammarec(-al+0.5d0)/gammarec(-al+1.d0)*
     &        al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
      endif
c
c     Interpolate to get the sigma_r(t)
c
      rjl(1)=dlog(t)
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      if (ier .eq. 34) then
         print *,' ICSEVU in GWEN: max xris=', xris(ninterp), ' r=',rjl
      endif
      sigr=dexp(yrjl(1))
c
c     beta(r)
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
      if (kani.eq.0) then
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c     Osipkov Merritt
c
      elseif (kani.eq.2) then
         bec=t*t/(t*t+cbe*cbe)
c
c     Mamon Lokas
c     
      elseif (kani.eq.1) then
         bec=0.5*t/(t+cbe)
c
c     Simplified Wojtak
c
      elseif (kani.eq.3) then
         bec=t**1.5/(t**1.5+cbe**1.5)
c
c     Simplified Tiret (modified ML)
c
      elseif (kani.eq.4) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=(1.-1./(cbe*cbe))*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c     modified Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
c
c     Tiret with Tied Anisotropy Light (TAL)
c     
      elseif (kani.eq.6) then
         bec=(1.-1./(cbe*cbe))*t/(t+rc)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with inner tangetial orbits, sigmar/sigmat=0.7
c     
      elseif (kani.eq.7) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=-5.25 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with inner radial orbits, sigmar/sigmat=1.3
c     
      elseif (kani.eq.8) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=0.6 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Hansen+Moore relation for NFW mass profile
c
      else
         rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
         bec=ahm-bhm*(rhm+3.*t)/(rhm+t)
      endif
c
c     sigma_z(R,r)
c
      sigmaz=sigr*dsqrt(1.-bec*xmin*xmin/(t*t))

c     add velocity error

      sigmaz=sqrt(sigmaz*sigmaz+ej*ej)
c
c     Integrand function (note that sigmaz is in units of v200, just
c     as vj)

c     we have a Gaussian MAMPOSSt and a Tsallis MAMPOSSt

      if (ktsa.lt.0.5) then
      
         gwen=2.d0*xnu/dsqrt(2.*pig)/(sigmaz)*
     *        t/dsqrt(t*t-xmin*xmin)*
     *        dexp(-vj*vj/(2.*sigmaz*sigmaz))

      else
         gweni=(1.d0-batsa*vj*vj/(2.d0*sigmaz*sigmaz))
         if (gweni.gt.0) then 
            gwen=2.d0*xnu/dsqrt(2.*pig)/(sigmaz)*
     *           t/dsqrt(t*t-xmin*xmin)*
     *           atatsa*gweni**(atsa+1.d0)
         else
            gwen=1.d-30
         endif
      endif

      return
      end
c
c SIGMARBYINTERP evaluates sigma_r at set of radii given in rvec 
c by cubic-spline interpolation
c
      subroutine sigmarbyinterp(rvec,sigrinterp,n)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer n
      real*8 rvec(n), sigrinterp(n)
      real*8 alsigma(100), alrvec(100)
      real*8 rjl(1), yrjl(1)
      integer ier
      include 'sr.i'
      if (n .gt. 100) then
         print *,' SIGMARBYINTERP: dimension of rvec = ', n, 
     &    ' is too large'
         stop
      endif
      do i = 1, n
         alrvec(i) = dlog(rvec(i))
      enddo
      call icsevu(xris,yris,ninterp,csr,icsr,alrvec,alsigma,n,ier)
      do i = 1, n
         sigrinterp(i) = dexp(alsigma(i))
      enddo
c
      return
      end

c
c GWENU returns integrand of g(R,vz) integrated over u = cosh-1 (r/R)
c
c     Integrand to derive the f(R,vlos)
c     distribution function for given velocity
c     vlos, given anisotropy and mass profile
c     (from Gwenael Boue')
c
      double precision function gwenu(u)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 rjl(1), yrjl(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'tsallis.i'
      external gammarec
c
      t = xmin*dcosh(u)

cc      write(*,*) ' vel and error is',vj,ej

c
c     tracer density nu
c
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
         xnu=rc/(2.*pig*t)/(t+rc)**3.
      else
         xnu=-1.d0/dsqrt(pig*1.d0)*
     &        gammarec(-al+0.5d0)/gammarec(-al+1.d0)*
     &        al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
      endif
c
c     Interpolate to get the sigma_r(t)
c
      rjl(1)=dlog(t)
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      if (ier .eq. 34) then
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
     &    rjl(1)
      endif
      sigr=dexp(yrjl(1))

c
c     beta(r)
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
      if (kani.eq.0) then
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c     Osipkov Merritt
c
      elseif (kani.eq.2) then
         bec=t*t/(t*t+cbe*cbe)
c
c     Mamon Lokas
c     
      elseif (kani.eq.1) then
         bec=0.5*t/(t+cbe)
c
c     Simplified Wojtak
c
      elseif (kani.eq.3) then
         bec=t**1.5/(t**1.5+cbe**1.5)
c
c     Simplified Tiret or modified ML
c
      elseif (kani.eq.4) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=(1.-1./(cbe*cbe))*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c     modified Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
c
c     Tiret with inner tangetial orbits, sigmar/sigmat=0.7
c     
      elseif (kani.eq.7) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=-5.25 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with inner radial orbits, sigmar/sigmat=1.3
c     
      elseif (kani.eq.8) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=0.6 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with Tied Anisotropy Light (TAL)
c     
      elseif (kani.eq.6) then
         bec=(1.-1./(cbe*cbe))*t/(t+rc)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Hansen+Moore relation for NFW mass profile
c
      else
         rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
         bec=ahm-bhm*(rhm+3.*t)/(rhm+t)
      endif
c
c     sigma_z(R,r)
c
      sigmaz=sigr*dsqrt(1.-bec*xmin*xmin/(t*t))

c     add velocity error

      sigmaz=sqrt(sigmaz*sigmaz+ej*ej)

c
c     Integrand function (note that sigmaz is in units of v200, just
c     as vj)

c     we have a Gaussian MAMPOSSt and a Tsallis MAMPOSSt

      if (ktsa.lt.0.5) then
      
         gwenu=2.d0*xnu/dsqrt(2.*pig)/sigmaz*dcosh(u)
     *        *dexp(-vj*vj/(2.*sigmaz*sigmaz))

      else
         gweni=(1.d0-batsa*vj*vj/(2.d0*sigmaz*sigmaz))
         if (gweni.gt.0) then 
            gwenu=2.d0*xnu/dsqrt(2.*pig)/(sigmaz)*dcosh(u)
     *         *atatsa*gweni**(atsa+1.d0)
         else
            gwenu=1.d-30
         endif
      endif

      return
      end
c
c     Subroutine for the computation of Max Lik
c     
      subroutine vmaxlik(nfv,xfv,f)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension xfv(nfv),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),wmbm(25000),embm(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'probs.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      if (dabs(al).lt.0.0001) n1=1

cc      write(*,*) ' Prima: rs, cbe = ',rs,cbe

      rs=xfv(1)
      cbe=xfv(2)

cc      write(*,*) ' Dopo: rs, cbe = ',rs,cbe

      ngambm=nga
      do j=1,nga
         rmbm(j)=r(j)
         vmbm(j)=v(j)
         embm(j)=e(j)
         wmbm(j)=w(j)
      enddo

      irule=2
      errabs=0.d0
      errrel=0.001d0
      rismin=1.d-190
      rismax=1.d190
      rinfinity=20.0d0

      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin
         risl = dcadre(sr2int,xris(i),dlog(2.*rinfinity),errabs,errrel,
     &    errest,ier)
         if (ier .gt. 128) then
            print *,'VMAXLIK-SR2INT: rmin=',xx1,' rmax=',xx2
            stop
         endif
         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)
      enddo

c
c compute spline coeffs for later interpolation of sigma_r
c
      call icsccu(xris,yris,ninterp,csr,icsr,ier)
      if (ier .gt. 128) then
         print *,' in VMAXLIK ...'
         do i = 1, ninterp
            print *,' i xris yris = ', i, xris(i), yris(i)
         enddo
      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      if (kintd.gt.0) then
         rinfinity=0.99999*r200
      endif

c     Now compute the distribution function at the radial distance of
c     each galaxy by interpolating the sigma_r(r) and then
c     evaluate the Max Lik, using weights if required

c     If required, make the estimate only on a grid of values
c     and then bispline interpolate in the R,v positions of the
c     galaxies (kbsp=1)
c
c     uses an average velocity error

      if (kbsp.eq.1) then

         rmin=1.e12
         vmin=1.e12
         rmax=-1.e12
         vmax=-1.e12
         eveave=0.
         do j=1,ngambm
            if (rmbm(j).lt.rmin) rmin=rmbm(j)
            if (abs(vmbm(j)).lt.vmin) vmin=abs(vmbm(j))
            if (rmbm(j).gt.rmax) rmax=rmbm(j)
            if (abs(vmbm(j)).gt.vmax) vmax=abs(vmbm(j))
            eveave=embm(j)+eveave
         enddo
         eveave=eveave/ngambm

         nr=rknots
         nv=vzknots
         kspo=3
         npts=nr*nv
         stepr=(rmax-rmin)/(nr-1)
         stepv=(vmax-vmin)/(nv-1)
         rmi1=rmin*0.9
         rma1=rmax*1.1
         if (rma1.ge.rinfinity) rma1=0.9999d0*rinfinity
         do j=1,nr
            rp(j)=10.**(dlog10(rmi1)+(j-1.)*dlog10(rma1/rmi1)/(nr-1.))
         enddo
         do j=1,nv
            vp(j)=vmin+stepv*(j-1)
         enddo

         do j=1,nr
            rj=rp(j)
            do i=1,nv
               vj=vp(i)
               ej=eveave
               xmin=rj
               k=k+1
c
c     re-define rinfinity as the radius where v_z/sigma_r<0.1
c     to avoid a problem in v_z=0
c
               if (kintd.eq.0) call rinf(rinfinity)
c
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c
               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)

               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN: rmin=',xx1,' rmax=',xx2
                  stop
               endif
c     
c     Add interloper contribution if requested
c
               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
     &                 errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
cc                  write(*,*) ' Number of obsd galaxies ',ngambm
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
cc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
              endif
               
               if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm


               else

                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                     sigmaden=sigmar1(rj)
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                     sigmaden=sigmar2(rj)
                  else
                     sigmaden=sigmar3(rj)
                  endif
                  p=g/sigmaden

               endif

               pv(k)=p
               rpv(k)=rj
               vpv(k)=vj

c     The log(probability) is stored in z(j,i)

               z(j,i)=dlog(p)
            enddo
         enddo
         npv=k

c     Bi-spline 2-d interpolation to obtain the p values
c     of all data points from the p values of the sparse grid
     
         call ibcccu(z,rp,nr,vp,nv,c,ic,wk,ier)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp,ier)
            psum=psum-dlogp*wmbm(l)
 
         enddo

      else

c     Use all the original data, no interpolation

         do j=1,ngambm
            rj=rmbm(j)
            vj=vmbm(j)
            ej=embm(j)
            xmin=rj

c     re-define rinfinity as the radius where v_z/sigma_r leq 0.1
c     to avoid a problem in v_z=0; if we require to consider the 
c     interloper contribution, rinfinity is set to 1 already

            if (kintd.eq.0) call rinf(rinfinity)

c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)

               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
                  stop
               endif


c     Add interloper contribution if requested

               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4.+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  n2=2
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,
     &                 errrel,errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
                  g=g+gi
               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally

               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm
               
               else

                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                     sigmaden=sigmar1(rj)
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                     sigmaden=sigmar2(rj)
                  else
                     sigmaden=sigmar3(rj)
                  endif
                  p=g/sigmaden

               endif

               pv(j)=p
               rpv(j)=rj
               vpv(j)=vj

               call temp(rj,sigrt1,sigrt2,sigzt1t2,bec)


               psum=psum-dlog(p)*wmbm(j)
               

            endif

         enddo

         npv=ngambm

      endif

c
c     The factor nga/wsum is =1 if no weights, but is not 1
c     if there are weights, and is needed to take into account
c     the real number of points (otherwise, we would estimate
c     a Max Lik too low by an average factor of 1/Ngr, where
c     Ngr is the average number of galaxies per group, since
c     the likelihood we want to estimate is the Sum of the
c     ngambm galaxy probabilities)
c     
      psum=psum/wsum*ngambm

      f=psum


      return
      end
c
c     Max Lik for beta-profile
c
      function fcn3(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      external sigmarnorm
      data pig /3.1415926535897932d0/

      fml=0.

      call sigmarnorm(x,fnorm)

      wsum=0.
      do j=1,nga
         fr=(1.+(r(j)/x)**2)**al
         fml=fml-dlog(fr/fnorm)*w(j)
         wsum=wsum+w(j)
      enddo

      fcn3=fml/w(j)

      return
      end
c
c     Max Lik function for projected NFW-profile
c
      function fcn1(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      external sigmarnorm
      data pig /3.1415926535897932d0/

      fml=0.

      call sigmarnorm(x,fnorm)

      c=r200/x
      gc=1./(dlog(c+1.)-c/(c+1.))
      facnorm=gc/(2.*pig)/(r200*r200)
      wsum=0.d0

      do j=1,nga
         cx=c*r(j)/r200
         uu=1./cx
         if (uu.lt.1) then
            cm1=dacos(uu)
            fr=c*c*(1.-1./dsqrt(dabs(cx*cx-1.))*cm1)/(cx*cx-1.)
         elseif (uu.gt.1) then
            cm1=dacosh(uu)
            fr=c*c*(1.-1./dsqrt(dabs(cx*cx-1.))*cm1)/(cx*cx-1.)
         else
            fr=c*c/3.
         endif
         fr=fr*facnorm
         fml=fml-dlog(fr/fnorm)*w(j)
         wsum=wsum+w(j)
      enddo

      fcn1=fml/wsum

      return
      end
c
c     Max Lik function for projected Hernquist-profile
c
      function fcn2(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension x2(2)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      external sigmarnorm
      data pig /3.1415926535897932d0/

      fml=0.

      call sigmarnorm(x,fnorm)

      wsum=0.d0
      do j=1,nga
         s=r(j)/x
         if (s.gt.1.d0) then
            xs=dacos(1./s)/dsqrt(s*s-1.)
            fr=((2.+s*s)*xs-3.)/(2.*pig*x*x*(1.-s*s)**2)
         elseif (s.lt.1.d0) then
            xs=dlog((1.+dsqrt(1.-s*s))/s)/dsqrt(1.-s*s)
            fr=((2.+s*s)*xs-3.)/(2.*pig*x*x*(1.-s*s)**2)
         endif
         if (abs(s-1.d0).lt.0.001d0) then
            fr=2./(15.*pig*x*x)
         endif
         if (r(j).le.0.or.fr.le.0) write(*,*) 'WARNING!',r(j),x,s,xs,fr
         fml=fml-dlog(fr/fnorm)*w(j)
         wsum=wsum+w(j)
      enddo

      fcn2=fml/wsum

      return
      end

c
c     determine infinity as the radius where sigma_z/v_z leq 0.1
c
      subroutine rinf(t)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 rjl(1), yrjl(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      external gammarec

      do j=1,1000
         t=0.01d0*dfloat(j)+10.d0
c
c     nu(t)
c
         if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
            c=r200/rc
            gc=1./(dlog(c+1)-c/(c+1))
            xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/
     ,           (r200*r200*r200)
         elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
            xnu=rc/(2.*pig*t)/(t+rc)**3.
         else
            xnu=-1.d0/dsqrt(pig*1.d0)*gammarec(-al+0.5d0)/
     &           gammarec(-al+1.d0)*
     &           al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
         endif
c     
c     Interpolate to get the sigma_r(t)
c
         rjl(1)=dlog(t)
         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      if (ier .eq. 34) then
         print *,' ICSEVU in RINF: max xris=', xris(ninterp), ' r=',rjl
      endif
         sigr=dexp(yrjl(1))
c
c     beta(r)
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
         if (kani.eq.0) then
            bec=1.-1./(cbe*cbe)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Osipkov Merritt
c
         elseif (kani.eq.2) then
            bec=t*t/(t*t+cbe*cbe)
c
c     Mamon Lokas
c     
         elseif (kani.eq.1) then
            bec=0.5*t/(t+cbe)
c
c     Simplified Wojtak
c
         elseif (kani.eq.3) then
            bec=t**1.5/(t**1.5+cbe**1.5)
c
c     Simplified Tiret or modified ML
c
         elseif (kani.eq.4) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            bec=(1.-1./(cbe*cbe))*t/(t+rm2)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     modified Tiret (non-zero central anisotropy)
c
         elseif (kani.eq.5) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
            if (bec.le.-1.) bec=-0.999d0 ! avoid unphysical values
c
c     Tiret with Tied Anisotropy Light (TAL)
c     
         elseif (kani.eq.6) then
            bec=(1.-1./(cbe*cbe))*t/(t+rc)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with inner tangetial orbits, sigmar/sigmat=0.7
c     
      elseif (kani.eq.7) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=-5.25 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Tiret with inner radial orbits, sigmar/sigmat=1.3
c     
      elseif (kani.eq.8) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4) rm2=1.521*rs ! Burkert
         bec0=0.6 
         bec=bec0+(1.-1./(cbe*cbe)-bec0)*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     Hansen+Moore relation for NFW mass profile
c
         else
            rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
            bec=ahm-bhm*(rhm+3.*t)/(rhm+t)
         endif
c     
c     sigma_z(R,r)
c
         sigmaz=sigr*dsqrt(1.-bec*xmin*xmin/(t*t))
c
c     Ratio sigma_z/v_z
c
         rat=dabs(sigmaz/vj)
         if (rat.le.0.1) return
      enddo

      return
      end


c


c
c     Integrand (eq. A15 in Mamon & Lokas 2005b)
c     for constant or null beta and NFW or beta-function
c     nu(r) - using eq. A16 in Mamon & Lokas but with
c     a modification introduced after complaints by 
c     Oleg Gnedin for cst-beta
c
      function fa(tlog)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 rjl(1), yrjl(1)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      dimension rsvalues(28),r100values(28)
      include 'paramsoptbcg.i'
      include 'sr.i'
      external betairec
      external gammarec
      
      data rsvalues/0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,
     ,     0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,
     ,     0.30,0.35,0.40,0.45,0.50,1.00/
      data r100values/1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,
     ,     1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,
     ,     1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476/

      t=dexp(tlog)

c     check on value of t and xmin (t must always be > xmin)

      if (xmin.ge.t) then
         fa=0.d0
         return
cc         t=xmin*1.001d0
      endif

c
c     physical units
c
      gm200=r200*v200*v200

c
      if (kmp.eq.1) then
c
c     DM mass profile is NFW
c
         if (kgas.eq.1) then
c
c     If requested, I consider the total mass as the sum of
c     the gas and DM mass profiles; beyond r100, the fraction
c     Mgas/Mtot is assumed to be constant, hence the shape of
c     Mtot becomes that of the DM mass profile, rescaled by
c     the appropriate factor (1.-0.167). The DM profile is
c     normalised in such a way as to give 0.857 (1-0.123-0.02) of the 
c     total mass at r200, which is what is generally found in groups 
c     and clusters. 0.02 is provided by the central BGG. The gas mass fraction 
c     becomes ~0.167 (the Universal baryon ratio) at ~r100.
c     Note that the radii have to be changed accordingly, r200 is
c     not r200 of the total mass: r200mod=r200*(1-0.143)^(1/3)
c
            ficm=0.123
            fbgg=0.020
            fac200=1.*(1.d0-ficm-fbgg)/
     &           (dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
            tm=t*(1.d0-ficm-fbgg)**(-1./3.)
            xmdm=(dlog(1.d0+tm/rs)-tm/rs/(1.d0+tm/rs))*fac200
c
c     Determine r100 by interpolation
c
            j1=0
            do j=1,28
               if (rsvalues(j).le.rs) j1=j
            enddo
            if (j1.eq.0) r100=r100values(1)
            if (j1.eq.28) r100=r100values(28)
            if (j1.gt.0.and.j1.lt.28) r100=(r100values(j1+1)-
     -           r100values(j1))/
     /           (rsvalues(j1+1)-rsvalues(j1))*(rs-rsvalues(j1))+
     +           r100values(j1)
c
c     IC gas mass profile: it is computed from integration
c     of a beta profile with r_c/r_200=0.01 and beta=0.5
c     i.e. close to the average values for the Osmond &
c     Ponman sample. The normalization from Sanderson
c     et al. (2003) but it is not very accurate. It indicates
c     Mgas/Mtot=0.02 at 0.3*r_200, a value representative
c     of groups. Such a value makes Mgas/Mtot=0.123 at r200,
c     and 0.167 (the universal baryon fraction, see Zhang et al. 06)
c     at 1.23*r_200, which is close to r_100.
c     I then adopt Mgas/Mtot=0.123 at r_200, for simplicity.
c     I use Mamon+Lokas 05b expression for Mgas(r) which proves 
c     accurate to within +-2.5% over the radial range 0.01-20 r_200.
c
            ga=-1.09051d0
            z=t/0.01d0
            xmg=1.87277d-4*((z**3.d0/3.d0)**ga+
     +           (2.d0/3.d0*z**(1.5d0))**ga)**(1.d0/ga)
            if (t.le.r100) then
               xm=xmdm+xmg+fbgg
            else
               xm=xmdm/(1.d0-ficm-fbgg)
            endif
         else
c
c     If not requested, the mass profile becomes
c     identical to NFW with the normalisation
c     set to 1 at r/r200=1
c

            bhyp1 = 3.d0-gam
            bhyp2 = 3.d0-gam
            ahyp = 4.d0-gam

            chyp = -r200/rs
            thyp = -t/rs

            hyptmp1=HYGFX(bhyp1,bhyp2,ahyp,chyp,hypgeoc)
            hyptmp2=HYGFX(bhyp1,bhyp2,ahyp,thyp,hypgeot)
      
            fac200g=(r200/rs)**(3.d0-gam) * hypgeoc
            gnfw   =(t/rs)**(3.d0-gam) * hypgeot

            xm=gnfw/fac200g

c$$$            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
c$$$            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200
         endif
      elseif (kmp.eq.2) then
c
c     M(r) is Hernquist; no gas contribution allowed for the time being
c     Normalisation set to 1 at r/r200=1
c
         fac200=(r200+rs)*(r200+rs)/(r200*r200)
         xm=t*t/(t+rs)/(t+rs)*fac200
      elseif (kmp.eq.3) then
c
c     M(r) is PIEMD; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         fac200=1./(rcut*datan(r200/rcut)-rs*datan(r200/rs))
         xm=fac200*(rcut*datan(t/rcut)-rs*datan(t/rs))
      elseif (kmp.eq.4) then
c
c     M(r) is Burkert; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         trs=t/rs
         rvrs=r200/rs
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))
         xm=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))
      elseif (kmp.eq.5) then
c
c     M(r) is Soft Isoth Sph; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         fac200=1./(r200-rs*datan(r200/rs))
         xm=fac200*(t-rs*datan(t/rs))
      else
c
c     M(r) is Einasto m=5; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm=fac200*gip
      endif

c     add the mass of the BCG in physical units
c     (set the free param xmstarl to zero if this is not required)

      xmbcg=xmstarl*xlumbcg*t/rjaf/(1.+t/rjaf)

      xm=xm*gm200 ! bring to physical units the cluster mass

      xm=xm+xmbcg


c
c     nu(r) from NFW, Hernquist or beta-model
c
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
         c=r200/rc
         gc=1./(dlog(c+1)-c/(c+1))
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
         xnu=rc/(2.*pig*t)/(t+rc)**3.
      else
         xnu=-1.d0/dsqrt(pig*1.d0)*
     &        gammarec(-al+0.5d0)/gammarec(-al+1.d0)*
     &        al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
      endif

      u=t/xmin
c
c     Constant anisotropy
c
      if (kani.eq.0) then
c
c     Convert from cbe=beta' to bec=beta 
c     in the case of constant anisotropy
c
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         app=1.d-2
c
c     Radial
c
         if(dabs(bec-1.d0).lt.app) then
            xk=pig/4.d0*u-0.5d0*dsqrt(1.d0-1.d0/(u*u))-0.5d0*u*
     *           dasin(1.d0/u)
c
c     Isotropic
c
         elseif (dabs(bec).lt.app) then
            xk=dsqrt(1.d0-1.d0/(u*u))
c
c     Negative integer
c
         elseif (dabs((bec-nint(bec))/bec).lt.app.and.bec.lt.0) then
            xk1=dsqrt(u*u-1.d0)/u**(1.d0-2.d0*bec)
            xk2=0.d0
            xk3=0.d0
            do k=0,nint(-bec)
               xk2 = xk2 + bico(nint(-bec),k)*
     &              (u*u-1.d0)**k/(2.d0*k+1.d0)
            enddo
            do k=0,nint(-bec-1)
               xk3=xk3+bico(nint(-bec-1),k)*
     &              (u*u-1.d0)**k/(2.d0*k+1.d0)
            enddo
            xk=xk1*(xk2-bec*xk3)
c
c     +1/2 or -1/2
c
         elseif ((dabs(bec-0.5d0).lt.app).or.
     &           (dabs(bec+0.5d0).lt.app)) then
            xk=u**(2.d0*bec-1.d0)*dacosh(u)-
     &           bec*dsqrt(1.d0-1.d0/(u*u))
c
c     Generic constant anisotropy, using
c     expression in Mamon+Lokas (2006, MNRAS, 370, 1582)
c
         else
            xxx=1.d0/(u*u)
            xk1=dsqrt(1.d0-xxx)/(1.d0-2.d0*bec)
            xk2=sqrt(pig)/2.d0*gammarec(bec-0.5d0)/gammarec(bec)
            xk3=(1.5d0-bec)*u**(2.d0*bec-1.d0)
            xk4=1.d0-betairec(xxx)
            xk=xk1+xk2*xk3*xk4
         endif
         fa=2.*xk*xm*xnu/t
c
c     cbe is the radius in eq.(61) of Mamon & Lokas 2005 
c     (Osipkov-Merritt anisotropy)
c
      elseif (kani.eq.2) then
         ua=cbe/xmin
         xk=(ua*ua+0.5d0)/(ua*ua+1.d0)**(1.5d0)*(u*u+ua*ua)/u*
     &        datan(dsqrt((u*u-1.d0)/(ua*ua+1.d0)))-
     &        0.5d0/(ua*ua+1.d0)*dsqrt(1.d0-1.d0/(u*u))
         fa=2.*xk*xm*xnu/t
c
c     cbe is the radius ra in eq.(60) of Mamon & Lokas 2005
c     in this case

      elseif (kani.eq.1) then
         ua=cbe/xmin
         if (dabs(ua-1.d0).lt.1.d-5) then
            xk=(1.d0+1.d0/u)*dcosh(u)-
     &           (8.d0/u+7.d0)/6.d0*dsqrt((u-1.d0)/(u+1.d0))
         else
            if (ua.gt.1.d0) then
               ccc=dacosh((ua*u+1.d0)/(u+ua))
               sss=1.d0
            else
               ccc=dacos((ua*u+1.d0)/(u+ua))
               sss=-1.d0
            endif
            fff=ua*(ua**2d0-0.5d0)/(dabs(ua**2.d0-1.d0))**1.5d0*
     &           (1.d0+ua/u)*ccc
            xk=0.5d0/(ua*ua-1.d0)*dsqrt(1.d0-1.d0/(u*u))+
     &           (1.d0+ua/u)*dacosh(u)-sss*fff
         endif
         fa=2.*xk*xm*xnu/t
c
c     if the anisotropy profile is a simplified Wojtak
c     or a simplified Tiret (modified ML) there are no 
c     analytical solution to provide K(r,ra) in eq.(A8) 
c     of ML05b, so we use the expression that includes
c     sigma_r for which we have a spline interpolation
c
      else

c     choose between simplified Wojtak, simpl. Tiret,
c     a model similar to Tiret (modified Tiret)
c     with non-zero central anisotropy, TAL Tiret and Hansen+Moore

         if (kani.eq.3) then
            b=t**1.5/(t**1.5+cbe**1.5)
         elseif (kani.eq.4) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            b=(1.-1./(cbe*cbe))*t/(t+rm2)
         elseif (kani.eq.5) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            b=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         elseif (kani.eq.6) then
            b=(1.-1./(cbe*cbe))*t/(t+rc)
         elseif (kani.eq.7) then
            bec0=-5.25
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            b=bec0+(1.-1./(cbe*cbe))*t/(t+rm2)
         elseif (kani.eq.8) then
            bec0=0.6
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4) rm2=1.521*rs ! Burkert
            b=bec0+(1.-1./(cbe*cbe))*t/(t+rm2)
         else   ! Hansen+Moore
            rhm=rc ! use nu(r)
cc            rhm=rs ! use rho(r)
            b=ahm-bhm*(rhm+3.*t)/(rhm+t)            
         endif
         if (b.ge.1.) b=0.999d0 ! avoid unphysical values
         if (b.le.-1.) b=-0.999d0 ! avoid unphysical values

c     interpolate sigma_r

         rjl(1)=tlog
         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
         sr=dexp(yrjl(1))
         fa=2.*xnu*sr**2.*(1.-b*(xmin/t)**2)*t/dsqrt(t*t-xmin*xmin)
      endif

ccc checking constant beta solution using double integral rather than kernel

ccc      fa1=fa
ccc      if (kani.eq.0) then
ccc         rjl(1)=tlog
ccc         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
ccc         sr=dexp(yrjl(1))
ccc         fa=2.*xnu*sr**2.*(1.-bec*(xmin/t)**2)*t/dsqrt(t*t-xmin*xmin)
ccc         if (abs(1.-fa/fa1).gt.0.1)  write(*,915) xmin,t,fa1/fa,bec
ccc 915     format(5(2x,f6.3))
ccc      endif

ccc remove the above part after check done

      fa=t*fa   ! integral in dlog

      return
      end
c
c     bootstrap y -> yboo
c
      subroutine boot(y,w,yboo,wboo,ny)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension y(ny),w(ny),yboo(ny),wboo(ny),r(25000)
      idum = 12234
      do j=1,ny
         r(j) = ran2(idum)
         jj=nint((ny-1)*r(j)+1)
         yboo(j)=y(jj)
         wboo(j)=w(jj)
      enddo
      return
      end

c
c
      SUBROUTINE REORDR(N,ARRAY,INDX)
      implicit real*8 (a-h,l-m,o-z)
      DIMENSION ARRAY(N), INDX(N)
      PARAMETER (NMAX=20000)
      REAL*8    HLP(NMAX)
      DO I=1,N
         HLP(I)=ARRAY(INDX(I))
      END DO         
      DO I=1,N
         ARRAY(I)=HLP(I)
      END DO         
      RETURN
      END 
c
c
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit real*8 (a-h,o-z)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

c     The MAMPOSSt procedure

      subroutine mamposst(di,ve,eve,rso,vso,npg1)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,clight=2.99792458d5)
      dimension di(npg1),ve(npg1),eve(npg1),rso(npg1),vso(npg1),
     &     xfv(2)
      dimension freepar(6),freelow(6),freeup(6),wpar(6000),xi(6,6)

      include 'paramsoptbcg.i'
      include 'datarv.i'
      include 'units.i'
      include 'free.i'
      include 'probs.i'
      include 'omegas.i'

      external fcn1,fcn2,fcn3,fa
      external sr2int,sr2out,sigmar1,sigmar2,sigmar3,gwenu

      print *,' Entering MAMPOSSt subroutine'

c      grra=1.4  ! this adjust the grid range, larger values => wider grid
      grra=1.2  ! this adjust the grid range, larger values => wider grid

c     Select only galaxy members within the requested radial range

      rup=-1.e10
      rlow=1.e10
      jsel=0

c     The radial selection range is in units of Mpc

      rlowinmpc=rlowin
      rupinmpc=rupin

      do j=1,npg1
         if (di(j).ge.rlowinmpc.and.
     &        di(j).le.rupinmpc) then
            jsel=jsel+1
            r(jsel)=di(j)
            v(jsel)=(ve(j)-va)/(1.+va/clight)
            e(jsel)=eve(j)
            w(jsel)=w(j)
            if (r(jsel).le.rlow) rlow=r(jsel)
            if (r(jsel).ge.rup) rup=r(jsel)
         endif
      enddo
      nga=jsel
cc      write(*,*) ' Datum 10 is ',r(10),v(10),e(10)
      write(*,932) nga,rlow,rup
 932  format(/' Using ',i5,' galaxies in the fits',/,
     +     ' in the radial range: ',f6.3,'--',f6.3,' Mpc')

c     Best-fit to N(R) external to MAMPOSSt, if required

      if (nrc.eq.-2) then

c     NFW, Hernquist or beta-model (fixed alpha)

         rf=1.e20
         do j=1,2000
            xtry=j*0.002
            if (knfit.eq.1) then
               res=fcn1(xtry)
            elseif (knfit.eq.2) then
               res=fcn2(xtry)
            else
               res=fcn3(xtry)
            endif                
            if (res.lt.rf) then
               xf=xtry
               rf=res
            endif
         enddo

         rc=xf
         rcg=xf

         write(*,195) rc
 195     format(' Best-fit to N(R) outside MAMPOSSt: r_tr=',f6.3)

      endif

c     Loop over r200, if requested, searching in the range
c     of 0.5 to 1.8 the real r200

      if (nr200.ge.2) then 
         dd=(grra-1.0)*(500./nga)**(0.2)
         deltarv=dd/nr200
         kr200=-nr200/2-1
      else
         kr200=0
      endif

      fmlmingrid=1.e20

      r200g=r200
      gamg=gam

      write(*,*) ' '
      write(*,*) ' ... now running MAMPOSSt procedure ... '

c     Now run MAMPOSSt

      kboby=0   ! used to keep track of minimization algorithms already tried
      knewu=0 
      kpowe=0

 174  continue

      ktried=kboby+knewu+kpowe

c     Using NEWUOA or BOBYQA or POWELL minimization 
c     The free parameters are 4 in general but can be less than
c     4 according to the choice of the user

      ipar(1)=1  ! identify parameters that are to be minimized
      ipar(2)=1
      ipar(3)=1
      ipar(4)=1
      ipar(5)=1
      ipar(6)=1

      facguess=1.0   ! factor to move away from best guess

c     define free parameters

      ip=0
      if (nr200.gt.0) then   ! virial radius
         ip=ip+1
         freepar(ip)=dlog10(facguess*r200g)
      else
         ipar(1)=0
      endif
      if (nrc.gt.0) then     ! scale radius of tracers
         ip=ip+1
         freepar(ip)=dlog10(facguess*rcg)
      else
         ipar(2)=0
      endif
      if (nrs.gt.0) then     ! scale radius of total matter
         ip=ip+1
         freepar(ip)=dlog10(facguess*rsg)
      else
         ipar(3)=0
      endif
      if (nbs.gt.0) then     ! tracer velocity anisotropy parameter
         ip=ip+1
         freepar(ip)=dlog10(facguess*cbeg)
      else
         ipar(4)=0
      endif
      if (ngam.gt.0) then     ! inner slope of gNFW
         ip=ip+1
         freepar(ip)=dlog10(facguess*gamg)
      else
         ipar(5)=0
      endif
      if (nml.gt.0) then     ! stellar mass to light ratio of BCG
         ip=ip+1
         freepar(ip)=dlog10(facguess*xmstarlg)
      else
         ipar(6)=0
      endif

      nfreepar=ip

      write(*,622) nfreepar,ipar,r200g,rcg,rsg,cbeg,gamg,xmstarlg,
     &     kmp,kani
 622  format(/,' Number of free parameters = ',i1,/
     &         ' r200, r_tr, r_s, anis, gamma, M*/L = ',6(i1,1x),/,
     &         ' Initial guess values: ',6(f6.3,1x),/,
     &     ' mass model= ',i2,/,
     &     ' anisotropy model= ',i2,/)

      npoints=(nfreepar+1)*(nfreepar+2)/2


c     Setting switches

      if (nrc.eq.-1) klfm=1   ! Light follows mass
      if (nrs.eq.-1) kmfl=1   ! Mass follows Light
      if (nrs.eq.-2) klcdm=1 ! LCDM c=c(M)
      if (nbs.eq.-1) kaml=1 ! a_ML = r_s
      if (nbs.eq.-1) khm=1  ! HM beta(r)

c     start optmization unless not required

      if (kopt.lt.0.) then
         r200new=r200g
         rcnew=rcg
         rsnew=rsg
         cbnew=cbeg
         gamnew=gamg
         xmstarlgnew=xmstarlg
         goto 732
      endif

c     newuoa and bobyqa only work with at least 2 free params

      if (nfreepar.gt.1..and.ktried.lt.2.9) then

         if (kopt.eq.0.and.kboby.eq.0) then

            write(*,*) ' '
            write(*,*) ' Running optimization algorithm BOBYQA'
            write(*,*) ' '

            rhoend=1.d-4 
            iprint=0
            maxfun=5000
            rhobeg=0.05         ! should be 10% greatest expected change to a variable

            addlowup=0.7        ! constant to subtract or add to get lower and upper interval
            do ilu=1,nfreepar
               freelow(ilu)=freepar(ilu)-addlowup
               freeup(ilu)=freepar(ilu)+addlowup
            enddo

            CALL BOBYQA (nfreepar,npoints,freepar,freelow,freeup,
     &           rhobeg,rhoend,iprint,maxfun,wpar)

            if (freepar(1).ne.freepar(1).and.inew.lt.10) then
               write(*,*) ' '
               write(*,*) ' BOBYQA failure!'
               write(*,*) ' Try another optimization'
               kboby=1
               kopt=1
               goto 174
            endif

         elseif (kopt.eq.1.and.knewu.eq.0) then

            inew=-1
            rhoend=1.d-4 
            iprint=0
            maxfun=5000
            rhobeg=0.05         ! should be 10% greatest expected change to a variable

            write(*,*) ' '
            write(*,*) ' Running optimization algorithm NEWUOA'
            write(*,*) ' '

cc 294        inew=inew+1


            CALL NEWUOA (nfreepar,npoints,freepar,rhobeg,rhoend,
     &           iprint,maxfun,wpar)


            if (freepar(1).ne.freepar(1).and.inew.lt.10) then
               write(*,*) ' '
               write(*,*) ' NEWUOA failure!'
               write(*,*) ' Try another optimization'
               knewu=1
               kopt=0
               goto 174
            endif

         else

            ftol=1.d-4          ! desired fractional tolerance in function to be minimized
            npars=6  ! max number of free params
            nfreefunc=nfreepar  ! passed through a common to function func

            do ixi=1,npars
               do jxi=1,npars
                  xi(ixi,jxi)=0.d0
                  if (ixi.eq.jxi) xi(ixi,jxi)=1.d0
               enddo
            enddo

            write(*,*) ' '
            write(*,*) ' Running optimization algorithm POWELL'
            write(*,*) ' '

            call POWELL(freepar,xi,nfreepar,npars,ftol,iter,fret)

            if (freepar(1).ne.freepar(1).and.inew.lt.10) then
               write(*,*) ' '
               write(*,*) ' POWELL failure!'
               write(*,*) ' Try another optimization'
               kpowe=1
               kopt=1
               goto 174
            endif

         endif            
      else                      ! when there is only 1 free param call POWELL

         ftol=1.d-4             ! desired fractional tolerance in function to be minimized
         npars=6                ! max number of free params
         nfreefunc=nfreepar     ! passed through a common to function func
         
         do ixi=1,npars
            do jxi=1,npars
               xi(ixi,jxi)=0.d0
               if (ixi.eq.jxi) xi(ixi,jxi)=1.d0
            enddo
         enddo
         
         write(*,*) ' '
         write(*,*) ' Running optimization algorithm POWELL'
         write(*,*) ' '
         
         call POWELL(freepar,xi,nfreepar,npars,ftol,iter,fret)
         
         if (freepar(1).ne.freepar(1).and.inew.lt.10) then
            write(*,*) ' '
            write(*,*) ' POWELL failure!'
         endif
         
      endif

      call calfun(nfreepar,freepar,f)
      
      call freepareva(nfreepar,freepar,r200new,rcnew,rsnew,cbnew,gamnew,
     &     xmstarlnew)

 732  continue

      do j=1,nga
         rso(j)=r(j)
         vso(j)=v(j)
      enddo

c
c     Grid search around the minimum found
c

      klfm=0
      if (nrc.ge.1) then
         rcg=rcnew
         dd=grra*(500./nga)**(0.2)  ! the term (500./nga)**(0.2) adjust the grid width to the number of galaxies
         deltarc=dd/nrc
         nrc1=-nrc/2
         nrc2=nrc/2
      else
         if (nrc.eq.-1) klfm=1
         deltarc=0.
         nrc1=1
         nrc2=1
      endif

      if (nrs.ge.1) then
         rsg=rsnew
c$$$         dd=grra*(500./nga)**(0.2) 
         dd=grra*(500./nga)**(0.2)   
         deltars=dd/nrs
         nrs1=-nrs/2
         nrs2=nrs/2
         kmfl=0
         klcdm=0
      else
         if (nrs.eq.-1) kmfl=1
         if (nrs.eq.-2) klcdm=1
         deltars=0.
         nrs1=1
         nrs2=1
      endif

      if (nbs.ge.1) then
         cbeg=cbnew
         if (kani.eq.0) then
            dd=(grra-0.5)*(500./nga)**(0.2)
            deltab=dd/nbs
         elseif (kani.eq.1) then
            dd=(grra+1.1)*(500./nga)**(0.2)
            deltab=dd/nbs
         elseif (kani.eq.2.or.kani.eq.3) then
            dd=(grra-0.1)*(500./nga)**(0.2)
            deltab=dd/nbs
         else
            dd=(grra-0.5)*(500./nga)**(0.2)
            deltab=dd/nbs
         endif
         nbs1=-nbs/2
         nbs2=nbs/2
      elseif (nbs.eq.-1) then
         kaml=1
         nbs1=1
         nbs2=1
      elseif (nbs.eq.-2) then
         khm=1
         nbs1=1
         nbs2=1
      else
         deltab=0.
         nbs1=1
         nbs2=1
      endif

      if (ngam.ge.1) then
         gamg=gamnew
c$$$         dd=grra*(500./nga)**(0.2)  
c$$$         deltagam=dd/ngam
c$$$         ngam1=-ngam/2
c$$$         ngam2=ngam/2
         ngam1=1
         ngam2=ngam
         gammax=3.
         gammin=0.
         if (gammax.lt.2.*gamnew) gammax=2.*gamnew
         deltagam=(gammax-gammin)/ngam
      else
         deltagam=0.
         ngam1=1
         ngam2=1
      endif

      if (nml.ge.1) then
         xmstarlg=xmstarlnew
         dd=grra*(500./nga)**(0.2)  
         deltaml=dd/nml
         nml1=-nml/2
         nml2=nml/2
      else
         deltaml=0.
         nml1=1
         nml2=1
      endif

      nfv=2
      r200g=r200new

      if (nr200.le.1.and.nrc.le.1.and.nrs.le.1.and.
     & nbs.le.1.and.nml.le.1) goto 891

 727  continue

      kr200=kr200+1
      if (nr200.ge.2) then
         r200=10.**(dlog10(r200g)+kr200*deltarv)
      endif

c     Determine m200, v200 and c appropriate for the m200 found
c     adopting the relation of Maccio, Dutton, van den Bosch 08
c     for relaxed halos

      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
      rm200=100.*hz*hz/grav*r200**3
cc      cduffy=5.78*(rm200/2.e12)**(-0.089)*1.1**(-0.52)
cc      cgao=10.**(-0.138*dlog10(rm200*0.7)+2.646)
cc      cmean=(cduffy+cgao)/2.
      cmean=6.76*(rm200/1.e12)**(-0.098)
      v200=10.*hz*r200

      write(*,239) za,va,r200,v200
 239  format(//' New run with  <z>,<v>,r200,v200 = ',f7.3,2x,f7.0,2x,
     &     2x,f6.3,2x,f5.0)

c     Best-fit to N(R) external to MAMPOSSt, if required

      if (nrc.eq.-2) then

         write(*,*) ' Fit to N(R) up to R/r200=',rup/r200

         write(*,293) rc,r200/rc
 293     format(' Best-fit N(R) scale parameter= ',f7.3,' (c=',f5.2,')')

      endif

c     Search on a grid

      do kgam=ngam1,ngam2
c$$$         gam=10.**(dlog10(gamg)+kgam*deltagam)
         gam=(kgam-1)*deltagam
         if (ngam1.eq.ngam2) gam=gamg

         do kml=nml1,nml2
            xmstarl=10.**(dlog10(xmstarlg)+kml*deltaml)

            do kb=nbs1,nbs2
               cbe=10.**(dlog10(cbeg)+kb*deltab)
            
               do ictr=nrc1,nrc2
                  rc=10.**(dlog10(rcg)+ictr*deltarc)
               
                  do k=nrs1,nrs2
c     
c     Mass follows light
c
                     if (kmfl.eq.1) then
                        rs=rc
                     elseif (klcdm.eq.1) then
c     
c     concentration from c=c(M)
c
                        rs=r200/cmean
                     else
c     
c     r_s as a free parameter
c     
                        rs=10.**(dlog10(rsg)+k*deltars)
                     endif
c
c     a_ML forced = r_s
c
                     if (kaml.eq.1) cbe=rs

c     
c     Light follows Mass (TLM)
c
                     if (klfm.eq.1) then
                        rc=rs   ! assume N(R) and M(r) chosen with same model
                        if (kmp.eq.1.and.knfit.eq.2) rc=2.*rs ! NFW M(r), Her N(R)
                        if (kmp.eq.2.and.knfit.eq.1) rc=rs/2. ! Her M(r), NFW N(R)
                     endif
c     
cc     If considering the universal distrib of interlopers
cc     the surface density profile cannot be as concentrated
cc     as the mass density profile; we use the factor found
cc     in MBM for R<r200
c
cc               if (kintd.eq.1) rc=0.85*rs
                  
                     xfv(1)=rs
                     xfv(2)=cbe

                     call vmaxlik(nfv,xfv,fml)

cc     If r200 changes, the number of objects selected changes
cc     as well, and one must scale the likelihood accordingly
c     [Obsolete, since we select galaxies within a radial range in Mpc
c      and not in units of r200 that may change]
cc               fml=fml*float(ngaref)/float(nga)
               
                     write(iu60,673) r200,rc,rs,cbe,gam,xmstarl,fml,kani
 673                 format(7(f13.3,2x),i2)
                     
                     if (fml.lt.fmlmingrid) then
                        fmlmingrid=fml
                        r200mingrid=r200
                        rsmingrid=rs
                        rcmingrid=rc
                        cbmingrid=cbe
                        gamgrid=gam
                        xmstarlgrid=xmstarl
                        xfn1best=rc
                        ngamin=nga
                        rlowmin=rlow
                        rupmin=rup
                        rcmin=rc
                        do j=1,nga
                           rso(j)=r(j)
                           vso(j)=v(j)
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      if (kr200.lt.nr200/2) goto 727

      write(*,*) ' number of galaxies is ',nga

      fmlmin=fmlmingrid
      r200min=r200mingrid
      rsmin=rsmingrid
      rcmin=rcmingrid
      cbmin=cbmingrid
      gammin=gamgrid
      xmstarlmin=xmstarlgrid
      rc=xfn1best
      nga=ngamin
      rlow=rlowmin
      rup=rupmin
      write(iu60,673) r200min,rcmin,rsmin,cbmin,gammin,
     &     xmstarlmin,fmlmin,kani

 891  continue

      if (kopt.lt.0) then          ! optimization skipped
         r200new=r200min
         rcnew=r200min
         rsnw=rsmin
         cbnew=cbmin
         gamnew=gammin
         xmstarlnew=xmstarlmin
         f=fmlmin
      endif
      write(iu60,673) r200new,rcnew,rsnew,cbnew,gamnew,xmstarlnew,f,kani

      r200=r200new
      cbe=cbnew
      rs=rsnew
      rc=rcnew
      gam=gamnew
      xmstarl=xmstarlnew

      write(*,292) r200,rc,r200/rc,rs,r200/rs,gam,cbe,xmstarl
 292  format(/' Best-fit ',
     &       /'   r_200    = ',f6.3,
     &       /'   r_tracer = ',f6.3,' (c_tracer= ',f7.2,')'
     &       /'   r_mass   = ',f6.3,' (c_mass=   ',f7.2,')'
     &       /'   gamma    = ',f6.3,
     &       /'   Anisotropy parameter = ',f8.4,
     &       /'   M*/L = ',f6.3,/)


c     output results for N(R) in a file

 737  continue
      sigma0=0.
      write(iu30,330) knfit,sigma0,rc,al
 330  format(i2,1x,e11.4,1x,e13.5,1x,f6.3)

c     output file of probs
      
      xfv(1)=rs
      xfv(2)=cbe

      call vmaxlik(nfv,xfv,fml)

      do jp=1,npv
         write(iu20,*) rpv(jp),vpv(jp),pv(jp)
      enddo
      close(iu20)

      r200=r200new
      rc=rcnew
      rs=rsnew
      cbe=cbnew
      gam=gamnew
      xmstarl=xmstarlnew

 777  continue
      return
      end

c
c     compute the difference in 'projected mass'
c     at the extreme of the radial range, for
c     3 different Sigma(R) profiles
c     
      subroutine sigmarnorm(x,fnorm)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      data pig /3.1415926535897932d0/


      if (knfit.eq.3) then

c     beta model

      if (al.eq.-1.) then
         fnup=dlog(x*x+rup*rup)
         if (rlow.gt.0) then
               fnlow=dlog(x*x+rlow*rlow)
            else
               fnlow=0.
            endif
         else
            fnup=(1.+(rup/x)**2)**(1.+al)/(1.+al)
            if (rlow.gt.0) then
               fnlow=(1.+(rlow/x)**2)**(1.+al)/(1.+al)
            else
               fnlow=0.
            endif
         endif
         fnorm=(fnup-fnlow)*pig*x*x

c     NFW

      elseif ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then

         c=r200/x
         gc=1./(dlog(c+1.)-c/(c+1.))
         cx=c*rup/r200
         uu=1./cx

cc         write(*,*) ' c_tracer is ',c

         if (uu.lt.1) then
            cm1=dacos(uu)
            fnup=cm1/dsqrt(dabs(cx*cx-1.))+dlog(cx/2.)
         elseif (uu.gt.1) then
            cm1=dacosh(uu)
            fnup=cm1/dsqrt(dabs(cx*cx-1.))+dlog(cx/2.)
         else
            fnup=1.+dlog(cx/2.)
         endif

         if (rlow.gt.0) then
            cx=c*rlow/r200
            uu=1./cx
            if (uu.lt.1) then
               cm1=dacos(uu)
               fnlow=cm1/dsqrt(dabs(cx*cx-1.))+dlog(cx/2.)
            elseif (uu.gt.1) then
               cm1=dacosh(uu)
               fnlow=cm1/dsqrt(dabs(cx*cx-1.))+dlog(cx/2.)
            else
               fnlow=1.+dlog(cx/2.)
            endif
         else
            fnlow=0.
         endif
cccc         fnorm=(fnup-fnlow)*gc/(r200*r200)
         fnorm=(fnup-fnlow)*gc
         
      else

c     Hernquist

         s=rup/x
         if (s.gt.1) then
            xs=dacos(1./s)/dsqrt(s*s-1.)
            fnup=s*s*(xs-1.)/(1.-s*s)
         elseif (s.lt.1) then
            xs=dlog((1.+dsqrt(1.-s*s))/s)/dsqrt(1.-s*s)
            fnup=s*s*(xs-1.)/(1.-s*s)
         else
            fnup=1./3.
         endif

         if (rlow.gt.0) then
            s=rlow/x
            if (s.gt.1) then
               xs=dacos(1./s)/dsqrt(s*s-1.)
               fnlow=s*s*(xs-1.)/(1.-s*s)
            elseif (s.lt.1) then
               xs=dlog((1.+dsqrt(1.-s*s))/s)/dsqrt(1.-s*s)
               fnlow=s*s*(xs-1.)/(1.-s*s)
            else
               fnlow=1./3.
            endif
         else
            fnlow=0.
         endif

         fnorm=fnup-fnlow

      endif

cc      write(*,*) ' rlow, rup ',rlow,rup

      return
      end


c
c     temp
c
      subroutine temp(t1,sigrt1,sigrt2,sigzt1t2,bec)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      real*8 rjl(1), yrjl(1)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'tsallis.i'
c
      t2 = 2.*t1
c
c     Interpolate to get the sigma_r(t)
c
      rjl(1)=dlog(t1)
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      if (ier .eq. 34) then
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
     &    rjl(1)
      endif
      sigrt1=dexp(yrjl(1))

      rjl(1)=dlog(t2)
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      if (ier .eq. 34) then
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
     &     rjl(1)
      endif
      sigrt2=dexp(yrjl(1))

c
c     beta(r)
c
c     Convert from cbe=beta' to bec=beta
c     in the case of constant anisotropy
c
      bec=1.-1./(cbe*cbe)
      if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
c
c     sigma_z(R,r)
c
      sigzt1t2=sigrt1*dsqrt(1.-bec*t1*t1/(t2*t2))

      return
      end

c
c     compute the projected mass within the virial sphere
c     for NFW model (App. B2 in MBM)
c     
      subroutine virsphmassp(rrr,vsmp)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      data pig /3.1415926535897932d0/

      c=r200/rc
      x=rrr/rc

      arglog=(c+1.)*(c-sqrt(c*c-x*x))/x
      argcos=(c+x*x)/((c+1.)*x)

      if (x.le.0.) then

         vsmp=0.

      elseif (x.gt.0.and.x.lt.1.and.x.lt.c) then

         vsmp=(dsqrt(c*c-x*x)-c)/(c+1.)+dlog(arglog)+
     +        dacosh(argcos)/dsqrt(1.-x*x)

      elseif (x.gt.1.and.x.lt.c) then

         vsmp=(dsqrt(c*c-x*x)-c)/(c+1.)+dlog(arglog)+
     +        dacos(argcos)/sqrt(x*x-1.)

      elseif (x.eq.1.and.x.lt.c) then

         vsmp=dlog((c+1.)*(c-dsqrt(c*c-1.)))-c/(c+1.)+
     +        2.*dsqrt((c-1.)/(c+1.))

      else

         vsmp=dlog(c+1.)-c/(c+1.)

      endif

      vsmp=vsmp/(alog(2.)-0.5)
      
      return
      end

c
c
c
      function sdint(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      include 'datarv.i'
      data pig /3.1415926535897932d0/

      aint=10.**(-1.061+0.364*x*x-0.580*x**4+0.533*x**6)
      sigmaint=0.612-0.0653*x*x
      bint=0.0075
      xkappa=4.

      argerf=xkappa/dsqrt(2.d0)/sigmaint

      sdint=x*(dsqrt(pig/2.d0)*aint*sigmaint*erf(argerf)+xkappa*bint)

      return
      end
c
c     This function evaluates the incomplete beta function
c     that enters the 1st of eqs. A16 in Mamon & Lokas,
c     using a recursive relation, see addendum in
c     ML 2006 MNRAS, 370, 1582
c
      function betairec(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptbcg.i'
      external gammarec

      bec=1.-1./(cbe*cbe)
      if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
      a=bec+0.5d0
      b=0.5d0

      if (a.gt.0.d0) then
         bir=dbetai(x,a,b)
      else
         nit=nint(0.5d0-a)
         birgj=0.d0
         do j=0,nit-1
            birgj=gammarec(a+b+j)/gammarec(a+1.d0+j)/gammarec(b)*
     *           x**(a+j)*(1.d0-x)**b+birgj
         enddo
         bir=birgj+dbetai(x,a+nit,b)
      endif

      betairec=bir
      return
      end

c
c     Function gammarec: uses recursive relation
c                Gamma(x)=Gamma(x+1)/x
c     when x<0 to call dgamma(arg) with arg>0
c

      function gammarec(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      nit=-nint(x-0.5)

      if (nit.le.0) then

         gammarec=dgamma(x)
         return

      else

         prod=1.
         do i=1,nit
            prod=prod*(x+i-1.)
         enddo
         gammarec=dgamma(x+nit)/prod
         return

      endif

      end


c
c     Subroutine for the computation of Max Lik
c     with BOBYQA or NEWUOA
c     
      subroutine calfun(n,x,f)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension x(n),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),wmbm(25000),embm(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'omegas.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      call freepareva(n,x,r200,rc,rs,cbe,gam,xmstarl)

ccc      write(*,429) r200,rc,rs,cbe,al
ccc 429  format(' In calfun: r200,rc,rs,cbe,al = ',5(f6.3,2x))

      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
      v200=10.*hz*r200

      if (dabs(al).lt.0.0001) n1=1
      
      ngambm=nga
      do j=1,nga
         rmbm(j)=r(j)
         vmbm(j)=v(j)
         embm(j)=e(j)
         wmbm(j)=w(j)
      enddo

      irule=2
      errabs=0.d0
      errrel=0.001d0
      rismin=1.d-190
      rismax=1.d190
      rinfinity=20.0d0

      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin
         risl = dcadre(sr2int,xris(i),dlog(2.*rinfinity),errabs,errrel,
     &    errest,ier)
         if (ier .gt. 128) then
            print *,'VMAXLIK-SR2INT: rmin=',xx1,' rmax=',xx2
            stop
         endif
         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)
      enddo

c
c compute spline coeffs for later interpolation of sigma_r
c
      call icsccu(xris,yris,ninterp,csr,icsr,ier)
      if (ier .gt. 128) then
         print *,' in VMAXLIK ...'
         do i = 1, ninterp
            print *,' i xris yris = ', i, xris(i), yris(i)
         enddo
      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.d0
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      if (kintd.eq.1) then
         rinfinity=0.99999*r200
      endif

c     Now compute the distribution function at the radial distance of
c     each galaxy by interpolating the sigma_r(r) and then
c     evaluate the Max Lik, using weights if required

c     If required, make the estimate only on a grid of values
c     and then bispline interpolate in the R,v positions of the
c     galaxies (kbsp=1)


      if (kbsp.eq.1) then

         rmin=1.e12
         vmin=1.e12
         rmax=-1.e12
         vmax=-1.e12
         eveave=0.
         do j=1,ngambm
            if (rmbm(j).lt.rmin) rmin=rmbm(j)
            if (abs(vmbm(j)).lt.vmin) vmin=abs(vmbm(j))
            if (rmbm(j).gt.rmax) rmax=rmbm(j)
            if (abs(vmbm(j)).gt.vmax) vmax=abs(vmbm(j))
            eveave=eveave+embm(j)
         enddo
         eveave=eveave/ngambm

         nr=rknots
         nv=vzknots
         kspo=3
         npts=nr*nv
         stepr=(rmax-rmin)/(nr-1)
         stepv=(vmax-vmin)/(nv-1)
         rmi1=rmin*0.9
         rma1=rmax*1.1
         if (rma1.ge.rinfinity) rma1=0.9999d0*rinfinity
         do j=1,nr
            rp(j)=10.**(dlog10(rmi1)+(j-1.)*dlog10(rma1/rmi1)/(nr-1.))
         enddo
         do j=1,nv
            vp(j)=vmin+stepv*(j-1)
         enddo

         do j=1,nr
            rj=rp(j)
            do i=1,nv
               vj=vp(i)
               ej=eveave
               xmin=rj
               k=k+1
c
c     re-define rinfinity as the radius where v_z/sigma_r<0.1
c     to avoid a problem in v_z=0
c
               if (kintd.eq.0) call rinf(rinfinity)
c
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c
               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)

               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN: rmin=',xx1,' rmax=',xx2
                  stop
               endif
c     
c     Add interloper contribution if requested
c
               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
     &                 errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
cc                  write(*,*) ' Number of obsd galaxies ',ngambm
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
cc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
              endif
               
              if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


              if (nrc.gt.0.5) then
                 
                 call sigmarnorm(rc,fnorm)
                 if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                 p=2.*pig*rj*g/fnorm
                 
                 
              else
                 
                 if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                    sigmaden=sigmar1(rj)
                 elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                    sigmaden=sigmar2(rj)
                 else
                    sigmaden=sigmar3(rj)
                 endif
                 p=g/sigmaden
                 
              endif

c     The log(probability) is stored in z(j,i)

              z(j,i)=dlog(p)
            enddo
         enddo

c     Bi-spline 2-d interpolation to obtain the p values
c     of all data points from the p values of the sparse grid
     
         call ibcccu(z,rp,nr,vp,nv,c,ic,wk,ier)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp,ier)
            psum=psum-dlogp*wmbm(l)
 
         enddo

      else

c     Use all the original data, no interpolation

         do j=1,ngambm
            rj=rmbm(j)
            vj=vmbm(j)
            ej=embm(j)
            xmin=rj

c     re-define rinfinity as the radius where v_z/sigma_r leq 0.1
c     to avoid a problem in v_z=0; if we require to consider the 
c     interloper contribution, rinfinity is set to 1 already

            if (kintd.eq.0) call rinf(rinfinity)

c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)

               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
                  stop
               endif


c     Add interloper contribution if requested

               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4.+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  n2=2
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,
     &                 errrel,errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
                  g=g+gi
               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally
               
               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm
               
               else

                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                     sigmaden=sigmar1(rj)
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                     sigmaden=sigmar2(rj)
                  else
                     sigmaden=sigmar3(rj)
                  endif
                  p=g/sigmaden

               endif

               call temp(rj,sigrt1,sigrt2,sigzt1t2,bec)


               psum=psum-dlog(p)*wmbm(j)
               

            endif

         enddo

      endif
c
c     The factor nga/wsum is =1 if no weights, but is not 1
c     if there are weights, and is needed to take into account
c     the real number of points (otherwise, we would estimate
c     a Max Lik too low by an average factor of 1/Ngr, where
c     Ngr is the average number of galaxies per group, since
c     the likelihood we want to estimate is the Sum of the
c     ngambm galaxy probabilities)
c     

      f=psum/wsum*ngambm

      return
      end

c
c     function for the computation of Max Lik
c     with POWELL
c     
      function func(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension x(6),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),embm(25000),wmbm(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptbcg.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'free.i'
      include 'omegas.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      n=nfreefunc

      call freepareva(n,x,r200,rc,rs,cbe,gam,xmstarl)

cc      if (cbe.le.0) then
cc         write(*,*) x
cc         write(*,*) cbe
cc         stop
cc      endif

      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
      v200=10.*hz*r200

      if (dabs(al).lt.0.0001) n1=1
      
      ngambm=nga
      do j=1,nga
         rmbm(j)=r(j)
         vmbm(j)=v(j)
         embm(j)=e(j)
         wmbm(j)=w(j)
      enddo

      irule=2
      errabs=0.d0
      errrel=0.001d0
      rismin=1.d-190
      rismax=1.d190
      rinfinity=20.0d0

      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin
         risl = dcadre(sr2int,xris(i),dlog(2.*rinfinity),errabs,errrel,
     &    errest,ier)
         if (ier .gt. 128) then
            print *,'VMAXLIK-SR2INT: rmin=',xx1,' rmax=',xx2
            stop
         endif
         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)
      enddo

c
c compute spline coeffs for later interpolation of sigma_r
c
      call icsccu(xris,yris,ninterp,csr,icsr,ier)
      if (ier .gt. 128) then
         print *,' in VMAXLIK ...'
         do i = 1, ninterp
            print *,' i xris yris = ', i, xris(i), yris(i)
         enddo
      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.d0
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      if (kintd.eq.1) then
         rinfinity=0.99999*r200
      endif

c     Now compute the distribution function at the radial distance of
c     each galaxy by interpolating the sigma_r(r) and then
c     evaluate the Max Lik, using weights if required

c     If required, make the estimate only on a grid of values
c     and then bispline interpolate in the R,v positions of the
c     galaxies (kbsp=1)


      if (kbsp.eq.1) then

         rmin=1.e12
         vmin=1.e12
         rmax=-1.e12
         vmax=-1.e12
         do j=1,ngambm
            if (rmbm(j).lt.rmin) rmin=rmbm(j)
            if (abs(vmbm(j)).lt.vmin) vmin=abs(vmbm(j))
            if (rmbm(j).gt.rmax) rmax=rmbm(j)
            if (abs(vmbm(j)).gt.vmax) vmax=abs(vmbm(j))
         enddo

         nr=rknots
         nv=vzknots
         kspo=3
         npts=nr*nv
         stepr=(rmax-rmin)/(nr-1)
         stepv=(vmax-vmin)/(nv-1)
         rmi1=rmin*0.9
         rma1=rmax*1.1
         if (rma1.ge.rinfinity) rma1=0.9999d0*rinfinity
         do j=1,nr
            rp(j)=10.**(dlog10(rmi1)+(j-1.)*dlog10(rma1/rmi1)/(nr-1.))
         enddo
         do j=1,nv
            vp(j)=vmin+stepv*(j-1)
         enddo

         do j=1,nr
            rj=rp(j)
            do i=1,nv
               vj=vp(i)
               xmin=rj
               k=k+1
c
c     re-define rinfinity as the radius where v_z/sigma_r<0.1
c     to avoid a problem in v_z=0
c
               if (kintd.eq.0) call rinf(rinfinity)
c
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c
               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)

               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN: rmin=',xx1,' rmax=',xx2
                  stop
               endif
c     
c     Add interloper contribution if requested
c
               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
     &                 errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
cc                  write(*,*) ' Number of obsd galaxies ',ngambm
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
cc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
              endif
               
               if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm


               else

                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                     sigmaden=sigmar1(rj)
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                     sigmaden=sigmar2(rj)
                  else
                     sigmaden=sigmar3(rj)
                  endif
                  p=g/sigmaden

               endif

c     The log(probability) is stored in z(j,i)

               z(j,i)=dlog(p)
            enddo
         enddo

c     Bi-spline 2-d interpolation to obtain the p values
c     of all data points from the p values of the sparse grid
     
         call ibcccu(z,rp,nr,vp,nv,c,ic,wk,ier)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp,ier)
            psum=psum-dlogp*wmbm(l)
 
         enddo

      else

c     Use all the original data, no interpolation

         do j=1,ngambm
            rj=rmbm(j)
            vj=vmbm(j)
            xmin=rj

c     re-define rinfinity as the radius where v_z/sigma_r leq 0.1
c     to avoid a problem in v_z=0; if we require to consider the 
c     interloper contribution, rinfinity is set to 1 already

            if (kintd.eq.0) call rinf(rinfinity)

c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)

               umax = dacosh(rinfinity/rj)
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
               if (ier .gt. 128) then
                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
                  stop
               endif


c     Add interloper contribution if requested

               if (kintd.eq.1) then
                  rjrv=rj/r200
                  vzvv=vj/v200
                  ctr=r200/rc
                  aint=10.**(-1.061+0.364*rjrv*rjrv-
     -                 0.580*rjrv**4.+0.533*rjrv**6)
                  sigmaint=0.612-0.0653*rjrv*rjrv
                  bint=0.0075
                  n2=2
                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
                  call virsphmassp(rup,xmprmax)
                  call virsphmassp(rlow,xmprmin)
                  xmin=rlow/r200
                  xmax=rup/r200
                  errabs=0.d0
                  errrel=0.001d0
                  xinteg=dcadre(sdint,xmin,xmax,errabs,
     &                 errrel,errest,ier)
                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
                  xnh=xnv*(xmprmax-xmprmin)/xmctr
                  xni=xnv*2.*pig*xinteg
                  gi=xnv/(v200*r200*r200)*gihat
                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
                  g=g+gi
               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally

               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm
               
               else

                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then
                     sigmaden=sigmar1(rj)
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
                     sigmaden=sigmar2(rj)
                  else
                     sigmaden=sigmar3(rj)
                  endif
                  p=g/sigmaden

               endif

               call temp(rj,sigrt1,sigrt2,sigzt1t2,bec)


               psum=psum-dlog(p)*wmbm(j)
               

            endif

         enddo

      endif
c
c     The factor nga/wsum is =1 if no weights, but is not 1
c     if there are weights, and is needed to take into account
c     the real number of points (otherwise, we would estimate
c     a Max Lik too low by an average factor of 1/Ngr, where
c     Ngr is the average number of galaxies per group, since
c     the likelihood we want to estimate is the Sum of the
c     ngambm galaxy probabilities)
c     

      func=psum/wsum*ngambm

      return
      end


c
c

      subroutine freepareva(nfreepar,freepar,r200new,rcnew,rsnew,cbnew,
     &    gamnew,xmstarlnew)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
       parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,clight=2.99792458d5)
      dimension freepar(nfreepar)
      include 'paramsoptbcg.i'
      include 'omegas.i'
c     These two cosmo params are now included with omegas.i
cc      omegal=0.7
cc      omega0=0.3
      hz=h0*sqrt(omega0*(1.+za)**3+omegal)

      r200g=r200

      ifp=0
      if (ipar(1).eq.0) then
         r200new=r200g
      else
         ifp=ifp+1
         r200new=10.**freepar(ifp)
      endif
      if (ipar(2).eq.0) then
         rcnew=rcg
      else
         ifp=ifp+1
         rcnew=10.**freepar(ifp)
      endif
      if (ipar(3).eq.0) then
         rsnew=rsg
      else
         ifp=ifp+1
         rsnew=10.**freepar(ifp)
      endif
      if (ipar(4).eq.0) then
         cbnew=cbeg
      else
         ifp=ifp+1
         cbnew=10.**freepar(ifp)
      endif
      if (ipar(5).eq.0) then
         gamnew=gamg
      else
         ifp=ifp+1
         gamnew=10.**freepar(ifp)
      endif
      if (ipar(6).eq.0) then
         xmstarlnew=xmstarlg
      else
         ifp=ifp+1
         xmstarlnew=10.**freepar(ifp)
      endif

      if (klcdm.eq.1) then      ! LCDM c=c(M)
         rm200=100.*hz*hz/grav*r200new**3
         cmean=6.76*(rm200/1.e12)**(-0.098)
         rsnew=r200new/cmean
      endif

      if (klfm.eq.1) then ! LfM (TLM in Mamon et al. 2013)
         rcnew=rsnew   ! assume N(R) and M(r) chosen with same model
         if (kmp.eq.1.and.knfit.eq.2) rcnew=2.*rsnew ! NFW M(r), Her N(R)
         if (kmp.eq.2.and.knfit.eq.1) rcnew=rsnew/2. ! Her M(r), NFW N(R)
      endif
      if (kmfl.eq.1) rsnew=rcnew  ! MfL
      if (kaml.eq.1) cbnew=rsnew  ! a_ML = r_s

      return
      end

