!******************** MG-MAMPOSSt routines and functions ***********************   
! made available in the libMAM.a library 
! when compiling, link all the libraries, stored in the build folder 
! e.g.:
! f95 -o testMAM.e testMAM.f -L build/ -lMAM -L build/GamI/ -lGamI -L build/Newuoa/ -lNewuoa -L build/JJin/ -lJJin -L build/Utili/ -lUtili -L build/Powell/ -lPowell
! **********************************************************************
! **********************************************************************

!>   \addtogroup MAM MG-MAMPOSSt main ingredients
!>   @details List of the main functions and subroutines needed for the execution
!!   of MG-MAMPOSSt. 
!>   This block includes all the fundamental routines and functions on which the 
!!   MG-MAMPOSSt method is based. The central core is represented by the
!!   mamposst() subroutine, which computes the tabulated \f$-\ln \mathcal{L}\f$.
!>
!>   All the routines in this block need external inputs as described below.
!!   In particular, the COMMON block of parameters "paramsoptS.i" should always be included.


!>  \ingroup MAM
      subroutine mamposst(di,ve,eve,rso,vso,npg1)
!> @authors A. Biviano, G. Mamon, G. Boue', L. Pizzuti
     
!>    @details MAMPOSSt subroutine computes the Log-likelihood in the projected
!!     phase space. 
!>     This is the core of the MG-MAMPOSSt code. Along with
!!     the input data, mamposst() requires as external
!!     input the values of all the input parameters/Options of 
!!     MG-MAMPOSSt 
!!     along with the file names. In particular:
!>     - name and id number of the output likelihood file: `iu60`
!>     - name and id number of the output file of projected position  
!!       and velocity with the associated probability p(R_i,v_i): `iu20`
!>     - name and id number of the output file containing the best fit of
!!      projected number density profile: `iu30`
!>     - input parameters in `test_pars.txt` described in 
!!       the documentation and named as in the COMMON block `paramsoptS.i`
!>     - The `Options.txt` file, read within the subroutine 
!>
!>  In the case of `Nclust>1` 
!! in `Options.txt`, the 
!! mamposst() subroutine requires the additional data-files:
!> `data/datphys_<i>.dat`,
!! where "i" is an integer number larger than 1.
 
!>    
!>     @param[in] di(npg1)  dimension REAL*8   projected distance in Mpc
!>     @param[in] ve(npg1)  dimension REAL*8   line-of-sight velocity [km/s]
!>     @param[in] eve(npg1) dimension REAL*8   velocity error [km/s]
!>     @param[in] npg1      INTEGER*4          number of galaxies in the 
!!     phase space    
!> 
!>     @param[out] rso(npg1) dimension REAL*8   projected distance in Mpc of the 
!!     galaxies considered in the fit  
!>     @param[out] vso(npg1)  dimension REAL*8  line-of-sight velocity [km/s] of
!!     the galaxies considered in the fit 

!>     Example of usage, assuming to have read data and all the input parameters
!!     in pars_test.txt and Options.txt:
!>
!>   \code{.f}
!>          
!>          open(20,file="rnvn.dat",status='unknown')
!>          open(30,file="Npar.dat",status='unknown')
!>          open(60,file="MaxLik.dat",status='unknown')
!>          iu20=20
!>          iu30=30
!>          iu60=60
!>          call mamposst(di,ve,eve,rso,vso,npg1)
!>          close(30)
!>          close(60) 
!>    \endcode
!>    Note that the file "iu20" is closed within the mamposst subroutine.

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      logical*4 eval
      parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,q0=-0.64,clight=2.99792458d5, npar=11)
      dimension di(npg1),ve(npg1),eve(npg1),rso(npg1),vso(npg1),
     &     xfv(2), sigma(npar) 
      dimension dpp(25000),vpp(25000),epp(25000),wpp(25000), 
     &     did(25000),viv(25000),eie(25000),wiw(25000)
      dimension freepar(npar),freelow(npar),freeup(npar),wpar(5000),
     &     xi(npar,npar)
      CHARACTER(len= 1) ARR1
      character(len= 2) arr2
      character(len = 13) :: filename
      character(len = 16) :: filebcg
      character(len=4) :: dat
      character(len=300) :: buffer, label, label2
      integer :: pos, pos1, pos2, posl1, posl2, posl, poscom
      integer, parameter :: fh = 15
      integer :: ios = 0
      integer :: line = 0
      logical is_cham, use_u

      include 'paramsoptS.i'
      include 'datavdp.i'
      include 'datarv.i'
      include 'units.i'
      include 'free.i'
      include 'probs.i'
      include 'barions.i'

      external fcn1,fcn2,fcn3,fa
      external Plens
      external sr2int,sr2out,sigmar1,sigmar2,sigmar3,gwenu

      print *,' Entering MAMPOSSt subroutine'
      
c     Opening log file for keeping track of all the options
      open(11,file="Run_info.txt",status='unknown')
      write(11,*) 'Information about the current run'
      write(11,*) ''      
      
      
c     Open files containing the gas and galaxy mass parameters      
      if (kbary.eq.1) then
       open(92, file="GasGal_mass.txt", status='unknown')
       write(92,'(4(a13,2x))') 'rhogas', 'rsgas', 'rhostar', 'rstar'
      endif
    
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
      write(*,932) nga,rlow,rup
      write(11,932) nga,rlow,rup
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
         write(11,195) rc
         
 195     format(' Best-fit to N(R) outside MAMPOSSt: r_tr=',f6.3)

      endif

C******************
      if (kbcg.eq.1.and.kbdt.eq.1) then
            !read data file of BCG
       !filebcg='data/vel_bcg.dat'
       !open (12, FILE=filebcg, STATUS='OLD')
        jj=0
        j0=-1
        dimin=1.e12
 406    continue
                
        read(81,*,end=176) dmpc,dvr,dev !12
        jj=jj+1
        rbcg(jj)=dmpc
        sbcg(jj)=dvr
        esbcg(jj)=dev
       write(*,*) rbcg(jj), sbcg(jj), esbcg(jj)
        goto 406
 176    continue
        close(81) !12
        nbcg=jj
        
       endif
       
       
      if(kbary.eq.1) then
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
!      do k = 1, 10000
!       x = 0.001+0.01*k
!       write(*,*) fbars(x), fbars2(x)
!      enddo
      endif

c     switches to control the grid search *****************	            
      open(fh, file='Options.txt')

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

      write(*,*) ''
      write(*,*) '********************************'
      write(*,*) ' Now reading Options'
      write(*,*) '********************************'
      write(*,*) ''
      do while (ios == 0)
         read(fh, '(A)', iostat=ios) buffer
         if (ios == 0) then
            line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        
        !devo togliere gli spazi all'inizio
            ic=1
            do while(scan(buffer(ic:ic),' ').ne.0) 
            !aumenta la posizione se trovi uno spazio
             ic=ic+1
  
            enddo
            buffer=buffer(ic:)

        
            pos1 = scan(buffer, '    ')
            pos2 = scan(buffer, '=')
            if (pos2.ge.pos1) then
               pos=pos2
            else
               pos=pos1
            endif         
         ! split equality sign and space.
            posl1=scan(buffer, '    ')
            posl2=scan(buffer, '=')
            if (posl1.ne.0.and.posl2.ne.0) then 
             posl=min(posl1,posl2)
             label = buffer(1:posl-1)
            else if (posl1.ne.0.and.posl2.eq.0) then 
             label = buffer(1:posl1)
             
            else if (posl1.eq.0.and.posl2.ne.0) then 
             label = buffer(1:posl2)
            endif 
            
            buffer = buffer(pos+1:)
            poscom= scan(buffer, '!')
            if (poscom.ne.0) buffer= buffer(:poscom-1)
            
            !devo togliere gli spazi all'inizio
c            ic=1
c            do while(scan(label(ic:ic),' ').ne.0) 
            !aumenta la posizione se trovi uno spazio
c             ic=ic+1
  
c            enddo
    
c            label2=label(ic:)

            select case (label)
                
                !write(*,"('  Read rst:     ',f10.2)")  rst
            case ('delta1')
                read(buffer, *, iostat=ios) wr2
                if (ios.eq.-1) then
                    ios=0
                    wr2=0.1
                    if (kmp.eq.8) wr2=0.3 
                endif
                !write(*,"('  Read delta1:   ',f10.2)")  wr2
            case ('delta2')
                read(buffer, *, iostat=ios) wr
                if (ios.eq.-1) then
                    ios=0
                    wr=0.3
                    if (kmp.eq.8) wr=0.005 
                endif
                !write(*,"('  Read delta2:   ',f10.2)")  wr
            case ('delta3')
                read(buffer, *, iostat=ios) wx
                if (ios.eq.-1) then
                    ios=0
                    wx=0.5
                    if (kmp.eq.8) wx=30 
                endif
                !write(*,"('  Read delta3:   ',f10.2)")  wx
                
!         parameter space limits                
            case ('kpro')
                read(buffer, *, iostat=ios) kpro
                if (ios.eq.-1) then
                    ios=0
                    kpro=1
                endif    
                !write(*,"('  Read kpro:   ',i10)")  kpro
                if (kpro.gt.1.or.kpro.lt.0) kpro=1  
                
!         joint analysis 
            case ('Nclust')
                read(buffer, *, iostat=ios) Nclust
                if (ios.eq.-1) then
                    ios=0
                    Nclust=1
                endif                 
                !write(*,"('  Read Nclust:   ',i10)") Nclust
                if(Nclust.gt.1.and.kbcg.eq.1) then
                 write(*,*) '**** WARNING *********'
                 write(*,*) 'no BCG data for more than one cluster'
                 write(*,*) ' setting kbcg=0'
                 kbcg=0
                endif
                
                
            case ('nsingle')
                read(buffer, *, iostat=ios) nsingle
                if (ios.eq.-1) then
                    ios=0
                    nsingle=0
                endif                 
                !write(*,"('  Read nsingle:   ',i10)") nsingle
                if (nsingle.gt.1.or.nsingle.lt.0) nsingle=0
            
!         additional options                
            case ('istop')            
                read(buffer, *, iostat=ios) istop
                if (ios.eq.-1) then
                    ios=0
                    istop=0
                endif                
                !write(*,"('  Read istop:   ',i10)") istop
                if (istop.gt.1.or.istop.lt.0) istop=0
                
                
            case ('teps')            
                read(buffer, *, iostat=ios) teps
                if (ios.eq.-1) then
                    ios=0
                    teps(1)=0.1
                    teps(2)=0.1
                    teps(3)=0.1
                    teps(4)=0.1
                    teps(5)=0.1
                    teps(6)=0.1
                    teps(7)=0.1
                    teps(8)=0.1
                    teps(9)=0.1
                    teps(10)=0.1
                    teps(11)=0.1
                endif 
                !write(*,"('  Read teps:   ',7(f10.2))") teps
            case ('delik')            
                read(buffer, *, iostat=ios) delik
                if (ios.eq.-1) then
                    ios=0
                    delik=0.2
                endif 
                !write(*,"('  Read delik:   ',f10.3)") delik
            case ('nskip')
                read(buffer, *, iostat=ios) nskip
                if (ios.eq.-1) then
                    ios=0
                    nskip=0
                endif
                !write(*,"('  Read nskip:   ',i10)") nskip
                if (nskip.gt.1.or.nskip.lt.0) nskip=0
            case ('nres')
                read(buffer, *, iostat=ios) nres
                if (ios.eq.-1) then
                    ios=0
                    nres=0
                endif
                !write(*,"('  Read nres:   ',i10)") nres
                if (nres.gt.1.or.nres.lt.0) nres=0                
            case ('nequ')
                read(buffer, *, iostat=ios) nequ
                if (ios.eq.-1) then
                    ios=0
                    nequ=0
                endif
                !write(*,"('  Read nequ:   ',i10)") nequ
                if (nequ.gt.1.or.nequ.lt.0) nequ=0  
            case ('nsame')
                read(buffer, *, iostat=ios) nsame
                if (ios.eq.-1) then
                    ios=0
                    nsame=0
                endif
                !write(*,"('  Read nsame:   ',i10)") nsame 
                if (nsame.gt.1.or.nsame.lt.0) nsame=0   
                if (nmcmc.eq.0) then
                 nsame=0
                 write(*,*) 'nsame option not valid in grid mode' 
                endif            
            case default
               ! print *, 'Skipping invalid label at line', line
            end select
         end if
      end do      

      nplens=10 !number of points in the lensing likelihood
c     ******************************************************

      write(*,601) nmcmc, nres,nequ,istop,nskip,nlens,Nclust,Nsample
      write(11,601) nmcmc,nres,nequ,istop,nskip,nlens,Nclust,Nsample
 601  format(/,' Switches selected for parameter exploration: ',/,
     &         ' grid or MCMC (0/1)=',i1,/,
     &         ' finer grid= ',i1,/,
     &         ' equal grid values = ',i1,/,
     &     ' stop if desired = ',i1,/,
     &     ' skip the optimization results =  ',i1,/,
     &     ' add lensing probability distribution =  ',i1,/,
     &     ' Number of halos in the analysis =  ',i2,/,
     &     ' Number of sampling points in MCMC = ',i7,/)

c     Loop over r200, if requested, searching in the range
c     of 0.5 to 1.8 the real r200
      
      if (nr200.ge.2) then 
        dd=0.21*(500./nga)**(0.2) 
        if(nres.eq.1) dd=0.1*(500./nga)**(0.2) !finer grid
        if (nequ.eq.1) then !equal spaced grid
          dd=0.23
          if(nres.eq.1) dd=0.08
        endif
        deltarv=dd/nr200
        kr200=-nr200/2
 
      else
         kr200=0
      endif

      fmlmingrid=1.e20

      r200g=r200

      write(*,*) ' '
      write(*,*) ' ... now running MAMPOSSt procedure ... '

c     Now run MAMPOSSt

      kboby=0   ! used to keep track of minimization algorithms already tried
      knewu=0 
      kpowe=0

 174  continue

      ktried=kboby+knewu+kpowe

c     Using NEWUOA or BOBYQA or POWELL minimization 
c     The free parameters are 6 in this version but can be less than
c     6 according to the choice of the user

      ipar(1)=1  ! identify parameters that are to be minimized
      ipar(2)=1
      ipar(3)=1
      ipar(4)=1
      ipar(5)=1
      ipar(6)=1
      ipar(7)=1
      ipar(8)=1
      ipar(9)=1
      ipar(10)=1
      ipar(11)=1

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
         freepar(ip)=(facguess*cbeg)
      else
         ipar(4)=0
      endif
      if (ntmass.gt.0) then   ! modified grav parameter
         ip=ip+1   
         if (kmp.ne.8) freepar(ip)=dlog10(facguess*tmassg)
         if (kmody.eq.1.and.kmp.ne.8.or.kmp.eq.17) then
         ! in the case of CS 
            freepar(ip)=(1-dexp(-facguess*tmassg*0.1))
         endif
	     if (kmp.eq.8) then
           freepar(ip)=(facguess*(tmassg))
           if (tmassg.eq.0.0) freepar(ip)=(facguess*(tmassg+1.e-4)) 
         endif !For BH gravity Y1 can be negative
      else
         ipar(5)=0
      endif
      
      if (nb2.gt.0) then     ! tracer velocity anisotropy parameter
         ip=ip+1
         freepar(ip)=(facguess*cbe0g)
      else
         ipar(7)=0
      endif
      
      if (njaf.gt.0) then     ! scale radius jaffe profile BCG
         ip=ip+1
         freepar(ip)=dlog10(facguess*rjafg)
      else
         ipar(8)=0
      endif
      
      if (nxmas.gt.0) then     ! m/L jaff profile BCG
         ip=ip+1
         freepar(ip)=dlog10(facguess*xmstarg)
         
      else
         ipar(9)=0
      endif
      
      if (nA3.gt.0) then    !dA3 MG parameter
       ip=ip+1
       freepar(ip)=dlog10(facguess*dA3g)
      else
         ipar(10)=0
      endif    
      
      if (nrbeta.gt.0) then    !rbeta parameter
       ip=ip+1
       freepar(ip)=dlog10(facguess*rbetag)
      else
         ipar(11)=0
      endif            
       
      if (nhone.gt.0) then   ! screening parameter
         ip=ip+1
         freepar(ip)=dlog10(facguess*screeg)
         if (kmody.eq.1.and.kmp.ne.8.or.kmp.eq.17) then
         !if(kmp.eq.9.or.kmp.eq.11.or.kmp.eq.12.or.kmp.eq.15) then
          freepar(ip)=facguess*screeg/(facguess*screeg+1)
         endif
!         if (kmp.eq.8) then
!           freepar(ip)=(facguess*(screeg))
!           if (screeg.eq.0.0) freepar(ip)=(facguess*(screeg+1.e-2)) !TEST2
!         endif !For BH gravity Y2 can be negative         
         
      else
         ipar(6)=0
      endif
     
      
      nfreepar=ip
      write(11,*) 'Running MG-MAMPOSSt with the selected mass profile:'
      select case(kmp)
       case(1)
	     write(11,*) 'NFW mass model'
       case(2)
	     write(11,*) 'Hernquist mass model'
       case(3)
	     write(11,*) 'PIEMD mass model'
       case(4)
	     write(11,*) 'Burkert mass model'
       case(5)
	     write(11,*) 'Soft Isothermal mass model'
       case(6)
	     write(11,*) 'Einasto mass model with n=5'
!       case(7)
!	     write(11,*) 'mNFW (linear Hordenski) mass model'
       case(8)
	     write(11,*) 'mNFW DHOST mass model'
       case(9)
	     write(11,*) 'mNFW chameleon screening mass model'
       case(10)
	     write(11,*) 'gNFW mass model'
       case(11)
         write(11,*) 'mBurkert chameleon screening mass model'
       case(12)
         write(11,*) 'mEinasto chamelon screening mass model'
       case(13)
         write(11,*) 'Refracted garvity mass model'
       case(15)
         write(11,*) 'mIsothermal chameleon screening mass model'
       case(17)
         write(11,*) 'multi-component chameleon screening mass model' 
       case(18)
         write(11,*) 'multi-component clustering dark energy model' 
      end select
      write(11,*) ' '
      if(kbcg.eq.1.and.kbdt.eq.1) then
       write(11,*) 'BCG data added to the fit'
      endif

      write(11,622) nfreepar,ipar,r200g,rcg,rsg,cbeg,tmassg,
     &             screeg,cbe0g,rjafg,xmstarg,dA3g,rbetag,kmp,kani
      write(*,622) nfreepar,ipar,r200g,rcg,rsg,cbeg,tmassg,
     &             screeg,cbe0g,rjafg,xmstarg,dA3g,rbetag,kmp,kani
 622  format(/,' Number of free parameters = ',i1,/
     & 'r200, r_tr, r_s, anis1, A1, A2, anis2, rjaf, M/L, A3, rbeta = ',
     &  11(i1,1x),/,
     &         ' Initial guess values: ',11(f6.3,1x),/,
     &     ' mass model= ',i2,/,
     &     ' anisotropy model= ',i2,/)
      
      !For DHOST, the optimization algorithm cannot work with Y2
      if(kmp.eq.8.and.nhone.gt.1) then
       ipar(6)=0 
       nfreepar=nfreepar-1
       write(*,*) "MG second parameter not optimized in kinematic BH"
       write(*,*) " "
      endif
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
         cb0new=cbe0g
         rbetanew=rbetag
         scrnew=screeg
         tmassnew=tmassg
         dA3new = dA3g
         rjafnew=rjafg
         xmsnew=xmstarg
         goto 732
      endif
      if (kopt.eq.2.and.kani.gt.10) then
c     POWELL doesn't work for gT, gOM with MG parametrizations
        if (nfreepar.gt.1) then
          write(*,*) 'WARNING: POWELL does not work efficiently with'
          write(*,*) 'gOM, gT: switching to NEWUOA algorithm OPT=1'
          write(*,*) ''
          kopt=1
        endif
      endif
c     newuoa and bobyqa only work with at least 2 free params
      temp=r200g
      if (nfreepar.gt.1.and.ktried.lt.2.9) then

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

            ftol=1.d-5          ! desired fractional tolerance in function to be minimized
            npars=11  ! max number of free params
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
         npars=11                ! max number of free params
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
      call calfun(nfreepar,freepar,f) !compute the likelihood at
				      !the minimum found
      
      call freepareva(nfreepar,freepar,r200new,rcnew,rsnew,cbnew,
     & tmassnew,scrnew,dA3new,cb0new,rbetanew,rjafnew,xmsnew)
      !At the end of the optimization, if Y2 was a free parameter,
      !put it back to free
      if(kmp.eq.8.and.nhone.gt.1) ipar(6)=1 
      call sigmaeva(npar,sigma) 

      
      
      write(*,*) 
      write(*,*) 'opt.results: r200,rc,rs,beta,A1,A2,rj,xmb,cbe0,-logL'


      nfv=2 !dimension parameter for vmaxlik
       
      rs=rsnew
      cbe=cbnew
      cbe0=cb0new
      rbeta=rbetanew
      xfv(1)=rs
      xfv(2)=cbe
      r200=r200new
      tmass=tmassnew
      screen=scrnew
      dA3 = dA3new
      rjaf=rjafnew
      xmstar=xmsnew
      
      if (kmp.eq.11.or.kmp.eq.12) then
        xscCH=alimi
        alimit=alimi
        alimbc=0.0001
        call find_xsc(xscCH,alimit)
        
      endif        
      
      !Environmental screening  - correct advanced modelling
      if (kmp.eq.17) then
       xlim = 1.0e-3
       xscCH = ScreenGen(xlim)
      endif
      
      
      call vmaxlik(nfv,xfv,fml2)


       
      write(*,429) r200new, rcnew, rsnew, cbnew, tmassnew, scrnew, 
     & dA3new, cb0new,rbetanew,rjafnew,xmsnew, fml2
     
 429  format(11(f6.3,2x),f10.3)
      rs=rsg
      rc=rcg
      cbe=cbeg
      cbe0=cbe0g
      rbeta=rbetag
      xmstar=xmstarg
      rjaf=rjafg
      xfv(1)=rs
      xfv(2)=cbe
      r200=temp
      
      dA3 = dA3g
      tmass=tmassg
      screen=screeg
      if (kmp.eq.11.or.kmp.eq.12) then
        xscCH=alimi
        alimit=alimi
        alimbc=0.0001
        call find_xsc(xscCH,alimit)

      endif
      
      
      !Environmental screening  - correct advanced modelling
      if (kmp.eq.17) then
       xlim = 1.0e-3
       xscCH = ScreenGen(xlim)
      endif
        
      call vmaxlik(nfv,xfv,fml3)
      write(*,*) 'likelihood of the initial value and Delta chi square:'
      write(*,*) 'parametrs, -logL, Delta chi square:'
      write(*,439) r200,rc, rs, cbe,tmass,screen,cbe0,rjaf,xmstar,dA3,
     &   fml3, 2*(fml3-fml2)
 439  format(10(f6.3,2x),(f10.3,2x),f6.3)


          
      
      rcdf=dabs((rc-rcnew)/rcnew)
      rsdf=dabs((rs-rsnew)/rsnew)
      r2df=dabs((r200-r200new)/r200new)
      cbdf=dabs((cbe-cbnew)/cbnew)
      cb2f=dabs((cbe0-cb0new)/cb0new)
      rbdf=dabs((rbeta-rbetanew)/rbetanew)
      tmdf=dabs((tmass-tmassnew)/tmassnew)
      dA3df=dabs((dA3-dA3new)/dA3new)
      scdf=dabs((screen-scrnew)/scrnew)
      rjdf=dabs((rjaf-rjafnew)/rjafnew)
      xmdf=dabs((xmstar-xmsnew)/xmsnew)
      dfm=2*(fml3-fml2)
      
      if (istop.eq.1) then  !stops if required
       if(dfm.gt.(2*delik)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
        endif
       if (r2df.gt.teps(1).or.rcdf.gt.teps(2).or.scdf.gt.teps(6)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
        endif
       if (rsdf.gt.teps(3).or.cbdf.gt.teps(4).or.tmdf.gt.teps(5)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
       endif
       if (cb2f.gt.teps(7).or.rbdf.gt.teps(11)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
       endif
       if (rjdf.gt.teps(8)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
       endif
       if (xmdf.gt.teps(9)) then
         write(*,*) ''
	     write(*,*) 'Initial guess  not ',
     &     'fulfilling the requirements'
         stop
       endif
       if (dA3df.gt.teps(10)) then
         write(*,*) ''
	     write(*,*) 'Initial guess not fulfilling the requirements'
         stop
       endif       
       
       
      else
       kreq=0
       if(dfm.gt.(2*delik)) kreq=1
        
       if(r2df.gt.teps(1).or.rcdf.gt.teps(2).or.scdf.gt.teps(6)) kreq=1
        
       if (rsdf.gt.teps(3).or.cbdf.gt.teps(4).or.tmdf.gt.teps(5)) kreq=1
       if (cb2f.gt.teps(7).or.rbdf.gt.teps(11)) kreq=1
       if (rjdf.gt.teps(8)) kreq=1
       if (xmdf.gt.teps(9)) kreq=1
       if (dA3df.gt.teps(10)) kreq=1
       
       
       if (kreq.eq.1) then 
        write(*,*) ''
        write(*,*) 'Initial guess not fulfilling the requirements'
      write(*,*) 'As requested, the parameter exploration starts anyway'
       endif
       
       
      endif

c   if the best fit are out of the range, stop the code
      if(r200new.gt.r2up.or.r200new.lt.r2low) then
       write(*,*) 'Warning: best fit outside the prior range'
       write(*,*) 'redefine the best fit as the average value'
       r200new=(r2up+r2low)/2
      endif
      if(rsnew.gt.rsup.or.rsnew.lt.rslow) then
       write(*,*) 'Warning: best fit "rs" outside the prior range'
       write(*,*) 'redefine the best "rs" fit as the average value'
       rsnew=(rsup+rslow)/2
      endif 
      if(rcnew.gt.rcup.or.rcnew.lt.rclow) then
         write(*,*) 'Warning: best fit "rc" outside the prior range'
        write(*,*) 'redefine the best fit "rc" as the average value'
        rcnew=(rcup+rclow)/2
        
      endif 
      if(cbnew.gt.bup.or.cbnew.lt.blow) then
        write(*,*) 'Warning: best fit "beta" outside the prior range'
        write(*,*) 'redefine the best fit "beta" as the average value'
        cbnew=(bup+blow)/2
      endif 
      if(tmassnew.gt.tmup.or.tmassnew.lt.tmlow) then
        write(*,*) 'Warning: best fit "tmass" outside the prior range'
        write(*,*) 'redefine the best fit "tmass" as the average value'
        tmassnew=(tmup+tmlow)/2
      endif
      if(scrnew.gt.scrup.or.scrnew.lt.scrlow) then
       
       write(*,*) 'Warning: best fit "screen" outside the prior range'
       write(*,*) 'redefine the best fit "screen"  as the average value'
       scrnew=(scrup+scrlow)/2
       
      endif
      if(cb0new.gt.cb0up.or.cb0new.lt.cb0low) then
        write(*,*) 'Warning: best fit "beta2" outside the prior range'
        write(*,*) 'redefine the best fit "beta2" as the average value'
        cb0new=(cb0up+cb0low)/2
      endif 

      if(rbetanew.gt.rbetaup.or.rbetanew.lt.rbetalow) then
        write(*,*) 'Warning: best fit "rbeta" outside the prior range'
        write(*,*) 'redefine the best fit "rbeta" as the average value'
        rbetanew=(rbetaup+rbetanew)/2
      endif
      
      if(rjafnew.gt.rjaup.or.rjafnew.lt.rjalow) then
        write(*,*) 'Warning: best fit "rjaffe" outside the prior range'
        write(*,*) 'redefine the best fit "rjaffe" as the average value'
        rjafnew=(rjaup+rjalow)/2
      endif 
      
      if(xmsnew.gt.xmaup.or.xmsnew.lt.xmalow) then
        write(*,*) 'Warning: best fit "M*/L" outside the prior range'
        write(*,*) 'redefine the best fit "M*/L" as the average value'
        xmsnew=(xmaup-xmalow)/2
      endif 
      
      if(dA3new.gt.dA3up.or.dA3new.lt.dA3low) then
        write(*,*) 'Warning: best fit "dA3" outside the prior range'
        write(*,*) 'redefine the best fit "dA3" as the average value'
        dA3new=(dA3up-dA3low)/2
      endif 
      
      
 732  continue  
      
      write(*,*)  ''   
      write(*,*)  '**********************************'
      if(kbsp.eq.1) write(*,*) 'Running in Fast mode'
      if(kbsp.eq.0) write(*,*) 'Running in Normal mode'
      write(*,*)  '**********************************'
      write(*,*)  ''
      itotal=0
      

      call sigmaeva(npar,sigma)
      
     

      if(kbcg.eq.0.and.kmp.ne.13) then
       if(nlens.eq.0) then
         WRITE(iu60,'(11(a13,2x))') 'r200','rc','rs','cbe','A1','A2',
     & 'A3','cbe0','rbeta','-log(P)', 'anisotropy'
         WRITE(iu60,'(11(a13,2x))') '[Mpc]','[Mpc]','[Mpc]','-','-','-',
     & '-','-','[Mpc]', '-' ,'-' 
       else
         WRITE(iu60,'(11(a13,2x))') 'r200','rc','rs','cbe','A1','A2',
     & 'A3','cbe0','rbeta','-log(P*Plens)', 'anisotropy'
        WRITE(iu60,'(11(a13,2x))') '[Mpc]','[Mpc]','[Mpc]','-','-','-',
     & '-','','[Mpc]', '-', '-'      
       endif
       write(iu60, *) ''  
      else  
       if(nlens.eq.0) then
         WRITE(iu60,'(13(a13,2x))') 'r200','rc','rs','cbe','A1','A2',
     &   'A3','cbe0','rbeta','rjaffe','M/L','-log(P)', 'anisotropy'
         WRITE(iu60,'(13(a13,2x))') '[Mpc]','[Mpc]','[Mpc]','-','-','-',
     & '-','-','[Mpc]','[Mpc]','-', '-','-' 
       else
         WRITE(iu60,'(13(a13,2x))') 'r200','rc','rs','cbe','A1','A2',
     & 'A3','cbe0','rbeta','rjaffe','M/L','-log(P*Plens)', 'anisotropy'
        WRITE(iu60,'(13(a13,2x))') '[Mpc]','[Mpc]','[Mpc]','-','-','-',
     & '-','-', '[Mpc]','[Mpc]','-', '-','-'      
       endif
       write(iu60, *) ''  
      endif
      
      
      nfv = 2
      if (nmcmc.eq.1) then  !Start the MCMC run 
       if (nsame.eq.0) then
        icount=0 !Steps in the MCMC
        call CPU_TIME(start) !uses cpu time for the random generator seed
        nstar=int(10304*start, kind=8)
        nseed = 123434789+(-1)**nstar*nstar
        
        if (kopt.lt.0) then
          rs=rsnew
          cbe=cbnew
          cbe0=cb0new
          rbeta=rbetanew
          rc=rcnew
          xfv(1)=rs
          xfv(2)=cbe
          r200=r200new
          tmass=tmassnew
          screen=scrnew
          dA3 = dA3new
          rjaf=rjafnew
          xmstar=xmsnew
        

         if (kmp.eq.11.or.kmp.eq.12) then
            xscCH=alimi
            alimit=alimi
            alimbc=0.0001
            call find_xsc(xscCH,alimit)
         endif            


      !Environmental screening  - correct advanced modelling
         if (kmp.eq.17) then
          xlim = 1.0e-3
          xscCH = ScreenGen(xlim)
        endif
        
        
         call vmaxlik(nfv,xfv,f)
         !write(*,*) r200,rs,cbe, f
         
        endif
        fmlb=f
        
!********* add the BCG analysis ****************************************
        if (kbcg.eq.1.and.kbdt.eq.1) then
         rtbcg=20. !Max radius for the integration of the BCG integral
         call vdplik(chibcg)
         fmlb=fmlb+chibcg
        endif

        nlonly=0    !do only lensing analysis 
! **********************************************************************

        if (nlens.eq.1.and.naddlens.eq.0) then
        !If you have to add the lensing distribution, do it only 
        !when you don't have it in an additional file

         if (kmody.eq.1.and.kmp.ne.8) then
           fmlb=fmlb-Plens(r200,rs)
           if (nlonly.eq.1) then
            fmlb=-Nclust*Plens(r200,rs)
           endif
         endif
         !TO ADD LENSING BCG (constraints on rjaffe and gamma)
         !First attempt - only rs and r200
         if (kmp.eq.17) then 
           fmlb=fmlb-Plens(r200,rs)

           if (nlonly.eq.1) then
            fmlb=-Nclust*Plens(r200,rs)
           endif
         endif
         
         if(kmp.eq.8) then
           call Likelens_bh(nplens,plen)
           fmlb=fmlb+plen
           if (nlonly.eq.1) then
            fmlb=plen*Nclust
            
           endif
         endif
         
         
        endif
        
c============ voglio trasformare i bounds  =============================
c             trasformare tmass
      
      is_cham = (kmp.eq.9 .or. kmp.eq.11 .or. kmp.eq.12
     &           .or. kmp.eq.15 .or. kmp.eq.16 .or. kmp.eq.17)

      use_u = (is_cham .and. nLogP.eq.0 .and. nhone.le.0)
      
      if (is_cham .and. nhone.le.0) then
       if (dabs(scrnew - 0.40d0) .le. 0.02d0) then
          write(*,*) "Running in f(R) gravity: exploring log parameter"
          nLogP = 0
        endif
      endif
      
            
!       if (use_u .and. tmassnew.le.0.d0) then
!         write(*,*) 'ERROR: tmassnew<=0 while using log10(tmass)'
!         stop
!      endif
!      if (use_u) then
!        dtm_log = dlog10(tmassnew)
!        
!      else 
!        dtm_log = tmassnew
!        
!      endif
      

c****** defines temporary variables ************************************
        rsst=rsnew       
        cbt=cbnew
        cb0t=cb0new
        rbt=rbetanew
        r2t=r200new
        tmt= tmassnew
        dA3t = dA3new
        rct=rcnew
        scrt=scrnew
        rjat=rjafnew
        xmst=xmsnew
        
        !Baryons
        rhgta = rhog1
        rhosta = rhost1
        rsgnta = rsga1
        rstnta = rst1
    
        !add the initial prior to the likelihood
        if (ipriorGauss.ne.1) then
          call priord(r2t,rct,rsst,cbt,tmt,scrt,dA3t,
     &    cb0t,rbt,rjat,xmst,poini)
        else
          call priorG(r2t,rct,rsst,cbt,tmt,scrt,dA3t,
     &    cb0t,rbt,rjat,xmst,poini) 
        endif 
        
        !Add the distribution of the baryons
        !If no baryons included, this is just a constant term
        !call propose_baryons_rw_log(nseed, 0.05d0)
        !call baryon_prior_penalty(pbar)
        

        fmlb = fmlb - poini !(at the beginning, the baryons posterior is always zero)
        fml2 = fmlb !Check
        
        
        
c********************************************************************** 
        write(*,*) 'MCMC over ', Nsample, 'trials'   
        cdc=0 !Dhost model with Y1=Y2=y

        write (*,428) r2t, rct, rsst, cbt,
     &       tmt, scrt, dA3t, cb0t,rbt, rjat,xmst,fml2, plen
     !always accept the first trial
c          if (kbcg.eq.1.or.kmp.eq.13) then
c           write(iu60,274) r2t,rct,rsst,cbt,tmt,scrt,dA3t,cb0t,rbt,rjat,
c     &      xmst,fmlb,kani
c          else 
c            write(iu60,674) r2t,rct,rsst,cbt,tmt,scrt,dA3t,
c     &       cb0t,rbt,fmlb,kani 
c          endif

        do while (icount<Nsample)

         !if (mod(i,100).eq.0) then 
         !Change the seed for random generator
         !   nseed=abs(nseed+(-1)**(nseed)*floor(nseed/4.))
         !endif
         
          r200n=r2t       
          rcn=rct
          rsn=rsst
          cbn=cbt
          cb0n=cb0t
          rben=rbt
          tmassn=tmt
          dA3n = dA3t
          scrn=scrt
          rjan=rjat
          xmsn=xmst
          
          if (cdc.eq.1) scrn=tmt 
c          do i=1,9
c           if (float(itotal/500).lt.3) then
c            sigma(i)=(1+float(itotal/500))*sigma(i)
c           endif
c          enddo 
c ======== Check =======================================================
c ====== DEBUG =====================================
c         write(*,*) 'before'
c          write (*,428) r200n, rcn, rsn, cbn,
c     &       tmassn, scrn,dA3n,rjan,xmsn,fml2, fmlb !plen fml2
c
c           if(kbdt.eq.1) then
c            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
c     &      'rjaf,xmstar, -log(Pmam), plen, -log(Lbcg)'
c            write (*,418) r200n, rcn, rsn, cbn,
c     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, plen, chibcg2
c           else
c            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
c     &      'rjaf,xmstar, -log(Pmam), plen'
c            write (*,428) r200n, rcn, rsn, cbn,
c     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, plen           
c           endif     
c           write(*,*) rgmass(3.2d0)
c =================================================
         call tranmod(r200n,rcn,rsn,cbn,
     &           tmassn,scrn,dA3n,cb0n,rben,rjan,xmsn,sigma,npar,nseed)
         
c         write(*,*) 'after'
c ==================================================
         dlp = 0.0 ! lensing distribution
         if (naddlens.eq.1) then
          
          if (nlens.eq.0) then
           stop("ERROR! You selected nlens = 0, no lensing here")
          endif
          !Read the lensing data from the file 
          select case(kmp)
           case(9,11,12,15,16)
              read(82,*, iostat=ios) r200_f, rs_f, dlp 
               if(ios < 0) then
                close(82)
                exit
              endif
              rsn = rs_f
              r200n = r200_f
              
           case(10,17) !when you have gamma as free parameter    
              read(82,*, iostat=ios) r200_f, rs_f, dA3_f, dlp  
              
              if(ios < 0) then
                close(82)
                exit
              endif
              rsn = rs_f
              r200n = r200_f
              dA3n = dA3_f
                                                                                             
           case(8) !For DHOST, you read also the MG parameters
              read(82,*, iostat=ios) r200_f, rs_f, dA1_f, dA2_f, dlp 
              if(ios < 0) then
                close(82)
                exit
              endif
              rsn = rs_f
              r200n = r200_f
              tmassn = dA1_f
              scrn  =  dA2_f
              
           case default
              read(82,*, iostat=ios) r200_f, rs_f, dlp 
              if(ios < 0) then
                close(82)
                exit
              endif
              rsn = rs_f
              r200n = r200_f
           end select
          
         endif 

        !Now add the baryons
        rhgn = rhgta
        rhosn = rhosta
        rsgn = rsgnta
        rstn = rstnta
        
        
        call propose_baryons_rw_lin(rhgn, rhosn,rsgn,rstn,
     &           nseed, 0.1d0)
    
        rs=rsn
        cbe=cbn
        cbe0=cb0n
        rbeta=rben
        r200=r200n
        tmass=tmassn
        dA3 = dA3n
        rjaf=rjan
        xmstar=xmsn
        screen=scrn
        !Redefine the baryons for the likelihood
        rhog1 = rhgn
        rhost1 = rhosn
        rsga1 = rsgn
        rst1 = rstn
        if(cdc.eq.1) screen=tmassn
        rc=rcn
                  
        omegal=Olam
        omega0=Omegam
        hz=h0*sqrt(omega0*(1.+za)**3+omegal)
        rm200=100.*hz*hz/grav*r200**3
        cmean=6.76*(rm200/1.e12)**(-0.098)
        v200=10.*hz*r200
        !Checkpoint 2
        !write(*,*) r200, rs, rbeta, cbe, cbe0

c     Mass follows light
c
         if (kmfl.eq.1) then
          rs=rc
         elseif (klcdm.eq.1) then
c
c     concentration from c=c(M)
c
          rs=r200/cmean

c     Analysis of the Euclid sample. Ha-line galaxies, should have rc = 4*rs          
          if (kfs2.eq.1) rc = 4 * rs
         endif
         
c     a_ML forced = r_s
c
         if (kaml.eq.1) cbe=rs

c     
c     Light follows Mass (TLM)
c
         if (klfm.eq.1) then
            rc=rs  ! assume N(R) and M(r) chosen with same model
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
         
c  Screening radius for Einasto or Burkert profile 
         if (kmp.eq.11.or.kmp.eq.12) then
            xscCH=alimi
            alimit=alimi
             alimbc=0.0001
        
            call find_xsc(xscCH,alimit) 

         endif         
         
      !Environmental screening  - correct advanced modelling
        if (kmp.eq.17) then !.or.kmp.eq.9) then
         xlim = 1.0e-3
         xscCH = ScreenGen(xlim)
       
        endif         

        
    !        write (*,428) r200n, rcn, rsn, cbn,
    ! &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, plen  
        

        call vmaxlik(nfv,xfv,fml2)
       
        
        Nrun=Nclust

!********* add the BCG analysis ****************************************

         if (Nclust.eq.1.and.kbcg.eq.1.and.kbdt.eq.1) then
          rtbcg=20. !Max radius for the integration of the BCG integral
          call vdplik(chibcg2)
          fml2=fml2+chibcg2
          
         endif
         if (nsingle.eq.1) then
          fml2=fml2*Nclust
          Nrun=1 !if nsingle=1 abort the cycle over Nrun different clusters
         endif
!        read all the files for 
 9        format(a)
          filename='data/datphys_'
          dat='.dat'
          do irun=2,Nrun
           if(irun.lt.10) then
            ARR1 = ''
            WRITE(ARR1,"(I1)") irun

            open (12, FILE=filename//ARR1//dat, STATUS='OLD')
           else
            ARR2 = ' '
            WRITE(ARR2,"(I2)") irun
            open (12, FILE=filename//arr2//dat, STATUS='OLD')
           endif
               read(12,9) line
               read(12,9) line
               j=0
               j0=-1
               dimin=1.e12
 442           continue
                
                read(12,*,end=111) dkpc,vkms,evkms 
                j=j+1
                dpp(j)=dkpc/1.e3
                vpp(j)=vkms
                epp(j)=evkms
                wpp(j)=1.0
                
                if (dpp(j).lt.dimin) then
                 dimin=dpp(j)
                 j0=j
               endif
               goto 442
 111          continue
            close(12)
              npp=j
              jgiu=0
              do j=1,npp
               if (dpp(j).ge.rlowinmpc.and.
     &          dpp(j).le.rupinmpc) then
                jgiu=jgiu+1
                did(jgiu)=dpp(j)
                viv(jgiu)=(vpp(j)-va)/(1.+va/clight)
                eie(jgiu)=epp(j)
                wiw(jgiu)=wpp(j)
                if (did(jgiu).le.rlow) rlow=did(jgiu)
                if (did(jgiu).ge.rup) rup=did(jgiu)
               endif               
              enddo
              ngam=jgiu
              
             call vmaxSUM(nfv,xfv,fprova,did,viv,eie,wiw,ngam)
             fml2=fml2+fprova 
          enddo

         plen=0.0d0
         if (nlens.eq.1) then  
          if (kmody.eq.1.and.kmp.ne.8) then !Chameleon modified gravity
            plen=-Plens(r200n,rsn) 
            if (naddlens.eq.1) plen = 0.0d0 !dlp do not resample twice!!
            fml2=fml2 + Nclust*plen
            
            if (nlonly.eq.1) fml2= Nclust*plen
          endif
          !TO ADD BCG LENSING
         !First attempt - only rs and r200
         if (kmp.eq.17) then 
            fml2=fml2-Nclust*Plens(r200n,rsn)
            plen=-Plens(r200n,rsn)
           if (nlonly.eq.1) then
            fml2=-Nclust*Plens(r200n,rsn)
           endif
         endif
          if(kmp.eq.8) then
           if (naddlens.eq.0) then
            call Likelens_bh(nplens,plen)
           else
            plen = 0.0d0 !dlp do not resample twice!!
           endif
           fml2=fml2+Nclust*plen
           if (nlonly.eq.1) fml2=plen*Nclust
          endif
         endif
c**************** compute the priors (flat prior) **********************   
         !Baryons
         call baryon_prior_penalty(pbar) 

!MCMC WP        
         if (ipriorGauss.ne.1) then
          call priord(r200n,rcn,rsn,cbn,tmassn,scrn,dA3n,
     &    cb0n,rben,rjan,xmsn,pout1)
         else
          call priorG(r200n,rcn,rsn,cbn,tmassn,scrn,dA3n,
     &    cb0n,rben,rjan,xmsn,pout1) 
         endif        
         !acceptance criterium:
         fml2 = fml2 - pout1 + pbar
         call acceptance(-fmlb,(-fml2),eval,nseed)
c ======== Check =======================================================
C ========= DEBUG ==================================
        
c        write(*,*) rhosta/1e13,tmassn,scrn,dA3n, fml2, fmlb
c           if(kbdt.eq.1) then
c            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
c     &      'rjaf,xmstar, -log(Pmam), plen, -log(Lbcg)'
c            write (*,418) r200n, rcn, rsn, cbn,
c     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, plen, chibcg2
c           else
c            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
c     &      'rjaf,xmstar, -log(Pmam), plen'
c            write (*,428) r200n, rcn, rsn, cbn,
c     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, plen           
c           endif     
c           write(*,*) rsga1, rhog1, rst1, rhost1
c        if (kmp.eq.13) then
c            write (*,987) tmassn, scrn, dA3n, fml2, pout1
c        endif
c 987     format(3(f7.5,2x),f15.3,2x, f15.3)           
c =================================================
         itotal=itotal+1
            
         ! make sigma smaller if you are not finding a new point
         if (mod(itotal,500).eq.0) then

           call CPU_TIME(start) !uses cpu time for the random generator seed
           nstar=int(10040*start, kind=8)
           nseed = nseed+(-1)**nstar*nstar
           !do i=1,npar !test
           ! sigma(i)=0.2*sigma(i)
           !enddo
         endif
         
         if(eval) then
          icount=icount+1
          itotal=0 !number of total trials

          call sigmaeva(npar,sigma)
     
          if (mod(icount,100).eq.0) then
           write(*,*) icount, 'trials accepted'
           write(*,*) ' '
           if(kbdt.eq.1) then
            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
     &      'rjaf,xmstar, -log(Pmam), plen, -log(Lbcg)'
            write (*,418) r200n, rcn, rsn, cbn,
     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, dlp, chibcg2
           else
            write(*,*) 'r200,rc,rs,cbe,tmass,screen,A3,cbe0,rbeta,',
     &      'rjaf,xmstar, -log(Pmam), plen'
            write (*,428) r200n, rcn, rsn, cbn,
     &       tmassn, scrn, dA3n, cb0n,rben,rjan,xmsn,fml2, dlp           
           endif     
          endif
          
c  trial accepted: update the temporary parameters *********************      
          r2t=r200n     
          rct=rcn
          rsst=rsn
          cbt=cbn
          cb0t=cb0n
          rbt = rben
          tmt=tmassn
          scrt=scrn
          dA3t = dA3n
          rjat=rjan
          xmst=xmsn
          scprova=scrt
          
          !Baryons
          rhgta = rhgn
          rhosta  = rhosn
          rsgnta = rsgn 
          rstnta = rstn 
          
          fmlb = fml2
c******************* write on the file *********************************          
          if (kbcg.eq.1.or.kmp.eq.13) then
           write(iu60,274) r2t,rct,rsst,cbt,tmt,scprova,dA3t,cb0t,rbt,
     &      rjat,xmst,fml2+dlp,kani
          else 
            write(iu60,674) r2t,rct,rsst,cbt,tmt,scprova,dA3t,
     &       cb0t,rbt,fml2+dlp,kani 
          endif
         
c***************** Baryons ********************************************* 
          if (kbary.eq.1) then
           write(92,271) rhgta/1e13,rsgnta,rhosta/1e13,rstnta
271        format(4(f15.5,2x))
          endif
c***********************************************************************              
          if (fml2.lt.fmlmingrid) then
                  plenmingrid=plen
                  fmlmingrid=fml2
                  r200mingrid=r2t
                  rsmingrid=rsst
                  rcmingrid=rct
                  cbmingrid=cbt
                  cb0mingrid=cb0t
                  rbmingrid=rbt
	  	          tmingrid=tmt
                  smingrid=scprova
                  dAmingrid = dA3t
                  rjmingrid=rjat
                  xmmingrid=xmst
                  
                  xfn1best=rct
                  ngamin=nga
                  rlowmin=rlow
                  rupmin=rup
                  rcmin=rc
                  do j=1,nga
                     rso(j)=r(j)
                     vso(j)=v(j)
                  enddo        
            endif
          
         endif
         
 428     format(11(f6.3,2x),f15.3,2x,f13.3)
 418     format(11(f6.3,2x),f15.3,2x,f10.3,2x,f10.3)
       enddo
       
       write(*,*) icount, 'trials accepted'
       plenmin=plenmingrid
       fmlmin=fmlmingrid
       r200min=r200mingrid
       rsmin=rsmingrid
       rcmin=rcmingrid
       cbmin=cbmingrid
       cb0min=cb0mingrid
       rbmin=rbmingrid
       tmmin=tmingrid
       scmin=smingrid
       dA3min = dAmingrid
       rjmin=rjmingrid
       xmmin=xmmingrid
       
       rc=xfn1best
       nga=ngamin
       rlow=rlowmin
       rup=rupmin
       write(iu60,*) '.................................................'
       
       if (kbcg.eq.1.or.kmp.eq.13) then
           write(iu60,274) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,dA3min,cb0min,rbmin,rjmin,xmmin,fmlmin,kani
          else 
           write(iu60,674) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,dA3min,cb0min,rbmin,fmlmin,kani
       endif
       
       write(iu60,*) '.................................................'
      
      
       if (kopt.lt.0) then          ! optimization skipped
         r200new=r200min
         rcnew=rcmin
         rsnew=rsmin
         cbnew=cbmin
         cb0new=cb0min
         rbetanew=rbmin
         scrnew=scmin
         tmassnew=tmmin
         dA3new = dA3min
         rjafnew=rjmin
         xmsnew=xmmin
         plenew=plenmin
         f=fmlmin
       endif
       
       if (kbcg.eq.1.or.kmp.eq.13) then
           write(iu60,274) r200new,rcnew,rsnew,cbnew,tmassnew,scrnew,
     &                 dA3new,cb0new,rbetanew,rjafnew,xmsnew,f,kani
          else 
           write(iu60,674) r200new,rcnew,rsnew,cbnew,tmassnew,scrnew,
     &                 dA3new, cb0new,rbetanew,f,kani
       endif       
       

 674           format(10(f15.5,2x),i2)
 274           format(12(f15.5,2x),i2)
c    
c
c    if requested, it compute the likelihood on fixed values of the
c    parameters which are the results of a previous MCMC
c    This is a feature introduce to obtain a combined likelihood
c    in the analysis of syntetic clusters produced with the aim of 
c    ClusterGEN code   
      elseif(nsame.eq.1) then 
      open(1,file='MaxLik_input.dat',iostat=niost,status='old')

      if (niost.ne.0) then
       write(*,*) 'ERROR: in nsame=1'
       write(*,*) 'ERROR: INPUT FILE "MaxLik_input.dat" NOT FOUND'
       stop
      endif
       
  22  continue
      
      READ (1,fmt="(11(f13.5,2x))",end=23) r200,rc,rs,cbe,
     & tmass, screen, dA3, cbe0,rbeta,rjaf,xmstar
       

         
         if (kmp.eq.11.or.kmp.eq.12) then
            xscCH=alimi
            alimit=alimi
            alimbc=0.0001
            call find_xsc(xscCH,alimit)

         endif              
       
        !Environmental screening  - correct advanced modelling
         if (kmp.eq.17) then
          xlim = 1.0e-3
          xscCH = ScreenGen(xlim)
         endif
       
       xfv(1)=rs
       xfv(2)=cbe
       call vmaxlik(nfv,xfv,fml2)
       if (kbcg.eq.1.and.kbdt.eq.1) then
          rtbcg=20. !Max radius for the integration of the BCG integral
          call vdplik(chibcg2)
          fml2=fml2+chibcg2
       endif 
       write(*,fmt="(12(f13.5,2x))") r200,rc,rs,cbe,tmass,
     & screen,dA3,cbe0,rbeta,rjaf,xmstar,fml2

          if (kbcg.eq.1.or.kmp.eq.13) then
            write(iu60,274) r200,rc,rs,cbe,tmass,screen,dA3,cbe0,rbeta,
     &          rjaf,xmstar,fml2,kani
          else 
            write(iu60,674) r200,rc,rs,cbe,tmass,screen,dA3,cbe0,rbeta,
     &         fml2,kani
          endif
       
          
           if (fml2.lt.fmlmingrid) then
                  fmlmingrid=fml2
                  r200mingrid=r200
                  rsmingrid=rs
                  rcmingrid=rc
                  cbmingrid=cbe
                  cb0mingrid=cbe0
                  rbmingrid=rbeta
	  	          tmingrid=tmass
                  smingrid=screen
                  dAmingrid = dA3  
                  rjmingrid=rjat
                  xmmingrid=xmst
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

       
      goto 22
  23  continue
       
      close(1)
       fmlmin=fmlmingrid
       r200min=r200mingrid
       rsmin=rsmingrid
       rcmin=rcmingrid
       cbmin=cbmingrid
       cb0min=cb0mingrid
       rbmin=rbmingrid
       tmmin=tmingrid
       scmin=smingrid
       dA3min = dAmingrid
       rjmin=rjmingrid
       xmmin=xmmingrid
       rc=xfn1best
       nga=ngamin
       rlow=rlowmin
       rup=rupmin
       write(iu60,*) '.................................................'
       if (kbcg.eq.1.or.kmp.eq.13) then
           write(iu60,274) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,dA3min,cb0min,rbmin,rjmin,xmmin,fmlmin,kani
          else 
           write(iu60,674) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,dA3min,cb0min,rbmin,fmlmin,kani
       endif
       

       write(iu60,*) '.................................................'
      
      
       if (kopt.lt.0) then          ! optimization skipped
         r200new=r200min
         rcnew=rcmin
         rsnew=rsmin
         cbnew=cbmin
         cb0new=cb0min
         rbetanew=rbmin
         scrnew=scmin
         tmassnew=tmmin
         dA3new = dA3min
         rjafnew=rjmin
         xmsnew=xmmin
         f=fmlmin
       endif
       if (kbcg.eq.1.or.kmp.eq.13) then
           write(iu60,274) r200new,rcnew,rsnew,cbnew,tmassnew,scrnew,
     &                 dA3new,cb0new,rbetanew,rjafnew,xmsnew,f,kani
          else 
           write(iu60,674) r200new,rcnew,rsnew,cbnew,tmassnew,scrnew,
     &                 dA3new,cb0new,rbetanew,f,kani
       endif  
      endif
      if (naddlens.eq.1) then
       close(82) !Close the lensing file if it was open
      endif
******************************** end MCMC run **************************        
      else if(nmcmc.eq.0) then
      
       if (kmp.eq.13.or.kmp.eq.14.or.kmp.eq.17.or.kani.eq.8) then
        stop('ERORR: RG,BS,BP and mTotal_GC setups work only with mcmc')
       endif
      fmlb=f

      do j=1,nga
         rso(j)=r(j)
         vso(j)=v(j)
      enddo

c
c     Grid search around the minimum found
c
c     Store the likelihood value for the initial guess parameters (true r200)

      if (kfs2.eq.1) then      
       nfv=2 !dimension parameter for vmaxlik
       
       r200=r200new
       
       omegal=Olam!6
       omega0=Omegam!2.*(q0+omegal)
       hz=h0*sqrt(omega0*(1.+za)**3+omegal)
       rm200=100.*hz*hz/grav*r200**3
cc      cduffy=5.78*(rm200/2.e12)**(-0.089)*1.1**(-0.52)
cc      cgao=10.**(-0.138*dlog10(rm200*0.7)+2.646)
cc      cmean=(cduffy+cgao)/2.
       cmean=6.76*(rm200/1.e12)**(-0.098)
       v200=10.*hz*r200
       
       rs=r200/cmean
       cbe=cbnew
       cbe0=cb0new
       xfv(1)=rs
       xfv(2)=cbe
       rc = 4 * rs
       tmass=tmassnew
       screen=scrnew
       dA3 = dA3new
       rjaf=rjafnew
       xmstar=xmsnew      
       call vmaxlik(nfv,xfv,fml3)
      endif
      
      
      klfm=0
      
      
         if (kpro.eq.1.and.nr200.gt.1) then
          call test_bound(r200g,r2low,r2up)
          kr200=-nr200/2
          deltr2m=dlog10(r2low/r200g)/(-nr200/2)
          deltr2p=dlog10(r2up/r200g)/(nr200/2)
         else
          deltr2m=0.0
          deltr2p=0.0
         endif
      

      if (nrc.ge.1) then
         sigma(2)=0.5d0   !TEST
         if (nskip.eq.0) rcg=rcnew
         dd=1.4*(500./nga)**(0.2)  ! this adjust the grid width to the number of galaxies
         deltarc=dd/nrc
         nrc1=-nrc/2
         nrc2=nrc/2
         if (kpro.eq.1.and.nrc.gt.1) then
          call test_bound(rcg,rclow,rcup)
          deltRcm=dlog10(rclow/rcg)/nrc1
          deltRcp=dlog10(rcup/rcg)/nrc2
         endif
      else
         sigma(2)=0.0d0
         if (nrc.eq.-1) klfm=1
         deltarc=0.
         deltRcm=0.0
         deltRcp=0.0
         nrc1=1
         nrc2=1
      endif
      
      klin=1 !linear spacing for horndeski coupling
      if (nhone.ge.1) then 
         sigma(6)=0.5d0   !TEST
         if (nskip.eq.0) screeg=scrnew !use the optimization-found value
         dd=2.3
         deltaSc=dd/nhone
         nhone1=-nhone/2
         nhone2=nhone/2
c***********************************************************************        
c*********************************************************************** 
         if (kpro.eq.1.and.nhone.gt.1) then
          if (kmp.ne.8.or.klin.ne.1) then 
         
            call test_bound(screeg,scrlow,scrup)
            if (scrlow.le.0) then
             write(*,*) 'inconsistent negative parameter!'
             scrlow=1e-2
            endif
            deltaSm=dlog10(scrlow/screeg)/nhone1
            deltaSp=dlog10(scrup/screeg)/nhone2

          else
           call test_bound(screeg,scrlow,scrup)
           ssafe=(dabs(scrlow)+0.05*dabs(scrlow))
           sl1=scrlow+ssafe
           smt=screeg+ssafe
           sl2=scrup+ssafe
                     
            deltaSm=dlog10(sl1/smt)/(-nhone/2)
            deltaSp=dlog10(sl2/smt)/(nhone/2)
            
          endif
         endif

c**********************************************************************
          
          
      else
         sigma(6)=0.0d0
         deltaSc=0.
         deltaSm=0.0
         deltaSp=0.0
         nhone1=1
         nhone2=1
      endif

      if (nrs.ge.1) then
         sigma(3)=0.5   !TEST
         if (nskip.eq.0) rsg=rsnew !use the optimization-found value
         dd=1.3*(500./nga)**(0.2)  !1.1 this adjust the grid width to the number of galaxies
	 
         if (nres.eq.1) dd=0.4*(500./nga)**(0.2) !finer grid
	    
	     if (nequ.eq.1) then
	      dd=1.1 !equal grid points 
          if (nres.eq.1) dd=0.4 !equal grid points and finer grid
         endif

         deltars=dd/nrs
         nrs1=-nrs/2
         nrs2=nrs/2
        if (kmp.eq.8) then !for beyond Horndeski gravity, shift towards large rs values
	      nrs1=-nrs/2+4
          nrs2=nrs/2+4
        endif
         kmfl=0
         klcdm=0
        if (kpro.eq.1.and.nrs.gt.1) then
          call test_bound(rsg,rslow,rsup)
          deltrsm=dlog10(rslow/rsg)/nrs1
          deltrsp=dlog10(rsup/rsg)/nrs2
        endif
      else
         if (nrs.eq.-1) kmfl=1
         if (nrs.eq.-2) klcdm=1
          deltars=0.
          deltrsm=0.0
          deltrsp=0.0
         nrs1=1
         nrs2=1
         sigma(3)=0.0d0
      endif

      if (nbs.ge.1) then
         sigma(4)=0.9   !TEST
	     if (nskip.eq.0) cbeg=cbnew
         if (kani.eq.0) then
            dd=0.9*(500./nga)**(0.2)
            if (nres.eq.1) dd=0.4*(500./nga)**(0.2)
	    
            if (nequ.eq.1) then
	          dd=0.9 !equal grid points 
              if (nres.eq.1) dd=0.4 !equal grid points and finer grid
            endif
            deltab=dd/nbs
         elseif (kani.eq.1) then
            dd=3.7*(500./nga)**(0.2)
            if (nres.eq.1) dd=1.4*(500./nga)**(0.2)
            if (nequ.eq.1) then
             dd=3.7 !equal grid points 
             if (nres.eq.1) dd=1.4 !equal grid points and finer grid
            endif

            deltab=dd/nbs

         elseif (kani.eq.2.or.kani.eq.3.or.kani.eq.21) then
            dd=1.3*(500./nga)**(0.2)
            deltab=dd/nbs
         else
            dd=0.9*(500./nga)**(0.2)
            
            if (nres.eq.1) dd=0.4*(500./nga)**(0.2)
            if (nequ.eq.1) then
             dd=0.9 !equal grid points 
             if (nres.eq.1) dd=0.4d0 !equal grid points and finer grid
            endif
	      deltab=dd/nbs
         endif
         nbs1=-nbs/2
         nbs2=nbs/2
         
         if (kpro.eq.1.and.nbs.gt.1) then
          call test_bound(cbeg,blow,bup)
          deltbm=dlog10(blow/cbeg)/nbs1
          deltbp=dlog10(bup/cbeg)/nbs2
         endif
      elseif (nbs.eq.-1) then 
         kaml=1
         nbs1=1
         nbs2=1
      elseif (nbs.eq.-2) then 
         khm=1
         nbs1=1
         nbs2=1
      else
         sigma(4)=0.0d0
         deltab=0.
         deltbm=0.
         deltbp=0.
         nbs1=1
         nbs2=1
      endif
c******************************************************************+      
            if (nrc.ge.1) then
         sigma(2)=0.5d0   !TEST
         rcg=rcnew
         dd=1.4*(500./nga)**(0.2)  ! this adjust the grid width to the number of galaxies
         deltarc=dd/nrc
         nrc1=-nrc/2
         nrc2=nrc/2
         if (kpro.eq.1.and.nrc.gt.1) then
          call test_bound(rcg,rclow,rcup)
          deltRcm=dlog10(rclow/rcg)/nrc1
          deltRcp=dlog10(rcup/rcg)/nrc2
         endif
      else
         sigma(2)=0.0d0
         if (nrc.eq.-1) klfm=1
         deltarc=0.
         deltRcm=0.0
         deltRcp=0.0
         nrc1=1
         nrc2=1
      endif
            
      
      if (nb2.ge.1.and.kani.gt.20) then

         sigma(7)=0.9   !TEST
	     if (nskip.eq.0) cbe0g=cb0new
         dd=0.9*(500./nga)**(0.2)
         if (nres.eq.1) dd=0.4*(500./nga)**(0.2)
         if (nequ.eq.1) then
           dd=0.9 !equal grid points 
           if (nres.eq.1) dd=0.4d0 !equal grid points and finer grid
         endif
         delta0b=dd/nb2
         nb21=-nb2/2
         nb22=nb2/2
         
         if (kpro.eq.1.and.nb2.gt.1) then
          call test_bound(cbe0g,cb0low,cb0up)
          deltb0m=dlog10(cb0low/cbe0g)/nb21
          deltb0p=dlog10(cb0up/cbe0g)/nb22
         endif
      else
         sigma(7)=0.0d0
         delta0b=0.
         deltb0m=0.0
         deltb0p=0.0
         nb21=1
         nb22=1
      endif
      
      
      kxb=0
      if (kbcg.eq.1.or.kmp.eq.13.or.kmp.eq.17) kxb=1
      
      if (njaf.ge.1.and.kxb.gt.0) then

         sigma(8)=0.5   
	     if (nskip.eq.0) rjaf=rjafg
         dd=0.9*(500./nga)**(0.2)
         if (nres.eq.1) dd=0.4*(500./nga)**(0.2)
         if (nequ.eq.1) then
           dd=0.9 !equal grid points 
           if (nres.eq.1) dd=0.4d0
         endif
         deljaf=dd/njaf
         njaf1=-njaf/2
         njaf2=njaf/2
         
         if (kpro.eq.1.and.njaf.gt.1) then
          call test_bound(rjafg,rjalow,rjaup)
          deltjam=dlog10(rjalow/rjafg)/njaf1
          deltjap=dlog10(rjaup/rjafg)/njaf2
         endif
      else
         sigma(8)=0.0d0
         deljaf=0.
         deltjam=0.0
         deltjap=0.0
         njaf1=1
         njaf2=1
      endif
      
      if (nxmas.ge.1.and.kxb.gt.0) then

         sigma(9)=0.2   
	     if (nskip.eq.0) xmstar=xmstarg
         dd=0.6*(500./nga)**(0.2)
         if (nres.eq.1) dd=0.4*(500./nga)**(0.2)
         if (nequ.eq.1) then
           dd=0.7 !equal grid points 
           if (nres.eq.1) dd=0.4d0
         endif
         delxm=dd/nxmas
         nxm1=-nxmas/2
         nxm2=nxmas/2

         if (kpro.eq.1.and.nxmas.gt.1) then
          call test_bound(xmstarg,xmalow,xmaup)
          delxmm=dlog10(xmalow/xmstarg)/nxm1
          delxmp=dlog10(xmaup/xmstarg)/nxm2
         endif
      else
         sigma(9)=0.0d0
         delxm=0.
         delxmm=0.0
         delxmp=0.0
         nxm1=1
         nxm2=1
      endif      
      
      
      
      if (ntmass.ge.1) then
        sigma(5)=0.9   
        if (nskip.eq.0) tmassg=tmassnew 
          
        if(kmp.eq.8) then !case Beyond Horndeski
           !in this case the grid is by default equally spaced
           if(tmassg.lt.1) then      
	        deltam=14.4/ntmass
           else
            deltam=10.0/ntmass
           endif 
           
           if (kpro.eq.1.and.ntmass.gt.1) then
             call test_bound(tmassg,tmlow,tmup)
             tsafe=(dabs(tmlow)+0.05*dabs(tmlow))
             tl1=tmlow+tsafe
             tmt=tmassg+tsafe
             tl2=tmup+tsafe
             deltmm=dlog10(tl1/tmt)/(-ntmass/2)
             deltmp=dlog10(tl2/tmt)/(ntmass/2)
           endif
       
        else !case general Chameleon 
          dd=3.9*(500./nga)**(0.2)  ! this adjust the grid width to the number of galaxies
          if (nres.eq.1) dd=3.8*(500./nga)**(0.2)
          if (nequ.eq.1) then
	       dd=3.9d0 !equal grid points 
           if (nres.eq.1) dd=3.8d0 !equal grid points and finer grid
          endif
          deltam=dd/ntmass
          if (kpro.eq.1.and.ntmass.gt.1) then
           call test_bound(tmassg,tmlow,tmup)
           deltmm=dlog10(tmlow/tmassg)/(-ntmass/2)
           deltmp=dlog10(tmup/tmassg)/(ntmass/2)
          endif
          
      endif
      
         ntm1=-ntmass/2
         ntm2=ntmass/2
         
      else
         sigma(5)=0.0d0
         deltam=0.
         deltmm=0.0
         deltmp=0.0
         ntm1=1
         ntm2=1
      endif

      nfv=2
      
      if (nskip.eq.0) r200g=r200new !use the optimization-found value

                     
      
      
      if (nr200.le.1.and.nrc.le.1.and.nrs.le.1.and.
     &   nbs.le.1.and.ntmass.le.1.and.nhone.le.1.and.nb2.le.1.and.
     &   njaf.le.1.and.nxmas.le.1) goto 891
      
 727  continue

      kr200=kr200+1
      if (nr200.ge.2) then
         r200=10.**(dlog10(r200g)+kr200*deltarv)
         if (kpro.eq.1) then
            if(kr200.lt.0) r200=10.**(dlog10(r200g)+kr200*deltr2m)
            if(kr200.ge.0) r200=10.**(dlog10(r200g)+kr200*deltr2p)
         endif
         
      endif

c     Determine m200, v200 and c appropriate for the m200 found
c     adopting the relation of Maccio, Dutton, van den Bosch 08
c     for relaxed halos

      omegal=Olam!6
      omega0=Omegam!2.*(q0+omegal)
      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
      rm200=100.*hz*hz/grav*r200**3
cc      cduffy=5.78*(rm200/2.e12)**(-0.089)*1.1**(-0.52)
cc      cgao=10.**(-0.138*dlog10(rm200*0.7)+2.646)
cc      cmean=(cduffy+cgao)/2.
      cmean=6.76*(rm200/1.e12)**(-0.098)
      v200=10.*hz*r200

c     Search on a grid
      do ms=ntm1,ntm2

       if(kmp.ne.8) then 
          tmass=10**(dlog10(tmassg)+ms*deltam)                    

          if (kpro.eq.1.and.ntmass.gt.1) then
           if(ms.lt.0) tmass=10.**(dlog10(tmassg)+ms*deltmm)
           if(ms.ge.0) tmass=10.**(dlog10(tmassg)+ms*deltmp)
          endif
          
       else
          !linear spacing if MG parameter is Y

          tmass=tmassg+ms*deltam
          if (kpro.eq.1.and.ntmass.gt.1) then
           tsafe=(dabs(tmlow)+0.05*dabs(tmlow))
           tl1=tmlow+tsafe
           tmt=tmassg+tsafe
           tl2=tmup+tsafe
           if(ms.lt.0) tmass=10.**(dlog10(tmt)+ms*deltmm)-tsafe
           if(ms.ge.0) tmass=10.**(dlog10(tmt)+ms*deltmp)-tsafe
          endif
            
       endif
       

       write(*,239) za,va,r200,v200,tmass
 239   format(//' New run with  <z>,<v>,r200,v200, MG param = ',
     &     f7.3,2x,f7.0,2x,2x,f6.3,2x,f5.0,2x,f10.3)

c     Best-fit to N(R) external to MAMPOSSt, if required

       if (nrc.eq.-2) then

c        write(*,*) ' Fit to N(R) up to R/r200=',rup/r200

c        write(*,293) rc,r200/rc
 293    format(' Best-fit N(R) scale parameter= ',f7.3,' (c=',f5.2,')')

       endif
       do nhc=nhone1,nhone2
        
        if(kmp.ne.8.or.klin.ne.1) then
         screen=10.**(dlog10(screeg)+nhc*deltaSc)
                
        
         if (kpro.eq.1.and.nhone.gt.1) then
           if(nhc.lt.0) screen=10.**(dlog10(screeg)+nhc*deltaSm)
           if(nhc.ge.0) screen=10.**(dlog10(screeg)+nhc*deltaSp)
         endif
        else !linear spacing
         screen=screeg+nhc*deltaSc
         if (kpro.eq.1.and.nhone.gt.1) then
           ssafe=(dabs(scrlow)+0.05*dabs(scrlow))
           sl1=scrlow+ssafe
           smt=screeg+ssafe
           sl2=scrup+ssafe
           if(nhc.lt.0) screen=10.**(dlog10(smt)+nhc*deltaSm)-ssafe
           if(nhc.ge.0) screen=10.**(dlog10(smt)+nhc*deltaSp)-ssafe
          endif
        endif 
        do kb=nbs1,nbs2         
         cbe=10.**(dlog10(cbeg)+kb*deltab)
         if (kpro.eq.1.and.nbs.gt.1) then
          if (kb.lt.0)  cbe=10.**(dlog10(cbeg)+kb*deltbm)
          if (kb.ge.0)  cbe=10.**(dlog10(cbeg)+kb*deltbp)
         endif
         
         do kb2=nb21,nb22         
         cbe0=10.**(dlog10(cbe0g)+kb2*delta0b)
         if (kpro.eq.1.and.nb2.gt.1) then
          if (kb2.lt.0)  cbe0=10.**(dlog10(cbe0g)+kb2*deltb0m)
          if (kb2.ge.0)  cbe0=10.**(dlog10(cbe0g)+kb2*deltb0p)
         endif
         
         do kjaf=njaf1,njaf2         
         rjaf=10.**(dlog10(rjafg)+kjaf*deljaf)
         if (kpro.eq.1.and.njaf.gt.1) then
          if (kjaf.lt.0)  rjaf=10.**(dlog10(rjafg)+kjaf*deltjam)
          if (kjaf.ge.0)  rjaf=10.**(dlog10(rjafg)+kjaf*deltjap)
         endif
         

         do kxm=nxm1,nxm2         
         xmstar=10.**(dlog10(xmstarg)+kxm*delxm)

         if (kpro.eq.1.and.nxmas.gt.1) then
          if (kxm.lt.0)  xmstar=10.**(dlog10(xmstarg)+kxm*delxmm)
          if (kxm.ge.0)  xmstar=10.**(dlog10(xmstarg)+kxm*delxmp)

         endif         

c         write(*,*) 'M/L of the BCG:', xmstar
         
 249     format(//'lambda, anisotropy =', f8.3,2x,f6.3)
 248     format(//'Y_1, anisotropy =', f8.3,2x,f6.3)
 247     format(//'Y_1, Y_2, anisotropy =', f8.3,2x,f8.3,2x,f6.3)
 246     format(//'Phi_infty, Q, anisotropy =', f8.3,2x,f8.3,2x,f6.3)

 259      format(//'anisotropy =', f6.3)

         do ictr=nrc1,nrc2
            rc=10.**(dlog10(rcg)+ictr*deltarc)
            if (kpro.eq.1.and.nrc.gt.1) then
             if(ictr.lt.0) rc=10.**(dlog10(rcg)+ictr*deltRcm)
             if(ictr.ge.0) rc=10.**(dlog10(rcg)+ictr*deltRcp)
            endif
            
            
            
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
                  !Flagship: assumes a specific value for rc
                  if (kfs2.eq.1) rc = 4 * rs
                  
                  write(*,*) 'with rc, rs: ', rc, rs
                  
               else
c
c     r_s as a free parameter
c     
                  rs=10.**(dlog10(rsg)+k*deltars)
                  if (kpro.eq.1.and.nrs.gt.1) then
                    if(k.lt.0) rs=10.**(dlog10(rsg)+k*deltrsm)
                    if(k.ge.0) rs=10.**(dlog10(rsg)+k*deltrsp)
                  endif
                  
               endif
c
c     a_ML forced = r_s
c
               if (kaml.eq.1) cbe=rs

c     
c     Light follows Mass (TLM)
c
               if (klfm.eq.1) then
                  rc=rs  ! assume N(R) and M(r) chosen with same model
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

               if (kmp.eq.11.or.kmp.eq.12) then
                xscCH=alimi
                alimit=alimi
                alimbc=0.0001
                call find_xsc(xscCH,alimit)

               endif
               
              !Environmental screening  - correct advanced modelling
               if (kmp.eq.17) then
                xlim = 1.0e-3
                xscCH = ScreenGen(xlim)
               endif
                                
               xfv(1)=rs
               xfv(2)=cbe

               call vmaxlik(nfv,xfv,fml)
                             
               nstar=0
               Nrun=Nclust
               if (nsingle.eq.1) then
                fml=fml*Nclust
   !if nsingle=1 abort the cycle over Nrun different clusters             
                Nrun=1 
              endif
!        read all the files for joint likelihood computation 
 901            format(a)
              filename='data/datphys_'
              dat='.dat'
              do irun=2+nstar,Nrun+nstar
               if(irun.lt.10) then
                ARR1 = ''
                WRITE(ARR1,"(I1)") irun
                open (12, FILE=filename//arr1//dat, STATUS='OLD')
               else
                ARR2 = ' '
                WRITE(ARR2,"(I2)") irun
                open (12, FILE=filename//arr2//dat, STATUS='OLD')
               endif
               read(12,901) line
               read(12,901) line
               j=0
               j0=-1
               dimin=1.e12
 404           continue
                
               read(12,*,end=171) dkpc,vkms,evkms 
               j=j+1
               dpp(j)=dkpc/1.e3
               vpp(j)=vkms
               epp(j)=evkms
               wpp(j)=1.0
                 
               if (dpp(j).lt.dimin) then
                dimin=dpp(j)
                j0=j
               endif
               goto 404
 171           continue
               close(12)
               npp=j
               jgiu=0
               do j=1,npp
                if (dpp(j).ge.rlowinmpc.and.
     &           dpp(j).le.rupinmpc) then
                 jgiu=jgiu+1
                 did(jgiu)=dpp(j)
                 viv(jgiu)=(vpp(j)-va)/(1.+va/clight)
                 eie(jgiu)=epp(j)
                 wiw(jgiu)=wpp(j)
                 if (did(jgiu).le.rlow) rlow=did(jgiu)
                 if (did(jgiu).ge.rup) rup=did(jgiu)
                endif
               
               enddo
               ngam=jgiu
              
               call vmaxSUM(nfv,xfv,fprova,did,viv,eie,wiw,ngam)
               fml=fml+fprova 
              enddo
   
               
            if (Nclust.eq.0) Nclust=1 !avoid problems in cluster number
                
            if (nlens.eq.1) then
             if(kmody.eq.1.and.kmp.ne.8) then
              fml=fml-Nclust*Plens(r200,rs)
              if (nlonly.eq.1) fml=-Plens(r200,rs)*Nclust
             endif
             !TO ADD BCG LENSING
             !First attempt - only rs and r200
             if (kmp.eq.17) then 
              fml2=fml2-Nclust*Plens(r200n,rsn)
              plen=-Plens(r200n,rsn)
              if (nlonly.eq.1) then
               fml2=-Nclust*Plens(r200n,rsn)
              endif
             endif
             if(kmp.eq.8) then
              if (naddlens.eq.0) then
               call Likelens_bh(nplens,plen)
              else
                plen = dlp
              endif
              fml=fml+Nclust*plen
              if (nlonly.eq.1) fml=plen*Nclust
             endif
            endif
            
            if (kbcg.eq.1.and.kbdt.eq.1) then
                rtbcg=20. !Max radius for the integration of the BCG integral
                call vdplik(chibcg2)
                fml=fml+chibcg2
            endif                
               
               
cc     If r200 changes, the number of objects selected changes
cc     as well, and one must scale the likelihood accordingly
c     [Obsolete, since we select galaxies within a radial range in Mpc
c      and not in units of r200 that may change]
cc               fml=fml*float(ngaref)/float(nga)
               scprova=screen
c               if(kmp.eq.9.and.nhone.gt.0) scprova=screen/(1+screen)

            if(kbcg.eq.1.or.kmp.eq.13) then
            
            write(iu60,173) r200,rc,rs,cbe,tmass,scprova,cbe0,rjaf,
     &      xmstar,fml,kani
           
            else
             write(iu60,673) r200,rc,rs,cbe,tmass,scprova,cbe0,fml,kani 
            endif
            
            
 673           format(8(f13.5,2x),i2)            
 173           format(10(f13.5,2x),i2)            
            
               if (fml.lt.fmlmingrid) then
                  fmlmingrid=fml
                  r200mingrid=r200
                  rsmingrid=rs
                  rcmingrid=rc
                  cbmingrid=cbe
                  cb0mingrid=cbe0
	  	          tmingrid=tmass
                  smingrid=screen
                  rjmingrid=rjaf
                  xmmingrid=xmstar
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
       enddo !rjaffe!!
       enddo !xmstar!!
      enddo
      if (kr200.lt.nr200/2) goto 727
      
      write(*,*) ' number of galaxies is ',nga
      
      fmlmin=fmlmingrid
      r200min=r200mingrid
      rsmin=rsmingrid
      rcmin=rcmingrid
      cbmin=cbmingrid
      cb0min=cb0mingrid
      tmmin=tmingrid
      scmin=smingrid
      rjmin=rjmingrid
      xmmin=xmmingrid
      rc=xfn1best
      nga=ngamin
      rlow=rlowmin
      rup=rupmin
      write(iu60,*) '.................................................'
      
      if(kbcg.eq.1.or.kmp.eq.13) then
       write(iu60,173) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,cb0min,rjmin,xmmin,fmlmin,kani
      else
       write(iu60,673) r200min,rcmin,rsmin,cbmin,tmmin,
     &                scmin,cb0min,fmlmin,kani
      endif
      write(iu60,*) '.................................................'
      
      
      
 891  continue

      if (kopt.lt.0) then          ! optimization skipped
         r200new=r200min
         rcnew=rcmin
         rsnew=rsmin
         cbnew=cbmin
         cb0new=cb0min
         scrnew=scmin
         rjafnew=rjmin
         xmsnew=xmmin
	     tmassnew=tmmin
         fmlb=fmlmin
      endif
      

      
      
      if(kbcg.eq.1.or.kmp.eq.13) then
       write(iu60,173) r200new,rcnew,rsnew,cbnew,
     & tmassnew,scrnew,cb0new,rjafnew,xmsnew,fmlb,kani
      else
       if(kfs2.eq.0) then
        write(iu60,673) r200new,rcnew,rsnew,cbnew,
     &  tmassnew,scrnew,cb0new,fmlb,kani
       else !this is for fs2 to save the likelihood of true r200
         write(iu60,673) r200g,rcnew,rsnew,cbnew,
     &  tmassnew,scrnew,cb0new,fml3,kani
       endif
      
      endif
      
      endif !END OF MCMC SAMPLING
      
      if (kbary.eq.1) close(92)
      
      call CPU_TIME(fine)
      

c      write(*,*) 'time of execution: ', fine-astar
      
      r200=r200new
      cbe=cbnew
      rs=rsnew
      rc=rcnew
      cbe0=cb0new
      screen=scrnew
      tmass=tmassnew
      dA3 = dA3new
      xmstar=xmsnew
      rjaf=rjafnew
      
      f=fmlb
      if (nlonly.eq.0) then
       write(*,292) r200,rc,r200/rc,rs,r200/rs,cbe,tmass,screen,dA3,
     &  cbe0,rjaf,xmstar,f 
       write(11,292) r200,rc,r200/rc,rs,r200/rs,cbe,tmass,screen,dA3,
     &  cbe0,rjaf,xmstar,f      
      else
       write(*,290) r200,rc,r200/rc,rs,r200/rs,cbe,tmass,screen,dA3,
     &  cbe0,rjaf,xmstar,f  
       write(11,290) r200,rc,r200/rc,rs,r200/rs,cbe,tmass,screen,dA3,
     &  cbe0,rjaf,xmstar,f        
      endif
 292  format(/' Best-fit from optimization ',
     &       /'   r_200    = ',f6.3,
     &       /'   r_tracer = ',f6.3,' (c_tracer= ',f7.2,')'
     &       /'   r_mass   = ',f6.3,' (c_mass=   ',f7.2,')'
     &       /'   Anisotropy parameter = ',f8.4,
     &       /'   First MG parameter =',f8.4,
     &       /'   Second MG parameter =',f8.4,
     &       /'   Third MG parameter = ',f8.4,
     &       /'   Second Anisotropy parameter =',f8.4,
     &       /'   BCG Jaffe radius =',f8.4,
     &       /'   BCG stellar M/L =',f8.4,
     &       /' Likelihood =',f13.5,/)
     
 290  format(/' Best-fit from optimization (kinematic only) ',
     &       /'   r_200    = ',f6.3,
     &       /'   r_tracer = ',f6.3,' (c_tracer= ',f7.2,')'
     &       /'   r_mass   = ',f6.3,' (c_mass=   ',f7.2,')'
     &       /'   Anisotropy parameter = ',f8.4,
     &       /'   First MG parameter =',f8.4,
     &       /'   Second MG parameter =',f8.4,
     &       /'   Third MG parameter = ',f8.4,     
     &       /'   Second Anisotropy parameter =',f8.4,   
     &       /'   BCG Jaffe radius =',f8.4,
     &       /'   BCG stellar M/L =',f8.4,
     &       /' Likelihood =',f13.5,/)
     
c     output results for N(R) in a file

 737  continue
      sigma0=0.
            
      write(iu30,330) knfit,rc,al !sigma0,
      
 330  format(i2,1x,e13.5,1x,f6.3) !e11.4,1x,

c     output file of probs

         if (kmp.eq.11.or.kmp.eq.12) then
            xscCH=alimi
            alimit=alimi
            alimbc=0.0001
            call find_xsc(xscCH,alimit)

         endif        
      
      
        !Environmental screening  - correct advanced modelling
         if (kmp.eq.17) then
            xlim = 1.0e-3
            xscCH = ScreenGen(xlim)
        endif      
      
        xfv(1)=rs
        xfv(2)=cbe

      call vmaxlik(nfv,xfv,fml)

      do jp=1,npv
         write(iu20,*) rpv(jp),vpv(jp),pv(jp)
      enddo
      close(iu20)
      close(11) !file of log
      
      r200=r200min
      rc=rcmin
      rs=rsmin
      cbe=cbmin
      cbe0=cb0min
      tmass=tmmin
      xmstar=xmmin
      rjaf=rjmin
      screen=scmin
      dA3 = dA3min
      f=fmlmin
      
 777  continue
 
      return
      end


      subroutine baryon_prior_penalty(pbar)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       include 'paramsoptS.i'
       ! Returns pbar = prior "penalty" (positive), in log-space
       ! We use lognormal-style prior: log(param) ~ N(log(mu), sigma_log)
       ! Here we omit additive constants (fine for MCMC).
       real*8 pbar, siglog, z, mu0, eps
       eps = 1.d-30
       pbar = 0.d0
       
       if (rho0gas.gt.0.d0.and.srhoga.gt.0.d0) then
        siglog = dsqrt(dlog(1.d0 + (srhoga/rho0gas)**2))
        mu0    = dlog(rho0gas)
        z      = (dlog(max(rhog1,eps)) - mu0)/siglog
        pbar   = pbar + 0.5d0*z*z
        
       endif

       if (rsgas.gt.0.d0.and.srsga.gt.0.d0) then
        siglog = dsqrt(dlog(1.d0 + (srsga/rsgas)**2))
        mu0    = dlog(rsgas)
        z      = (dlog(max(rsga1,eps)) - mu0)/siglog
        pbar   = pbar + 0.5d0*z*z
       endif

       if (rhostar.gt.0.d0.and.srhost.gt.0.d0) then
        siglog = dsqrt(dlog(1.d0 + (srhost/rhostar)**2))
        mu0    = dlog(rhostar)
        z      = (dlog(max(rhost1,eps)) - mu0)/siglog
        pbar   = pbar + 0.5d0*z*z
       endif

       if (rstar.gt.0.d0.and.srstar.gt.0.d0) then
        siglog = dsqrt(dlog(1.d0 + (srstar/rstar)**2))
        mu0    = dlog(rstar)
        z      = (dlog(max(rst1,eps)) - mu0)/siglog
        pbar   = pbar + 0.5d0*z*z
       endif
      end subroutine baryon_prior_penalty


!=======================================================================
!> @brief Propose new baryonic parameters using a random walk in log-space
!!
!> @details
!> This routine performs a symmetric Gaussian random walk in logarithmic
!> space for the baryonic nuisance parameters:
!>  - gas density normalization (rhog1)
!>  - gas scale radius (rsga1)
!>  - stellar density normalization (rhost1)
!>  - stellar scale radius (rst1)
!>
!> The proposal has the form:
!>   log(x_new) = log(x_old) + N(0, s_rw)
!> ensuring strict positivity of all parameters after exponentiation.
!>
!> This proposal is appropriate when the prior on the baryonic parameters
!> is defined in log-space (i.e. lognormal prior). If the prior is Gaussian
!> in linear space, this routine may lead to large prior penalties and low
!> acceptance rates.
!>
!> The proposed values are written directly to the global baryonic state
!> variables defined in 'paramsoptS.i'.
!>
!> @param[in,out] nseed Random number generator seed (must be updated externally)
!> @param[in]     s_rw  Standard deviation of the Gaussian step in log-space
!>
!> @note
!> Typical values for s_rw are in the range 0.01–0.10. Larger values may
!> cause runaway behaviour due to multiplicative jumps.
!>
!> @warning
!> This routine assumes that the current baryonic parameters are strictly
!> positive. A small numerical floor (eps) is enforced to avoid log(0).
!=======================================================================
        subroutine propose_baryons_rw_log(nseed, s_rw)
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
          include 'paramsoptS.i'
          ! Symmetric RW in log-space for each positive baryon param.
          ! s_rw is a *fractional* step scale (e.g. 0.05-0.20 typical).
          ! Use your r_normal_ab(mu, sigma, seed).
          real*8 s_rw, lg, eps
          integer*4 nseed
          eps = 1.d-30

          if (rho0gas.gt.0.d0) then
            lg   = dlog(max(rhog1,eps))
            lg   = r_normal_ab(lg, s_rw, nseed)
            rhog1 = dexp(lg)
          endif

          if (rsgas.gt.0.d0) then
            lg   = dlog(max(rsga1,eps))
            lg   = r_normal_ab(lg, s_rw, nseed)
            rsga1 = dexp(lg)
          endif

          if (rhostar.gt.0.d0) then
            lg    = dlog(max(rhost1,eps))
            lg    = r_normal_ab(lg, s_rw, nseed)
            rhost1 = dexp(lg)
          endif

          if (rstar.gt.0.d0) then
            lg   = dlog(max(rst1,eps))
            lg   = r_normal_ab(lg, s_rw, nseed)
            rst1 = dexp(lg)
          endif
        end subroutine propose_baryons_rw_log



!=======================================================================
!> @brief Propose new baryonic parameters using a random walk in linear space
!!
!> @details
!> This routine performs a Gaussian random walk directly in linear space
!> for the baryonic nuisance parameters:
!>  - gas density normalization (rhog1)
!>  - gas scale radius (rsga1)
!>  - stellar density normalization (rhost1)
!>  - stellar scale radius (rst1)
!>
!> The step size is expressed in units of the corresponding observational
!> uncertainty (sigma), i.e.:
!>   x_new = x_old + N(0, s_abs * sigma_x)
!>
!> Positivity is enforced by reflecting negative proposals using abs(x)
!> and applying a small numerical floor.
!>
!> This proposal is particularly suitable when the prior on the baryonic
!> parameters is Gaussian in linear space.
!>
!> The proposed values are written directly to the global baryonic state
!> variables defined in 'paramsoptS.i'.
!>
!> @param[in,out] nseed Random number generator seed (must be updated externally)
!> @param[in]     s_abs Step size in units of the prior standard deviation
!!
!> @note
!> Typical values for s_abs are in the range 0.1–1.0, corresponding to
!> 0.1–1 sigma jumps.
!>
!> @warning
!> The reflection via abs(x) slightly distorts the proposal near zero.
!> This effect is negligible when the probability mass is far from zero.
!=======================================================================
        subroutine propose_baryons_rw_lin(rhgn, rhosn,rsgn,rstn,
     &           nseed, s_abs)
          implicit real*8 (a-h,o-z)
          include 'paramsoptS.i'
          integer*4 nseed
          real*8 s_abs, xnew, eps
          eps = 0.01d0

          if (rho0gas.gt.0.d0 .and. srhoga.gt.0.d0) then
            xnew = r_normal_ab(rhgn, s_abs*srhoga, nseed)   ! step in unità di sigma
            rhgn = max(abs(xnew), eps)
          endif

          if (rsgas.gt.0.d0 .and. srsga.gt.0.d0) then
            xnew = r_normal_ab(rsgn, s_abs*srsga, nseed)
            rsgn = max(abs(xnew), eps)
          endif

          if (rhostar.gt.0.d0 .and. srhost.gt.0.d0) then
            xnew = r_normal_ab(rhosn, s_abs*srhost, nseed)
            rhosn = max(abs(xnew), eps)
          endif

          if (rstar.gt.0.d0 .and. srstar.gt.0.d0) then
            xnew = r_normal_ab(rstn, s_abs*srstar, nseed)
            rstn = max(abs(xnew), eps)
          endif
        end subroutine


C
C      ________________________________________________________
C     |                                                        |
C     |            SORT AN ARRAY IN INCREASING ORDER           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY OF NUMBERS                       |
C     |                                                        |
C     |         Y     --WORKING ARRAY (LENGTH  AT LEAST N)     |
C     |                                                        |
C     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --SORTED ARRAY                           |
C     |________________________________________________________|
C      From napack of Netlib
      SUBROUTINE SORTP(X,N)
      real*8 X(N),Y(N),S,T
      INTEGER I,J,K,L,M,N
      I = 1
10    K = I
20    J = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      Y(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      Y(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      S = X(I)
      J = Y(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      T = X(J)
      L = Y(J)
      X(I) = L
60    IF ( S .GT. T ) GOTO 70
      Y(M) = S
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      S = X(I)
      GOTO 60
70    Y(M)= T
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      T = X(J)
      GOTO 60
80    Y(M) = T
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      Y(M) = X(M+I)
      GOTO 90
100   X(I) = J
      L = J
110   Y(M) = S
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = X(I)
140   X(I) = Y(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      Y(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .EQ. 1 ) RETURN
      GOTO 40
      END
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


c
c     Subroutine for the computation of mean velocity,
c     velocity dispersion and r200 (in Mpc, h=0.7)
c
      subroutine kine(d,v,n,va,sv,rv,eva,esv,jack)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (cl=2.99792458e5)
      dimension d(n),v(n),x(25000),y(25000),ww(25000)
c
c     Assume velocities are rest-frame (cosmologically corrected already)
c
      nm=n
      do j=1,n
         x(j)=d(j)
         y(j)=v(j)
      enddo
c
c     Use Gapper (ibwt=0; set ibwt=1 for biweight)
c
      ibwt=0
      call robusti(y,nm,ibwt,ave,sig)
      va=va+ave
      sv=sig
c
c     Jacknife error estimate, if requested
c
      if (jack.gt.0.5) then
         rn=float(nm)
         n1=nm-1
         rn1=float(n1)
         zall=sv*rn
         zstar=0.0
         zstar2=0.0
         do k=1,nm
            l=0
            do j=1,nm
               if (j.ne.k) then
                  l=l+1
                  ww(l)=y(j)
               endif
            enddo
            call robusti(ww,n1,ibwt,avej,sigj)
            xstar=xstar+xall-rn1*avej
            zstar=zstar+zall-rn1*sigj
            zstar2=zstar2+(zall-rn1*sigj)*(zall-rn1*sigj)
         enddo
         sstarz=(zstar2-zstar*zstar/rn)/(rn*rn1)
         dof=float(nm-1)
         tprob=.68
         print *,' jack = 0 is unavailable'
         stop
         esv=sqrt(sstarz)*t68
         eva=t68*sv/sqrt(rn)
      else
         eva=0.
         esv=0.
      endif
      return
      end
c
c
      
      subroutine robusti(xin,n,ibwt,c,s)

C     Routine that uses robust techniques to estimate the central location C
C     and the spread S, for a distribution of N sorted values X.
C     Based on the work of Beers, Flynn and Gebhardt, AJ 100, 32, 1990.

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

      call sortp(x,n)
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
      call sortp(xred,n)
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
c     N(R) - projected NFW; use eq. (41), (42) in Lokas & Mamon (2001)
!>   \ingroup MAM
      function sigmar1(tt)
!>    @authors A. Biviano, G. Mamon, G. Boue'
!>    @details N(R) - NFW-model. This function computes the projected
!!    cumulative number profile for a NFW model which is used to calculate
!!    the radial velocity dispersion profile (VDP). It is based on eq. (41)
!!    and (42) of Lokas & Mamon (2001)
!>    It requires the inclusion of the external file paramsoptS.i
!!    and the values of the number density scale radius rc and of the  
!!    virial radius of the mass model, r200.

!>   Note that the normalization is not considered as it simplifies in
!!   the equation of the VDP.

!>   @param[in] tt  (projected) radial distance at which the profile is 
!!   computed, in unit of Mpc. 


!>  Example of usage:

!>  \code{.f}
!>  rc=0.3
!>  r200=1.41
!>  tt=0.5 
!>  tnNFW=sigmar1(tt)  
!>  write(*,*) tnNFW
!>  \endcode

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
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
!>   \ingroup MAM
      function sigmar2(tt)
!>    @authors A. Biviano, G. Mamon, G. Boue'
!>    @details N(R) - Hernquist model. This function computes the projected
!!    cumulative number profile for a Hernquist model which is used to calculate
!!    the radial velocity dispersion profile (VDP). See Hernquist (1990).

!>    It requires the inclusion of the external file paramsoptS.i
!!    and the values of the number density scale radius rc.

!>   Note that the normalization is not considered as it simplifies in
!!   the equation of the VDP.

!>   @param[in] tt  (projected) radial distance at which the profile is 
!!   computed, in unit of Mpc.

!>  Example of usage:

!>  \code{.f}
!>  rc=0.3
!>  tt=0.5 
!>  tnHer=sigmar2(tt)  
!>  write(*,*) tnHer
!>  \endcode
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
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


!>    \ingroup MAM
      function sigmar3(tt)
!>    @authors A. Biviano, G. Mamon, G. Boue'
     
!>    @details N(R) - beta-function. This function computes the projected
!!    cumulative number profile for a beta-model which is used to calculate
!!    the radial velocity dispersion profile (VDP).
!>    It requires the inclusion of the external file paramsoptS.i
!!    and the values of the number density scale radius rc and of the exponent
!!    of the beta-model al.

!>   Note that the normalization is not considered as it simplifies in
!!   the equation of the VDP.

!>   @param[in] tt  (projected) radial distance at which the profile is 
!!   computed.   

!>  Example of usage:

!>  \code{.f}
!>  rc=0.3
!>  tt=0.5 
!>  tnbeta=sigmar3(tt)  
!>  write(*,*) tnbeta
!>  \endcode
      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptS.i'
      t=tt/rc
      sigmar3=(1.d0+t*t)**al
      return
      end


!>   \ingroup MAM
      function sr2int(alr)
!>     @authors A. Biviano, G. Mamon, G. Boue', L. Pizzuti 
!>     @details Integrand function for the determination of
!!     \f$ \nu(r)*\sigma_r^2(r) \f$, given beta(r) and M(r)
!!     from eqs. A3, A4, A6 in Mamon & Lokas (2005)
!!     or from eq.(A3) in Mamon, Biviano & Boue' (2013). It includes the
!!     modified gravity parametrizations of the mass profile described in
!!     Pizzuti et al., 2021

!>    It requires the inclusion of the COMMON parameters defined in paramsoptS.i,
!! the values of all the paramters r200, rs, rc, cbe, cbe0, 
!! tmass, screen, rcut, al, kani, kmp, knfit, v200=10.*hz*r200.
!> A special case corresponds to the number of steps in the anisotropy
!! parameter nbs=-1, where the Hansen & Moore model is selected
!! (the anisotropy depends on the mass density profile \f$ \beta(r)=
!! a+b d\rho(r)/d\ln r \f$). In this case,
!! the parameter ahm=a and bhm=b should be specified externally.

!>   @param[in] alr  natural logarithm of the radial distance at which the profile is 
!!   computed.  


!>  Example of usage:
!>
!>
!>  \code{.f}
!>       nbs=1
!> !     Values of the parameters of Hansen & Moore's (2006)
!> !     beta = a + b dln(rho)/dln(r) relation
!>       ahm=-0.15
!>       bhm=-0.192
!> 
!>       r200=1.41
!>       rc=0.3
!>       rs=0.3
!>       cbe=1.2
!>       cbe0=1.0 !not used unless kani=21, 41
!>       tmass=0.10
!>       screen=0.4
!>       kmp=8    !mass profile model. kmp=7 corresponds to mNFW_BH
!>       v200=r200*10*85
!>       knfit=1  !Number density profile model. knfit=1 Corresponds to pNFW
!>       kani=4   !Anisotropy model. kani=4 corresponds to Tiret profile
!>       kscr=1   !screening option
!>       rcut=1.0
!>       alr=dlog(2.5d0) 
!>       sint=sr2int(alr)  
!>       write(*,*)  sint
!>  \endcode
!>
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      complex gamma,z,hypgeo,z200m, hfz
      dimension rsvalues(28),r100values(28)
      include 'paramsoptS.i'
      include 'barions.i'
      data rsvalues/0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,
     ,     0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,
     ,     0.30,0.35,0.40,0.45,0.50,1.00/
      data r100values/1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,
     ,     1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,
     ,     1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476/
      external gammarec
      external frLin
      external hypgeo
      external fbars
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
            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200
         endif


c     attempt to insert generalized NFW model         
      elseif (kmp.eq.10) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)

       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
       
       xm =  hfx/hfy*(t/r200)**(3-dA3)
c     Clustering Dark Energy 
c     dA3 = exponent of gNFW
c     dA2 = f(cs,w)
c     dA1 = cs2 (in unit of c**2)        
      elseif (kmp.eq.18) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)

       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
                                                 
       xm=hfx/hfy*(t/r200)**(3-dA3)*(1.0d0+3*10**(screen-3.5)*(1+tmass))

      
      elseif (kmp.eq.8)  then	
c     M(r) is the modified NFW profile with a beyond Horndeski model 
c     See Saltas+16, Sakestein+16
c     the free parameter tmass is now the coupling constant Y 
	
       fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
       abey=tmass*t*t*(rs-t)/(rs+t)**3
       xm=((dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))+abey/4)*fac200
        
      elseif (kmp.eq.9) then
c     M(r) is the modified NFW profile in chameleon screening 
c      (Pizzuti et al.,2021)    
       
       !fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
       !xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200+dphidr(t/rs)
       xm = tM_bnfw_CS(t)
       
      
      elseif (kmp.eq.15) then
c     M(r) is the modified gNFW profile in chameleon screening 
c     gamma for now is fixed to eim
        xm = tM_gnfw_CS(t)
        
      elseif (kmp.eq.16) then
c     M(r) is the modified Isothermal profile in chameleon screening 
c     the exponent is 3/2       
       xm = tM_iso_CS(t)

      elseif (kmp.eq.17) then 
      !M(r) is multicomponent general Chameleon. gNFW+Jaffe+beta+NFW
        xm = grav*tmass_Gen(t)/gm200 ! To divide by gm200 as this is already
                                ! including the normalization factor
        
        
       !Check the dA3
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
c
c
      elseif (kmp.eq.11) then
c        mass is modified Burkert in Chameleon Gravity   
         trs=t/rs
         rvrs=r200/rs
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))
         xm1=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))
         xm=xm1+dphib(trs)
         
      elseif (kmp.eq.12) then
c     mass is modified Einasto in Chameleon Gravity
      !   eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm1=fac200*gip
         
         xm=xm1+dphiE(t/rs)            
                  
         
      elseif (kmp.eq.5) then
c
c     M(r) is Soft Isoth Sph; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         fac200=1./(r200-rs*datan(r200/rs))
         xm=fac200*(t-rs*datan(t/rs))
      
      elseif (kmp.eq.13) then
c     M(r) is Refracted Gravity; no DM contribution allowed for the time begin
c     only BCG parameters kept free
c
        !refracted gravity with the function of Valeria
        xm= grav*rgmass(t)/gm200
        !fbarr=0.0
        !if (kbary.eq.1) fbarr=fbars(t)
        
        !Covariant formulation. It is not good for the moment
        !xm=grav*(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)/gm200
      
      elseif (kmp.eq.14) then
      
        !M(r) is boson star + gas and galaxies
        fbarr=0.0
        if (kbary.eq.1) fbarr=fbars(t)
         
        xt = t/rs
        if (xt.le.2.0d0) then 
         xm = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
        else
         xm =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
        endif 
        xm = xm + grav*fbarr/gm200
      
      else
c
c     M(r) is Einasto ; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
       !  eim=5.

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

C******************* BCG ***********************************************      
 
c  test

      fbarr=0.0
      if (kbary.eq.1) fbarr=fbars(t)
      
      if (kbcg.eq.1) then

       !GR setup
       if(kmp.eq.10.or.kmp.eq.18.or.kmp.le.6) then !attempt to add clustering DE
        fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
        xm=xm+fjaffe !inclusion of gas,galaxies and bcg
        
        
        !Chameleon setup
       !elseif (kmody.eq.1) then 
        !dpj=dphiJ(t/rjaf) !Include the chameleon field of the BCG
      !fjaffe=grav*(xmstar*xlumbcg*(t/rjaf/(1.+t/rjaf)+dpj)+fbarr)/gm200
        !xm=xm+fjaffe !inclusion of gas,galaxies and bcg
        
        !Refracted gravity setup
       elseif (kmp.eq.13) then
       
       !refracted gravity with the function of Valeria
        xm= grav*rgmass(t)/gm200
        
        !Covariant formulation. It is not good for the moment
        !the mass is all given by the BCG and gas+baryons, no DM
        !xm=grav*(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)/gm200
       elseif (kmp.eq.14) then
        !Boson star + BCG
        fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
        
        xt = t/rs
        if (xt.le.2.0d0) then 
         xm = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
        else
         xm =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
        endif 
        
        xm = fjaffe + xm
        
        
       endif 
       
        
        
      endif
c***********************************************************************
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
c     Generalized Osipkov Merritt
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.21) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         bec=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         if (dabs((2*bec0-2)*dlog10(t)).gt.60) then
          write(*,*) 'problem in cbe0 = ', cbe0
          stop('ABORTED: value of anisotropy makes overflow')
         endif   
         if (dabs(((bec-bec0))*dlog10(t+rm2)).gt.60) then
          write(*,*) 'problem in cbe = ', cbe, 'and cbe0=', cbe0
          stop('ABORTED: value of anisotropy makes overflow')
         endif
         
         
         sr2int=xnu*xm*t**(2.*bec0-2.)*(t*t+rm2*rm2)**(bec-bec0)
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
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2-dA3)*rs
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(2.*bec)/(t*t)
c
c     Generalized Tiret (with additional parameter cbe0)         
c         

      elseif (kani.eq.41) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         if(ifreeBeta.gt.0) rm2 = rbeta
         
         bec=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values

         if (dabs((2*bec0-2)*dlog10(t)).gt.60) then
          write(*,*) 'problem in cbe0 = ', cbe0
          stop('ABORTED: value of anisotropy makes overflow')
          
         endif 
         if (dabs((2.*(bec-bec0))*dlog10(t+rm2)).gt.60) then
          write(*,*) 'problem in cbe=', cbe, 'and cbe0=', cbe0
          stop('ABORTED: value of anisotropy makes overflow')
         endif
         !if ((bec-bec0).lt.-6.) bec=-5.999d0+bec0
         
         sr2int=xnu*xm*t**(2*bec0-2.)*(t+rm2)**(2.*(bec-bec0))

      elseif (kani.eq.8) then
c     WORK IN PROGRESS
c     Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values

         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free 

        
        db1 =dexp(-(bi-b0)*dexp(-(t**2)/(rm2**2)))
        db2 = t**(2.0*b0-2)*(t+rm2)**(2.0*bi-2.0*b0)
        
        if (dabs((2.*(bec-bec0))*dlog10(t+rm2)).gt.60) then
          write(*,*) 'problem in cbe=', cbe, 'and cbe0=', cbe0
          stop('ABORTED: value of anisotropy makes overflow')
        endif
        
        sr2int=xnu*xm*db1*db2

        
c fare poi un check!
         
c
c     modified Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2int=xnu*xm*(t+rm2)**(4.*bec)*t**(-2.-2.*bec)
c
c     Hansen+Moore relation for NFW mass profile
c
      else
         rhm=rc ! use nu(r)
cc         rhm=rs ! use rho(r)
         sr2int=xnu*xm*t**(2.*(ahm-bhm-1.))*(rhm+t)**(-4.*bhm)
      endif
      sr2int = t*sr2int*gm200    ! integral in dlog and G*M scaling
      
      
      !write(*,*) sr2int
c     
      return
      end

c **********************************************************************
!>   \ingroup MAM
      function sr2out(t)
!>     @authors A. Biviano, G. Mamon, G. Boue',
!>     @details This function computes the factor outside the integral for the
!!     \f$ \sigma_r^2 \f$ formula. 

!>    It requires the inclusion of the COMMON parameters in file paramsoptS.i,
!! the values of r200, rc, cbe, cbe0,
!! al, kani, knfit.
!> 
!> A special case corresponds to the number of steps in the anisotropy
!! parameter nbs=-1, where the Hansen & Moore model is selected
!! (the anisotropy depends on the mass density profile \f$ \beta(r)=
!! a+b d\rho(r)/d\ln r \f$). In this case,
!! the parameter ahm=a and bhm=b should be specified externally.
!>
!>   @param[in] t radial distance at which the profile is 
!!   computed.  
!>
!>  Example of usage:
!>
!>
!>  \code{.f}
!>      nbs=1
!>      !Values of the parameters of Hansen & Moore's (2006)
!>      !beta = a + b dln(rho)/dln(r) relation
!>      ahm=-0.15
!>      bhm=-0.192 
!>      r200=1.41
!>      rc=0.3 
!>      cbe=1.2
!>      cbe0=1.0 !not used unless kani=21, 41
!>      knfit=2  !number density profile model. knfit=2 corresponds to pHer
!>      kani=4  !velocity anisotropy model. kani=4 corresponds to Tiret profile
!>  
!>      alr=dlog(0.5) 
!>      sout=sr2out(t)  
!>      write(*,*) sout
!>  \endcode
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
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

c     Generalized Osipkov Merritt
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.21) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         bec=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         sr2out=t**(-2.*bec0)*(t*t+rm2*rm2)**(-1.*(bec-bec0))


         
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
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-2.*bec)

c
c     generalized Tiret (non-zero central anisotropy  as cbe0 parameter)
c
         
         
      elseif (kani.eq.41) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta
         
         bec=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         sr2out=t**(-2.*bec0)*(t+rm2)**(-2.*(bec-bec0))
         
      elseif (kani.eq.8) then
c     WORK IN PROGRESS
c     Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free 
        db1 =dexp((bi-b0)*dexp(-(t**2)/(rm2**2)))
        db2 = t**(-2.0*b0)*(t+rm2)**(-2.0*bi+2.0*b0)
        sr2out = db1*db2
        
c
c     modfied Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=1.-1./(cbe*cbe)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
         sr2out=(t+rm2)**(-4.*bec)*t**(2.*bec)
      else
c
c     Hansen+Moore relation for NFW mass profile
c
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
      real*8 rjl(1), yrjl(1), ydev1(1), ydev2(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
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
      call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  ydev1,ydev2)
      
c      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)

c      if (ier .eq. 34) then
c         print *,' ICSEVU in GWEN: max xris=', xris(ninterp), ' r=',rjl
c      endif
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

c     Generalized Osipkov Merritt
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.21) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t*t/(t*t+rm2*rm2)        
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         
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
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=(1.-1./(cbe*cbe))*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values

c
c     Generalized Tiret 

      elseif (kani.eq.41) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t/(t+rm2)   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values


      elseif (kani.eq.8) then
c     WORK IN PROGRESS
c     Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free 
         bec = b0+(bi-b0)*(t/(t+rm2)+t**2/rm2**2*dexp(-(t/rm2)**2) )
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c
      elseif (kani.eq.5) then
c     modified Tiret (non-zero central anisotropy)
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
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
      real*8 alsigma(100), alrvec(100), ydev1(100),ydev2(100)
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
c      call icsevu(xris,yris,ninterp,csr,icsr,alrvec,alsigma,n,ier)
      call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,alrvec,alsigma,
     &  ydev1,ydev2)
      
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
      real*8 rjl(1), yrjl(1), ydev1(1), ydev2(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
      include 'sr.i'
      include 'vlos.i'
      include 'tsallis.i'
      external gammarec
c
      t = xmin*dcosh(u)

cc       write(*,*) ' vel and error is',vj,ej

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
c      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  ydev1,ydev2)
c      if (ier .eq. 34) then
c         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
c     &    rjl(1)
c      endif
      sigr=dexp(yrjl(1)) !check

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

c     Generalized Osipkov Merritt
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.21) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t*t/(t*t+rm2*rm2)        
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values         
         
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
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=(1.-1./(cbe*cbe))*t/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         
c     Generalized Tiret
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.41) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t/(t+rm2)        
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values         

      elseif (kani.eq.8) then
c     WORK IN PROGRESS
c     Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values
                
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free 
         bec = b0+(bi-b0)*(t/(t+rm2)+t**2/rm2**2*dexp(-(t/rm2)**2) )
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
c
c     modified Tiret (non-zero central anisotropy)
c
      elseif (kani.eq.5) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         
         bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values
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

      sigmaz=dsqrt(sigmaz*sigmaz+ej*ej)

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


!>    \ingroup MAM
c     Subroutine for the computation of Max Lik
c     
      subroutine vmaxlik(nfv,xfv,f)
!>    @authors A. Biviano, G. Mamon, G. Boue', L. Pizzuti
!>    @details this subroutine compute the value of the likelihood \f$ -\ln \mathcal{L}(\theta) \f$,
!!    where \f$\theta \f$ is the vector of values for the model parameters.
!>    It requires all the values of the input parameters in `pars_test.txt`,
!!    and in `Options.txt`, as well as the COMMON blocks:
!>     `paramsoptS.i`,
!>     `sr.i`,
!>     `vlos.i`,
!>     `datarv.i`
!>     `probs.i`.
!>
!!    Data of projected positions (*in Mpc*), velocities (*in km/s*) and 
!!   and velocity uncertainties (*in km/s*) should be passed as a common set of  vectors 
!!    `r(ndata), v(ndata) e(ndata), w(ndata)`, where `ndata` is the total number of data points
!!   and `w(ndata)` represents the weights (assumed to be 1 in general)
!>
!>    @param[in] nfv=2 dimension of the vector xfv
!>    @param[in] xfv=(rs,cbe) two-dimensional vector with the values of rs and cbe
!>    @param[out] f           value of the (-)log likelihood    



      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension xfv(nfv),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension y2a(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),wmbm(25000),embm(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      dimension dlogp2(1), xi(1), yi(1)
      include 'paramsoptS.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'probs.i'
      include 'barions.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
       
       
      !generate the value or rho0gas and rsgas
      !According to their uncertainties
      
c      if (rho0gas.gt.0) then
c        call CPU_TIME(start) !uses cpu time for the random generator seed
c        nst=int(134023213*start, kind=8)
c        nseed = 123434789+((-1)**nst)*nst
c        rhog1 = dabs(r_normal_ab (rho0gas,srhoga,nseed) )
      
c        coo = 0
c        do while (rhog1.le.0)
c            coo = coo + 1
c            nseed = 123434789+((-1)**nst)*(nst+coo)
c            rhog1 = dabs(r_normal_ab(rho0gas,srhoga,nseed)) 
c        end do
c      endif      
      
c      if (rsgas.gt.0) then
c        call CPU_TIME(start) !uses cpu time for the random generator seed
c        nst=int(14451231*start, kind=8)
c        nseed = 123462714+(-1)**nst*nst
       
c        rsga1 =  dabs(r_normal_ab (rsgas,srsga,nseed) )
c        coo = 0
c        do while (rsga1.le.0)
c            coo = coo + 1
c            nseed = 123434789+((-1)**nst)*(nst+coo)
c            rsga1 =  dabs(r_normal_ab (rsgas,srsga,nseed) )
c        end do
      
c      endif
c      !check
c      !write(*,*) rhog1, rsga1, coo
      
c      if(rhostar.gt.0) then
c        call CPU_TIME(start) !uses cpu time for the random generator seed
c        nst=int(1343000*start, kind=8)
c        nseed = 123434789+((-1)**nst)*nst
c        rhost1 = dabs(r_normal_ab(rhostar,srhost,nseed))
    
c        coo = 0
c        do while (rhost1.le.0)
c            coo = coo + 1
c            nseed = 123434789+((-1)**nst)*(nst+coo)
c            rhost1 = r_normal_ab(rhostar,srhost,nseed) 
c        end do
c      endif
      
c      if (rstar.gt.0) then
c        call CPU_TIME(start) !uses cpu time for the random generator seed
c        nst=int(144231325*start, kind=8)
c        nseed = 123462714+(-1)**nst*nst
       
c        rst1 =  dabs(r_normal_ab ((rstar),srstar,nseed) )
c        coo = 0
c        do while (rst1.le.0)
c            coo = coo + 1
c            nseed = 123434789+((-1)**nst)*(nst+coo)
c            rst1 =  dabs(r_normal_ab ((rstar),srstar,nseed) ) 
c        end do
c      endif

      !check
      !write(*,*) rhost1, rst1, coo
c      do i = 1, 100
c       xt = i*0.025
c       write(*,*) xt, fbars(xt)
c      enddo
c      stop('check')

      

c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      if (dabs(al).lt.0.0001) n1=1
      rs=xfv(1)
      cbe=xfv(2)

      ngambm=nga
      do j=1,nga
         rmbm(j)=r(j)
         vmbm(j)=v(j)
         embm(j)=e(j)
         wmbm(j)=w(j)
      enddo
      irule=2
      errabs=0.0d0
      errrel=0.001d0
      rismin=1.d-190
      rismax=1.d190
      rinfinity=25.0d0
c      rlow=0.051

      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
         xx1=xmin 

         call dgaus8 (sr2int,xris(i),dlog(2.*rinfinity), errrel, 
     &    risl, IERR)

         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)
         !write(*,*) xris(i),yris(i)

      enddo

c
c compute spline coeffs for later interpolation of sigma_r
c
      call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
c      call icsccu(xris,yris,ninterp,csr,icsr,ier)
      
c      if (ier .gt. 128) then
c         print *,' in VMAXLIK ...'
c         do i = 1, ninterp
c            print *,' i xris yris = ', i, xris(i), yris(i)
c         enddo
c      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.
      psum2=0.
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      !if (kintd.eq.1) then
       !  rinfinity=0.99999*r200 
      !endif

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
         k=0
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
               !if (kintd.eq.0) call rinf(rinfinity)
               
               !NOTE: to speed up the code, I keep rinfinity at 25 Mpc
               !otherwise it will be way slower
               !call rinf(rinfinity)
c
c
               umax = dacosh(rinfinity/rj)
               errrel=0.001d0
               
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
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

                  call dgaus8 (sdint,xmin,xmax, errrel, 
     &                xinteg, IERR)  
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
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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

c     

c         call ibcccu(z,rp,nr,vp,nv,c,ic,wk,ier)
         call splicoff2(rp,vp,z,nr,nv,y2a)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
            xi(1)=xx
            yi(1)=yy
c            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp4,ier)
            call spli2d(rp,vp,z,y2a,nr,nv,xx,yy,dlogp3)
            psum=psum-dlogp3*wmbm(l)
c            psum2=psum2-dlogp4*wmbm(l)
         enddo

      else
      
c     Use all the original data, no interpolation

         do j=1,ngambm
            rj=rmbm(j) !dati
            vj=vmbm(j)
            ej=embm(j)
            xmin=rj

c     re-define rinfinity as the radius where v_z/sigma_r leq 0.1
c     to avoid a problem in v_z=0; if we require to consider the 
c     interloper contribution, rinfinity is set to 1 already


            !for now no recomputation of rinfinity
            !if (kintd.eq.0) call rinf(rinfinity)

            
            
c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)
               umax = dacosh(rinfinity/rj)
	       errrel=0.001d0
            call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif


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
c                  xinteg=dcadre(sdint,xmin,xmax,errabs,
c     &                 errrel,errest,ier)
                  call dgaus8 (sdint,xmin,xmax, errrel, 
     &                xinteg, IERR)  
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
               
               else !sono qui

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

c      write(*,*) psum
      return
      end
      
      
      subroutine vmaxSUM(nfv,xfv,f,rr,vv,ww,ee,ngam)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension xfv(nfv),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension y2a(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),wmbm(25000),embm(25000)
      dimension rr(25000),vv(25000),ww(25000),ee(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptS.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'probs.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
c      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      if (dabs(al).lt.0.0001) n1=1
      rs=xfv(1)
      cbe=xfv(2)

      ngambm=ngam
      do j=1,ngam
         rmbm(j)=rr(j)
         vmbm(j)=vv(j)
         embm(j)=ee(j)
         wmbm(j)=ww(j)
      enddo

      irule=2
      errabs=0.d0
      errrel=0.001d0
      rismin=1.d-190
      rismax=1.d190
      rinfinity=25.0d0
c      rlow=0.051
      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin
        call dgaus8 (sr2int,xris(i),dlog(2.*rinfinity),errrel,
     &   risl, IERR) 

c         risl = dcadre(sr2int,xris(i),dlog(2.*rinfinity),errabs,errrel,
c     &    errest,ier)
c         if (ier .gt. 128) then
c            print *,'VMAXLIK-SR2INT: rmin=',xx1,' rmax=',xx2
c            stop
c         endif
         risok=dsqrt(risl*sr2out(xmin))
c	     write(*,*) (xris(i)), risok
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)

      enddo
c      write(*,*) 'END'

c
c compute spline coeffs for later interpolation of sigma_r
c
      call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
c      call icsccu(xris,yris,ninterp,csr,icsr,ier)
c      if (ier .gt. 128) then
c         print *,' in VMAXLIK ...'
c         do i = 1, ninterp
c            print *,' i xris yris = ', i, xris(i), yris(i)
c         enddo
c      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      !if (kintd.eq.1) then
       !  rinfinity=0.99999*r200
      !endif

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
         k=0
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
               !if (kintd.eq.0) call rinf(rinfinity)
               call rinf(rinfinity)
c               
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c
               umax = dacosh(rinfinity/rj)
               errrel=0.001d0
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c     
c     Add interloper contribution if requested
c
c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.005d0
cc                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
cc     &                 errest,ier)
c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR)  
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
ccc                  write(*,*) ' Number of obsd galaxies ',ngambm
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
ccc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c              endif
               
               if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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

         call splicoff2(rp,vp,z,nr,nv,y2a)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
            call spli2d(rp,vp,z,y2a,nr,nv,xx,yy,dlogp3)
            psum=psum-dlogp3*wmbm(l)
 
         enddo

      else

c     Use all the original data, no interpolation

         do j=1,ngambm
            rj=rmbm(j) !dati
            vj=vmbm(j)
            ej=embm(j)
            xmin=rj

c     re-define rinfinity as the radius where v_z/sigma_r leq 0.1
c     to avoid a problem in v_z=0; if we require to consider the 
c     interloper contribution, rinfinity is set to 1 already

            !if (kintd.eq.0) call rinf(rinfinity)
            call rinf(rinfinity)
c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)
               umax = dacosh(rinfinity/rj)
	       errrel=0.001d0
           call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif


c     Add interloper contribution if requested

c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4.+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  n2=2
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.001d0
cc                  xinteg=dcadre(sdint,xmin,xmax,errabs,
cc     &                 errrel,errest,ier)
c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR) 
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c                  g=g+gi
c               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally

               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
                  p=2.*pig*rj*g/fnorm
               
               else !sono qui

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

c      write(*,*) psum
      return
      end
      
      
      
c
c     Max Lik for beta-profile
c
      function fcn3(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptS.i'
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
      include 'paramsoptS.i'
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
      include 'paramsoptS.i'
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
      real*8 rjl(1), yrjl(1),yy1(1), yy2(1)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
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
c         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
         call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  yy1,yy2)
         
         
         
c      if (ier .eq. 34) then
c         print *,' ICSEVU in RINF: max xris=', xris(ninterp), ' r=',rjl
c      endif
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
            
c     Generalized Osipkov Merritt
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.21) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t*t/(t*t+rm2*rm2)        
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values            
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
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !  avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            
            bec=(1.-1./(cbe*cbe))*t/(t+rm2)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values

c     Generalized Tiret
c     In this case there is an additional free parameter which is 
c     the inner anisotropy for r=0, cbe0
         
      elseif (kani.eq.41) then
         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
         !avoid unphyiscal values 
         if (rm2.lt.0) rm2=1e-4
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
         
         becin=1.-1./(cbe*cbe)
         bec0=1.-1./(cbe0*cbe0)
         if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
         if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
         bec=bec0+(becin-bec0)*t/(t+rm2)        
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values   

      elseif (kani.eq.8) then
c     WORK IN PROGRESS
c     Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values

         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free 
         bec = b0+(bi-b0)*(t/(t+rm2)+t**2/rm2**2*dexp(-(t/rm2)**2) )
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values


c
c     modified Tiret (non-zero central anisotropy)
c
         elseif (kani.eq.5) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
            !avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            
            bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
            if (bec.le.-1.) bec=-0.999d0 ! avoid unphysical values
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
!>    \ingroup MAM
      function fa(tlog)
!>     @authors A. Biviano, G. Mamon, G. Boue', L. Pizzuti 
!>     @details Integrand function for the determination of
!!     \f$ N(R)*\sigma_{LOS}^2(R) \f$, given beta(r), M(r) and \f$ \nu(r) $\f
!!     from eqs. A15, A16 in Mamon & Lokas (2005). It includes the
!!     modified gravity parametrisations of the mass profile described in
!!     Pizzuti et al., 2021

!>    It requires the inclusion of the external file paramsoptS.i,
!! the values of all the paramters r200, rs, rc, cbe, 
!! tmass, screen, rcut, al, kani, kmp, knfit, kscr, v200=10.*hz*r200.
!> A special case corresponds to the number of steps in the anisotropy
!! parameter nbs=-1, where the Hansen & Moore model is selected
!! (the anisotropy depends on the mass density profile \f$ \beta(r)=
!! a+b d\rho(r)/d\ln r \f$). In this case,
!! the parameter ahm=a and bhm=b should be specified externally.
!> NOTE: this function cannot be used as it is, unless a previous computation of spline coefficients for the radial 
!! velocity dispersion profile \f$ \sigma_r^2(r) \f$
!! is done. Indeed fa(tlog) include a cubic spline interpolation of \f$ \sigma_r^2(r) \f$. The spline
!! is computed by using the routines by John Burkardt, spline_cubic_set and spline_cubic_val (see spline_cubic.f90 in folder GamI 
!! for source codes, credits and usage).

!>   @param[in] tlog  natural logarithm of the projected radial distance at which the profile is 
!!   computed.  


!>  Example of usage:
!>
!>
!>  \code{.f}
!>
!> 
!>  !compute spline coeffs for later interpolation of sigma_r
!>  call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
!>  !xris(ninterp) and yris(ninterp) are global variables storing the 
!>  !log-radii and the values of sigma^2_r computed before.
!>  nbs=1
!> !  Values of the parameters of Hansen & Moore's (2006)
!> !     beta = a + b dln(rho)/dln(r) relation
!>  ahm=-0.15
!>  bhm=-0.192
!> 
!>  r200=1.41
!>  v200=r200*10*85
!>  rc=0.3
!>  rs=0.3
!>  cbe=1.2
!>  tmass=1.0
!>  screen=0.4
!>  kmp=7
!>  knfit=1
!>  kani=4
!>  rcut=1.0
!>  
!>
!>  tlog=dlog(0.5d0) 
!>  sVDP=fa(tlog)  
!>  write(*,*) sVDP
!>  \endcode      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 rjl(1), yrjl(1), ydev1(1), ydev2(1)
      parameter (pig=3.1415926535897932d0)
      dimension rsvalues(28),r100values(28)
      include 'paramsoptS.i'
      include 'sr.i'
      include 'datavdp.i'
      external betairec
      external gammarec
      external hypgeo
      external fbars
      
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
            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200
         endif

c     attempt to insert generalized NFW model         
      elseif (kmp.eq.10) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)  
    
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
       
       xm =  hfx/hfy*(t/r200)**(3-dA3)     
c     Clustering Dark Energy 
c     dA3 = exponent of gNFW
c     dA2 = f(cs,w)
c     dA1 = cs2 (in unit of c**2)        
      elseif (kmp.eq.18) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)

       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
       
       xm=hfx/hfy*(t/r200)**(3-dA3)*(1.0d0+3*10**(screen-3.5)*(1+tmass))
             
      elseif (kmp.eq.8)  then	
c     M(r) is the modified NFW profile with a beyond Horndeski model 
c     See Saltas+16, Sakestein+16
c     the free parameter tmass is now the coupling constant Y 
	
       fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
       abey=tmass*t*t*(rs-t)/(rs+t)**3
       xm=((dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))+abey/4)*fac200
       
      elseif (kmp.eq.9) then

c       fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
c       xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200+dphidr(t/rs)
       
       xm = tM_bnfw_CS(t)

      elseif (kmp.eq.15) then
c     M(r) is the modified gNFW profile in chameleon screening 
c     gamma is fixed to dA3
        xm = tM_gnfw_CS(t)

      elseif (kmp.eq.16) then
c     M(r) is the modified Isothermal profile in chameleon screening 
c     the exponent is 3/2       
       xm = tM_iso_CS(t)
       
      elseif (kmp.eq.17) then 
      !M(r) is multicomponent general Chameleon. gNFW+Jaffe+beta+NFW
        xm = grav*tmass_Gen(t)/gm200 ! To divide by gm200 as this is already
                                ! including the normalization factor       
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
         
      elseif (kmp.eq.11) then
c        mass is modified Burkert in Chameleon Gravity   
         trs=t/rs
         rvrs=r200/rs
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))
         xm1=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))
         xm=xm1+dphib(trs)         


      elseif (kmp.eq.12) then
c     mass is modified Einasto in Chameleon Gravity
       !  eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm1=fac200*gip
         
         xm=xm1+dphiE(t/rs) 
         
          
      elseif (kmp.eq.13) then
c     M(r) is Refracted Gravity; no DM contribution allowed for the time begin
c     only BCG parameters kept free
c
        !refracted gravity with the function of Valeria
        xm= grav*rgmass(t)/gm200
        !fbarr=0.0
        !if (kbary.eq.1) fbarr=fbars(t)
        
        !Covariant formulation. It is not good for the moment
        !xm=grav*(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)/gm200
      elseif (kmp.eq.14) then
      
        !M(r) is boson star + gas and galaxies
        fbarr=0.0
        if (kbary.eq.1) fbarr=fbars(t)
        xt = t/rs
        if (xt.le.2) then 
         xm = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
        else
         xm =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
        endif 
         xm = xm + grav*fbarr/gm200
      
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
       !  eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm=fac200*gip
      endif


!*********************** BCG *******************************************
      fbarr=0.0
      if (kbary.eq.1) fbarr=fbars(t)      
      
      if (kbcg.eq.1) then

       !GR setup
       if(kmp.eq.10.or.kmp.le.6.or.kmp.eq.18) then 
        
        fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
        xm=xm+fjaffe !inclusion of gas,galaxies and bcg
        
        !INCLUDED IN THE GENERAL CHAMELEON
        !Chameleon setup
        !elseif (kmody.eq.1) then 
        !dpj=dphiJ(t/rjaf) !Include the chameleon field of the BCG
        !fjaffe=grav*(xmstar*xlumbcg*(t/rjaf/(1.+t/rjaf)+dpj)+fbarr)/gm200
        !xm=xm+fjaffe !inclusion of gas,galaxies and bcg
        
        !Refracted gravity setup
        !Hernquist mass for the BCG
       elseif (kmp.eq.13) then
        !the mass is all given by the BCG and gas+baryons, no DM
        
        !refracted gravity with the function of Valeria
        xm= grav*rgmass(t)/gm200

        !Covariant formulation. It is not good for the moment
        
        !xm=grav*(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)/gm200
        
       elseif (kmp.eq.14) then
        !Boson star + BCG
        fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
        
        xt = t/rs
        if (xt.le.2) then 
         xmer = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
        else
         xmer =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
        endif 
        
         xm = xmer+fjaffe
         
       endif  
        
      endif     
!***********************************************************************
      xm=xm*gm200 ! physical units


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
               xk2 = xk2 + dbinom(nint(-bec),k)*
     &              (u*u-1.d0)**k/(2.d0*k+1.d0)
            enddo
            do k=0,nint(-bec-1)
               xk3=xk3+dbinom(nint(-bec-1),k)*
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
c     if the anisotropy profile is a simplified Wojtak,
c     simplified Tiret (modified ML), generalized Tiret,
c     generalized OM or modifed Tiret (Opposite), there are no 
c     analytical solutions to provide K(r,ra) in eq.(A8) 
c     of ML05b, so we use the expression that includes
c     sigma_r for which we have a spline interpolation
c
      else

c     choose between simplified Wojtak, simpl. Tiret,
c     a model similar to Tiret (modified Tiret)
c     with non-zero central anisotropy, and Hansen+Moore

         if (kani.eq.3) then
            b=t**1.5/(t**1.5+cbe**1.5)
            
         elseif (kani.eq.21) then
            rm2=rs                     ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
            !avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
            
            becin=1.-1./(cbe*cbe)
            bec0=1.-1./(cbe0*cbe0)
            if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
            if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
            b=bec0+(becin-bec0)*t*t/(t*t+rm2*rm2)        
            if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values            
         
         elseif (kani.eq.4) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
            !avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            
            b=(1.-1./(cbe*cbe))*t/(t+rm2)
            
         elseif (kani.eq.41) then
            rm2=rs                     ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
            !avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free
            
            becin=1.-1./(cbe*cbe)
            bec0=1.-1./(cbe0*cbe0)
            if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
            if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
            b=bec0+(becin-bec0)*t/(t+rm2)        
            if (b.ge.1.) bec=0.999d0  ! avoid unphysical values 

        elseif (kani.eq.8) then
c       WORK IN PROGRESS
c       Biviano&Pizzuti model
         bi=1.-1./(cbe*cbe)
         b0=1.-1./(cbe0*cbe0)
         if (bi.ge.1.) bec=0.999d0  ! avoid unphysical values
         if (b0.ge.1.) bec0=0.999d0  ! avoid unphysical values

         rm2=rs                     ! NFW or other mass profiles
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
         if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
         if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
         if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
         if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
        
         if(ifreeBeta.gt.0) rm2 = rbeta !rbeta free  

         b = b0+(bi-b0)*(t/(t+rm2)+t**2/rm2**2*dexp(-(t/rm2)**2) )
         if (b.ge.1.) bec=0.999d0  ! avoid unphysical values

        
         elseif (kani.eq.5) then
            rm2=rs              ! NFW or other mass profiles
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist
            if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
            if (kmp.eq.10.or.kmp.eq.15.or.kmp.eq.17) rm2=(2-dA3)*rs !gNFW
            if (kmp.eq.18) rm2 = (2 - dA3) *rs !clustering DE gNFW
            if (kmp.eq.16) rm2=dsqrt(2.0d0)*rs !Isothermal
            !avoid unphyiscal values 
            if (rm2.lt.0) rm2=1e-4
            
            b=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
         else   ! Hansen+Moore
            rhm=rc ! use nu(r)
cc            rhm=rs ! use rho(r)
            b=ahm-bhm*(rhm+3.*t)/(rhm+t)            
         endif

c     interpolate sigma_r

         rjl(1)=tlog
c         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
         call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  ydev1,ydev2)
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
         r(j) = r_uniform_01(idum)
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

c
c     compute the difference in 'projected mass'
c     at the extreme of the radial range, for
c     3 different Sigma(R) profiles
c     
      subroutine sigmarnorm(x,fnorm)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptS.i'
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
      real*8 rjl(1), yrjl(1), ydev1(1),ydev2(1)
      include 'paramsoptS.i'
      include 'sr.i'
      include 'vlos.i'
      include 'tsallis.i'
c
      t2 = 2.*t1  !raggio proiettato
c
c     Interpolate to get the sigma_r(t)
c
      rjl(1)=dlog(t1)
c      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  ydev1,ydev2)
c      if (ier .eq. 34) then
c         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
c     &    rjl(1)
c      endif
      sigrt1=dexp(yrjl(1))

      rjl(1)=dlog(t2)
c      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
      call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
     &  ydev1,ydev2)
c      if (ier .eq. 34) then
c         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',
c     &     rjl(1)
c      endif
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
      include 'paramsoptS.i'
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
      include 'paramsoptS.i'
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
      include 'paramsoptS.i'
      external gammarec

      bec=1.-1./(cbe*cbe)
      if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values
      a=bec+0.5d0
      b=0.5d0

      if (a.gt.0.d0) then
         call incob(a,b,x,bir)
c         bir=dbetai(x,a,b)
      else
         nit=nint(0.5d0-a)
         birgj=0.d0
         do j=0,nit-1
            birgj=gammarec(a+b+j)/gammarec(a+1.d0+j)/gammarec(b)*
     *           x**(a+j)*(1.d0-x)**b+birgj
         enddo
         call incob(a+nit,b,x,bic)
         bir=birgj+bic !dbetai(x,a+nit,b)
      endif

      betairec=bir
      return
      end

c
c     Subroutine for the computation of Max Lik
c     with BOBYQA or NEWUOA
c     
      subroutine calfun(n,x,f)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)!,
c     &     omegal=0.7,omega0=0.3)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension x(n),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension y2a(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),wmbm(25000),embm(25000)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptS.i'
      include 'datavdp.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
c      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20

      ic = rknots
      ninterp=21*2

      call freepareva(n,x,r200,rc,rs,cbe,tmass,screen,
     & dA3,cbe0,rbeta,rjaf,xmstar)
ccc     write(*,429) r200,rc,rs,cbe,tmass,f
ccc 429  format(' In calfun: r200,rc,rs,cbe,tmass,L = ',5(f6.3,2x),f10.3)
      omegal=Olam
      omega0=Omegam
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
      rinfinity=25.0d0

      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin
         call dgaus8 (sr2int,xris(i),dlog(2.*rinfinity),errrel,
     &   risl, IERR) 


         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)

      enddo

c
c compute spline coeffs for later interpolation of sigma_r
c
       call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
c      call icsccu(xris,yris,ninterp,csr,icsr,ier)
c      if (ier .gt. 128) then
c         print *,' in VMAXLIK ...'
c         do i = 1, ninterp
c            print *,' i xris yris = ', i, xris(i), yris(i)
c         enddo
c      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.d0
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      !if (kintd.eq.1) then
      !   rinfinity=0.99999*r200
      !endif

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
               !if (kintd.eq.0) call rinf(rinfinity)
               call rinf(rinfinity)
c
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c
               umax = dacosh(rinfinity/rj)
               errrel=0.001d0
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)

c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif
c     
c     Add interloper contribution if requested
c
c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.001d0
cc                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
cc     &                 errest,ier)
c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR)  
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
ccc                  write(*,*) ' Number of obsd galaxies ',ngambm
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
ccc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c              endif
               
              if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


              if (nrc.gt.0.5) then
                 
                 call sigmarnorm(rc,fnorm)
                 !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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
c         call pwl_interp_2d (nr, nv, rp, vp, z, xx, yy,dlogp)
         
c        spline interpolation
         call splicoff2(rp,vp,z,nr,nv,y2a)
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
c            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp,ier)
            call spli2d(rp,vp,z,y2a,nr,nv,xx,yy,dlogp3)
            psum=psum-dlogp3*wmbm(l)
 
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

            !if (kintd.eq.0) call rinf(rinfinity)
            call rinf(rinfinity)

c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)

               umax = dacosh(rinfinity/rj)
               errrel=0.001d0 
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif


c     Add interloper contribution if requested

c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4.+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  n2=2
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.001d0
cc                  xinteg=dcadre(sdint,xmin,xmax,errabs,
cc     &                 errrel,errest,ier)
c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR) 
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c                  g=g+gi
c               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally
               
               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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
      if(kbcg.eq.1.and.kbdt.eq.1) then 
         rtbcg=20. !Max radius for the integration of the BCG integral
         
         call vdplik(chibcg)
         f=f+chibcg
c         write(*,*) r200, rs, xmstar, cbe, chibcg
      endif
      return
      end

c
c     function for the computation of Max Lik
c     with POWELL
c     
      function func(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)!,
c     &     omegal=0.7,omega0=0.3)
      integer RKNOTS, VZKNOTS
      parameter (RKNOTS=10,VZKNOTS=6)
      dimension x(9),rp(RKNOTS),vp(VZKNOTS),z(RKNOTS,VZKNOTS)
      dimension y2a(RKNOTS,VZKNOTS)
      dimension bsc(RKNOTS,VZKNOTS),xknot(13),yknot(9)
      dimension rmbm(25000),vmbm(25000),embm(25000),wmbm(25000)
c      dimension rmbcg(2500),vmbcg(2500),embcg(2500)
      real*8 c(2,RKNOTS,2,VZKNOTS)
      real*8 wk(2*RKNOTS*VZKNOTS+2*10)
      include 'paramsoptS.i'
      include 'datavdp.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'free.i'
      external sr2int,sr2out,sigmar1,sigmar2,gwenu
      external sigmarnorm
      external sdint
c      call uerset(1,levold)


c     Max Lik fit: start by computing the sigma_r(r) function in 
c     Ninterp points logarithmically spaced between 0.001 and 20
      
      omegal=Olam
      omega0=Omegam
      
      ic = rknots
      ninterp=21*2
      
      n=nfreefunc

      call freepareva(n,x,r200,rc,rs,cbe,tmass,screen,
     & dA3,cbe0,rbeta,rjaf,xmstar)
      
      call CPU_TIME(start) !uses cpu time for the random generator seed
      nst=int(1330*start, kind=8)
      nseed = 123434789+(-1)**nst*nst
      
      
      if(r200.gt.r2up.or.r200.lt.r2low) then
c       write(*,*) 'Warning: best fit outside the prior range'
c       write(*,*) 'redefine the best fit as the average value'
       r200=(r2up-r2low)/2*(1+r_uniform_01 (nseed))
      endif
      if(rs.gt.rsup.or.rs.lt.rslow) then
c       write(*,*) 'Warning: best fit "rs" outside the prior range'
c       write(*,*) 'redefine the best "rs" fit as the average value'
       rs=(rsup-rslow)/2*(1+r_uniform_01 (nseed))
      endif 
      if(rc.gt.rcup.or.rc.lt.rclow) then
c        write(*,*) 'Warning: best fit "rc" outside the prior range'
c        write(*,*) 'redefine the best fit "rc" as the average value'
        rc=(rcup-rclow)/2*(1+r_uniform_01 (nseed))
      endif 
      if(cbe.gt.bup.or.cbe.lt.blow) then
c        write(*,*) 'Warning: best fit "beta" outside the prior range'
c        write(*,*) 'redefine the best fit "beta" as the average value'
        cbe=(bup-blow)/2*(1+r_uniform_01 (nseed))
      endif 
      if(tmass.gt.tmup.or.tmass.lt.tmlow) then
c        write(*,*) 'Warning: best fit "tmass" outside the prior range'
c        write(*,*) 'redefine the best fit "tmass" as the average value'
        tmass=(tmup-tmlow)/2*(1+r_uniform_01 (nseed))
      endif
      if(screen.gt.scrup.or.screen.lt.scrlow) then
       
c       write(*,*) 'Warning: best fit "screen" outside the prior range'
c       write(*,*) 'redefine the best fit "screen"  as the average value'
       screen=(scrup-scrlow)/2*(1+r_uniform_01 (nseed))
       write(*,*) scrnew
      endif
      if(cbe0.gt.cb0up.or.cbe0.lt.cb0low) then
c        write(*,*) 'Warning: best fit "beta2" outside the prior range'
c        write(*,*) 'redefine the best fit "beta2" as the average value'
        cbe0=(cb0up-cb0low)/2*(1+r_uniform_01 (nseed))
      endif 
      if(rbeta.gt.rbetaup.or.rbeta.lt.rbetalow) then
c        write(*,*) 'Warning: best fit "beta2" outside the prior range'
c        write(*,*) 'redefine the best fit "beta2" as the average value'
        rbeta=(rbetaup-rbetalow)/2*(1+r_uniform_01 (nseed))
      endif       
      if(rjaf.gt.rjaup.or.rjaf.lt.rjalow) then
c        write(*,*) 'Warning: best fit "rjaffe" outside the prior range'
c        write(*,*) 'redefine the best fit "rjaffe" as the average value'
        rjaf=(cb0up-cb0low)/2*(1+r_uniform_01 (nseed))
      endif 
      
      if(xmstar.gt.xmaup.or.xmstar.lt.xmalow) then
c        write(*,*) 'Warning: best fit "M*/L" outside the prior range'
c        write(*,*) 'redefine the best fit "M*/L" as the average value'
        xmstar=(xmaup-xmalow)/2*(1.0d0+r_uniform_01 (nseed))
      endif      
      
      
      
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
      rinfinity=25.0d0
      
      do i=1,ninterp
         xx2=2.d0*rinfinity
         xris(i) = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)
         xmin=dexp(xris(i))
     
c     dqdagi integrates from a lower bound (xmin) to +infinity
c     while dqdag does not

         xx1=xmin

         call dgaus8 (sr2int,xris(i),dlog(2.*rinfinity),errrel,
     &   risl, IERR) 

         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         yris(i)=dlog(risok)

      enddo
c
c compute spline coeffs for later interpolation of sigma_r
c

      call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
      
c      call icsccu(xris,yris,ninterp,csr,icsr,ier)
c      if (ier .gt. 128) then
c         print *,' in VMAXLIK ...'
c         do i = 1, ninterp
c            print *,' i xris yris = ', i, xris(i), yris(i)
c         enddo
c      endif

c     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]

      psum=0.d0
      wsum=0.d0


c     If interlopers are required, the g function is composed
c     of 2 parts; one is the usual integral, truncated to r200

      !if (kintd.eq.1) then
        ! rinfinity=0.99999*r200
      !endif

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
               !if (kintd.eq.0) call rinf(rinfinity)
               call rinf(rinfinity)
c
c     here I uses dqdag to handle peak singularities 
c     that dqdagi does not handle 
c              
               errrel=0.001d0
               umax = dacosh(rinfinity/rj)
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
c
c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif
c     
c     Add interloper contribution if requested
c
c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.001d0
cc                  xinteg=dcadre(sdint,xmin,xmax,errabs,errrel,
cc     &                 errest,ier)
c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR) 
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
ccc                  write(*,*) ' Number of obsd galaxies ',ngambm
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
ccc                  print *,' Nh, Ni, Nv, Nobs = ',xnh,xni,xnv,ngambm
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c              endif
               
               if (g.ne.g) write(*,*) 'BEWARE! ',rj,vj,g

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally


               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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
     
c         call ibcccu(z,rp,nr,vp,nv,c,ic,wk,ier)
         call splicoff2(rp,vp,z,nr,nv,y2a)
         
         do l=1,ngambm
            xx=rmbm(l)
            yy=abs(vmbm(l))
            wsum=wsum+wmbm(l)
c            call ibcevl(rp,nr,vp,nv,c,ic,xx,yy,dlogp,ier)
            call spli2d(rp,vp,z,y2a,nr,nv,xx,yy,dlogp3)
            psum=psum-dlogp3*wmbm(l)
 
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

            !if (kintd.eq.0) call rinf(rinfinity)
            call rinf(rinfinity)

c     check the limits of the integral

            if (rj.lt.rinfinity) then
               wsum=wsum+wmbm(j)
     
c     compute the integral of gwen, yielding g(R,vz)

               umax = dacosh(rinfinity/rj)
               errrel=0.001d0
               call dgaus8 (gwenu,0.d0,umax,errrel,gdgau, IERR) 
               g= rj*gdgau
               
c               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)
c               if (ier .gt. 128) then
c                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2
c                  stop
c               endif


c     Add interloper contribution if requested

c               if (kintd.eq.1) then
c                  rjrv=rj/r200
c                  vzvv=vj/v200
c                  ctr=r200/rc
c                  aint=10.**(-1.061+0.364*rjrv*rjrv-
c     -                 0.580*rjrv**4.+0.533*rjrv**6)
c                  sigmaint=0.612-0.0653*rjrv*rjrv
c                  bint=0.0075
c                  n2=2
c                  gihat=aint*dexp(-0.5*(vzvv/sigmaint)**2)+bint
c                  call virsphmassp(rup,xmprmax)
c                  call virsphmassp(rlow,xmprmin)
c                  xmin=rlow/r200
c                  xmax=rup/r200
c                  errabs=0.d0
c                  errrel=0.001d0

c                  call dgaus8 (sdint,xmin,xmax, errrel, 
c     &                xinteg, IERR) 
c                  xmctr=(dlog(1.+ctr)-ctr/(1.+ctr))/(dlog(2.d0)-0.5d0)
c                  xnv=ngambm/((xmprmax-xmprmin)/xmctr+2.*pig*xinteg)
c                  xnh=xnv*(xmprmax-xmprmin)/xmctr
c                  xni=xnv*2.*pig*xinteg
c                  gi=xnv/(v200*r200*r200)*gihat
c                  gi=gi/xnv    ! g is in units of Nv, so I use the same units for gi
c                  g=g+gi
c               endif


c     From g, one then has the probability p

c     There are two expressions: one if for the case in which c_nu
c     is part of the fitting parameters of MAMPOSSt
c     the other is for the case in which c_nu is fitted externally

               if (nrc.gt.0.5) then

                  call sigmarnorm(rc,fnorm)
                  !if (kintd.eq.1) fnorm=fnorm*(xnh+xni)/xnh
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
      
      if(kbcg.eq.1.and.kbdt.eq.1) then 
         rtbcg=20. !Max radius for the integration of the BCG integral
         
         call vdplik(chibcg)
         func=func+chibcg
      endif

      return
      end
    
c
c

      subroutine freepareva(nfreepar,freepar,r200new,rcnew,rsnew,cbnew,
     &           tmassnew,scrnew,dA3new,cb0new,rbetanew,rjafnew,xmsnew)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
       parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,q0=-0.64,clight=2.99792458d5)
      dimension freepar(nfreepar)
      include 'paramsoptS.i'
      include 'datavdp.i'
      omegal=Olam !0.7
      omega0=Omegam!2.*(q0+omegal)

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
         cbnew=freepar(ifp)
      endif
      if (ipar(5).eq.0) then
         tmassnew=tmassg
      else
         ifp=ifp+1
         if(kmp.ne.8) tmassnew=10.**freepar(ifp)
         if (kmody.eq.1.and.kmp.ne.8.or.kmp.eq.17) then !explore from 0 to 1 the range of parameter
         ! in the case of CS 
          
          if (freepar(ifp).lt.0) freepar(ifp)=dabs(freepar(ifp)) !avoid < 0
          if (freepar(ifp).ge.1) freepar(ifp)=1/freepar(ifp) !avoid >1
          tmassnew=-dlog(1-freepar(ifp))/0.1
         endif
         if (kmp.eq.8) then
           tmassnew=freepar(ifp)
           if (tmassnew.eq.0.0d0) tmassnew=1.e-4
         endif
      endif
      
      if (ipar(7).eq.0) then
         cb0new=cbe0g
      else
         ifp=ifp+1
         cb0new=freepar(ifp)
      endif

      if (ipar(8).eq.0) then
         rjafnew=rjafg
      else
         ifp=ifp+1
         rjafnew=10.**freepar(ifp)
      endif      
      if (ipar(9).eq.0) then
         xmsnew=xmstarg
      else
         ifp=ifp+1
         xmsnew=10.**freepar(ifp)
      endif
      if (ipar(10).eq.0) then
         dA3new=dA3g
      else
         ifp=ifp+1
         dA3new=10.**freepar(ifp)
      endif      
      if (ipar(11).eq.0) then
         rbetanew=rbetag
      else
         ifp=ifp+1
         rbetanew=10.**freepar(ifp)
      endif             
      if (ipar(6).eq.0) then  !TEST2
         scrnew=screeg
      else
         ifp=ifp+1
         scrnew=10.**freepar(ifp)
         if (kmody.eq.1.and.kmp.ne.8) then !explore from 0 to 1 the range of parameter
         ! in the case of CS 
          if (freepar(ifp).lt.0) freepar(ifp)=dabs(freepar(ifp)) !avoid < 0
          if (freepar(ifp).ge.1) freepar(ifp)=1/freepar(ifp) !avoid >1
          scrnew=freepar(ifp)/(1-freepar(ifp))
         endif

      endif



      if (klcdm.eq.1) then      ! LCDM c=c(M)
         rm200=100.*hz*hz/grav*r200new**3
         cmean=6.76*(rm200/1.e12)**(-0.098)
         rsnew=r200new/cmean
         if (kfs2.eq.1) rcnew = 4 * rsnew
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
  
c************************************************************************
!>   \addtogroup MOD Modified Gravity
!>   \details This set of functions are used to compute some quantities needed 
!! for the kinematic analysis in modified gravity

!> \ingroup MOD 
      
      function dphidr(x) 
!>   @author L. Pizzuti
!>   @details  compute the effective mass due to chameleon field `Meff=Q/G *r**2 * dphidr`
!!   assuming a NFW density profile. Note that it has been generalized by dphibn.
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen`. 
!!  The background field value `tmass` is given by tmass (in unit of 1e-5)  
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za`                                 
!>  - mass model `kmp` 
!>  - exponential cutoff `ncutof`: if greater than 0, it sets a cutoff to the 
!!  fifth force.
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       h0=70
!>       Omegam=0.3
!>       Omegal=0.7
!>       kmp=9
!>       r200=1.2
!>       rs=0.5
!>       tmass=100.23
!>       nhone=10
!>       screen=0.3
!>       write(*,*) dphidr(2.12d0)
!>  \endcode
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, ampl=2.43e27, ah0=1.48e-33,
     & ampc=1.57e29,rhoc=3.88e-11, clight=2.99792458d5)
       
       parameter (grav=4.302e-9, alimit=0.001)
       include 'paramsoptS.i'

       hz=h0*dsqrt(Omegam*(1+za)**3+Olam)
       c200=r200/rs
       phinf=tmass*1e-5 !tmass is phinf/Mpl in units of 1e-5clight**2
       fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))

       if(kmp.eq.9.and.nhone.ge.0) then
         bcoup=screen !the modified parameter of the screening becomes
                      !the coupling constant in the genneric chameleon run

       else
         bcoup=1./dsqrt(6.0d0)
       endif
       !this factor is rho_s*Q*rs**2/Mpl**2. In unit of clight**2
       B=bcoup*200*c200**3*fac200*hz**2*rs**2/clight**2
       
       xczero=(B/phinf-1)
       
       Czero=-B*dlog(1+xczero)+phinf*xczero
        if (xczero.lt.0.001) then
          xczero=0.0   
          Czero=0.0  
        endif
        
        !Exponential cutoff
        dminf=1
        if (ncutof.gt.0) then
            dminf=dexp(-x*rs*dsqrt(1e-7*bcoup/phinf))
        endif
    
        if (x.le.xczero) then
            dphidr=0
        else 
            cdav=rs*bcoup*clight**2/(100*hz*hz*r200**3)
            dphidr=cdav*(Czero-B*(x/(1.0d0+x)-dlog(1.0d0+x)))*dminf

        endif
        
        return
       end 

!modified gravity functions for other mass profiles ********************
!> \ingroup MOD        
      function  dphibn(x,rhos,rss,b,Q,phinf) 
!>   @author L. Pizzuti
!>   @details  this function computes `Meff=Q/G *r**2 * dphidr` 
!!   (where phi is Phi/Mpl) for a chameleon field `Phi`, assuming a
!!   b-NFW mass density profile.
!!  It requires the parameters `ncutof` in the COMMON block paramsoptS,
!!  which set an exponential cutoff to the fifth force if greater than 0.
!>   @param[in] x REAL*8, value at which the function is computed
!>   @param[in] rhos REAL*8, characteristic density [M_sun/Mpc**3]
!>   @param[in] rss  REAL*8  scale radius of the density profile [Mpc]
!>   @param[in] b  REAL*8    exponent of the b-NFW profile
!>   @param[in] Q  REAL*8    Coupling of the field with matter
!>   @param[in] phinf  REAL*8    background value of the scalar field (units of c**2)
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       rhos=1.2e15
!>       rss=0.5
!>       Q=1.23
!>       Phinf=4.3e-5
!>       ncutof = 1
!>       b = 3
!>       write(*,*) dphibn(2.12d0,rhos,rss,b,Q,phinf)
!>  \endcode       

      
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, ampl=2.43e27, ah0=1.48e-33,
     & ampc=1.57e29,rhoc=3.88e-11, clight=2.99792458e5, cl=3e5)
       
       parameter (grav=4.302e-9, alimit=0.0005) !limiting value for the screening radius
       include 'paramsoptS.i'    
       
       !Facda=Q*rss**2*rhos*8*pig*grav/clight**2 !factor containing the planck mass
       Facda=Q*rss**2*8*pig*grav*rhos/cl**2
       !Exponential cutoff
        dminf=1.0
        if (ncutof.gt.0) then
            dminf=dexp(-x*rss*dsqrt(1e-7*Q/phinf))
        endif
        
        if (b.gt.2) then !bNFW case
         cdav=(Q**2)*(rss**3)*rhos*8*pig/(clight**2) !factor in front of everything
         
         xcs=(phinf*(b-1.0d0)/Facda)**(1.0d0/(1.0d0-b))-1.0d0
         
         if (xcs.ge.alimit) then
            Cs=(1+xcs)**(1.0d0-b)*(1.0d0-xcs+b*xcs)/((b-1)*(b-2))
         else
          xcs=0.0d0            
          Cs=1/((b-1)*(b-2))
         endif 
         
         if (x.gt.xcs) then
          dout1 = ((x+1)**(1-b)*(x-b*x-1)/((b-1)*(b-2))+Cs)
          dout = dout1*cdav*clight**2*dminf
         else
          dout = 0.0d0
         endif
        
        else if (b.eq.2) then !NFW case
c          write(*,*) 'check'
          xcs=Facda/phinf-1.0d0
          if (xcs.ge.alimit) then
            !write(*,*) xcs, 'here', Facda
            Cs=-phinf-Facda*dlog(1+xcs)
          else
            xcs=0.0d0
            Cs=-Facda
          endif
          if (x.gt.xcs) then
          
           dout1=(Facda*(1/(x+1.0d0)+dlog(1.0d0+x))+Cs)
           dout = dout1*rss*Q*clight**2/grav*dminf
          else
           dout = 0.0d0
          endif
        
        else
         dout = 0.0d0
            
        endif
        
        dphibn = dout 
        
        return  
       end
       
!> \ingroup MOD        
      function tM_bnfw_CS(x) 
!>   @author L. Pizzuti
!>   @details  compute the effective mass due to chameleon field for a
!!   b-NFW mass density profile, which exponent is defined by `eim >=2`.
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen, eim`. 
!!  The background field value `tmass` is given by tmass (in unit of 1e-5)  
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za`                                 
!>  - mass model `kmp` 
!>  - number of steps in the second MG parameter `nhone`. If equal to -1
!!  it forces f(R) chameleon case.
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       h0=70
!>       Omegam=0.3
!>       Omegal=0.7
!>       kmp=9
!>       r200=1.2
!>       rs=0.5
!>       tmass=100.23
!>       nhone=10
!>       screen=4.3
!>       ncutof = 1
!>       eim = 3
!>       write(*,*) tM_bnfw_CS(2.12d0)
!>  \endcode       
        
       
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, ampl=2.43e27, ah0=1.48e-33,
     & ampc=1.57e29,rhoc=3.88e-11, clight=2.99792458d5)
       
       parameter (grav=4.302e-9, alimit=0.001) !limiting value for the screening radius
       include 'paramsoptS.i'
       
       cx = x/rs
       hz=h0*dsqrt(Omegam*(1+za)**3+Olam)
       c200=r200/rs
       phinf=tmass*1e-5 !in realtà è phinf/Mpl in unità di 1e-5clight**2
       
       bcoup=screen !the modified parameter of the screening becomes
                      !the coupling constant in the genneric chameleon run
       
       dm200=100*hz**2/grav*r200**3
       !bnfw
       if (eim.gt.2) then
        b = int(eim)
        f200=(c200+1)**(-b)*((c200+1.0)**b-b*c200*(c200+1)+c200**2-1.0)
        dma=dm200*(1.0d0 + (1 + cx)**(1 - b)*(-1.0d0 + cx - b*cx))/f200
        rhos=dm200*(b-1)*(b-2)/(4*pig*rs**3*f200)
        
        !Sono arrivato qui
         if (tmass.gt.0.001.and.bcoup.gt.0.001) then
         
          dphir = dphibn(cx,rhos,rs,b,bcoup,phinf)
         else
          dphir = 0.0d0
          
         endif
         dmtes =(dphir+dma)/dm200
        !nfw
        else if (eim.eq.2.0d0) then
         b = 2.0d0 !eim
         fac200=dlog(1+c200)-c200/(1+c200) 
         dma=(dlog(1+cx)-cx/(1+cx))/fac200
         rhos=dm200/(4*pig*rs**3*fac200)
         
         if (tmass.ne.0.and.bcoup.ne.0) then
          dphir = dphibn(cx,rhos,rs,b,bcoup,phinf)/dm200
         else
          dphir = 0.0d0
          
         endif
         dmtes = dphir+dma  
               
        else if (eim.lt.2) then
         dmtes = 0.0d0
        
        endif      
        !In MG-MAMPOSSt values should be in units of M200

        tM_bnfw_CS = dmtes
       return
       end
         
         
! ================ TEST FOR THE SCREENING VALUES ======================= 
!> \ingroup MOD         
      subroutine ScreenNFW(xcs)   
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, ampl=2.43e27, ah0=1.48e-33,
     & ampc=1.57e29,rhoc=3.88e-11, clight=2.99792458d5)
       
       parameter (grav=4.302e-9, alimit=0.001) !limiting value for the screening radius
       include 'paramsoptS.i'
         hz=h0*dsqrt(Omegam*(1+za)**3+Olam)
         c200=r200/rs
         phinf=tmass*1e-5
         fac200=dlog(1+c200)-c200/(1+c200) 
         dm200=100*hz**2/grav*r200**3
         rhos=dm200/(4*pig*rs**3*fac200)
         
         
         
         Facda=screen*rs**2*rhos*8*pig*grav/clight**2 !factor containing the planck mass
         
         xcs=(Facda/phinf-1.0d0)*rs
         if (xcs.lt.0) xcs = 0.0d0  
         
        return 
        end
c***********************************************************************      
      function dphigNFW(x,rhos,rss,bgam,Q,phinf) 
!>   @author L. Pizzuti
!>   @details  this function computes `Meff=Q/G *r**2 * dphidr` 
!!   (where phi is Phi/Mpl) for a chameleon field `Phi`, assuming a
!!   gNFW mass density profile.
!!  It requires the parameters `ncutof` in the COMMON block paramsoptS,
!!  which set an exponential cutoff to the fifth force if greater than 0.
!>   @param[in] x REAL*8, value at which the function is computed
!>   @param[in] rhos REAL*8, characteristic density [M_sun/Mpc**3]
!>   @param[in] rss  REAL*8  scale radius of the density profile [Mpc]
!>   @param[in] bgam  REAL*8    exponent of the b-NFW profile
!>   @param[in] Q  REAL*8    Coupling of the field with matter
!>   @param[in] phinf  REAL*8    background value of the scalar field (units of c**2)
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       rhos=1.2e15
!>       rss=0.5
!>       Q=1.23
!>       Phinf=4.3e-5
!>       ncutof = 1
!>       bgam = 3
!>       write(*,*) dphibn(2.12d0,rhos,rss,bgam,Q,phinf)
!>  \endcode  
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, clight= 2.99792458d5)
       parameter (grav=4.302e-9, alimit=0.001) !limiting value for the screening radius
       include 'paramsoptS.i'    
       
       Facda=Q*rss**2*rhos*8*pig*grav/clight**2 !factor containing the planck mass
       !Exponential cutoff
       dminf=1.0d0
       
       if (ncutof.gt.0) then
            dminf=dexp(-x*rss*dsqrt(1e-7*Q/phinf))
       endif
      
       call hygfx ( 3-bgam, 3-bgam, 4-bgam, -x, hfx ) 
       
       
       facsc = 1 + (phinf/Facda)*(bgam - 2.0d0)

       if (facsc.lt.0) then
        xsc = alimit
       else
        xsc = 1/(1 - ( facsc )**(1/(2-bgam)) ) - 1 
       endif
       
       if (isnan(xsc)) xsc = alimit
       
       
       !Constant factor containing rs
       Clogg = 1 + dlog(rss)
    
       if (xsc.gt.alimit) then
        
        call hygfx ( 3-bgam, 3-bgam, 4-bgam, -xsc, hxsc)
        !Compute the constant of integration
        
        Cgnfw = - Facda  * ( Clogg + (xsc**(3-bgam))/(3 - bgam)*hxsc)
        
        if (x.lt.xsc) then
        
         dout=0.0d0
        else            
         dout= Cgnfw + Facda * ( Clogg  + x**(3 - bgam)/ (3-bgam)*hfx) 
        endif
        
       else
        Ccg = - Facda  * (Clogg)
        dout = Ccg + Facda  * ( Clogg + x**(3- bgam)/ (3-bgam)*hfx)
       endif 
       
       dphigNFW =  dout*Q*rss*clight**2/grav*dminf
       
       return 
      end

c######=================================================================
      function tM_gnfw_CS(x) 
!>   @author L. Pizzuti
!>   @details  compute the effective mass due to chameleon field for a
!!   gNFW mass density profile, which exponent is defined by `dA3 <=2`.
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen, dA3`. 
!!  The background field value `tmass` is given by tmass (in unit of 1e-5)  
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za`                                 
!>  - mass model `kmp` 
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       h0=70
!>       Omegam=0.3
!>       Omegal=0.7
!>       kmp=9
!>       r200=1.2
!>       rs=0.5
!>       tmass=100.23
!>       nhone=10
!>       screen=4.3
!>       ncutof = 1
!>       dA3 = 1.2
!>       write(*,*) tM_gnfw_CS(2.12d0)
!>  \endcode       
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, clight=2.99792458d5)
       
       parameter (grav=4.302e-9, alimit=0.001) !limiting value for the screening radius
       include 'paramsoptS.i'
       
       cx = x/rs
       hz=h0*dsqrt(Omegam*(1+za)**3+Olam)
       c200=r200/rs
       phinf=tmass*1e-5
       
       bcoup=screen !the modified parameter of the screening becomes
                      !the coupling constant in the genneric chameleon run
       
       dm200=100*hz**2/grav*r200**3
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -c200, h2x)
       
       rhos=(dm200*r200**(-3.0+dA3)*rs**(-dA3)*(3.0-dA3))/(4*pig*h2x)
       
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -cx, h2y)
       dmass = (x/r200)**(3-dA3)*dm200*h2y/h2x
       
       if (bcoup.ne.0.and.tmass.ne.0) then 
         dphidr= dphigNFW(cx,rhos,rs,dA3,bcoup,phinf) 
      else
         dphidr=0.0d0
      endif
        tM_gnfw_CS = (dphidr + dmass)/dm200 !given in units of m200
      return
        
      end


c***********************************************************************     
     
     

      function r_normal_ab ( a, b, nseed )

!***********************************************************************

!! R8_NORMAL_AB returns a scaled pseudonormal R8.

!  Discussion:

!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    06 August 2013

!  Author:

!    John Burkardt

!  Parameters:

!    Input, real ( kind = 8 ) A, the mean of the PDF.

!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.

!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.

!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.

         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         
         parameter (r_pi = 3.141592653589793D+00)
         real ( kind = 8 ) r_uniform_01
         r1 = r_uniform_01 ( nseed )
         r2 = r_uniform_01 ( nseed )
         x = dsqrt ( - 2.0D+00 * dlog(r1) )*dcos(2.0D+00 *r_pi*r2)

         r_normal_ab = a + b * x

         return
      end

      function r_uniform_01 ( seed )

!***********************************************************************

!! R8_UNIFORM_01 returns a unit pseudorandom R8.

!  Discussion:

!    This routine implements the recursion

!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )

!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.

!    If the initial seed is 12345, then the first three computations are

!      Input     Output      R8_UNIFORM_01
!      SEED      SEED

!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    31 May 2007

!  Author:

!    John Burkardt

!  Reference:

!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.

!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.

!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.

!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.

!  Parameters:

!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.

!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.

       implicit none

       integer ( kind = 4 ) k
       real ( kind = 8 ) r_uniform_01
       integer ( kind = 4 ) seed

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
        seed = seed + 2147483647
       end if

!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!

       r_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

       return
      end
c***********************************************************************
      subroutine tranmod(r200new,rcnew,rsnew,cbnew,tmassnew,scrnew, !checkpoint
     &         dA3new,cb0new,rbetanew,rjafnew,xmsnew,sigma,npar,neseed)
     
c     This routine compute the new values of the free parameters in 
c     a MCMC using a gaussian distribution with std sigma(nfreepar)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       dimension sigma(npar)
       include 'paramsoptS.i'  
       external r_normal_ab 
    
       if (sigma(1).gt.0.0d0) then
        
        r200_old = r200new
        if (kmp.le.6) then
          almu = dlog(r200_old)
 11       zr2l=r_normal_ab (almu,sigma(1),neseed)
          r200new=dexp(zr2l)
          
          !To avoid very small, unphysical values when running assumig
          !Gaussian Priors, a lower limit of r200min = 0.1 r200g is set
          if (ipriorGauss.eq.1 .and. r200new.le.0.1d0*r200g) goto 11
          
        else
 10       r200new = dabs(r_normal_ab (r200_old,sigma(1),neseed))
        
          !To avoid very small, unphysical values when running assumig
          !Gaussian Priors, a lower limit of r200min = 0.1 r200g is set
          if (ipriorGauss.eq.1 .and. r200new.le.0.1d0*r200g) goto 10
        endif
       endif

       if (sigma(2).gt.0.0d0) then
         rc_old = rcnew
         alrc = dlog(rc_old)
         zrcl = r_normal_ab (alrc,sigma(2),neseed)
         rcnew=dexp(zrcl)
       endif
       
       if (sigma(3).gt.0.0d0) then
        rs_old = rsnew
        if (kmp.le.6) then
          alrs = dlog(rs_old)
 13       zrsl=r_normal_ab (alrs,sigma(3),neseed)
          rsnew=dexp(zrsl)
          !To avoid very small, unphysical values when running assumig
          !Gaussian Priors, a lower limit of rsmin = 0.05 Mpc is set
          if (ipriorGauss.eq.1 .and. rsnew.le.0.05d0) goto 13
        else
 12       rsnew=dabs(r_normal_ab (rs_old,sigma(3),neseed))
         !To avoid very small, unphysical values when running assumig
         !Gaussian Priors, a lower limit of rsmin = 0.05 Mpc is set
          if (ipriorGauss.eq.1 .and. rsnew.le.0.05d0) goto 12
        endif
       endif
        
       !Anisotropy parameter at infinity
       if (sigma(4).gt.0.0d0) then
        cb_old = cbnew
        alcb = dlog(cb_old)
        zcbl =r_normal_ab (alcb,sigma(4),neseed)
        cbnew=dexp(zcbl)
        !case of gT or gOM:
        if (kani == 41 .or. kani == 21 .or. kani == 8) then
         if (cbnew <= 0.4d0) then
            do itry = 1, 100
                zcbl  = r_normal_ab(log(max(cbnew,1.0d-12)), sigma(4), 
     &           neseed)
                cbnew = exp(zcbl)
                if (cbnew > 0.4d0) exit
             end do
            if (cbnew <= 0.4d0) cbnew = 0.5d0   ! fallback finale
         endif
        endif
       
       endif
       
       !anisotropy parameter a r=0
       if (sigma(7).gt.0.0d0) then
        cb0_old = cb0new
        alcb0 = dlog(cb0_old)
        zcb0l =r_normal_ab (alcb0,sigma(4),neseed)
        cb0new=dexp(zcb0l)
        !case of gT or gOM:
        !avoid very small values of beta which makes the integral overflow
        if (cb0new <= 0.4d0) then
            do itry = 1, 100
                zcb0l  = r_normal_ab(log(max(cb0new,1.0d-12)), sigma(4), 
     &           neseed)
                cb0new = exp(zcb0l)
                if (cb0new > 0.4d0) exit
             end do
            if (cb0new <= 0.4d0) cb0new = 0.5d0   ! fallback finale
        endif
       endif

       
       !BCG ***************************************************
       if (sigma(8).gt.0.0d0) then
        rjafnew=r_normal_ab (dlog(rjafnew),sigma(8),neseed)
        rjafnew=dexp(rjafnew)
       endif
       if (sigma(9).gt.0.0d0) then
        xmsnew=r_normal_ab (dlog(xmsnew),sigma(9),neseed)
        xmsnew=dexp(xmsnew)
       endif
       !****************************************************
    
       if (sigma(10) > 0.0d0) then
        ! special case if A3 = gamma gNFW
        !=============================================================
        if (kmp.eq.10 .or. kmp.eq.17 .or. kmp.eq.18) then
         
            imaxtry = 50   ! massimo numero di tentativi

            do itry = 1, imaxtry
                dA3new = r_normal_ab(dA3new, sigma(10), neseed)

                ! only physical values 0 < gamma < 2
                if (dA3new > 0.0d0 .and. dA3new < 2.0d0) exit
            end do

            ! after imaxtry trials, fallback
            if (dA3new <= 0.0d0 .or. dA3new >= 2.0d0) then
                dA3new = min(1.5d0, max(0.5d0, dA3g))  ! safety value
            endif
        
        else if (kmp.eq.13) then
            dA3new = dabs(r_normal_ab((dA3new), sigma(10), neseed))
        else
            dla3 = r_normal_ab(dlog(dA3new), sigma(10), neseed)
            dA3new = dexp(dla3)
        endif
       endif

c  ======================OLD ==========================================         
c       if (sigma(5).gt.0.0d0) then 
c       !DHOST case treated separately
c        if (kmp.eq.8) then
c        !!
c        !For mg Y parameter not log: it could be negative!
c        !!
c         tmassnew=r_normal_ab (tmassnew,sigma(5),neseed)
c        else

c         if(kmody.eq.1.and.kmp.ne.8.or.kmp.eq.17) then
         
c         !case of general chameleon. The parameter space is explored
c         !in terms of the rescaled variables Q_2 phi_2 (see e.g. Terukina
c         !et al., 2014, Pizzuti et al., 2021)
        
c          nLogP=1 !Set equal to 0 to explore the phi space using the logarithm
c         ! of the parameter. 
        
c         !Exploration is made using the phi_2 parameter
         
c          kcham=0
c          if(kmp.eq.9.or.kmp.eq.11.or.kmp.eq.12.or.kmp.eq.15) kcham=1
c          if (kmp.eq.16.or.kmp.eq.17) kcham=1
c          if (kcham.eq.1.and.nhone.ge.0.and.nLogP.eq.1) then
c           tmatemp=1-dexp(-tmassnew*0.1)

c           if (tmatemp.lt.0) then
c             tmatemp=dabs(r_normal_ab (0.d0,sigma(5),neseed))
c           else
c             tmatemp=dabs(r_normal_ab (tmatemp,sigma(5),neseed))
c           endif  
c             !be sure that tmassnew is smaller than 1
c           if (tmatemp.gt.1.0d0) then
c            tmatemp = r_uniform_01 ( neseed )
c           endif
                 
c           tmassnew=(-10*dlog(1-tmatemp))

c          !Exploration is made using the log(phi) parameter  
c          else 
c           tmassnew=r_normal_ab (dlog10(tmassnew),sigma(5),neseed)
c           tmassnew=10**(tmassnew)
c          endif
c         else if (kmp.eq.13.or.kmp.eq.14) then 
c         !Refracted gravity, BS: tmass is strictly positive
c          tmassnew=r_normal_ab (dlog10(tmassnew),sigma(5),neseed)
c          tmassnew=10**(tmassnew) 
c         else !all other cases for now
c          tmassnew=r_normal_ab (tmassnew,sigma(5),neseed)
c          !CDE: tmass striclty positive (is the cs^2 value)
c          if (kmp.eq.18.and.tmassnew.lt.0) tmassnew = dabs(tmassnew)
c         !tmassnew=dexp(tmassnew) I set as a trial a gaussian exploration with fixed step
c         !This is because gamma for gNFW can be close to zero, and we need to avoid problems
c         endif
c        endif
       
c       endif
c  ======================OLD ==========================================  

       
       if (sigma(5) > 0.0d0) then
        select case (kmp)

        !--------------------------------
         case (8)  ! DHOST: può essere negativo
         tmassnew = r_normal_ab(tmassnew, sigma(5), neseed)

        !--------------------------------
       
         case (9, 11, 12, 15, 16, 17)  ! Chameleon  !TOGLI 1            
            !kmody = 1
            
            if (kmody.eq.1) then
                
                !TO BE INCLUDED OUTSIDE
                !if I have only f(R), I have to explore in log
                !nLogP=1 !Set equal to 0 to explore the phi space using the logarithm
                ! of the parameter. 
                
                if (nhone.gt.0.or.nLogP.eq.1) then

                    eps = 1.d-12

                    t_low = 1.d0 - dexp(-0.1d0 * tmlow)
                    t_up  = 1.d0 - dexp(-0.1d0 * tmup)

                    tmatemp = 1.d0 - dexp(-0.1d0 * tmassnew)
                    tmatemp = min(max(tmatemp, t_low+eps), t_up-eps)

                    tmpdd = r_normal_ab(0.d0, sigma(5), neseed)
                    tmatnew = reflect_ab( tmatemp + tmpdd, t_low, t_up )

                    ! solo safety numerico estremo
                    tmatnew = min(max(tmatnew, eps), 1.d0-eps)

                    tmassnew = -10.d0 * dlog(1.d0 - tmatnew)

                    ! update current (if you accept later, do it there instead)
                    ! tmatemp = tmatnew

                else
                ! ---- log-phi space ----
                
             tmassnep = r_normal_ab(dlog10(tmassnew), sigma(5), neseed)
                    tmassnew = 10.0d0**tmassnep
                end if
            else
                ! If kmody ≠ 1, general case
                tmassnew = r_normal_ab(tmassnew, sigma(5), neseed)
            end if

            !--------------------------------
        case (13, 14)  ! Refracted gravity or Boson star: tmass > 0
            tmass_old = tmassnew
        
 15        tmassnew = dabs(r_normal_ab(tmass_old, sigma(5), neseed))    
          !To avoid very small, unphysical values 
c           if ( tmassnew.le.(0.002d0)) goto 13            
           !tmassnew = 10.0d0**tmassnew
            !--------------------------------
        case (18)  ! CDE: tmass positive
            tmassnew = r_normal_ab(tmassnew, sigma(5), neseed)
            if (tmassnew < 0.0d0) tmassnew = abs(tmassnew)

            !--------------------------------
        case default  ! other cases
            tmassnew = r_normal_ab(tmassnew, sigma(5), neseed)

        end select
       end if


       if (sigma(6) > 0.0d0) then
        select case (kmp)
        case(8,18,13)
         scrnew = r_normal_ab(scrnew, sigma(6), neseed)
        case default
        ! Chameleon 
            if (scrnew <= 0.0d0) then
            ! screening radius or coupling may be zero (not negative)
            ! Consider the absolute value
                scrnew = dabs(r_normal_ab(scrnew, sigma(6), neseed))
            else
            ! Trasform in log-space for positive values
                scrnew = r_normal_ab(dlog(scrnew), sigma(6), neseed)
                scrnew = dexp(scrnew)
            endif

            if (kmody.eq.1.or.kmp.eq.17) then
                eps = 1.d-12
                ! --- define scrtemp from current scrnew
                if (scrnew.le.0.d0) then
                   scrtemp = eps
                else
                   scrtemp = scrnew / (1.d0 + scrnew)
                endif

                ! --- bounds in scrtemp implied by scrlow/scrup
                s_low = scrlow / (1.d0 + scrlow)
                s_up  = scrup  / (1.d0 + scrup)

                ! keep current inside bounds (safety)
                scrtemp = min(max(scrtemp, s_low+eps), s_up-eps)

                ! ----- RW on scrtemp in [s_low, s_up] with reflection
                sctest = r_normal_ab(0.d0, sigma(6), neseed)
                scrtemp_new = reflect_ab(scrtemp + sctest, s_low, s_up)

                ! back to scrnew (guaranteed within [scrlow,scrup])
                scrtemp_new = min(max(scrtemp_new, eps), 1.d0-eps)  ! only numerical safety
                scrnew = scrtemp_new / (1.d0 - scrtemp_new)

            endif
        end select        
        
        !OLD. Keep for security reasons ================================
        
c        if (kmp.ne.8 .and. kmp.ne.18) then
c        ! Chameleon 
c            if (scrnew <= 0.0d0) then
c            ! screening radius or coupling may be zero (not negative)
c            ! Consider the absolute value
c                scrnew = dabs(r_normal_ab(scrnew, sigma(6), neseed))
c            else
c            ! Trasform in log-space for positive values
c                scrnew = r_normal_ab(dlog(scrnew), sigma(6), neseed)
c                scrnew = dexp(scrnew)
c            endif
c
c            if (kmody.eq.1 .or. kmp.eq.17) then
c            ! trasformazione in [0,1] per evitare accumulation ai bordi
c                scrtemp = scrnew / (1.0d0 + scrnew)
c            ! calcolo sigma logit-safe
c                sigma_logit = sigma(6) / (scrtemp * (1.0d0 - scrtemp))
c                ! passo Gaussiano in logit-space
c                y1 = dlog(scrtemp / (1.0d0 - scrtemp))
c                y1_new = r_normal_ab(y1, sigma_logit, neseed)
c                scrtemp = 1.0d0 / (1.0d0 + dexp(-y_new))
c                
c                scrnew = scrtemp / (1.0d0 - scrtemp)
c            endif
c
c        else
c        ! Vainshtein o CDE: coupling può essere negativo
c         scrnew = r_normal_ab(scrnew, sigma(6), neseed)
c        endif
       endif


       if (sigma(11).gt.0.0d0) then       
        rbeta_old = rbetanew
        alrbeta = dlog(rbeta_old)
 14     zlrbeta=r_normal_ab (alrbeta,sigma(11),neseed)
        rbetanew=dexp(zlrbeta)
        !To avoid very small, unphysical values when running assumig
        !Gaussian Priors, a lower limit of rsmin = 0.01 Mpc is set
        if (ipriorGauss.eq.1 .and. rbetanew.le.0.01) goto 14
       endif  

       
       return
      end
c
c


C       TEST REFLECTION
        double precision function reflect01(u)
          implicit none
          double precision, intent(in) :: u
          double precision :: x
          x = u
          do
             if (x .ge. 0.d0 .and. x .le. 1.d0) exit
             if (x .lt. 0.d0) x = -x
             if (x .gt. 1.d0) x = 2.d0 - x
          end do
          reflect01 = x
        end function reflect01
        
        double precision function reflect_ab(x, a, b)
          implicit none
          double precision, intent(in) :: x, a, b
          double precision :: y
          y = x
          do
             if (y .ge. a .and. y .le. b) exit
             if (y .lt. a) y = a + (a - y)
             if (y .gt. b) y = b - (y - b)
          end do
          reflect_ab = y
        end function reflect_ab





      subroutine sigmaeva(nparr,sigma)
!> evaluate the sigma values for the MCMC exploring based on the free 
!>  parameters in the MAMPOSSt subroutine      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
       parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,q0=-0.64,clight=2.99792458d5)
      dimension sigma(nparr)
      include 'paramsoptS.i'

      r200g=r200
      !r200
      if (ipar(1).eq.0) then
         sigma(1)=0.0
      else
         sigma(1)=0.1
         if(kmp.lt.6) sigma(1)=0.05 !if no MG, less degeneracy, smaller step
         
         if (kpro.eq.1) then
          sigma(1)=(dlog(r2up)-dlog(r2low))*0.05d0
          if (kbdt.eq.1) sigma(1)=(dlog(r2up)-dlog(r2low))*0.01d0
         endif
         
      endif

      !dA3 (Refracted gravity model, gamma exp for gNFW)
      if (ipar(10).eq.0) then
         sigma(10)=0.0
      else
         sigma(10)=0.05
         if (kpro.eq.1) then
           sigma(10)=(dlog(dA3up)-dlog(dA3low))*0.04d0
           if (kbdt.eq.1) sigma(10)=(dlog(dA3up)-dlog(dA3low))*0.01d0
           if (kmp.eq.13) sigma(10)=((dA3up)-(dA3low))*0.04d0
         endif
         
      endif      

      !rc 
      if (ipar(2).eq.0) then
         sigma(2)=0.0
      else
         sigma(2)=0.1
         if (kpro.eq.1) then
          sigma(2)=(dlog(rcup)-dlog(rclow))*0.05d0
          if (kbdt.eq.1) sigma(2)=(dlog(rcup)-dlog(rclow))*0.01d0
         endif
      endif
      !rs
      
      if (ipar(3).eq.0) then
         sigma(3)=0.0
      else
         sigma(3)=0.03
         if (kpro.eq.1) then
          sigma(3)=(dlog(rsup)-dlog(rslow))*0.05d0
          if (kbdt.eq.1) sigma(3)=(dlog(rsup)-dlog(rslow))*0.01d0
         endif
      endif
      !cbe
      if (ipar(4).eq.0) then
         sigma(4)=0.0
      else
         sigma(4)=0.04
         if (kpro.eq.1) then
          sigma(4)=(dlog(bup)-dlog(blow))*0.05d0
          if (kbdt.eq.1) sigma(4)=(dlog(bup)-dlog(blow))*0.01d0
         endif
       
      endif

      if (ipar(5).eq.0) then
         sigma(5)=0.0
      else
         sigma(5)=0.2
         
         !case chameleon fR
         if(kmp.eq.9.and.nhone.lt.0) sigma(5)=1
         !if(kmp.eq.10) sigma(5)=0.1 !not used anymore
         !case refracted grav
         !if (kmp.eq.13) sigma(5)=0.01
         if (kpro.eq.1) then
          dfact = 0.08d0
          if (kbdt.eq.1) then
            dfact = 0.04d0
            if (kmp.eq.18) dfact = 0.01d0 
            if (kmp.eq.13) dfact = 0.01d0
          endif
          sigma(5)=(dlog(tmup)-dlog(tmlow))*dfact
          !In the case of RG, explore in linear space
          if (kmp.eq.13) then
             dfact = 0.01d0
             sigma(5)=((tmup)-(tmlow))*dfact
          endif
          !NOW CHANGE
          if (kmody.eq.1.or.kmp.eq.17) then
            if(nLogP.eq.1.or.nhone.gt.0) then
             sigma(5) = dfact !In this case it shoud run from 0 to 1
            else
             
             sigma(5)=(dlog10(tmup)-dlog10(tmlow))*dfact
            endif
          endif

         endif

      endif
       !!!!
      !second MG parameter
      if (ipar(6).eq.0) then
         
         sigma(6)=0.0
      else
         sigma(6)=0.02
         if (kpro.eq.1) then
          dfact = 0.05d0
          if (kbdt.eq.1) then
            dfact = 0.02d0
            if (kmp.eq.18) dfact = 0.01d0
          endif
          if (kmp.eq.13) then
           dfact = 0.005d0
           sigma(6)=((scrup)-(scrlow)) * dfact 
          else
           sigma(6)=(dlog(scrup)-dlog(scrlow))* dfact 
          endif 
          if (kmody.eq.1.or.kmp.eq.17) sigma(6) = 0.08d0 !In this case it shoud run from
                                            ! 0 to 1
         endif
         
      endif
      
      !cbe0
      if (ipar(7).eq.0) then
         sigma(7)=0.0
      else
         sigma(7)=0.04
         if (kpro.eq.1) then
          sigma(7)=(dlog(cb0up)-dlog(cb0low))*0.05d0
          if (kbdt.eq.1) sigma(7)=(dlog(cb0up)-dlog(cb0low))*0.02d0
         endif
      endif
      !rjaffe
      if (ipar(8).eq.0) then
         sigma(8)=0.0
      else
         sigma(8)=0.04
         if (kpro.eq.1) then
          sigma(8)=(dlog(rjaup)-dlog(rjalow))*0.05d0
         endif
      endif
      !xmstar
      if (ipar(9).eq.0) then
         sigma(9)=0.0
      else
         sigma(9)=0.04
         if (kpro.eq.1) then
          sigma(9)=(dlog(xmaup)-dlog(xmalow))*0.01d0
          if (kmp.eq.13) sigma(9)=(dlog(xmaup)-dlog(xmalow))*0.01d0
         endif 
      endif
      !Rbeta parameter (Work in progress!!)
      if (ipar(11).eq.0) then
         sigma(11)=0.0
      else
         sigma(11)=0.03
         if (kpro.eq.1) then
          sigma(11)=(dlog(rbetaup)-dlog(rbetalow))*0.05d0
          if (kbdt.eq.1) sigma(11)=(dlog(rbetaup)-dlog(rbetalow))*0.01d0
         endif
      endif

      return
      end
c************* prior function for the mcmc run *************************
      subroutine priord(x,y,z,q,w,s,da,b2,rb,rj,xma,pout) 
      !(r200n,rcn,rsn,cbn,tmassn,scrn,beta2n,pout1)
      USE ieee_arithmetic
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
      external fmg
      include 'paramsoptS.i'

        IF (ieee_support_inf(rinf)) THEN
         rinf = ieee_value(rinf,  ieee_negative_inf)
        END IF

       if(kpro.eq.0) then
        if (x.lt.1.0d0.or.x.gt.3.9d0) then
          pout=rinf
          return
        endif
        
        if(y.lt.0.05.or.y.gt.3.9.or.z.lt.0.05.or.z.gt.3.9) then
         pout=rinf
         return 
        endif
        
        if (q.lt.0.2.or.q.gt.7.00) then
         pout=rinf
         return
        endif 
         
        if (da.lt.0.001.or.da.gt.3.00) then
          pout=rinf
          return
        endif         
         
        if (b2.lt.0.2.or.b2.gt.7.00) then
         pout=rinf
         return
        endif
        
        if(kbg.eq.1.or.kmp.eq.13) then
         if (rj.lt.0.01.or.rj.gt.1.00) then
          pout=rinf
          return
         endif
        
         if (xma.lt.0.1.or.xma.gt.20.00) then
          pout=rinf
          return
         endif
        
        endif
        
        if (kmp.eq.8) then
         if (w.gt.12.or.w.lt.-0.6) then
          pout=rinf
          return
         endif
         
        elseif (kmody.eq.1) then
         if (w.gt.150.or.w.lt.0.001) then
          
          pout=rinf
          return
         endif
         
         if(nhone.gt.1) then
          if(s.gt.29.or.s.lt.1e-2) then
            pout=rinf
            return
          endif
         endif
        endif
        pout= 0.0
        return
        
       else
c ======================================================================       
        !R200
        if (x.lt.r2low.or.x.gt.r2up) then
         pout=rinf
         return
        endif
        !rc, rs
        if(y.lt.rclow.or.y.gt.rcup.or.z.lt.rslow.or.z.gt.rsup) then
         pout=rinf
         return 
        endif
        !Beta 
        if (q.lt.blow.or.q.gt.bup) then
         pout=rinf
         return
        endif
        !beta0
        if (b2.lt.cb0low.or.b2.gt.cb0up) then
         pout=rinf
         return
        endif
        
        !tmass depending on the MG model   
        select case(kmp)
         case(9,11,12,15,16,17)
          if (nhone.gt.0.and.nLogP.eq.1) then
            if (w.gt.tmup.or.w.lt.tmlow) then
             pout=rinf
             return
            endif
          else 
            dlw = dlog10(w)
            if (dlw.gt.dlog10(tmup).or.dlw.lt.dlog10(tmlow)) then
             pout=rinf
             return
            endif
          endif   
        case default
          if (w.gt.tmup.or.w.lt.tmlow) then
             pout=rinf
             return
            endif
        end select
          
                
        if (rj.lt.rjalow.or.rj.gt.rjaup) then
          pout=rinf
          return
        endif
      
         if (xma.lt.xmalow.or.xma.gt.xmaup) then
          !write(*,*) 'here', xm,xmalow,xmaup
          pout=rinf
          return
         endif
         
         
        if (da.lt.dA3low.or.da.gt.dA3up) then
           pout=rinf
           return 
        endif 

        if (rb.lt.rbetalow.or.rb.gt.rbetaup) then
           pout=rinf
           return 
        endif 
         
        if(s.gt.scrup.or.s.lt.scrlow) then
           pout=rinf
           return
        endif

         
        pout = 0.0
        !Flat prior in log 
        !BE CAREFUL: When running with A1 = 0!!!
        !if (nhone.le.0.and.nLogP.eq.0.and.ntmass.gt.0) then
        !  pout = pout -dlog(w)
        !endif
        return
       endif
       icoun=0
       if (kmp.eq.8) then
       !in this model masses can be negative. In this case the sample
       !is automatically discharged
       fac200=1./(dlog(1.d0+x/z)-(x/z)/(1.d0+x/z))
         do i=1,4000
            t=0.01*i
            
            abey=w*t*t*(z-t)/(z+t)**3
            xm=((dlog(1.d0+t/rs)-t/z/(1.d0+t/rs))+abey/4)*fac200
            
            if (fmg(t,w,s).lt.0.or.xm.lt.0) icoun=icoun+1
         enddo
         if (icoun.gt.0) pout=rinf
        endif
       end


      

c =========== GAUSSIAN PRIOR =========================================== 
!> \ingroup MAM     
      subroutine priorG(x,y,z,q,w,s,da,b2,rb,rj,xma,pout) 
!>  @author L. Pizzuti      
!>  @details This subroutine computes the multivariate Gaussian prior 
!!  for the set of input free parameters.
!!  It requires the definition of kmp and of all the bounds in input_pars.
!>  For each parameter, the mean is given by 0.5*(xlow +xup) and the std
!!  by 0.5*(xup - xup).
!>   @param[in] x  value of r200
!>   @param[in] y value of rc
!>   @param[in] z value of rs
!>   @param[in] q value of cbe
!>   @param[in] w value of tmass
!>   @param[in] s value of screening
!>   @param[in] da value of A3
!>   @param[in] b2 value of cbe0
!>   @param[in] rb value of rbeta
!>   @param[in] rj value of rjaffe
!>   @param[in] xma value of xmstar

!>   @param[out] pout value of the log probability

      USE ieee_arithmetic
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      
      include 'paramsoptS.i'  
      
      !handle infinities
      IF (ieee_support_inf(rinf)) THEN
         rinf = ieee_value(rinf,  ieee_negative_inf)
      END IF
      
      !compute log prior for r200
      r2mean = 0.5 * (r2low + r2up)
      r2sig = 0.5 * (r2up - r2low)
      
      if (r2sig.eq.0) then
       dlogpr2 = 0
      else
       dlogpr2 = -0.5 * ((x - r2mean) / r2sig) ** 2
     &           - dlog(r2sig) - 0.5 * dlog(2.0 * pig)
      endif
      
      !compute log prior for rc
      rcmean = 0.5 * (rclow + rcup)
      rcsig = 0.5 * (rcup - rclow)
      if (rcsig.eq.0) then
        dlogprc = 0.0
      else
       dlogprc = -0.5 * ((y - rcmean) / rcsig) ** 2
     &           - dlog(rcsig) - 0.5 * dlog(2.0 * pig)
      endif
      
      !compute log prior for rs
           
      rsmean = 0.5 * (rslow + rsup)
      rssig = 0.5 * (rsup - rslow)
      if (rssig.eq.0) then 
       dlogprs = 0.0
      else 
       dlogprs = -0.5 * ((z - rsmean) / rssig) ** 2
     &           - dlog(rssig) - 0.5 * dlog(2.0 * pig)      
      endif

      !compute log prior for rbeta
      rbmean = 0.5 * (rbetalow + rbetaup)
      rbsig = 0.5 * (rbetalow - rbetaup)

      if (rbsig.eq.0) then 
       dlogprb = 0.0d0
      else
       dlogprb = -0.5 * ((rb - rbmean) / rbsig) ** 2
     &           - dlog(rbsig) - 0.5 * dlog(2.0 * pig) 
      endif

      !compute log prior for cbe
      cbmean = 0.5 * (blow + bup)
      cbsig = 0.5 * (bup - blow)
      
      if (cbsig.eq.0) then 
       dlogpbt = 0.0
      else 
       dlogpbt = -0.5 * ((q - cbmean) / cbsig) ** 2
     &           - dlog(cbsig) - 0.5 * dlog(2.0 * pig)       
      endif
      !if requested (1), Beta has still flat priors
      if (ipriorBeta.eq.1) then
       dlogpbt = 0.0
       if (q.lt.blow.or.q.gt.bup) then
         pout=rinf
         return
       endif
      endif 
      
      !compute log prior for cbe2
      cb2mean = 0.5 * (cb0low + cb0up)
      cb2sig = 0.5 * (cb0up - cb0low)
      
      if (cb2sig.eq.0) then 
       dlogpb2 = 0.0
      else       
       dlogpb2 = -0.5 * ((b2 - cb2mean) / cb2sig) ** 2
     &           - dlog(cb2sig) - 0.5 * dlog(2.0 * pig)       
      endif
      !if requested (1), Beta0 has flat priors 
      if (ipriorBeta0.eq.1) then 
        dlogpb2 = 0.0
        if (b2.lt.cb0low.or.b2.gt.cb0up) then
         pout=rinf
         return
        endif
      
      endif 
      
      !compute log prior for tmass
      tmean = 0.5 * (tmlow + tmup)
      tmsig = 0.5 * (tmup - tmlow)

      if (tmsig.eq.0) then
       dlogptm = 0.0
      else  
       dlogptm = -0.5 * ((w - tmean) / tmsig) ** 2
     &           - dlog(tmsig) - 0.5 * dlog(2.0 * pig)   
      
       if (w.lt.0.0.and.kmp.ne.8) then
       !tmass cannot be negative  
        pout = rinf
        return
       endif  
      endif
      
      !compute log prior for rj (BCG)
      rjmean = 0.5 * (rjalow + rjaup)
      rjsig = 0.5 * (rjaup - rjalow)
      if (rjsig.eq.0) then
       dlogprj = 0.0
      else  
       dlogprj = -0.5 * ((rj - rjmean) / rjsig) ** 2
     &           - dlog(rjsig) - 0.5 * dlog(2.0 * pig)  
      endif

      !compute log prior for xmstar (BCG)
      xmamean = 0.5 * (xmalow + xmaup)
      xmasig = 0.5 * (xmaup - xmalow)
      if (xmasig.eq.0) then
       dlogpxma = 0.0
      else  
       dlogpxma = -0.5 * ((xma - xmamean) / xmasig) ** 2
     &           - dlog(xmasig) - 0.5 * dlog(2.0 * pig) 
      endif
      !compute log prior for A3
      damean = 0.5 * (dA3low + dA3up)
      dasig = 0.5 * (dA3up - dA3low)

      if (dasig.eq.0) then
       dlogpda = 0.0
      else  
       dlogpda = -0.5 * ((da - damean) / dasig) ** 2
     &           - dlog(dasig) - 0.5 * dlog(2.0 * pig)

       if (kmp.eq.17) then
       !avoid unphysical values for gamma as exponent in gNFW models
        if (da.lt.1e-3.or.da.gt.1.999) then
         pout = rinf
         return
        endif
       endif
      endif
      
      !compute log prior for screen
      scmean = 0.5 * (scrlow + scrup)
      scsig = 0.5 * (scrup - scrlow)

      if (scsig.eq.0) then
       dlogpsc = 0.0
      else  

       dlogpsc = -0.5 * ((s - scmean) / scsig) ** 2
     &           - dlog(scsig) - 0.5 * dlog(2.0 * pig)

       if (s.lt.0.0.and.kmp.ne.8) then
        !screen cannot be negative  
        pout = rinf
        return
       endif 
      endif
      pout1 = dlogpr2 + dlogpbt + dlogpb2 + dlogprs + dlogprc 
      pout2 = dlogpsc + dlogpda + dlogpxma + dlogprj + dlogptm 
      
      !total log probability
      pout = pout1+pout2
      return
      end 
       
c     criteria for acceptance in the MCMC run
       subroutine acceptance(x, x_new,eval,nseed)
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        logical*4 eval
        if (x_new.gt.x) then
         eval=.True.
        else
          accept=r_uniform_01 (nseed)
          
        ! Since we did a log likelihood, 
        ! we need to exponentiate in order to compare to the random number
        ! less likely x_new are less likely to be accepted
          if (dlog(accept).lt.(x_new-x)) then
            eval=.True.
          else 
            eval=.False.
          endif

        endif
        return 
       end


c***********************************************************************
c Lensing Monte Carlo for Beyond Horndeski models **********************       
       


!>   \addtogroup WL Lensing Simulation
!>   @details This set of routines and functions are needed to
!!    compute the lensing
!!    likelihood in Vainsthein screening and Chameleon screening


!> \ingroup WL
        function fmg(x,ya,yb)
!>  @author L. Pizzuti      
!>  @details This function computes the analytic expression 
!!  for the radial dependence of the projected surface density
!!  of a NFW profile in DHOST
!!  gravity \f$ f_\text{mg}(x,Y_1,Y_2) \f$, defined by
!!  \f$ \Sigma_\text{mg}(R)=2\rho_\text{s}\,r_\text{s}\times f_\text{mg} \f$.
!>  The expression has been obtained by integrating eq. (36) with the surface
!!  density (41) in Pizzuti et al. (2021).
!>
!>   @param[in] x  dimensionless radius \f$ x=r/r_\text{s} \f$
!>   @param[in] ya first modified gravity parameter
!>   @param[in] yb second modifed gravity parameter
!>
!>   Example of usage:
!>
!>   \code{.f}
!>   ya=0.5d0
!>   yb=-2.0d0
!>   write(*,*) fmg(1.2d0,ya,yb)
!>   \endcode
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
    !effective polynomials (Pizzuti et al., 2021)
          apu=(8*x**4+x**2*(10*ya+15*yb-16)+5*ya-15*yb+8)/2
          apd=(x**4*(2*ya+5*yb+8)+x**2*(11*ya+5*yb-16)+2*(ya-5*yb+4))/2
          
          if (x.lt.1.0d0) then
           un=((-1+x**2)*(2*apu)+dsqrt(1-x**2)*(2*apd)*dacosh(1/x))
           afmg=un/(8*(x**2-1)**4)

          else if (x.gt.1.0d0) then
          
           uno=(2*apu)/(x**2-1.0d0)**3
           due=-((2*apd)*dacos(1/x))/dsqrt((x**2-1.0d0)**(7))
          
            afmg=(uno+due)/8 

          else if(x.eq.1.0d0) then
            afmg=-(ya + 7*(-2 + yb))/(21*2)
          endif
          fmg=afmg
          return 
         end

c***********************************************************************
ccontribution of the modified gravity average surface density profile
cfactor rs already included
        function ftoint(w)
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
          include 'paramsoptS.i'
          
          stron=1.0d0-dsqrt(1.0d0-(Rx*w)**2)
          
          ftoint=(4+Aa*w+Bb*w**2)*stron/(w*(1+w)**4)

           return
          end

!> \ingroup WL
        function gmg(x,ya,yb)
!>  @author L. Pizzuti      
!>  @details This function computes the expression 
!!  for the radial dependence of the projected mass density
!!  of a NFW profile in DHOST
!!  gravity \f$ g_\text{mg}(x,Y_1,Y_2) \f$, defined by
!!  \f$ M_\text{P,mg}(R)=2\rho_\text{s}\,r_\text{s}\times g_\text{mg} \f$.
!>  The expression has been obtained by integrating eq. (37) with the surface
!!  density (41) in Pizzuti et al. (2021).
!>
!>   @param[in] x  dimensionless radius \f$ x=r/r_\text{s} \f$
!>   @param[in] ya first modified gravity parameter
!>   @param[in] yb second modifed gravity parameter
!>
!>   Example of usage:
!>
!>   \code{.f}
!>   ya=0.5d0
!>   yb=-2.0d0
!>   write(*,*) gmg(1.2d0,ya,yb)
!>   \endcode
          
        
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
          include 'paramsoptS.i'
          external ftoint
          Aa=8-2*ya-5*yb
          Bb=4+ya-5*yb
          Rx=x
    
          amic=-x*(8+x*(16-ya+5*yb+x*(8+ya+5*yb)))
          am=amic/(2*(1+x)**3)+4*dlog(1+x)
          rmax=7.0e17
          errabs=0.0d0
          errrel=0.0001d0
          ux=1.0d0/Rx
          
          call dgaus8 (ftoint,0.00000000001d0,ux,errrel,risl, IERR) 

          gmg=(risl+am)/(2*x**2)
         return
         end



!> \ingroup WL
        function f_nfw(x)
!>  @author L. Pizzuti      
!>  @details This function computes the analytic expression 
!!  for the radial dependence of the projected surface density
!!  of a NFW profile in GR.
!!  \f$ f_\text{NFW}(x) \f$, defined by
!!  \f$ \Sigma_\text{NFW}(R)=2\rho_\text{s}\,r_\text{s}\times f_\text{NFW} \f$.
!>  See e.g. Bartelmann (1996).
!>   @param[in] x  dimensionless radius \f$ x=r/r_\text{s} \f$
!>   Example of usage:
!>   \code{.f}
!>   write(*,*) f_nfw(1.2d0)
!>   \endcode
         
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
          include 'paramsoptS.i'
         if (x.lt.1.0d0) then
          fn=(-1+2*datanh(dsqrt((1-x)/(1+x)))/(dsqrt(1-x**2)))/(1-x**2)
         else if (x.eq.1.0d0) then
          fn=1.0d0/3
         else
          fn=(1-2*datan(dsqrt((x-1)/(1+x)))/(dsqrt(x**2-1)))/(x**2-1)
         endif
        f_nfw=fn
        return
        end



!> \ingroup WL
        function g_nfw(x)
!>  @author L. Pizzuti      
!>  @details This function computes the analytic expression 
!!  for the radial dependence of the projected mass density
!!  of a NFW profile in GR.
!!  \f$ g_\text{NFW}(x) \f$, defined by
!!  \f$ \M_\text{P,NFW}(R)=2\rho_\text{s}\,r_\text{s}\times g_\text{NFW} \f$.
!>  See e.g. Bartelmann (1996).
!>   @param[in] x  dimensionless radius \f$ x=r/r_\text{s} \f$
!>   Example of usage:
!>   \code{.f}
!>   write(*,*) g_nfw(1.2d0)
!>   \endcode
          implicit real*8 (a-h,o-z)
          implicit integer*4 (i-n)
          include 'paramsoptS.i'        
        
         if (x.lt.1.0d0) then
           un=(2*datanh(dsqrt((1-x)/(1+x)))/(dsqrt(1-x**2))+dlog(x/2))
           gn=2*un/(x**2)
         else if (x.eq.1.0d0) then
          gn=2*(1+dlog(1.0d0/2))
         else
           un=(2/(dsqrt(x**2-1))*atan(dsqrt((x-1)/(1+x)))+dlog(x/2))
           gn=2*un/(x**2)
         endif
         g_nfw=gn
         return
        end 

!***********************************************************************
       subroutine test_bound(xmed,xlow,xup)
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
    
        include 'paramsoptS.i'
        if(xlow.ge.0.995*xmed) then
          xlow=xmed-dabs(xmed)/4.0
          
          write(*,199) xlow
 199      format(' WARNING:'/ 
     &     ' lower bound of parameter range too close to guess value',/
     &     'automatically reset to',f9.4 )
        endif
        if(xup.le.1.005*xmed) then
          xup=xmed+dabs(xmed)/4

          write(*,198) xup
 198      format(' WARNING:'/ 
     &     ' upper bound of parameter range too close to guess value',/
     &     'automatically reset to',f9.4 )
        endif
        return
       end

c************ Cosmology ************************************************
      function Ez(z) !E(z)^-1 time evolution of the Hubble factor
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)

      parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
      parameter (sigcrit_ave=464.0515847354301,flens=1.3544638077221953)
      parameter (sigmael=0.25, sigmalls=0.005, ngal=30)

      include 'paramsoptS.i'
       
       Ez=1/dsqrt(Oml*(1+z)**3+Oll)
       return
      end
      
c comoving distance (Mpc) **********************************************
       function Chi(z)
       
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0,grav=4.302e-9)

        parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
        
        include 'paramsoptS.i'
        
        external Ez
        
        errrel=0.0001d0

        call dgaus8 (Ez,0.0d0,z,errrel,ach, IERR)

        chi=ach*cl/h0lens
        return
        end  
        
        
        function Da(z)
        
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0,grav=4.302e-9)

        parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         
        external Chi
        
         Da=Chi(z)/(1+z)
        return 
        end
        

C************** this is built assuming a cluster at z=0.44 *************
c***********************************************************************
c         convergence maps
c         true GR averaged over all redshift


        function ak_true(R)
        
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)

         parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         parameter (sigcrit_ave=464.0515847354301, zc=0.44)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez,f_nfw
         
         !zc=za
         !if (za.eq.0.0d0) zc=0.44 
         ct=r200t/rst !concentration
         fac200t=dlog(1+ct)-ct/(1+ct)
        
         x=R/rst
         facfront=200*sigcrit_ave*r200t**3*h0lens**2/Ez(zc)**2
         facfront=facfront/(rst**2*cl**2*fac200t)
         ak_true=facfront*f_nfw(x)
         
         return
         
        end

c       MODIFIED averaged over all redshift    
        function ak_mg(R,r2,rss,ya,yb)
        
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)

         parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         parameter (sigcrit_ave=464.0515847354301, zc=0.44)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez,fmg 
         
         !zc=za
         !if (za.eq.0.0d0) zc=0.44 
         x=R/rss
         c2=r2/rss
         fac200=dlog(1+c2)-c2/(1+c2)
         facfront=200*sigcrit_ave*r2**3*h0lens**2/Ez(zc)**2
         
         facfront=facfront/(rss**2*cl**2*fac200)
         ak_mg=facfront*fmg(x,ya,yb)
         return
        end
c ************* shear maps maps
c*************** true GR averaged over all redshift
       function gamma_true(R)
       
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)

         parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         parameter (sigcrit_ave=464.0515847354301, zc=0.44)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez,f_nfw,g_nfw
         
         !zc=za
         !if (za.eq.0.0d0) zc=0.44 
         ct=r200t/rst !concentration
         fac200t=dlog(1+ct)-ct/(1+ct)
        
        
         x=R/rst
         facfront=200*sigcrit_ave*r200t**3*h0lens**2/Ez(zc)**2
         facfront=facfront/(rst**2*cl**2*fac200t)
         
         gamma_true=facfront*(g_nfw(x)-f_nfw(x))
         
         return
        end


c************* MODIFIED averaged over all redshift *********************
        function gamma_mg(R,r2,rss,ya,yb)
        
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)

         parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         parameter (sigcrit_ave=464.0515847354301, zc=0.44)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez,fmg, gmg 
         !zc=za
         !if (za.eq.0.0d0) zc=0.44 
         x=R/rss
         c2=r2/rss
         fac200=dlog(1+c2)-c2/(1+c2)
         facfront=200*sigcrit_ave*(r2**3)*h0lens**2/Ez(zc)**2
         
         facfront=facfront/(rss**2*cl**2*fac200)

         gamma_mg=facfront*(gmg(x,ya,yb)-fmg(x,ya,yb))    

         return
        end
!***********************************************************************


!> \ingroup WL
       function Plens(r2x,rsx)
!>  @author L. Pizzuti      
!>  @details This funcion computes the logarithm of a bivariate gaussian 
!!   probability distribution centered on the values r2x and rsx,
!!   with relative standard deviations and correlation
!!   given externally by wr2, wr, wx respectively. "Relative" means that 
!!   e.g. the variance of r2x is given by r2x*wr2. 
!>   @param[in] r2x REAL*8, central value of the first parameter
!>   @param[in] rsx REAL*8, central value of the second parameter

!>   Example of usage:
!>
!>   \code{.f}
!>   wr2=0.2
!>   wr=0.1
!>   wx=0.5
!>   write(*,*) Plens(r2x,rsx)
!>   \endcode 

        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)
         include 'paramsoptS.i'
         
         rssg=rst
         r22g=r200t
         
         
         if (rst.le.0.03d0) rssg=rsg
         if (r200t.le.0.1d0) r22g=r200g

         corr=wx
         sigmars=rssg*wr
         sigmar2=r22g*wr2
         p1= (rsx-rssg)*(rsx-rssg)/(sigmars*sigmars)
         p2=(r2x-r22g)*(r2x-r22g)/(sigmar2*sigmar2)
         p3=-2*corr*(rsx-rssg)*(r2x-r22g)/(sigmars*sigmar2)
         Plens=-(p1+p2+p3)/(2-2*corr*corr)
         return 
         
       end
       



!>     \ingroup WL
        function gt_true(R)
!>     @authors L. Pizzuti, G. Mamon
!>     @details This routine computes the average tangential shear profile,
!!     eq. (34) in Pizzuti et al., 2021, at the projected radius R for a NFW
!!     model described by the external parameters `r200t, rst`.
!>     @param[in] R REAL*8, projected radius from the center of the mass distribution

!>     Additional requirements are the mass profile parameters `r200t, rst`  
!!     defined in the COMMON block `paramsoptS.i`. 

!>    Example of usage:
!>
!>   \code{.f}
!>   r200t=1.41
!>   rst=0.3
!> 
!>   write(*,*) gt_true(1.2d0)
!>   \endcode 
  
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)
         parameter (cl=2.99792458e5)
         parameter (sigcrit_ave=464.0515847354301, zc=0.44)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez, gamma_true, ak_true 
         
         stro=(1-flens*ak_true(R))
         gt_true=gamma_true(R)/stro
         return
        end
                      
c************************ Modified: ************************************                      
!>    \ingroup WL
        function gt_mod(R,r2,rss,ya,yb)
!>     @authors L. Pizzuti, G. Mamon
!>     @details This routine computes the average tangential shear profile,
!!     eq. (34) in Pizzuti et al., 2021, at the projected radius R for a 
!!     modified NFW profile in Vainsthein screening (eq. (41) in Pizzuti et al., 2021).
!>     @param[in] R REAL*8, projected radius from the center of the mass distribution
!>     @param[in] r2 REAL*8, virial radius r200 of the NFW model  
!>     @param[in] rss REAL*8, scale radius rs of the NFW model
!>     @param[in] ya REAL*8, first MG coupling  
!>     @param[in] yb REAL*8, second MG coupling


!>     example of usage:
!>
!>   \code{.f}
!>   r2=1.41
!>   rss=0.3
!>   ya=0.1
!>   yb=0.5
!>   write(*,*) gt_mod(1.2d0,r2,rss,ya,yb)
!>   \endcode        
         implicit real*8 (a-h,o-z)
         implicit integer*4 (i-n)
         parameter (pig=3.1415926535897932d0,grav=4.302e-9)

         parameter (cl=2.99792458e5, h0lens=67.5, Oml=0.32, Oll=0.68)
         parameter (sigcrit_ave=464.0515847354301)
         parameter (flens=1.3544638077221953)
         
         include 'paramsoptS.i'
         
         external Ez, gamma_mg, ak_mg  
         
      
         stro=(1-flens*ak_mg(R,r2,rss,ya,yb)) 
         gt_mod=gamma_mg(R,r2,rss,ya,yb)/stro
         
         return 
         
        end
     
!>  \ingroup WL
      subroutine Likelens_bh(npoint,plen) 
!>    @details This subroutine computes the simulated lensing likelihood 
!!  for a galaxy cluster in Vainsthein screening gravity.
!>     the "true" mass profile is assumed to be a NFW model.
!>     @param[in] npoint INTEGER*4, number of points where the 
!!     likelihood is evaluated
!>     @param[out] plen REAL*8, value of the -log(likelihood)

!>  The redshift of the halo is fixed to z=0.44, as an average redshift
!!     for optimal lensing observations (see Pizzuti et al., 2022) 
!>  Note that in this version of the code,
!!     the lensing simulation is assumed to be at a fixed cosmology with
!!     \f$ H_0=67.5 \f$ Mpc, \f$\Omega_m=0.32, \Omega_\Lambda=0.68 \f$
!>           
!>    It requires additional external parameters from `pars_test.txt`, from `Options.txt` and the 
!!    COMMON block `paramsoptS.i`. In particular: 
!>  - the mass profile parameters 
!! `r200, rs, tmass, screen`
!>  - maximum radius of the kinematic analysis `rupin`
!>  - lensing "true" value of the mock profile `rst, r200t`
!>  - lensing uncertainties: sigma_ellipticity, sigma_lss,
!!    Number of galaxies per arcmin^2 (min 10)

!>
!>    Example of usage:

!>     \code{.f}
!>      r200=1.41
!>      rs=0.3 
!>      cbe=1.2
!>      tmass=1.2
!>      screen=0.4
!>      rupin=r200
!>      rst=0.3    !lensing "true" value of the scale radius
!>      r200t=1.41 !lensing "true" value of the virial raidius
!>      wr2=0.25  !intrinsic ellipticity
!>      wr=0.005  !lss noise
!>      wx=30     !number of galaxies per arcmin^2
!>      call  Likelens_bh(10,plen)
!>      write(*,*) plen
!>    \endcode
       
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)

      parameter (cl=2.99792458e5)
      parameter (sigcrit_ave=464.0515847354301,flens=1.3544638077221953)
      parameter (zc=0.44)
      dimension xth(npoint+1),rth(npoint+1),xbin(npoint),error(npoint)
      
      include 'paramsoptS.i'
      include 'sr.i'
      include 'vlos.i'
      include 'datarv.i'
      include 'tsallis.i'
      include 'probs.i'
      external ak_mg,ak_true,f_nfw,g_nfw,gt_true,gt_mod,Da
      external gamma_true, gamma_mg
      
      if (rst.le.0.03d0) rst=rsg
      
      if (r200t.le.0.1d0) r200t=r200g


      sigmael=wr2
      sigmalls=wr
      Oll=Olam
      Oml=Omegam
      h0lens=h0
      ngal=nint(wx)
      if(ngal.lt.10) ngal=10
      xr2min=rupin*0.12
      xr2max=rupin*2.9
      !Minimum and maximum angular radii for the lensing simulation
      xtmax=xr2max*(180*60.0d0)/(pig*Da(zc))     
      xtmin=xr2min*(180*60.0d0)/(pig*Da(zc))       

      do i=1,npoint+1
          
         xth(i) = dlog(xtmin)+dfloat(i-1)*
     &    dlog(1.001*xtmax/xtmin)/dfloat(npoint+1-1)
         xth(i)=dexp(xth(i))
         rth(i)=xth(i)*Da(zc)*pig/(180.0d0*60.0d0)
      enddo
      psum=0.0d0 !initialize the chisquare
      do k=1,npoint
       estat=(sigmael**2)/(pig*(xth(k+1)**2-xth(k)**2)*ngal)
       error(k)=dsqrt(estat+sigmalls)
       xbin(k)=(rth(k)+(rth(k)-rth(k))/2) 
       
       pp=gt_mod(xbin(k),r200,rs,tmass,screen)
       scarto=(gt_true(xbin(k))-pp)/error(k)
       psum=psum+scarto**2/2 !factor 1/2 
      enddo 
       plen=psum
      return
       
      end
C********************************* BCG ADDITION ************************


      function fintegrand(xl)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptS.i'
 !     external fmass,fbetaint,frho
      x=dexp(xl)
      fintegrand=frho(x)*fmass(x)/x*fbetaint(x) !integral in dlog!
      return
      end



!>  \ingroup MAM
      function fmass(t)
!>    @details This subroutine computes the total mass profile of the  
!!    cluster assuming spherical symmetry. 
!>     @param[in] t REAL*8, value where the function is evaluated
!>           
!>    It requires additional external parameters from `pars_test.txt`, from `Options.txt` and the 
!!    COMMON block `paramsoptS.i`. In particular: 
!>  - the mass profile switch `kmp`
!>  - the mass profile parameters 
!! `r200, rs, tmass, screen`
!>  - the value of velocity dispersion at r200, `v200`
!>  - the switch for the inclusion of the BCG `kbcg`
!>  - if `kbcg=1`, the BCG-mass profile parameters: `rjaf,xmstar,xlumbcg`

!>
!>    Example of usage:

!>     \code{.f}
!>      kbcg=1
!>      r200=1.41
!>      v200=r200*10*70
!>      kmp=10
!>      rs=0.3 
!>      cbe=1.2
!>      tmass=1.2
!>      screen=0.0d0
!>      xmstar=4.5d0
!>      xlumbcg=492000000000.0d0
!>      rjaf=0.039
!>      x=0.25d0
!>      write(*,*) fmass(x)
!>    \endcode     

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex gamma,z,hypgeo,z200m, hfz
      dimension rsvalues(28),r100values(28)
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)
      
      include 'paramsoptS.i'
      include 'barions.i'
      
      data rsvalues/0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,
     ,     0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,
     ,     0.30,0.35,0.40,0.45,0.50,1.00/
      data r100values/1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,
     ,     1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,
     ,     1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476/
      external gammarec
      external hypgeo
      external fbars
c    
      gm200=r200*v200*v200
      
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
            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200
         endif


c     attempt to insert generalized NFW model         
      elseif (kmp.eq.10) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)

       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
       
       xm =  hfx/hfy*(t/r200)**(3-dA3) 
c     Clustering Dark Energy 
c     dA3 = exponent of gNFW
c     dA2 = f(cs,w)
c     dA1 = cs2 (in unit of c**2)        
      elseif (kmp.eq.18) then

c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -t/rs, hfx ) 
c       call hygfx ( 3-tmass, 3-tmass, 4-tmass, -r200/rs, hfy)
       
c       xm =  hfx/hfy*(t/r200)**(3-tmass)

       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -t/rs, hfx ) 
       call hygfx ( 3-dA3, 3-dA3, 4-dA3, -r200/rs, hfy)
       
       xm=hfx/hfy*(t/r200)**(3-dA3)*(1.0d0+3*10**(screen-3.5)*(1+tmass))
               
      elseif (kmp.eq.8)  then	
c     M(r) is the modified NFW profile with a beyond Horndeski model 
c     See Saltas+16, Sakestein+16
c     the free parameter tmass is now the coupling constant Y 
	
       fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
       abey=tmass*t*t*(rs-t)/(rs+t)**3
       xm=((dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))+abey/4)*fac200
      
      elseif (kmp.eq.9) then
c     M(r) is the modified NFW profile in chameleon screening 
c      (Pizzuti et al.,2021)    
       !fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))
       !xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200+dphidr(t/rs)
       
       xm = tM_bnfw_CS(t)
       
      elseif (kmp.eq.15) then
c     M(r) is the modified gNFW profile in chameleon screening 
c     gamma is fixed to dA3
        xm = tM_gnfw_CS(t)       

      elseif (kmp.eq.16) then
c     M(r) is the modified Isothermal profile in chameleon screening 
c     the exponent is 3/2       
       xm = tM_iso_CS(t)
      
c      elseif (kmp.eq.17) then
c       !M(r) is multicomponent general Chameleon. gNFW+Jaffe+beta+NFW
c        xm = grav*tmass_Gen(t)/gm200 ! To divide by gm200 as this is already
c                                ! including the normalization factor  
      
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


      elseif (kmp.eq.11) then
c        mass is modified Burkert in Chameleon Gravity   
         trs=t/rs
         rvrs=r200/rs
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))
         xm1=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))
         xm=xm1+dphib(trs)           

      elseif (kmp.eq.12) then
c     mass is modified Einasto in Chameleon Gravity
       !  eim=5.
         agp=3.*eim
         xgp=2.*eim*(r200/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip200)
         fac200=1./gip200
         xgp=2.*eim*(t/rs)**(1./eim)
         call dincog(agp,xgp,gin,gim,gip)
         xm1=fac200*gip
         
         xm=xm1+dphiE(t/rs)  
      
      elseif (kmp.eq.14) then
      !Dark matter profile is Boson Star
        xt = t/rs
        if (xt.le.2.0d0) then 
         xm = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
        else
         xm =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
        endif 
        
         
      elseif (kmp.eq.5) then
c
c     M(r) is Soft Isoth Sph; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
         fac200=1./(r200-rs*datan(r200/rs))
         xm=fac200*(t-rs*datan(t/rs))
      else
c
c     M(r) is Einasto m=eim; no gas contribution allowed for the time begin
c     Normalisation set to 1 at r/r200=1
c
       !  eim=5.
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

!     Jaffe M(r) and baryonic mass in galaxies and gas
       
        fbarr=0.0
        if (kbary.eq.1) fbarr=fbars(t) !inclusion of baryons 
       
        fjaffe=xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr
      
      
        !NOT NEEDED. INCLUDED IN THE GENERAL CHAMELEON CASE
        !if (kmody.eq.1) then 
        ! dpj=dphiJ(t/rjaf) !Include the chameleon field of the BCG
        ! fjaffe= xmstar*xlumbcg*(t/rjaf/(1.+t/rjaf)+dpj)+fbarr
        !endif

       !else
        fmass=fjaffe+xm*gm200/grav   !+fmassgal+fmassgas+fjaffeicl

! Diaferio's model of refracted gravity: *******************************       
        if (kmp.eq.13) then 
        
        !refracted gravity with the function of Valeria
        fmass= rgmass(t)
        !fbarr=0.0
        !if (kbary.eq.1) fbarr=fbars(t)
        
        !Covariant formulation. It is not good for the moment
    
         !fmass=(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)
                  
        elseif (kmp.eq.17) then
        !M(r) is multicomponent general Chameleon. gNFW+Jaffe+beta+NFW
         fmass = tmass_Gen(t)  
          
        endif 
 
!
      return
      end



      function fbars(x)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       parameter (pig=3.1415926535897932d0)
       include 'paramsoptS.i'
       include 'barions.i'
       
       !call hunt(rgas, ngas,x,jlo) 
       !fbars=xmgas(jlo)
       !Aggiungi profilo beta con parametri LUCA
       !Test con RG
       if (srhoga.eq.0) rhog1 = rho0gas
       if (srsga.eq.0)  rsga1 = rsgas
       
       
       if (srhost.eq.0) rhost1 = rhostar
       if (srstar.eq.0) rst1 = rstar
       
       xg = x/rsga1
       if (kalpha.eq.1) then
        !New: Gas (this is a beta profile with beta = 1)
        fbb1 = (-1.0d0/(xg**2*dsqrt(1+xg**2))+ (dasinh(xg))/xg**3)                  
       else
        !Profile with beta = 2/3 
        fbb1 = rsga1**2 * (x - rsga1*datan(xg)) / x**3
        
       endif
       
       !New: galaxies: NFW profile
       xstar = x/rst1
       fbb2 = (-1.0d0/(xstar**2*(1+xstar))+(dlog(1+xstar))/xstar**3)
       
       fbars = 4*pig*x**3*(rhost1*fbb2 + rhog1*fbb1)


       !fbars = 4*pig*rhog1*rsga1**2*(x-rsga1*datan(x/rsga1))
       return
      end
    
      !OLD FUNCTION -> Interpolates
      function fbars2(x)
      !interpolate the baryons
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       include 'paramsoptS.i'
       include 'barions.i'
       if (x.le.(rgas(1)+0.01)) then
        fbars2 = xmgas(1)
       else if (x.ge.(rgas(ngas)-0.001)) then
        fbars2 = xmgas(ngas)
       else
        call hunt(rgas, ngas,x,jlo) 
        fbars2= xmgas(jlo)
       endif 
       
       return
      end

!************ FUNCTIONS FOR INTEGRAL IN REFRACTED GRAVITY **************
      function fintRG(x)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       include 'paramsoptS.i'
       include 'barions.i'
       
       external fbars
         fbarr=0.0
         if (kbary.eq.1) fbarr=fbars(x)
        fintRG=(fbarr)/x**2 !xmstar*xlumbcg*x**2/(x+rjaf)**2+
       return 
      end 
      
      function fmodRG(x)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       parameter (grav=4.302e-9, clight=2.99792458d5)
       
       include 'paramsoptS.i'
       include 'barions.i'
       
       external fintRG
       errrel=0.0001d0

        call dgaus8 (fintRG,0.0d0,x,errrel,ares, IERR)
        
        test=xmstar*xlumbcg*x/rjaf/(x+rjaf)
        fmodRG=(1-grav*(test+ares)/clight**2)
       
       return
       end

! **********************************************************************

!     exp(2*int(beta/x dx)) - OM model / constant model
!     !!! Note that the constant part of the integration is in the MAIN pgm

      function fbetaint(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramsoptS.i'
      if (kanbcg.eq.1) then !OM
       fbetaint=rabcg*rabcg+x*x
      elseif (kanbcg.eq.0) then 
       bcgbetac=rabcg
       fbetaint=x**(2.*bcgbetac) ! beta=bcgbetac=constant at all radii 
      endif
      
      return
      end


!     density profile (nu(r) or rho*(r))
!     Jaffe or Hernquist

      function frho(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0)

      include 'paramsoptS.i'
        
      external frhoHer
!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf

!      frho=xmstarl*xlumbcg*rjaf/(4.*pi*x*x*(rjaf+x)*(rjaf+x))
       if (kdbcg.ne.1) then
        frho=rjaf/(4.*pi*x*x*(rjaf+x)*(rjaf+x))
       else
        frho=frhoHer(x)
       endif
      return
      end
   
!     density profile  (nu(r) or rho*(r))
!     Hernquist  
!     This is to test the model of Diaferio
      function frhoHer(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0)

      include 'paramsoptS.i'

!      frhoHer=xmstarl*xlumbcg*rjaf/(2.*pi*x*(rjaf+x)**3)
       frhoHer=rjaf/(2.*pi*x*(rjaf+x)**3)

      return
      end   
      
      
      function fbeta(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      include 'paramsoptS.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
       if(kanbcg.eq.0) then
        fbeta=rabcg
       else if(kanbcg.eq.1) then
        fbeta=x*x/(rabcg*rabcg+x*x)
       endif
       
      return
      end

! ********************* BCG RUN ****************************************


!> \ingroup MAM
      function frpHer(tt) 
!>  @author L. Pizzuti      
!>  @details This function computes the analytic expression 
!!  for the radial dependence of the projected surface density
!!  of a Hernquist profile in GR excluding the normalization. 
!!  This is needed for the BCG, it is equivalent to sigmar2(tt).
!!  It requires the inclusion of 'paramsoptS.i' for `rjaf`
!>
!>   @param[in] R radius in Mpc
!>
!>   Example of usage:
!>
!>   \code{.f}
!>   rjaf=0.03d0

!>   write(*,*) frpHer(1.2d0)
!>   \endcode


      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      include 'paramsoptS.i'
      s=tt/rjaf

      if (s.gt.1) then
         xs=acos(1./s)/sqrt(s*s-1.)
         fr=((2.+s*s)*xs-3.)/(2.*pig*rjaf*rjaf*(1.-s*s)**2)
      elseif (s.lt.1) then
         xs=dlog((1.+sqrt(1.-s*s))/s)/sqrt(1.-s*s)
         fr=((2.+s*s)*xs-3.)/(2.*pig*rjaf*rjaf*(1.-s*s)**2)
      else
         xs=1.
         fr=2./(15.*pig*rjaf*rjaf)
      endif

      frpHer=fr

      return
      end




c!     projected density profile (N(R) or Sigma*(R) or I(R))
c!     Jaffe or Hernquist

      function frhoproj(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0, Dmj=0.184565d0)
      
      include 'paramsoptS.i'
      include 'datavdp.i'
    
      external frpHer
      
      
      if (kdbcg.ne.1) then
c!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf

       s=x/rjaf

       if (abs(s-1.).lt.0.001) then

         s1=0.999
         cs1=dacosh(1.d0/s1)
         s2=1.001
         cs2=dacos(1.d0/s2)
         
         ticaz1=(1.d0-(2.d0-s1*s1)*cs1/dsqrt(dabs(s1*s1-1.d0)))
         ticaz1=ticaz1/2.d0/(s1*s1-1.d0)  
         frhoproj1=pi/(4.d0*s1)-ticaz1
     
         ticaz2=(1.d0-(2.d0-s2*s2)*cs2/dsqrt(dabs(s2*s2-1.d0)))
         ticaz2=ticaz2/2.d0/(s2*s2-1.d0)               

         frhoproj2=pi/(4.d0*s2)-ticaz2

         frhoproj=(frhoproj2-frhoproj1)/(s2-s1)*(s-s1)+frhoproj1

c!     R<>rjaf

       else

         if (s.lt.1.) then
            cs=dacosh(1.d0/s)
         else
            cs=dacos(1.d0/s)
         endif
        pipp=(1.d0-(2.d0-s*s)*cs/dsqrt(dabs(s*s-1.d0)))/2.d0/(s*s-1.d0)
        frhoproj=pi/(4.d0*s)-pipp
       endif
      
       frhoproj=frhoproj/(pi*rjaf*rjaf)

      else
       !case of a Hernquist profile for the BCG
       frhoproj=frpHer(x)
      endif
      
      return
      end


      function fabbeta(xl)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1),ydev1(1),ydev2(1)

      include 'paramsoptS.i'
      include 'datavdp.i'


      ! common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
      ! common/spline/xx(500),yy(500),csc(500,3),icsc,nd
      ! common/limits/ax,rpsf

!      external fbeta

      u(1)=xl

      
      call SPLINE_CUBIC_VAL(nd,xbc,ybc,ypb2,u,s,
     &  ydev1,ydev2)
      
      s1=s(1)

      x=dexp(xl)
      b1=fbeta(x)

      a=dexp(ax)

      if (abs(1.-x/a).lt.0.001) then
         fabbeta=0.d0
      else
         fabbeta=dexp(s1)*x*x/dsqrt(x*x-a*a)*(1.-b1*a*a/(x*x))
      endif


      return
      end



c***********************************************************************
!     
!     Solving the integral in log space to get ln(nu*vr^2)
!     
      subroutine vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter (pi=3.1415927d0,grav=4.302d-9,iii=5000)
      parameter (nout=100) !200

      dimension xmnul(iii),vr2nu(iii),u(iii),s(iii)
      dimension rvdp(nvdp),vdp(nvdp),sig2i(iii)
      dimension bcgr(nvdp),bcgvdp(nvdp)

      include 'datavdp.i'      
      include 'paramsoptS.i'
      
      external fintegrand
      external fabbeta
      external frhoproj
      external frho
      external fmass
      !Try to see if the BCG should be integrated up the top or we can do less
      errabs=0.001d0
      errel=0.001d0
      bx=dlog(rtbcg/10.0d0) !rtbcg
      a1=dlog(1.d-4)
      deltalog=(bx-a1)/nout
      lout=0
      

      do i=1,nout-1
         ax=a1+deltalog*i
         res=0.d0
         call dgaus8 (fintegrand,ax,bx, errel, 
     &    res, IERR)
     
         
         if (kanbcg.eq.0) then
          rabcg=bcgbetac
          denbeta=dexp(2.*bcgbetac*ax) ! constant part of the integrated beta(r) for the model beta=bcgbetac=constant at all radii
         else if (kanbcg.eq.1) then
          denbeta=rabcg*rabcg+dexp(2.*ax) ! constant part of the integrated beta(r)
         endif
         res=res*grav/denbeta 
         
         if (IERR.gt.130) write(*,*) ' int_fvdm: ',dexp(ax),res
         if (res.eq.res) then 
            lout=lout+1
            xmnul(lout)=ax
            vr2nu(lout)=dlog(res)
            
         endif
        
      enddo
c!
c!     Spline fitting of ln(nu*vr^2)
c!
      nd=lout-2
      do k=1,nd
         xbc(k)=xmnul(k)
         ybc(k)=vr2nu(k)
      enddo
      
      call SPLINE_CUBIC_SET(nd,xbc,ybc,2,0.d0,2,0.d0,ypb2)

!     Abel projection of nu*vr^2 in nout points,
!     yielding I*VDP^2; then divide by I(R) and take the sqrt

      jout=0
      !call CPU_TIME(stbcg)
      do j=1,lout-1
         ax=a1+deltalog*j
         resp=0.0d0
         errel=0.001d0
         
         call dgaus8 (fabbeta,ax,bx,errel,resp,IERR)

         if (IERR.gt.130) write(*,*) ' int_fabbeta: ',dexp(ax),resp
         if (resp.lt.0.) resp=0.
         a=dexp(ax)
         fir=frhoproj(a)
         if (resp.eq.resp) then 
            jout=jout+1
            rvdp(jout)=a
            vdp(jout)=dsqrt(2.d0*resp/fir)

            sig2i(jout)=2.d0*resp 

            bcgr(jout)=rvdp(jout)
            bcgvdp(jout)=vdp(jout)            
c! END
         endif
      enddo
      
      !call CPU_TIME(enbcg)
      !write(*,*) lout, 'time of execution BCG = ', enbcg-stbcg
      
      iout=jout
      return
      end

      subroutine vdplik(chibcg)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0,grav=4.302d-9)
      parameter  (nout=200) !200
      
      dimension bcgr(nout),bcgvdp(nout),rvdp(nout),vdp(nout)


      include 'paramsoptS.i'
      include 'datavdp.i'

      nvdp=nout-1
         
      omegal=Olam
      omega0=Omegam
      hz=h0*sqrt(omega0*(1.+za)**3+omegal)
                                            
      xm200=100.*hz*hz/grav*r200**3


      call vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)
      chiq=0.d0
        
      rlim = 2.*rbcg(nbcg)
      do k=1,iout
        if (bcgr(k).le.rlim) imax=k
      enddo
      do j=1,nbcg
        rbcgj=rbcg(j)             
        call hunt(rvdp,imax,rbcgj,jlo)
        sbcge= vdp(jlo)+ ( vdp(jlo+1)-vdp(jlo) )/(rvdp(jlo+1) - 
     -                rvdp(jlo)) * (rbcgj - rvdp(jlo) )
  
        chiq=chiq+(sbcg(j)-sbcge)*(sbcg(j)-sbcge)/(esbcg(j)*esbcg(j))

      enddo

      chibcg = chiq/2 !The -log(Likelihood) is chiq/2

      return 
      end

      SUBROUTINE hunt(xx,n,x,jlo)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension xx(n)
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END


! Find roots of a function 
c      subroutine find_root( f, xinit, tol, maxiter, result, success )

c        interface
c             function f(x)
c              real*8, intent(in) :: x
c            end function f
c        end interface

c        real*8, intent(in)     :: xinit
c        real*8, intent(in)     :: tol
c        integer, intent(in)  :: maxiter
c    	real*8, intent(out)    :: result
c    	logical, intent(out) :: success

c    	real*8                 :: eps = 1.0e-4
c    	real*8                 :: fx1
c    	real*8                 :: fx2
c    	real*8                 :: fprime
c    	real*8                 :: x
c    	real*8                 :: xnew
c    	integer              :: i

c    	result  = 0.0
c    	success = .false.

c   	    x = xinit
c        do i = 1,max(1,maxiter)
c        	fx1    = f(x)
c        	fx2    = f(x+eps)
c        	write(*,*) i, fx1, fx2, eps
c        	fprime = (fx2 - fx1) / eps

c        	xnew   = x - fx1 / fprime

c        	if ( abs(xnew-x) <= tol ) then
c            	 success = .true.
c            	 result  = xnew
c            	 exit
c        	endif

c        	x = xnew
c        	write(*,*) i, x
c     	enddo

c      end subroutine find_root

       


!**************** MODIFIED GRAVITY CHAMELEON ***************************    
! This is the function where the root is the screening radius in the case
! of the Burkert chameleon field profile
       double precision function dfbur( x ) 
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
         
        phinf=tmass*1e-5 
        tm200=100/Grav*(r200**3)*h0**2*(Omegam*(1.0d0+za)**(3)+Olam)
        fac200=dlog(1+(r200/rs)**2)+2.*dlog(1.0d0+r200/rs)-2.0d0*
     &   datan(r200/rs)
         facda=tm200/fac200/(pig*rs**3)*  !rhos
     &   screen*(rs**2)*8*pig*Grav/(cl**2) !the rest of the factor

        
        fac1=(2*datan(x)-2*dlog(1.0d0+x) + dlog(1.0d0 + x**2))*facda
        dfbur=((1.0d0+x)*(4*phinf - facda*pig+fac1))/(4*x)

        return 
        end
         
!      It finds the root of the dfbur function     
       subroutine froot_Bur(x1,x2)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external dfbur
       include 'paramsoptS.i'
       limit=1000
       i=1
      !x1 = 0.05d0 ; x2 = 10.d0 ; 
       e = 1.0e-5
      
      DO 
        IF (i > limit) THEN
          WRITE(*,*) "Function not converging"
         EXIT
        END IF
        d = (x2 - x1) / (dfbur(x2) - dfbur(x1)) * dfbur(x2)

        IF (ABS(d) < e) THEN
         !WRITE(*,"(A,F18.15)") "Root found at x = ", x2
         EXIT    
        END IF
        x1 = x2
        x2 = x2 - d
        i = i + 1
      END DO
      end
      
      subroutine find_xsc(xsc,alimit)
      
      !alim=0.001 limiting value for screening radius to be considered 0
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external dfbur,froot_Bur, fein, froot_Ein
       include 'paramsoptS.i'
       
       xsc=alimit
       x0=alimit
       dftest=0.0
       
       if (kmp.eq.11) then
       ! see if the function is negative at the beginning
        dftest=dfbur(alimit)
        if (dftest.lt.0) then  
         do i=1,4000
            xx=x0+0.01d0

            dfx0=dfbur(x0)
            dfxx=dfbur(xx)
            dfprod=dfx0*dfxx
            !if you find the screening radius, refine the value
            if (dfprod.lt.0) then 

                call froot_Bur(x0,xx)
                xsc=xx           
                return
            endif
            x0=xx
         enddo
         !if not found, the screening radius is the maximum of the interval
         if (xsc.eq.alimit) xsc=xx
         return
        endif 
 
!      case Einasto                
       else if (kmp.eq.12) then 
         dftest=fein(alimit)
        if (dftest.lt.0) then  !monotonically increasing function
         do i=1,4000
c        !search for the screening radius in the interval 0-40 Mpc
            xx=x0+0.01d0
            dfx0=fein(x0)
            dfxx=fein(xx)
            dfprod=dfx0*dfxx            
!if you find the screening radius, refine the value
            if (dfprod.lt.0) then 

                call froot_Ein(x0,xx)
                xsc=xx           
                return
            endif
            x0=xx
         enddo
         !if not found, the screening radius is the maximum of the interval
         if (xsc.eq.alimit) xsc=xx
         return
        endif             

       endif
       
       return 
       end
 
!     Same for the BCG       
      subroutine find_xscB(xscB,alimbc)
      
      !alim=0.001 limiting value for screening radius to be considered 0
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external fjaf_r, froot_Jaffe
       include 'paramsoptS.i'

       x0=alimbc
       
       xscB=alimbc     
    
c       write(*,*) alimbc, fjaf_r(alimbc)
       if (kbcg.eq.1) then

    
         if (fjaf_r(alimbc).lt.0) then    
          do i=1,4000
            xp=x0+0.0005d0 
            if (fjaf_r(x0)*fjaf_r(xp).lt.0) then 
                call froot_Jaffe(x0,xp)
                xscB=xp           
                return
            endif
            x0=xp
          enddo
          
          if (xscB.eq.alimbc) xscB=xp
         endif       
       
       else 
         xscB=0.0d0 
       endif
       return
       end
       
       double precision function dphib(x)
!>   @author L. Pizzuti
!>   @details  compute the term in the effective mass  Q*r**2/G*dphi/dr 
!!   due to chameleon field for a Burkert mass model.  
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen`. 
!!  The background field value `tmass` is given by tmass (in unit of 1e-5)  
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za`                                 
!>  - mass model `kmp` 
!>  - number of steps in the second MG parameter `nhone`. If equal to -1
!!  it forces f(R) chameleon case.
!>  - the value of the screening radius `xsccH` previously computed by find_xsc 
!>  - the value `alimi` of the screening radius to be considered 0
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       h0=70
!>       Omegam=0.3
!>       Omegal=0.7
!>       kmp=9
!>       r200=1.2
!>       rs=0.5
!>       tmass=100.23
!>       nhone=10
!>       screen=0.3
!>       alimi=0.001
!>       xsccH=alimi
!>       alimit=alimi
!>       call find_xsc(xsccH,alimit)
!>       write(*,*) 'test dphib= ', dphib(1.0d0)
!>  \endcode       
       
       
       
       !Gives the value of Q*r**2/G*dphi/dr for burkert profile in unit
       !of M200
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
         
        phinf=tmass*1e-5 
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3
        !SONO QUI!!!!!
        !Force f(R) gravity (Q = 1/sqrt(6)) 
        if(kmp.eq.9.and.nhone.lt.0) then
         bcoup=1./dsqrt(6.0d0)         
       else
         bcoup=screen !the modified parameter of the screening becomes
                      !the coupling constant in the genneric chameleon run
       endif
        
        
        c200=r200/rs
        fac200=dlog(1+c200*c200)+2.*dlog(1.0+c200)-2.*datan(c200)
        rhos=tm200/fac200/(pig*rs**3)
        
        Facda=screen*rs**2*8*pig*Grav*rhos/cl**2
        
        Cb=Facda*(2.0d0*datan(xscCH)-2*dlog(1+xscCH)-dlog(1+xscCH**2))/4
        
        !(-4*phinf + Facda*pig - 2*Facda*dlog(1 + xscCH**2))/4. 
        
        !exponential cutoff
        dminf=1
        if (ncutof.gt.0) then
            
            dminf=dexp(-x*rs*dsqrt(1e-7*screen/phinf))
        endif
        
        if (xscCH.gt.alimi) then
           if (x.lt.xscCH) then
            dout=0.0d0
           else
           !this is x**2*dphi/dx
            dout=Cb+Facda/(4.)*(-2*datan(x)+2*dlog(1+x)+dlog(1+x**2))
           endif
        else
         dout=Facda/(4.)*(-2*datan(x)+2*dlog(1+x)+dlog(1+x**2))
        endif
        
        dphib=dout*rs*screen*cl**2/Grav/tm200*dminf !note the last factor!
                                            !masses are given in units
                                            ! of m200
        return
        end
        
        
        
       double precision function fein( x )
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
         
        phinf=tmass*1e-5 
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3
        
        c200=r200/rs
        !eim=5.
        call GAMMA(2*eim,tgam2)
        call GAMMA(3*eim,tgam3)
        agp=3.*eim
        bgp=2*eim
        xgp=2.*eim*(r200/rs)**(1./eim)
        
        xsc=2*eim*x**(1./eim)
        
        call dincog(agp,xgp,gin,gim,gip200) !(1-gammaincc)
        call dincog(bgp,xsc,gin,gim,gipx) !1-gammaincc -> 
        tgamc=1-gipx
        fac200=gip200*tgam3  
        rho1=tm200*(8**eim)*eim**(3*eim-1)
        rhos=rho1/(4*pig*rs**3*fac200*dexp(2*eim))
        Facda=screen*rs**2*rhos*8*pig*Grav/cl**2
        f1=4**eim*phinf*eim**(2*eim-1.0d0)/(Facda*dexp(2*eim)*tgam2)
        fein=f1-tgamc
        
        return 
        end        
        
       subroutine froot_Ein(x1,x2)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external fein
       include 'paramsoptS.i'
       limit=1000
       i=1
      !x1 = 0.05d0 ; x2 = 10.d0 ; 
       e = 1.0e-5
      
      DO 
        IF (i > limit) THEN
          WRITE(*,*) "Function not converging"
         EXIT
        END IF
        d = (x2 - x1) / (fein(x2) - fein(x1)) * fein(x2)

        IF (ABS(d) < e) THEN
         !WRITE(*,"(A,F18.15)") "Root found at x = ", x2
         EXIT    
        END IF
        x1 = x2
        x2 = x2 - d
        i = i + 1
      END DO
      end


!=======================================================================

       double precision function dphiso(x,rhos,rss,Q,phinf)
!>   @author L. Pizzuti
!>   @details  compute the term in the effective mass due to chameleon field
!!   for a Isothermal mass model with beta = 1 (exp = 3/2). 
!!  It requires the parameter `alimi` in the COMMON block paramsoptS.i:                     
!>  - the value `alimi` of the screening radius to be considered 0
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       alimi=0.001
!>       write(*,*) 'test dphiso= ', dphib(1.0d0)
!>  \endcode       
       
       
       
       !Gives the value of Q*r**2/G*dphi/dr for burkert profile in unit
       !of M200
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
         
        Facda=Q*rss**2*8*pig*Grav*rhos/cl**2
    
        
        !exponential cutoff
        dminf=1
        if (ncutof.gt.0) then
            
            dminf=dexp(-x*rss*dsqrt(1e-7*Q/phinf))
        endif
        
        if ((Facda/phinf).gt.1.0d0) then
        
    
         xsc = dsqrt((Facda/phinf)**2 - 1.0d0) 
         Ciso = phinf*xsc - Facda*dlog(dsqrt(1+xsc**2)+xsc)
       else
        xsc = alimi !the limiting value is defined in gomamposstoptS.f
        Ciso = 0.0d0
       endif
       if (x.ge.xsc) then
        dnum1 = Ciso-(Facda*x)/dsqrt(1.0d0+x**2)
        dnum2 = -Facda*dlog(1/(x+dsqrt(1 + x**2))) 
        dout =(dnum1+dnum2)*rss*Q*cl**2/Grav*dminf

       else
        dout = 0.0d0
       endif
        
    
        dphiso= dout 
        
        return
        end


c######=================================================================
      function tM_iso_CS(x) 
!>   @author L. Pizzuti
!>   @details  compute the effective mass due to chameleon field for a
!!   Isothermal mass density profile, which exponent is 3/2
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen`. 
!!  The background field value `tmass` is given by tmass (in unit of 1e-5)  
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za` 
!>  - cutoff of the fifth force (Y/N) `ncutof=0/1`.                                
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>        
!>       h0=70
!>       Omegam=0.3
!>       Omegal=0.7
!>       kmp=9
!>       r200=1.2
!>       rs=0.5
!>       tmass=100.23
!>       screen=4.3
!>       ncutof = 1
!>       eim = 1.2
!>       write(*,*) tM_iso_CS(2.12d0)
!>  \endcode       
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       parameter (pig=3.1415926535d0, clight=2.99792458d5)
       
       parameter (grav=4.302e-9, alimit=0.001) !limiting value for the screening radius
       include 'paramsoptS.i'
       
       cx = x/rs
       hz=h0*dsqrt(Omegam*(1+za)**3+Olam)
       c200=r200/rs
       phinf=tmass*1e-5
       bcoup = screen
       
       dm200=100*hz**2/grav*r200**3
       
       
       sqc = dsqrt(1+c200**2) 
       rhoss= - dm200*sqc/((c200+sqc*dlog(dsqrt(c200**2 +1) - c200) ) )
       
       rhos = rhoss/( 4*pig*rs**3 )
       dnim = -cx/(dsqrt(1+cx**2) ) + dlog( dsqrt( 1 + cx**2 ) + cx )
       dmass = dnim*rhoss 

       
       if (bcoup.ne.0.and.tmass.ne.0) then 
         dphid = dphiso(cx,rhos,rs,bcoup,phinf)
         
      else
         dphid = 0.0d0
      endif
        
        tM_iso_CS = (dphid + dmass)/dm200 !given in units of m200
       return
        
      end



!=======================================================================
        
        
        double precision function dphiE(x)
        
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'   
        
        phinf=tmass*1e-5 
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3
        
        c200=r200/rs
!        eim=5.
        call GAMMA(2*eim,tgam2)
        call GAMMA(3*eim,tgam3)
        agp=3.*eim
        bgp=2*eim
        xgp=2.*eim*(r200/rs)**(1./eim)
        
        xss=2*eim*x**(1./eim)
        
        xsc=2*eim*xscCH**(1./eim)
        
        
        !exponential cutoff
        dminf=1
        if (ncutof.gt.0) then
            
            dminf=dexp(-x*dsqrt(1e-7*screen/phinf))
        endif
        
!       Incomplete gamma function for  fac200
        call dincog(agp,xgp,gin,gim,gip200) !(1-gammaincc(c200))
        
        !Incomplete gamma function for  the x value
        call dincog(agp,xss,gin,gim,gipx) !(1-gammaincc(x))
        
        !Incomplete gamma function for the screening radius
        call dincog(agp,xsc,gin,gim,gipsc) !1-gammaincc(xscCH)
        
        
        tgamc=1.0-gipx
        fac200=gip200*tgam3  
        rho1=tm200*(8**eim)*eim**(3*eim-1)
        rhos=rho1/(4*pig*rs**3*fac200*dexp(2*eim))
        Facda=screen*rs**2*rhos*8*pig*Grav/cl**2
        
        
        Ce=8**(-eim)*dexp(2*eim)*Facda*eim**(1.-3*eim)*(1-gipsc)*tgam3

        if (xscCH.gt.alimi) then
          if (x.le.xscCH) then
            dout=0.0d0
          else
            dout1=(8)**(-eim)*dexp(2*eim)*Facda*eim**(1.-3*eim)
            dout=Ce-dout1*tgamc*tgam3
          endif
        else
         Cce=8**(-eim)*dexp(2*eim)*eim**(1 -3*eim)*tgam3*Facda
         dout1=(8)**(-eim)*dexp(2*eim)*Facda*eim**(1.-3*eim)
         dout=Cce-dout1*tgamc*tgam3
        endif 
        dphiE= dout*rs*screen*cl**2/Grav/tm200*dminf
        return 
        
        end


!xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)


       double precision function fjaf_r( x )
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
        
        phinf=tmass*1e-5 
        rhos=xmstar*xlumbcg/(4*pig*rjaf**3)
        
        Facda=screen*rjaf**2*8*pig*Grav*rhos/cl**2
        
        fjaf_r=Facda*(x+1)*(dlog(x/(x+1)))+Facda+phinf*(x+1)
        return
       end 

       subroutine froot_Jaffe(x1,x2)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external fjaf_r
       include 'paramsoptS.i'
       limit=1000
       i=1
       e = 1.0e-5
      
       DO 
        IF (i > limit) THEN
          WRITE(*,*) "Function not converging"
         EXIT
        END IF
        d = (x2 - x1) / (fjaf_r(x2) - fjaf_r(x1)) * fjaf_r(x2)

        IF (ABS(d) < e) THEN
         !WRITE(*,"(A,F18.15)") "Root found at x = ", x2
         EXIT    
        END IF
        x1 = x2
        x2 = x2 - d
        i = i + 1
       END DO
      end


      double precision function dphiJ(x)
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)
        
        include 'paramsoptS.i'
         
        phinf=tmass*1e-5 
        
        rhos=xmstar*xlumbcg/(4*pig*rjaf**3)
        
        Facda=screen*rjaf**2*8*pig*Grav*rhos/cl**2
        
        CJ= phinf*xsBCG+Facda*(xsBCG*dlog(xsBCG/(1+xsBCG))+1)
        
        !exponential cutoff
        dminf=1
        if (ncutof.gt.0) then
            
            dminf=dexp(-x*dsqrt(1e-7*screen/phinf))
        endif
        
        if (xsBCG.gt.0.0001) then
           if (x.lt.xsBCG) then
            dout=0.0d0
           else
            dout=CJ-Facda/(x+1)
           endif
        else
           dout=Facda*(1-1/(x+1))
        endif
       
        dphiJ=dout*rjaf*screen*cl**2/Grav/(xmstar*xlumbcg)*dminf !note the last factor!
                                            !masses are given in units
                                            ! of m200
        return
        end
        


!General screening implementation


      double precision function fscrtot(x)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       !complex hfy
       parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
    
       include 'paramsoptS.i'
       include 'barions.i'

       external hypgeo

       Facda = screen*(8*pig*Grav)/cl**2
       phinf=tmass*1e-5
       
       !Define the central jaffe density
       rhob = (xmstar*xlumbcg)/(4*pig*rjaf**3)

       !Use A3 as gamma
       dgami = dA3
       
       !Define the central gNFW profile
       !sono qui
        
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3

        call hygfx ( 3-dgami, 3-dgami, 4-dgami, -r200/rs, hfy)
        
    
        rhop1=(r200/rs)**(dgami-3.0d0)
        
        rhop2 = tm200/(4*pig*rs**3*hfy)*(3-dgami)
        rhos = rhop1*rhop2
       

       
       
        djaf = (rjaf**3*rhob)/(x+rjaf)+rjaf**2*rhob*dlog(x/(x + rjaf)) 
       
        dgas = -(rsga1**3 *rhog1)/dsqrt(x**2+rsga1**2)
        dgal =  -(rst1**3 * rhost1)/(x + rst1)
        ddm =rs**2*(1.0d0-(x/(x +rs))**(2-dgami))*rhos/(-2.0d0+dgami)
        !write(*,*) dgas, rsga1, rhog1
        fscrtot =  phinf + Facda*(dgas+dgal+ddm+djaf) 
        
      return
      end
      
      
      subroutine froot_GEN(x1,x2)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
      
       external fscrtot
       include 'paramsoptS.i'
       limit=10000
       i=1
      !x1 = 0.05d0 ; x2 = 10.d0 ; 
       e = 1.0e-6
      
       DO 
        IF (i > limit) THEN
          WRITE(*,*) "Function not converging"

         EXIT
        END IF
        d = (x2 - x1) / (fscrtot(x2) - fscrtot(x1)) * fscrtot(x2)

        IF (ABS(d) < e) THEN
         !WRITE(*,"(A,F18.15)") "Root found at x = ", x2
         EXIT    
        END IF
        x1 = x2
        x2 = x2 - d
        i = i + 1
      END DO
      end
      
      
! This computes the general multicomponent screening. Tested      
      double precision function ScreenGen(xlimit)
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
    
       include 'paramsoptS.i'
       include 'barions.i'

       external fscrtot
            
        xsc = xlimit
        x0 = xlimit   
        if(tmass.lt.0.005.or.screen.lt.0.005) then
         xsc = 0.0d0
        else
         if(fscrtot(x0).lt.0.0d0) then
          do i = 1,80000
            xx=x0+0.001
            if (fscrtot(x0)*fscrtot(xx).lt.0.0d0) then
                xtemp = xx !save option when the root in NaN

                !call froot_GEN(x0,xx)
                !if (isnan(xx)) xsc = xtemp
                xsc = xx
                ScreenGen = xsc        
                return
            endif
            x0=xx
          enddo
          if (xsc.eq.xlimit) xsc=xx
       
        else 
         xsc=0.0d0 
        endif
       endif
       ScreenGen = xsc

      return 
      end
      

!  fifth force in General multicomponent Chameleon
!  Validated against the python script
      double precision function dphiGen(x)
!>   @author L. Pizzuti
!>   @details  compute the term in the effective mass  Q*r**2/G*dphi/dr 
!!   due to chameleon field for a four multicomponent mass model.  
!!  It requires the parameters in the COMMON block paramsoptS.i, in particular:
!>  - (effective) mass profile parameters: `r200, rs, tmass, screen, dA3`.
!!  The background field value `tmass` is given in unit of 1e-5 c**2.  
!>  - the BCG parameters: `xlumbcg, xmstar, rjaf`.
!>  - The baryonic parameters: rho1g, rsga1, rhost1, rst1 
!>  - cosmological parameters: hubble constant `h0` (in km/s/Mpc),
!!  density parameters `Omegam, Omegal`, redshift `za`                                 
!>  - mass model `kmp` 
!>  - the value of the screening radius `xsccH` previously computed by ScreenGen 
!>   @param[in] x REAL*8, value at which the function is computed
!>
!>
!>   Example of usage:
!>   \code{.f}
!>       rs = 0.5
!>       r200 = 1.9
!>       tmass = 10.0d0
!>       screen = 0.5
!>       dA3 = 1.4
!>       xlumbcg = 4.9e11
!>       xmstar = 4.6
!>       rjaf = 0.039
!>       
!>       h0=70
!>       Omegam=0.3
!>       Olam=0.7
!>       za=0.44
!>       rhog1 = 1e14
!>       rsg1  = 0.14
!>       rhost1 = 2e13
!>       rst1 = 0.36
!>       xlim = 1.0e-6
!>       xsc = ScreenGen(xlim)
!>       write(*,*) "Value of the screeining:", xsc
!>       xscCH = xsc
!>       write(*,*) "Value of the effective mass:", dphiGen(2.0d0)
!>  \endcode       

       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       !complex hfy
       parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
    
       include 'paramsoptS.i'  
       
       xlimit = 1e-6 !limiting value for the screening
       
       Facda = screen*(8*pig*Grav)/cl**2    
       
       !Define the central jaffe density
       rhob = (xmstar*xlumbcg)/(4*pig*rjaf**3)

       !Use A3 as gamma
        dgami = dA3
       
       !Define the central gNFW density
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3

        call hygfx ( 3-dgami, 3-dgami, 4-dgami, -r200/rs, hfy)

        rhos=(r200/rs)**(dgami-3.0d0)*tm200/(4*pig*rs**3*hfy)*(3-dgami)
       
       
       xsc = xscCH !Define the screening using the common parameter
       
       call hygfx ( 3-dgami, 3-dgami, 4-dgami, -x/rs, hfx)
       call hygfx ( 3-dgami, 3-dgami, 4-dgami, -xsc/rs, h2sc)
        
        
       !Integration constant   
       duno = (rjaf**4*rhob)/(xsc+rjaf) 
       ddue = (xsc**3*(xsc/rs)**(-dgami)*rhos*h2sc)/(-3.0d0 + dgami)
       
       dtrep1 = xsc/dsqrt(xsc**2+ rsga1**2)
       dtrep2 = dlog(-xsc+ dsqrt(xsc**2 + rsga1**2))
       dtre = rsga1**3*rhog1*(dtrep1+dtrep2)  
       
       dquat=-rst1**3*rhost1*(rst1/(xsc + rst1)+dlog(xsc + rst1))
       
       Ctot = Facda*(duno + ddue + dtre + dquat)
       
       if (xsc.gt.xlimit) then
        if (x.lt.xsc) then
            dout=0.0d0
        else
            dun1=x**(3.0-dgami)*(1.0d0+x/rs)**(-dgami)*(x+rs)**(dgami)
            duno1 = dun1*rhos*hfx/(3.0d0-dgami) 
            ddue2 = -((rjaf**4 *rhob)/(x + rjaf))
            
            dtrp13 = -(x/dsqrt(x**2 + rsga1**2)) 
            dtrp23 = - dlog(-x + dsqrt(x**2 + rsga1**2))
            dtre3 = rsga1**3*rhog1* ( dtrp13 + dtrp23)
             
            dquat4 = rst1**3*rhost1*(rst1/(x+rst1)+dlog(x + rst1) )
            
            dout= Ctot + Facda*(duno1 + ddue2 + dtre3 + dquat4)
       
        endif
       else
             Ctot=Facda*(2*rjaf**3*rhob+rsga1**3*rhog1*dlog(rsga1**2) - 
     -         2*rst1**3*rhost1*(1.0d0+dlog(rst1)) )/2
     
            dun1=x**(3.0-dgami)*(1.0d0+x/rs)**(-dgami)*(x+rs)**(dgami)
            duno1 = dun1*rhos*hfx/(3.0d0-dgami) 
            ddue2 = -((rjaf**4 *rhob)/(x + rjaf))
            
            dtrp13 = -(x/dsqrt(x**2 + rsga1**2)) 
            dtrp23 = - dlog(-x + dsqrt(x**2 + rsga1**2))
            dtre3 = rsga1**3*rhog1* ( dtrp13 + dtrp23)
             
            dquat4 = rst1**3*rhost1*(rst1/(x+rst1)+dlog(x + rst1) )
            
            dout= Ctot + Facda*(duno1 + ddue2 + dtre3 + dquat4)            
       
       endif
       !validated
       dphiGen = dout*screen*cl**2/Grav 
        
      return 
      end

!     validated!!!  
      double precision function tmass_Gen(x)     
       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
       
       !complex hfy
       parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
    
       include 'paramsoptS.i'  
       
       xlimit = 1e-6 !limiting value for the screening
       
       Facda = screen*(8*pig*Grav)/cl**2    
       
       !Define the central jaffe density
       rhob = (xmstar*xlumbcg)/(4*pig*rjaf**3)

       !Use A3 as gamma
        dgami = dA3
       
       !Define the central gNFW density
        Hz2=h0**2*(Omegam*(1+za)**3+Olam)
        tm200=100*Hz2/Grav*r200**3

        call hygfx ( 3-dgami, 3-dgami, 4-dgami, -r200/rs, hfy)

        rhos=(r200/rs)**(dgami-3.0d0)*tm200/(4*pig*rs**3*hfy)*(3-dgami)
       

       xd = x/rs
       xg =  x/rsga1
       xstar = x/rst1
       xj = x/rjaf
        
       xsc = xscCH !Define the screening using the common parameter
       
       call hygfx ( 3-dgami, 3-dgami, 4-dgami, -x/rs, h2x)
       
       tmg1 = rhob/(xj**2 + xj**3)-rhog1/(xg**2 * dsqrt(1.0d0+xg**2))
       tmg2 = -rhost1/(xstar**2*(1.0d0+xstar))+(rhog1*dasinh(xg))/xg**3
       tmg3 = ((xd)**(3.0d0-dgami)*rhos*h2x)/(xd**3*(3-dgami))
       tmg4 = (rhost1*dlog( 1.0d0 + xstar) )/xstar**3
       
       tmgen = 4*pig*x**3*(tmg1 +tmg2 + tmg3 + tmg4)
       if (screen.ne.0.and.tmass.ne.0) then  
         dphir= dphiGen(x)
       else
         dphir=0.0
       endif
       tmass_Gen = tmgen + dphir
      return         
      end

        
!=======================================================================

!> \ingroup MAM
      double precision function rgmass(x)  
!>  @author V. Broccolato, L. Pizzuti, F. Fantoccoli      
!>  @details This function computes the analytic expression 
!!  for the radial gravitational potential in Refracted Gravity (RG)
!!  assuming a NFW profile for the gas and galaxy distribution with fixed 
!!  parameters. The other mass component is given by the BCG with the free
!!  parameters `rjaf`, `xmstar`.
!>  It requires the inclusion of 'paramsoptS.i' and a definition of `rjaf`,
!!  `xmstar`, `xlumbcg`, cosmology (`h0`, `Omegam`, `Olam`,`za`), as well as the
!!  modified gravity parameters, expressed by `tmass` (epsilon 0), `screen`
!!  (log10 of critical density) and `dA3` (steepness), as well as the  
!!  baryonic profile parameters `rho0gas`, `rsgas`. 

!>
!>   @param[in] R radius in Mpc
!>
!>   Example of usage:
!>
!>   \code{.f}
!>   rjaf=0.03d0
!>   xlumbcg=492000000000.0d0
!>   xmstar=4.3
!>
!>   !Cosmology
!>   h0=70
!>   Omegam=0.3
!>   Olam=0.7
!>   za=0.44
!>   !RG
!>   tmass=0.089d0
!>   dA3=0.47d0
!>   screen=14.6


!>   write(*,*) rgmass(1.2d0)
!>   \endcode
      
      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
      parameter (dln10 = 2.3025850929940459d0)
      include 'paramsoptS.i'
      include 'barions.i'
       
      external fbars 
      !r200gas=17.066 !Mpc
      !baryonic mass fitted with a beta profile with exponent 2/3
      
      !included in the paramsoptS.i file
      !rsgas= 0.12   !0.813 !Mpc
      !rho0gas=89.08e13  
      
       if (srhoga.eq.0) rhog1 = rho0gas
       if (srsga.eq.0)  rsga1 = rsgas
       
       
       if (srhost.eq.0) rhost1 = rhostar
       if (srstar.eq.0) rst1 = rstar
      
       !Test for using r200 as rhog1
       !if (srhoga.eq.0) rhog1 = r200 * 1e13
      
      
      !TO BE CHANGED IN THE NEW CONFIGURATION WITH THE TOTAL BARYIONIC MASS!
      !xmbar= 4*pig*rhog1*rsga1**2*(x-rsga1*datan(x/rsga1))
      xmbar = fbars(x)
      
      !Critical density
      !arhoc=10**(screen)
    
      ! NFW gas
      !c2gas=r200gas/rsgas
      !fac200=dlog(1.d0+c2gas)-c2gas/(1.d0+c2gas)
      !h2z=h0**2*((1+za)**3*Omegam+Olam)
      !dm200=h2z*100/Grav*r200gas**3
      !massa totale
      !rho0gas=dm200/(4*pig*(rsgas**3)*fac200)
      !xmbar=dm200/fac200*(dlog(1.d0+x/rsgas)-(x/rsgas)/(1.d0+x/rsgas))
      
      !total mass: bcg + gas + galaxies
      
      dmtot=xmstar*xlumbcg*x/rjaf/(1.+x/rjaf)+ xmbar 
      
      !now compute the effective permittivity
      rbc= x/rjaf 
      ! density of the BCG (jaffe profile)
      rhobcg=xmstar*xlumbcg/(4.d0*pig*rjaf)/((rbc)**2*(1.d0+rbc)**2)
      
      !OLD ========================
      !density baryons (beta model)
      !rhobar  = rhog1/(1.d0+(x/rsga1)**2)     !rho0gas/((x/rsgas)*(1.d0+x/rsgas)**2)
      ! ===========================
      
      !NEW! Density gas and galaxies:
      if (kalpha.eq.1) then
       rhogas = rhog1/(1.0d0 +(x/rsga1)**2)**(3.0d0/2)
      else
       rhogas = rhog1/(1.0d0 +(x/rsga1)**2)
      endif
      rhogal = rhost1/( x/rst1 * (1.0d0 + x/rst1)**2)
      
      rhotot=rhobcg + rhogas + rhogal  !rhobar
      
      dlrho = dlog(rhotot) - screen*dln10
      
      
      
      epsi=tmass+(1-tmass)*(1.d0/2)*(dtanh(dA3*dlrho)+1)
      
      
      
      !Effective mass 
      rgmass= (dmtot) /(epsi) ! (Grav*dmtot) /(epsi*x**2) 
      
      return 

      end
      
c      double precision function fgrt(x)
c      implicit real*8 (a-h,o-z)
c      implicit integer*4 (i-n)
c      parameter (pig=3.1415926535897932d0, Grav=4.302e-9, cl=3e5)       
c      external fbars
c      include 'paramsoptS.i'
c      include 'barions.i'

c      !parametri su cui MAMPOSSt fitta:
c      !r200, rs, rc, cbe, cbe0, tmass, screen, xmstar, rjaf

c      !r200gas = valore trovato dal fit di massa barionica
c      r200gas=17.066 !Mpc
c      rsgas= 0.12   !0.813 !Mpc
c      !rsgas = valore trovato dal fit di massa barionica
      
c      arhoc=10**(screen)
    
c      ! NFW gas
c      c2gas=r200gas/rsgas
c      fac200=dlog(1.d0+c2gas)-c2gas/(1.d0+c2gas)
c      h2z=h0**2*((1+za)**3*Omegam+Olam)
c      dm200=h2z*100/Grav*r200gas**3
c      !massa totale
c      rho0gas=89.08e13   !dm200/(4*pig*(rsgas**3)*fac200)
c      !xmbar=dm200/fac200*(dlog(1.d0+x/rsgas)-(x/rsgas)/(1.d0+x/rsgas))
      
c      xmbar= 4*pig*rho0gas*rsgas**2*(x-rsgas*datan(x/rsgas))
      
      
c      dmtot=xmstar*xlumbcg*x/rjaf/(1.+x/rjaf)+ xmbar!fbars(x)
      
c      !densità
c      rbc= x/rjaf 
c      ! density of the BCG
c      rhobcg=xmstar*xlumbcg/(4.d0*pig*rjaf)/((rbc)**2*(1.d0+rbc)**2)
c      !density baryons
c      rhobar  = rho0gas/(1.d0+(x/rsgas)**2)       !rho0gas/((x/rsgas)*(1.d0+x/rsgas)**2)
c      rhotot=rhobcg+rhobar
c      epsi= 1.0d0 !tmass+(1-tmass)*(1.d0/2)*(dtanh(rs*dlog(rhotot/screen))+1)

c      fgrt = (dmtot) /(epsi) ! (Grav*dmtot) /(epsi*x**2)   calcolo del gradiente del potenziale

c      return 

c      end



c#Field derivative:
cdef dphiJ(x,rs,rhos,Q,phinf):
c    Facda=Q*rhos*rs**2/Mpl
c    def fjaffe(y):
c        return Facda*(y+1)*(np.log(y/(y+1)))+Facda+phinf*(y+1)
c    xsc=limit
c    x0=limit
c    if fjaffe(limit)<0:    
c        for i in range(0,4000):
c            xx=x0+0.01*i
c            if fjaffe(x0)*fjaffe(xx) < 0:
c                sol=opt.root(fjaffe,xx)
c                if sol.success:
c                    xsc=sol.x[0]
c                else:
c                    xsc=xx                
c                break
c            x0=xx
c        if (xsc==limit):
c            xsc=xx

c    CJ= phinf*xsc+Facda*(xsc*np.log(xsc/(1+xsc))+1)   
c    if (xsc>limit):
c        if (x<xsc):
c            out=0
c        else:
c            out=CJ/(rs*x**2)-Facda/(rs*x**2*(x+1))
c    else:
c        out=Facda/(rs*x**2)-Facda/(rs*x**2*(x+1))
c    return out    
cdphiJ=np.vectorize(dphiJ,otypes=[float],excluded=[1,2,3,4])

c#Field profile:
cdef phiJ(x,rs,rhos,Q,phinf):
c    Facda=Q*rhos*rs**2/Mpl
c    def fjaffe(y):
c        return Facda*(y+1)*(np.log(y/(y+1)))+Facda+phinf*(y+1)
c    xsc=limit
c    x0=limit
c    if fjaffe(limit)<0:    
c        for i in range(0,4000):
c            xx=x0+0.01*i
c            if fjaffe(x0)*fjaffe(xx) < 0:
c                sol=opt.root(fjaffe,xx)
c                if sol.success:
c                    xsc=sol.x[0]
c                else:
c                    xsc=xx                
c                break
c            x0=xx
c        if (xsc==limit):
c            xsc=xx

c    CJ= phinf*xsc+Facda*(xsc*np.log(xsc/(1+xsc))+1)   
c    if (xsc>limit):
c        if (x<xsc):
c            out=0
c        else:
c            out=-CJ/(x)+Facda*(1/x+np.log(x/(1+x)))+phinf
c    else:
c        out=-Facda/(x)+Facda*(1/x+np.log(x/(1+x)))+phinf
c    return out




        

