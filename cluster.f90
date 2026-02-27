    !-----------------------------------------------------------!
    ! Module for cluster counts likelihood
    !-----------------------------------------------------------!


    module clcosmology
    ! module computing cosmological functions
    USE precision
    use Interpolation
    implicit none
    ! note that we switch to the generalized Linder parametrization now
    public
    TYPE cospar
        REAL(dl) :: rs,r200,rjaf,cbe,gam,xmstarl,rabcg
    END TYPE cospar
    Type (cospar), SAVE :: cosmopar

    end module clcosmology

!####Mio:modulo  interno a questo file. Qui si interfaccia con il resto di COSMOMC
    module clcounts
    use clcosmology !####Mio: set my cosmo par
    use settings  !####Mio:modulo definito in settings.f90, usato ovunque
    use CosmologyTypes !####Mio:modulo definito in CosmologyTypes.f90, usato ovunque
    use CosmoTheory !####Mio:modulo definito in  CosmoTheory.f90, usato ovunque CosmoTheory.f90
    use Calculator_Cosmology !####Mio:modulo definito in Calculator_Cosmology.f90, usato ovunque
    use Likelihood_Cosmology !####Mio:modulo definito in Likelihood_Cosmology.f90, usato ovunque
    use IniObjects !####Mio:modulo definito in IniObjects.f90, usato ovunque

    implicit none
    private
    logical :: do_cl_init = .true.  !####Mio: il primo loop inizializzo e poi non piu'
    logical :: use_cl = .false.     !####Mio: leggo dall'input file se far girare Cosmomc con i cluster o meno
    real(dl) :: extdata1,extdata2 ! switches for priors =0 no prior, =1 prior   !####Mio:nomi dei setting ti prior 
    real :: sz_kmax = 4.0  !controlla
    character(LEN=*), parameter :: CL_version =  'Apr_2017'

!####Mio: add in TCosmoCalcLikelihood the ClLikelihood
    type, extends(TCosmoCalcLikelihood) :: ClLikelihood
    contains
!####Mio: qui credo aggiunga il valore trovato della likelihood al loop di Cosmomc
    procedure :: LogLike => Cl_Cash
    end type ClLikelihood

    PUBLIC :: ClLikelihood_Add, Cl_Cash

    contains

    subroutine ClLikelihood_Add(LikeList, Ini)
    ! interface of the module with cosmomcplanck code
    ! user-defined settings are interfaced here
    implicit none
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    integer count
    Type(ClLikelihood), pointer :: this

!####Mio: usare i cluster. Si chima nell'input file
    if (Ini%Read_Logical('use_Cl',.false.)) then
        allocate(this)
        this%needs_background_functions = .false.
        this%needs_powerspectra = .false.
        this%kmax = sz_kmax   !controlla: e' il valore settato poche righe prima, non e' detto vada bene nel mio caso
        this%needs_sigmaR = .false.
        this%version = CL_version
        call this%loadParamNames(trim(DataDir)//'CL_mamp.paramnames')
        call LikeList%Add(this)
        this%LikelihoodType = 'CL'
        this%name='CL'

!####Mio: setto a zero il valore iniziale delle flagh delle prior per poi leggere se le ho chieste nell'input file nelle righe successive. La likelihood verra' moltiplicata per questo valore e poi per la prior, quindi se e' zero e' come non mettere prior su quel parametro
        extdata1=0.  !swithes for the priors (get multiplied to priors in Cl_Cash)
        extdata2=0.    !either 0 or 1

        if (Ini%Read_Logical('prior_extdata1',.false.)) then
            extdata1=1.
            print*,'External prior on clsuter mass'
        endif

       if (Ini%Read_Logical('prior_extdata2',.false.)) then
            extdata1=1.
            print*,'External prior on cluster mass'
        endif

        CALL Cl_init

    endif
    end subroutine ClLikelihood_Add

    subroutine Cl_init
    !implicit none

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535d0,iga=25000,grav=4.302e-9,clight=2.99792458d5)
      dimension pars(17),parsvdp(2)
!=== for MAMP      
      include 'datarv.i'
      include 'datagasgal.i'
      include 'paramcomm.i'
      include 'sr.i'
      include 'vkr.i'
!=== for VDP      
      include 'datavdp.i'
      include 'paramvdp.i'

!     Values of the parameters of Hansen & Moores (2006)
!     beta = a + b dln(rho)/dln(r) relation
!
      ahm=-0.15
      bhm=-0.192

!=== parameters for Mamposst       
      open(10,file='data_mamp/vel_gal.dat',status='old')
      open(29,file='data_mamp/mamp_par.dat',status='old')

      do i=1,17
         read(29,*) pars(i)
      enddo
      close(29)

      h0=pars(1)          ! Hubble constant at z=0
      omega0=pars(2)       ! omega matter
      omegal=pars(3)       ! omega lambda

      za=pars(4)     ! average redshift of the cluster (needed to evaluate Hz)
      rlowin=pars(5)        ! Inner radius for sample selection (Mpc) 
      rupin=pars(6)         ! Outer radius for sample selection (Mpc) 
      knfit=nint(pars(7))   ! N(R) model, projected NFW / projected Hernquist / beta-model (1/2/3)
                            !             
!      rcg=pars(8)          ! N(R) scale radius initial guess (Mpc)
      al=pars(8)            ! N(R) negative exponent (only if knfit=3)
      kmp=nint(pars(9))     ! rho(r) model: gNFW/Hernquist/PIEMD/Burkert/SoftIS/Einasto_m=5 (1/2/3/4/5/6)
      ! r200g=pars(11)        ! r200 initial guess (Mpc)
      ! r200=r200g             ! set r200 initial value to guess value
      ! rsg=pars(12)           ! rho(r) scale radius initial guess (Mpc)
      rcut=pars(10)  ! PIEMD model rcut in Mpc (only used if kmp=3)

      ! gamg=pars(14)             ! inner slope of gNFW M(r) >0
      ! gam=gamg

      kani=nint(pars(11))       ! Anisotropy model, beta'=constant, MamLok, OsiMer, simplified Wojtak, simplified Tiret, modified Tiret (0,1,2,3,4,5)
! use kani=-1 to force Hansen&Moore dependence of beta(r) on rho(r)

!      cbeg=pars(16)         ! Anisotropy initial guess, beta', a_ML, a_OM, a_W, beta'_inf
      rc=pars(12)         ! scale radius of BCG Jaffe lum dens. prof (Mpc)
 !     rjaf=pars(12)         ! scale radius of BCG Jaffe lum dens. prof (Mpc)
      xlumbcg=pars(13)      ! BCG total luminosity
 !     xmstarlg=pars(19)     ! free parameter: mass-to-light ratio of the
                            ! stellar population of te BCG (1st guess)
      !      kintd=nint(pars(14))   ! Use universal surface density of interlopers? N/Y=0/1
!      kbsp=nint(pars(15))   ! run MAMPOSSt in fast mode? N/Y=0/1

      xlabelgas=nint(pars(14))   ! gas mass:Y/N= 1/0 
      xlabelgal=nint(pars(15))   ! gal mass:Y/N= 1/0 
      rjicl=pars(16)          ! r_Jaffe ICL
      xlumicl=pars(17)        ! ICL total lum

!
!     Stop if H&M beta(r) required and M(r) difft from NFW
!
      if (kani.lt.0.and.knfit.ne.1) then
         write(*,*) ' Hansen & Moore beta(r) currently '
         write(*,*) ' implemented only for pNFW N(R) model '
         stop
      endif
         
!     read radial positions, velocities and velocity errors;
!     positions are in kpc, vels are in km/s, assumed rest-frame
!     errors are in km/s

      read(10,9) line
      read(10,9) line
  9   format(a)
      j=0
      j0=-1
      rlow=1.e12
      rup=-1.

 222  continue
      read(10,*,end=111) dkpc,vkms,evkms,wei
 !     radial selection
      if (dkpc/1.e3.ge.rlowin.and.dkpc/1.e3.le.rupin) then
         j=j+1
         r(j)=dkpc/1.e3
         v(j)=vkms
         e(j)=evkms
         w(j)=wei
         if (r(j).lt.rlow) rlow=r(j)
         if (r(j).gt.rup) rup=r(j)
      endif
      goto 222
 111  continue
      nga=j
      close(10)
      write(*,*) nga,' galaxies in the sample'
      write(*,*) 'in the radial range (Mpc): ',rlow,rup

!=== Read mass gas and mass star
!---  Read data

      open(12,file='data_mamp/massgas.dat',status='old')
      i=0
 777  continue
      read(12,*,end=999) rgasr,xmgasr,egasr
      i=i+1
      rgas(i)=rgasr
      xmgas(i)=xmgasr*1.e+13
      emgas(i)=egasr*1.e+13
      goto 777
 999  continue
      ngas=i

      close(12)
      
      open(14,file='data_mamp/massgal.dat',status='old')
      j=0
 732  continue
!      read(14,*,end=733) rgalr,xmgalr
      read(14,*,end=733) rgalr,xmgalr,emgalr
      j=j+1
      rgal(j)=rgalr
      xmgal(j)=xmgalr*1.e+13
      emgal(j)=egalr*1.e+13
      goto 732
 733  continue
      ngal=j

      close(14)
 

!===  VDP part
!---  Read parameters
      ! open(28,file='data_mamp/vdp_par.dat',status='old')

      ! do i=1,2
      !    read(28,*) parsvdp(i)
      ! enddo
      ! close(28)

!     BCG parameters:
!      ra=parsvdp(1)       ! anisotropy radius of beta(r) for the BCG (Mpc)
!     observational parameter:
!      spsf=parsvdp(1)     ! sigma of the PSF (Mpc)
!     numerical parameter:
!      rt=parsvdp(2)       ! limiting radius for the integration (Mpc)
      rt = 20.             ! limiting radius for the integration (Mpc)
!---  Read data


      write(*,*) ' File of BCG velocity disperion profile'
      write(*,*) '     (format: R [Mpc], sigma [km/s], sigma_error):'
      open(11,file='data_mamp/vel_bcg.dat',status='old')

      i=0
 734  continue
      read(11,*,end=735) xdata,ydata,edata
      i=i+1
      rbcg(i)=xdata
      sbcg(i)=ydata
      esbcg(i)=edata
      goto 734
 735  continue
      nbcg=i

      close(11)


!####Mio:fine

    do_cl_init = .false.

    print*,'PARAM IN: ',xlumbcg,za,h0,omega0,omegal,rt,rc
    print*,'End CL initialization'
    end subroutine Cl_init

!#######################################
! from vdmtovdpfuncmamp
      subroutine vdplik(nfreepar,freeparvdp,chibcg)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0,grav=4.302d-9)
      parameter  (nout=200)

!x      dimension pars(ipar),rbcg(iii),sbcg(iii),esbcg(iii)
      dimension bcgr(nout),bcgvdp(nout),rvdp(nout),vdp(nout)
      dimension freeparvdp(nfreepar)
!x      dimension r200a(imam),rsa(imam),gama(imam),xmla(imam),xlika(imam)

!x      character*60 fpar,fvdpbcg,fgridmam,fchiout

      include 'datagasgal.i'
      include 'datavdp.i'
      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'splinevdp.i'
      include 'limitsvdp.i'


!x      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
!x      common/spline/xx(500),yy(500),csc(500,3),icsc,nd
!x      common/limits/ax,rpsf

!x      external fintegrand,fmass,fbetaint,frho,frhoproj,fabbeta,fbeta
!x      external fpsfnum,fpsfden

      call uerset(1,levold)
      icsc=500


 !      write(*,*) ' File of grid values from MAMPOSSt: '
 !      read(*,9) fgridmam
 !      open(unit=98,file=fgridmam,status='old')

 !      i=0
 ! 634  continue
 !      read(98,*,end=635) x1,x2,x3,x4,x5,x6,x7,x8
 !      i=i+1
 !      r200a(i)=x1
 !      rsa(i)=x3
 !      gama(i)=x5
 !      xmla(i)=x6
 !      xlika(i)=x7
 !      goto 634
 ! 635  continue
 !      igrid=i

 !      close(98)

      ! write(*,*) ' Output file of chi^2 values: '
      ! read(*,9) fchiout
      ! open(unit=97,file=fchiout,status='unknown')



      nvdp=nout-1
!      do i=1,igrid

         r200=freeparvdp(1)
         hz=h0*sqrt(omega0*(1.+za)**3+omegal)                                      
         xm200=100.*hz*hz/grav*r200**3
         rs=freeparvdp(2)
         gam=freeparvdp(3)
         xmstarl=freeparvdp(4)
         ra=freeparvdp(5)
         rjaf=freeparvdp(6)

!         print*,'PARAM VDP: ',ra,rjaf,xlumbcg,spsf,za,h0,omega0,omegal,rt

         call vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)
         chiq=0.d0

         rlim = 2.*rbcg(nbcg)
         do k=1,iout
            if (bcgr(k).le.rlim) imax=k
         enddo
         do j=1,nbcg
            rbcgj=rbcg(j)
! ! Put these lines if you want the PSF convolution
!             if (rbcgj.gt.3.*spsf) then
! !     Use the VDP without convolution for PSF at large radii
! !     since at large radii the convolution does not work fine
! ! END               
               call hunt(rvdp,imax,rbcgj,jlo)         
!print*,' 2',rvdp(jlo),imax,rbcgj,jlo,vdp(jlo)
               sbcge= vdp(jlo)+ ( vdp(jlo+1)-vdp(jlo) )/(rvdp(jlo+1) - rvdp(jlo)) * (rbcgj - rvdp(jlo) )

! ! Put these lines if you want the PSF convolution
!             else

!                !     Use the VDP convoluted for PSF at small radii
               
!                call hunt(bcgr,imax,rbcgj,jlo)         
!                sbcge= bcgvdp(jlo)+ ( bcgvdp(jlo+1)-bcgvdp(jlo) )/(bcgr(jlo+1) - bcgr(jlo)) * (rbcgj - bcgr(jlo) )

!             endif
! ! END
            chiq=chiq+(sbcg(j)-sbcge)*(sbcg(j)-sbcge)/(esbcg(j)*esbcg(j))
         enddo

         ! if (chiq.lt.0.6) then
         !    do ii=1,imax
         !       write (34,*) bcgr(ii),bcgvdp(ii)
         !    enddo
         ! endif

 !         write(97,997) r200,rs,gam,xmstarl,xlika(i),chiq
 ! 997     format(4(F7.3,2x),F12.4,2x, D17.4)

!      enddo

!      close(97)                

!      stop

      chibcg = chiq

      return 
      end

!     
!     Solving the integral in log space to get ln(nu*vr^2)
!     
      subroutine vdpbcg(rvdp,vdp,bcgr,bcgvdp,nvdp,iout)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      ! parameter (pi=3.1415927d0,grav=4.302d-9,ipar=9,iii=5000)
      ! parameter (imam=500000,nout=200)

      parameter (pi=3.1415927d0,grav=4.302d-9,iii=5000)
      parameter (nout=200)

      dimension xmnul(iii),vr2nu(iii),u(iii),s(iii)
      dimension rvdp(nvdp),vdp(nvdp),sig2i(iii)
!      dimension rbcg(iii),sbcg(iii),esbcg(iii)
      dimension bcgr(nvdp),bcgvdp(nvdp)

      include 'datagasgal.i'
      include 'datavdp.i'
      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'splinevdp.i'
      include 'limitsvdp.i'

      ! common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
      ! common/spline/xx(500),yy(500),csc(500,3),icsc,nd
      ! common/limits/ax,rpsf

      ! external fintegrand,fmass,fbetaint,frho,frhoproj,fabbeta,fbeta
      ! external fpsfnum,fpsfden

      errabs=0.0d0
      errel=0.001d0
      bx=dlog(rt)
      a1=dlog(1.d-4)
      deltalog=(bx-a1)/nout
      lout=0

      do i=1,nout-1
         ax=a1+deltalog*i
         res=0.d0
         res=dcadre(fintegrand,ax,bx,errabs,errel,errest,ier)
         denbeta=ra*ra+dexp(2.*ax)  ! constant part of the integrated beta(r)
         res=res*grav/denbeta
         if (ier.gt.130) write(*,*) ' int_fvdm: ',dexp(ax),res
         if (res.eq.res) then 
            lout=lout+1
            xmnul(lout)=ax
            vr2nu(lout)=dlog(res)
         endif
         
      enddo
      
!
!     Spline fitting of ln(nu*vr^2)
!
      nd=lout-2
      do k=1,nd
         xx(k)=xmnul(k)
         yy(k)=vr2nu(k)
      enddo
      call icsccu(xx,yy,nd,csc,icsc,ier)


!     Abel projection of nu*vr^2 in nout points,
!     yielding I*VDP^2; then divide by I(R) and take the sqrt
      
      jout=0
      do j=1,lout-1
         ax=a1+deltalog*j
         res=0.d0
         res=dcadre(fabbeta,ax,bx,errabs,errel,errest,ier)
         if (ier.gt.130) write(*,*) ' int_fabbeta: ',dexp(ax),res
         if (res.lt.0.) res=0.
         a=dexp(ax)
         fir=frhoproj(a)
         if (res.eq.res) then 
            jout=jout+1
            rvdp(jout)=a
            vdp(jout)=dsqrt(2.d0*res/fir)
            sig2i(jout)=2.d0*res 
! Take away these lines if you want the PSF convolution
           bcgr(jout)=rvdp(jout)
           bcgvdp(jout)=vdp(jout)
! END
         endif
      enddo
! Take away these lines if you want the PSF convolution
     iout=jout
! END

! ! Put these lines if you want the PSF convolution
! !
! !     Spline fitting of ln(I*sigmap^2) vs. ln(R)
! !
!       nd=jout-2
!       do k=1,nd
!          xx(k)=dlog(rvdp(k))
!          yy(k)=dlog(sig2i(k))
!       enddo
!       call icsccu(xx,yy,nd,csc,icsc,ier)

! !     Convolution with PSF
! !

!       ax=a1

!       iout=0
!       do j=1,jout-1

!          rpsf=dexp(a1+deltalog*j)

!          rnum=0.d0
!          rnum=dcadre(fpsfnum,ax,bx,errabs,errel,errest,ier)
!          if (ier.gt.130) write(*,*) ' int_fpsf_num: ',rpsf,rnum
!          if (rnum.lt.0.) rnum=0.

!          den=0.d0
!          den=dcadre(fpsfden,ax,bx,errabs,errel,errest,ier)
!          if (ier.gt.130) write(*,*) ' int_fpsf_den: ',rpsf,den
!          if (den.lt.0.) den=0.
!          if (rnum.eq.rnum.and.den.eq.den.and.den.ne.0..and.rnum*den.ge.0.) then
!             iout=iout+1
!             bcgr(iout)=rpsf
!             bcgvdp(iout)=dsqrt(rnum/den)
! !            write(35,*) bcgr(iout),bcgvdp(iout)
!          endif
!       enddo
! ! END      

      return
      end


!     integrand function to get rho*sigma_r^2

      function fintegrand(xl)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
 !     external fmass,fbetaint,frho
      x=dexp(xl)
      fintegrand=frho(x)*fmass(x)/x*fbetaint(x)
      return
      end

!
!     mass profile (gNFW+Jaffe+Jaffe_ICL)
!     
      function fmass(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)


      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
 
      t=x/rs
      c=r200/rs

!     NFW M(r) - non serve, è il caso gNFW con gam=1.
!
!     fac200=1.d0/(dlog(1.d0+c)-c/(1.d0+c))
!
!     fnfw=(dlog(1.d0+t)-t/(1.d0+t))*fac200*xm200


!     Jaffe M(r)

      fjaffe=xmstarl*xlumbcg*x/rjaf/(1.+x/rjaf)
      fjaffeicl=xmstarl*xlumicl*x/rjicl/(1.+x/rjicl)

!     GAS the gas mass can be ignored with labelgas = 0 
      call hunt(rgas, ngas,x,jlo) 
      fmassgas =xlabelgas* xmgas(jlo)
!      print*,x,rgas(jlo),jlo,xmgas(jlo)," GAS"
!     GAL the gal mass can be ignored with labelgal = 0 
      call hunt(rgal, ngal,x,jlo1) 
      fmassgal =xlabelgal* xmgal(jlo1)

!     gNFW M(r), da eq.(6) in Walker et al. 2009, ApJ, 704, 1274

!      gam =  0.5d0
      bhyp1 = 3.d0-gam
      bhyp2 = 3.d0-gam
      ahyp = 4.d0-gam

      chyp = -c
      thyp = -t

      call HYGFX(bhyp1,bhyp2,ahyp,chyp,hypgeoc)
      call HYGFX(bhyp1,bhyp2,ahyp,thyp,hypgeot)
      
      fac200g=c**(3.d0-gam) * hypgeoc
      gnfw   =t**(3.d0-gam) * hypgeot

      gnfw=gnfw * xm200/fac200g

!      fmass=fnfw+fjaffe
!      write(*,*) x,gnfw/1.e+13,fjaffe/1.e+13,fmassgal/1.e+13,fmassgas/1.e+13," MASS"
      fmass=gnfw+fjaffe+fmassgal+fmassgas+fjaffeicl

      return
      end


!     exp(2*int(beta/x dx)) - OM model

!     !!! Note that the constant part of the integration is in the MAIN pgm

      function fbetaint(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
 
      fbetaint=ra*ra+x*x

      return
      end


!     density profile (nu(r) or rho*(r))
!     Jaffe

      function frho(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0)
      

      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf

!      frho=xmstarl*xlumbcg*rjaf/(4.*pi*x*x*(rjaf+x)*(rjaf+x))
      frho=rjaf/(4.*pi*x*x*(rjaf+x)*(rjaf+x))

      return
      end


!     projected density profile (N(R) or Sigma*(R) or I(R))
!     Jaffe

      function frhoproj(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pi=3.1415927d0)

      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf

      s=x/rjaf

!     R=rjaf (approximation by interplation)

      if (abs(s-1.).lt.0.001) then

         s1=0.999
         cs1=dacosh(1.d0/s1)
         s2=1.001
         cs2=dacos(1.d0/s2)
         
         frhoproj1=pi/(4.d0*s1)-(1.d0-(2.d0-s1*s1)*cs1/dsqrt(dabs(s1*s1-1.d0)))/2.d0/(s1*s1-1.d0)

         frhoproj2=pi/(4.d0*s2)-(1.d0-(2.d0-s2*s2)*cs2/dsqrt(dabs(s2*s2-1.d0)))/2.d0/(s2*s2-1.d0)

         frhoproj=(frhoproj2-frhoproj1)/(s2-s1)*(s-s1)+frhoproj1

!     R<>rjaf

      else

         if (s.lt.1.) then
            cs=dacosh(1.d0/s)
         else
            cs=dacos(1.d0/s)
         endif

         frhoproj=pi/(4.d0*s)-(1.d0-(2.d0-s*s)*cs/dsqrt(dabs(s*s-1.d0)))/2.d0/(s*s-1.d0)


!         write(*,*) ' normal: ',x,rjaf,s,frhoproj

      endif

!c      frhoproj=frhoproj*xmstarl*xlumbcg/(pi*rjaf*rjaf)
      frhoproj=frhoproj/(pi*rjaf*rjaf)

      return
      end


!     nu*<vr2>(r) Abel integrand in dln(x)
!     including beta(r)

      function fabbeta(xl)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1)

      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'
      include 'splinevdp.i'
      include 'limitsvdp.i'


      ! common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
      ! common/spline/xx(500),yy(500),csc(500,3),icsc,nd
      ! common/limits/ax,rpsf

!      external fbeta

      u(1)=xl
      call icsevu(xx,yy,nd,csc,icsc,u,s,1,ier)
      s1=s(1)

      x=dexp(xl)
      b1=fbeta(x)

      a=dexp(ax)

      if (abs(1.-x/a).lt.0.001) then
         fabbeta=0.d0
      else
         fabbeta=dexp(s1)*x*x/dsqrt(x*x-a*a)*(1.-b1*a*a/(x*x))
      endif

!      if (a.gt.0.1.and.a.lt.0.11) 
!     +     write(*,*) 'fabbeta: ',a,x,s1,b1,fabbeta

      return
      end

!     OM beta(r) model

      function fbeta(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'

!      common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
 
      fbeta=x*x/(ra*ra+x*x)

      return
      end

!     Modified Bessel function

      FUNCTION BESSI0(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7
      REAL*8 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (DABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=DABS(X)
        Y=3.75/AX
        BESSI0=(DEXP(AX)/DSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

!     Integrand of numerator for PSF convolution
!     Binney+Merrifield 1998 & de Zeeuw+Carollo 1996

      function fpsfnum(xl)      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1)


      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'
      include 'splinevdp.i'
      include 'limitsvdp.i'

      ! common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
      ! common/spline/xx(500),yy(500),csc(500,3),icsc,nd
      ! common/limits/ax,rpsf

!      external bessi0
      x=dexp(xl)

      u(1)=xl
      call icsevu(xx,yy,nd,csc,icsc,u,s,1,ier)
      s1=dexp(s(1)-50.d0)

      arg=rpsf*x/(spsf*spsf)
      CALL CALCI0(ARG,bessel,1)
!      besselold=bessi0(rpsf*x/(spsf*spsf))
      fexp=dexp(-(x*x+rpsf*rpsf)/(2.d0*spsf*spsf))

      fpsfnum=x*s1*bessel*fexp*x   ! *x for d(ln(x)) integration

!      write(*,*) arg,besselold,bessel

      return
      end


!     Integrand of denominator for PSF convolution

      function fpsfden(xl)      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1)

      include 'paramvdp.i'
      include 'paramcomm.i'
      include 'datagasgal.i'
      include 'splinevdp.i'
      include 'limitsvdp.i'

      ! common /params/rt,r200,rs,gam,ra,rjaf,xlumbcg,xmstarl,xm200,spsf
      ! common/spline/xx(500),yy(500),csc(500,3),icsc,nd
      ! common/limits/ax,rpsf

!      external bessi0,frhoproj

      x=dexp(xl)
      arg=rpsf*x/(spsf*spsf)
      CALL CALCI0(ARG,bessel,1)
!      besselold=bessi0(rpsf*x/(spsf*spsf))
      fexp=dexp(-(x*x+rpsf*rpsf)/(2.d0*spsf*spsf))

      fpsfden=x*frhoproj(x)*dexp(-50.d0)*bessel*fexp*x   ! *x for d(ln(x)) integration

!      write(*,*) x,frhoproj(x),bessel,fpsfden

      return
      end



      SUBROUTINE CALCI0(ARG,RESULT,JINT)
!--------------------------------------------------------------------

! This packet computes modified Bessel functions of the first kind
!   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!   arguments X.  It contains two function type subprograms, BESI0
!   and BESEI0, and one subroutine type subprogram, CALCI0.
!   The calling statements for the primary entries are 

!                   Y=BESI0(X)
!   and
!                   Y=BESEI0(X)

!   where the entry points correspond to the functions I0(X) and
!   EXP(-ABS(X))*I0(X), respectively.  The routine CALCI0 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCI0 with the statement
!          CALL CALCI0(ARG,RESULT,JINT)
!   where the parameter usage is as follows

!      Function                     Parameters for CALCI0
!       Call              ARG                  RESULT          JINT

!     BESI0(ARG)    ABS(ARG) .LE. XMAX        I0(ARG)           1
!     BESEI0(ARG)    any real ARG        EXP(-ABS(ARG))*I0(ARG) 2

!   The main computation evaluates slightly modified forms of
!   minimax approximations generated by Blair and Edwards, Chalk
!   River (Atomi! Energy of Canada Limited) Report AECL-4928,
!   October, 1974.  This transportable program is patterned after 
!   the machine-dependent FUNPACK packet NATSI0, but cannot match
!   that version for efficiency or accuracy.  This version uses
!   rational functions that theoretically approximate I-SUB-0(X)
!   to at least 18 significant decimal digits.  The accuracy
!   achieved depends on the arithmeti! system, the compiler, the
!   intrinsi! functions, and proper selection of the machine-
!   dependent constants.

!*******************************************************************
!*******************************************************************

! Explanation of machine-dependent constants

!   beta   = Radix for the floating-point system
!   maxexp = Smallest power of beta that overflows
!   XSMALL = Positive argument such that 1.0 - X = 1.0 to
!            machine precision for all ABS(X) .LE. XSMALL.
!   XINF =   Largest positive machine number; approximately
!            beta**maxexp
!   XMAX =   Largest argument acceptable to BESI0;  Solution to
!            equation:  
!               W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
!            where  W(X) = EXP(X)/SQRT(2*PI*X)


!     Approximate values for some important machines are:

!                          beta       maxexp       XSMALL

! CRAY-1        (S.P.)       2         8191       3.55E-15
! Cyber 180/855
!   under NOS   (S.P.)       2         1070       3.55E-15
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       2          128       2.98E-8
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       2         1024       5.55D-17
! IBM 3033      (D.P.)      16           63       6.95D-18
! VAX           (S.P.)       2          127       2.98E-8
! VAX D-Format  (D.P.)       2          127       6.95D-18
! VAX G-Format  (D.P.)       2         1023       5.55D-17


!                               XINF          XMAX

! CRAY-1        (S.P.)       5.45E+2465     5682.810
! Cyber 180/855
!   under NOS   (S.P.)       1.26E+322       745.893
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       3.40E+38         91.900
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       1.79D+308       713.986
! IBM 3033      (D.P.)       7.23D+75        178.182
! VAX           (S.P.)       1.70D+38         91.203
! VAX D-Format  (D.P.)       1.70D+38         91.203
! VAX G-Format  (D.P.)       8.98D+307       713.293

!*******************************************************************
!*******************************************************************
!
! Error returns

!  The program returns XINF for BESI0 for ABS(ARG) .GT. XMAX.


!  Intrinsic functions required are:

!     ABS, SQRT, EXP


!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!  Latest modification: June 7, 1988

!--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      INTEGER*4 I,JINT
!S    REAL
!$$$      DOUBLE PRECISION
!$$$     1       A,ARG,B,EXP40,FORTY,ONE,ONE5,P,PP,Q,QQ,RESULT,
!$$$     2       REC15,SUMP,SUMQ,TWO25,X,XINF,XMAX,XSMALL,XX
      DIMENSION P(15),PP(8),Q(5),QQ(7)
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
!S    DATA ONE/1.0E0/,ONE5/15.0E0/,EXP40/2.353852668370199854E17/,
!S   1     FORTY/40.0E0/,REC15/6.6666666666666666666E-2/,
!S   2     TWO25/225.0E0/
      DATA ONE/1.0D0/,ONE5/15.0D0/,EXP40/2.353852668370199854D17/,FORTY/40.0D0/,REC15/6.6666666666666666666D-2/,TWO25/225.0D0/
!--------------------------------------------------------------------
!  Machine-dependent constants
!--------------------------------------------------------------------
!S    DATA XSMALL/2.98E-8/,XINF/3.40E38/,XMAX/91.9E0/
      DATA XSMALL/5.55D-17/,XINF/1.79D308/,XMAX/713.986D0/
!--------------------------------------------------------------------
!  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
!--------------------------------------------------------------------
!S    DATA  P/-5.2487866627945699800E-18,-1.5982226675653184646E-14,
!S   1        -2.6843448573468483278E-11,-3.0517226450451067446E-08,
!S   2        -2.5172644670688975051E-05,-1.5453977791786851041E-02,
!S   3        -7.0935347449210549190E+00,-2.4125195876041896775E+03,
!S   4        -5.9545626019847898221E+05,-1.0313066708737980747E+08,
!S   5        -1.1912746104985237192E+10,-8.4925101247114157499E+11,
!S   6        -3.2940087627407749166E+13,-5.5050369673018427753E+14,
!S   7        -2.2335582639474375249E+15/
!S    DATA  Q/-3.7277560179962773046E+03, 6.5158506418655165707E+06,
!S   1        -6.5626560740833869295E+09, 3.7604188704092954661E+12,
!S   2        -9.7087946179594019126E+14/
      DATA  P/-5.2487866627945699800D-18,-1.5982226675653184646D-14,-2.6843448573468483278D-11,-3.0517226450451067446D-08,-2.5172644670688975051D-05,-1.5453977791786851041D-02,-7.0935347449210549190D+00,-2.4125195876041896775D+03,-5.9545626019847898221D+05,-1.0313066708737980747D+08,-1.1912746104985237192D+10,-8.4925101247114157499D+11,-3.2940087627407749166D+13,-5.5050369673018427753D+14,-2.2335582639474375249D+15/
      DATA  Q/-3.7277560179962773046D+03, 6.5158506418655165707D+06,-6.5626560740833869295D+09, 3.7604188704092954661D+12,-9.7087946179594019126D+14/
!--------------------------------------------------------------------
!  Coefficients for 15.0 .LE. ABS(ARG)
!--------------------------------------------------------------------
!S    DATA PP/-3.9843750000000000000E-01, 2.9205384596336793945E+00,
!S   1        -2.4708469169133954315E+00, 4.7914889422856814203E-01,
!S   2        -3.7384991926068969150E-03,-2.6801520353328635310E-03,
!S   3         9.9168777670983678974E-05,-2.1877128189032726730E-06/
!S    DATA QQ/-3.1446690275135491500E+01, 8.5539563258012929600E+01,
!S   1        -6.0228002066743340583E+01, 1.3982595353892851542E+01,
!S   2        -1.1151759188741312645E+00, 3.2547697594819615062E-02,
!S   3        -5.5194330231005480228E-04/
      DATA PP/-3.9843750000000000000D-01, 2.9205384596336793945D+00,-2.4708469169133954315D+00, 4.7914889422856814203D-01,-3.7384991926068969150D-03,-2.6801520353328635310D-03,9.9168777670983678974D-05,-2.1877128189032726730D-06/
      DATA QQ/-3.1446690275135491500D+01, 8.5539563258012929600D+01,-6.0228002066743340583D+01, 1.3982595353892851542D+01,-1.1151759188741312645D+00, 3.2547697594819615062D-02,-5.5194330231005480228D-04/
!--------------------------------------------------------------------
      X = ABS(ARG)
      IF (X .LT. XSMALL) THEN
            RESULT = ONE
         ELSE IF (X .LT. ONE5) THEN
!--------------------------------------------------------------------
!  XSMALL .LE.  ABS(ARG)  .LT. 15.0
!--------------------------------------------------------------------
            XX = X * X
            SUMP = P(1)
            DO 50 I = 2, 15
              SUMP = SUMP * XX + P(I)
   50       CONTINUE
            XX = XX - TWO25
            SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))*XX+Q(5)
            RESULT = SUMP / SUMQ
            IF (JINT .EQ. 2) RESULT = RESULT * EXP(-X)
         ELSE IF (X .GE. ONE5) THEN
            IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
                  RESULT = XINF
               ELSE
!--------------------------------------------------------------------
!  15.0  .LE.  ABS(ARG)
!--------------------------------------------------------------------
                  XX = ONE / X - REC15
                  SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+PP(4))*XX+PP(5))*XX+PP(6))*XX+PP(7))*XX+PP(8)
                  SUMQ = ((((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+QQ(4))*XX+QQ(5))*XX+QQ(6))*XX+QQ(7)
                  RESULT = SUMP / SUMQ
                  IF (JINT .EQ. 2) THEN
                        RESULT = (RESULT - PP(1)) / SQRT(X)
                     ELSE
!--------------------------------------------------------------------
!  Calculation reformulated to avoid premature overflow
!--------------------------------------------------------------------
                        IF (X .LE.(XMAX-ONE5)) THEN
                              A = EXP(X)
                              B = ONE
                           ELSE
                              A = EXP(X-FORTY)
                              B = EXP40
                        END IF
                        RESULT = ((RESULT*A-PP(1)*A)/SQRT(X))*B
                  END IF
            END IF
      END IF
!--------------------------------------------------------------------
!  Return for ABS(ARG) .LT. XSMALL
!--------------------------------------------------------------------
      RETURN
!----------- Last line of CALCI0 -----------
      END
!S    REAL 
      DOUBLE PRECISION FUNCTION BESI0(X)
!--------------------------------------------------------------------
!
! This long precision subprogram computes approximate values for
!   modified Bessel functions of the first kind of order zero for
!   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI0).
!
!--------------------------------------------------------------------
      INTEGER*4 JINT
!S    REAL  
      DOUBLE PRECISION X, RESULT
!--------------------------------------------------------------------
      JINT=1
      CALL CALCI0(X,RESULT,JINT)
      BESI0=RESULT
      RETURN
!---------- Last line of BESI0 ----------
      END
!S    REAL 
      DOUBLE PRECISION FUNCTION BESEI0(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   modified Bessel function of the first kind of order zero
!   multiplied by EXP(-ABS(X)), where EXP is the
!   exponential function, ABS is the absolute value, and X
!   is any argument.
!
!--------------------------------------------------------------------
      INTEGER JINT
!S    REAL  
      DOUBLE PRECISION X, RESULT
!--------------------------------------------------------------------
      JINT=2
      CALL CALCI0(X,RESULT,JINT)
      BESEI0=RESULT
      RETURN
!---------- Last line of BESEI0 ----------
      END




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

!#############################################


      subroutine mamlik(n,x,f)                                                 
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
      include 'datarv.i'
      include 'datagasgal.i'                                                       
      include 'paramcomm.i'                                                    
      include 'sr.i'                                                           
      include 'vkr.i'                                                          

      call uerset(0,levold)                                                    
      icsr = 100                                                               
      icsrk = 100                                                              
                                                                               
                                                                               
!     Max Lik fit: start by computing the sigma_r(r) function in               
!     Ninterp points logarithmically spaced between 0.001 and 20               
                                                                               
      ic = rknots                                                              
      ninterp=21*2                                                             
                                                                               
      hz=h0*sqrt(omega0*(1.+za)**3+omegal)                                     
                                                                               
      r200=10.**x(1)                                                           
      rjaf=10.**x(2)                                                             
      rs=10.**x(3)                                                             
      cbe=10.**x(4)                                                            
      gam=10.**x(5)                                                            
      xmstarl=10.**x(6)                                                        
                                                                               
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
         xris(i) = dlog(rlow)+dfloat(i-1)* dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)                         
         xmin=dexp(xris(i))                                                    
                                                                               
!     dqdagi integrates from a lower bound (xmin) to +infinity                 
!     while dqdag does not                                                     
                                                                               
         xx1=xmin                                                              
         risl = dcadre(sr2int,xris(i),dlog(2.*rinfinity),errabs,errrel,errest,ier)                                                          
         if (ier .gt. 128) then                                                
            print *,'VMAXLIK-SR2INT: rmin=',xx1,' rmax=',xx2                   
            stop                                                               
         endif                                                                 
         risok=dsqrt(risl*sr2out(xmin))                                        
         if (risl.gt.1.8d195) risok=rismax                                     
         if (risl.le.0.d0) risok=rismin                                        
         yris(i)=dlog(risok)                                                   
      enddo                                                                    
                                                                               
!                                                                              
! compute spline coeffs for later interpolation of sigma_r                     
!                                                                              
      ier=0                                                                    
      call icsccu(xris,yris,ninterp,csr,icsr,ier)                              
      if (ier .gt. 128) then                                                   
         print *,' in VMAXLIK ...'                                             
         do i = 1, ninterp                                                     
            print *,'ierr  i xris yris = ', ier, i, xris(i), yris(i)           
         enddo                                                                 
      endif                                                                    
                                                                               
!     initialize psum = Sum[-log(lik)] and wsum = Sum[galaxy_weights]          
                                                                               
      psum=0.d0                                                                
      wsum=0.d0                                                                
                                                                               
                                                                               
!     Now compute the distribution function at the radial distance of          
!     each galaxy by interpolating the sigma_r(r) and then                     
!     evaluate the Max Lik, using weights if required                          
                                                                               
                                                                               
!     Use all the original data, no interpolation                              
                                                                               
         do j=1,ngambm                                                         
            rj=rmbm(j)                                                         
            vj=vmbm(j)                                                         
            ej=embm(j)                                                         
            xmin=rj                                                            
                                                                               
!     re-define rinfinity as the radius where v_z/sigma_r leq 0.1              
!     to avoid a problem in v_z=0; if we require to consider the               
!     interloper contribution, rinfinity is set to 1 already                   
                                                                               
!            if (kintd.eq.0) call rinf(rinfinity)  
!            changed with                             
             call rinf(rinfinity)                               
                                                                               
!     check the limits of the integral                                         
                                                                               
            if (rj.lt.rinfinity) then                                          
               wsum=wsum+wmbm(j)                                               
                                                                               
!     compute the integral of gwen, yielding g(R,vz)                           
                                                                               
               umax = dacosh(rinfinity/rj)                                     
               g = rj*dcadre(gwenu,0.d0,umax,errabs,errrel,errest,ier)         
               if (ier .gt. 128) then                                          
                  print *,'VMAXLIK-GWEN2: rmin=',xx1,' rmax=',xx2              
                  print *,'ier &  g = ',ier,g                                  
!                  stop                                                        
               endif                                                           
                                                                               

                                                                                                                                                            
!     From g, one then has the probability p                                   
                                                                               
!     There are two expressions: one if for the case in which c_nu             
!     is part of the fitting parameters of MAMPOSSt                            
!     the other is for the case in which c_nu is fitted externally             
                                                                               
                                                                               
                  if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then           
                     sigmaden=sigmar1(rj)                                      
                  elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then       
                     sigmaden=sigmar2(rj)                                      
                  else                                                         
                     sigmaden=sigmar3(rj)                                      
                  endif                                                        
                  p=g/sigmaden                                                 
                                                                               
                                                                                
               call temp(rj,sigrt1,sigrt2,sigzt1t2,bec)                        
                                                                               
                                                                               
               psum=psum-dlog(p)*wmbm(j)                                       
                                                                               
                                                                               
            endif                                                              
                                                                               
         enddo                                                                 
                                                                               
!                                                                              
!     The factor nga/wsum is =1 if no weights, but is not 1                    
!     if there are weights, and is needed to take into account                 
!     the real number of points (otherwise, we would estimate                  
!     a Max Lik too low by an average factor of 1/Ngr, where                   
!     Ngr is the average number of galaxies per group, since                   
!     the likelihood we want to estimate is the Sum of the                     
!     ngambm galaxy probabilities)                                             
!                                                                              
                                                                               
      f=psum/wsum*ngambm                                                       
                                                                               
      return                                                                   
      end                                                                      
!                                                                              

!   imsl routine name   - mgamad=dgamma

!   computer            - vax/double

!   latest revision     - june 1, 1981

!   purpose             - evaluate the gamma function of a double
!                           precision argument

!   usage               - result = dgamma(x)

!   arguments    x      - input double precision argument.
!                           dgamma is set to machine infinity, with the
!                           proper sign, whenever
!                             x is zero,
!                             x is a negative integer,
!                             abs(x) .le. xmin, or
!                             abs(x) .ge. xmax.
!                             xmin is of the order of 10**(-39) and
!                             xmax is at least 34.8. the exact values
!                             of xmin and xmax may allow larger ranges
!                             for x on some computers.
!                             see the programming notes in the manual
!                             for the exact values.
!                dgamma - output double precision value of the gamma
!                           function. dgamma must be typed double
!                           precision in the calling program.

!   precision/hardware  - double/h32                                    
!                       - not available/h36,h48,h60                     
!                         note - dgamma may not be supplied by imsl if
!                           it resides in the mathematical subprogram
!                           library supplied by the manufacturer.

!   reqd. imsl routines - uertst,ugetio

!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp

!   remarks      an error message printed by uertst from dgamma should
!                be interpreted as follows
!                ier    - error indicator
!                         terminal error
!                           ier = 129 indicates that the absolute value
!                             of the input argument x was specified
!                             greater than or equal to xmax. dgamma
!                             is set to machine infinity.
!                           ier = 130 indicates that the input argument
!                             x was specified as zero or a negative
!                             integer or that the absolute value of the
!                             input argument was less than or equal to
!                             xmin. dgamma is set to machine infinity
!                             with the proper sign for the dgamma
!                             function. if x is zero or an even
!                             negative integer, gamma has a negative    
!                             sign. otherwise it has a positive sign.   

!   copyright           - 1978 by imsl, inc. all rights reserved.

!   warranty            - imsl warrants only that imsl testing has been 
!                           applied to this code. no other warranty,
!                           expressed or implied, is applicable.

      function dgamma (x)
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                              
      dimension          p(9),q(8),p4(7)
      logical            mflag
!                                  coefficients for minimax
!                                  approximation to gamma(x),
!                                  2.0 .le. x .le. 3.0
      data               p(1)/-5.966047488753637d01/,p(2)/5.864023793062003d01/,p(3)/-1.364106217165365d03/,p(4)/-8.117569271425580d02/,p(5)/-1.569414683149179d04/,p(6)/-1.525979925758372d04/,p(7)/-7.264059615964330d04/,p(8)/-8.972275718101010d-01/,p(9)/3.349618189847578d00/
      data               q(1)/4.103991474182904d02/,q(2)/-2.262590291514875d03/,q(3)/2.494325576714903d03/,q(4)/2.362106244383048d04/,q(5)/-5.741873227396418d04/,q(6)/-7.257239715408240d04/,q(7)/-9.491399521686949d00/,q(8)/-3.255006939455704d01/
!                                  coefficients for minimax
!                                  approximation to ln(gamma(x)),
!                                  12.0 .le. x
      data               p4(1)/8.40596949829d-04/,p4(2)/-5.9523334141881d-04/,p4(3)/7.9365078409227d-04/,p4(4)/-2.777777777769526d-03/,p4(5)/8.333333333333333d-02/,p4(6)/9.189385332046727d-01/,p4(7)/-1.7816839846d-03/
      data               pi/3.141592653589793d0/
      data               xinf/1.7d+38/
!                                  dgamma(xmin) .approx. xinf
!                                  dgamma(big1) .approx. xinf
      data               xmin/5.8775d-39/
      data               big1/34.844d0/
!                                  first executable statement
      ier = 0
      mflag = .false.
      t = x
      if (dabs(t).gt.xmin) go to 5
      ier = 130
      dgamma = xinf
      if (t.le.0.0d0) dgamma = -xinf
      go to 9000
    5 if (dabs(t).lt.big1) go to 10
      ier = 129
      dgamma = xinf
      go to 9000
   10 if (t.gt.0.0d0) go to 25
!                                  argument is negative
      mflag = .true.
      t = -t
      r = dint(t)
      sign = 1.0d0
      if (dmod(t,2.0d0).eq.0.0d0) sign = -1.0d0
      r = t-r
      if (r.ne.0.0d0) go to 20
      ier = 130
      dgamma = xinf
      if (sign.eq.-1.0d0) dgamma = -xinf
      go to 9000
!                                  argument is not a negative integer
   20 r = pi/dsin(r*pi)*sign
      t = t+1.0d0
!                                  evaluate approximation for dgamma(t)
!                                    t .gt. xmin
   25 if (t.gt.12.0d0) go to 60
      i = t
      a = 1.0d0
      if (i.gt.2) go to 40
      i = i+1
      go to (30,35,50),i
!                                  0.0 .lt. t .lt. 1.0
   30 a = a/(t*(t+1.0d0))
      t = t+2.0d0
      go to 50
!                                  1.0 .le. t .lt. 2.0
   35 a = a/t
      t = t+1.0d0
      go to 50
!                                  3.0 .le. t .le. 12.0
   40 do 45 j=3,i
         t = t-1.0d0
         a = a*t
   45 continue
!                                  2.0 .le. t .le. 3.0
   50 top = p(8)*t+p(9)
      den = t+q(8)
      do 55 j=1,7
         top = top*t+p(j)
         den = den*t+q(j)
   55 continue
      y = (top/den)*a
      if (mflag) y = r/y
      dgamma = y
      go to 9005
!                                  t .gt. 12.0
   60 top = dlog(t)
      top = t*(top-1.0d0)-.5d0*top
      t = 1.0d0/t
      b = t*t
      a = p4(7)
      do 65 j = 1,5
   65 a = a*b+p4(j)
      y = a*t+p4(6)+top
      y = dexp(y)
      if (mflag) y = r/y
      dgamma = y
      go to 9005
 9000 continue
      call uertst(-ier,6hmgamad)                                        
      call uertst(ier,6hdgamma)                                         
 9005 return
      end

!                                                                              
!     Function gammarec: uses recursive relation                               
!                Gamma(x)=Gamma(x+1)/x                                         
!     when x<0 to call dgamma(arg) with arg>0                                  
!     


      function gammarec(x)                                                     
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)
!      external dgamma
                                                                               
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
                                                                               
!                                                                              
!                                                                              
!                                                                              
      function sdint(x)                                                        
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      include 'paramcomm.i'                                                    
      include 'datarv.i'
      include 'datagasgal.i'                                                       
      data pig /3.1415926535897932d0/                                          
                                                                               
      aint=10.**(-1.061+0.364*x*x-0.580*x**4+0.533*x**6)                       
      sigmaint=0.612-0.0653*x*x                                                
      bint=0.0075                                                              
      xkappa=4.                                                                
                                                                               
      argerf=xkappa/dsqrt(2.d0)/sigmaint                                       
                                                                               
      sdint=x*(dsqrt(pig/2.d0)*aint*sigmaint*erf(argerf)+xkappa*bint)          
                                                                               
      return                                                                   
      end                                                                      
!                                                                              
!     compute the difference in 'projected mass'                               
!     at the extreme of the radial range, for                                  
!     3 different Sigma(R) profiles                                            
!                                                                              
      subroutine sigmarnorm(x,fnorm)                                           
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      include 'paramcomm.i'                                                    
      include 'datarv.i'
      include 'datagasgal.i'                                                       
      data pig /3.1415926535897932d0/                                          
                                                                               
                                                                               
      if (knfit.eq.3) then                                                     
                                                                               
!     beta model                                                               
                                                                               
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
                                                                               
!     NFW                                                                      
                                                                               
      elseif ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then                   
                                                                               
         c=r200/x                                                              
         gc=1./(dlog(c+1.)-c/(c+1.))                                           
         cx=c*rup/r200                                                         
         uu=1./cx                                                              
                                                                               
!c         write(*,*) ' c_tracer is ',c                                        
                                                                               
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
!ccc         fnorm=(fnup-fnlow)*gc/(r200*r200)                                 
         fnorm=(fnup-fnlow)*gc                                                 
                                                                               
      else                                                                     
                                                                               
!     Hernquist                                                                
                                                                               
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
                                                                               
!c      write(*,*) ' rlow, rup ',rlow,rup                                      
                                                                               
      return                                                                   
      end                                                                      
                                                                               
!                                                                              
!     N(R) - projected Hernquist - see Hernquist (1990)                        
!                                                                              
      function sigmar2(tt)                                                     
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      parameter (pig=3.1415926535897932d0)                                     
      include 'paramcomm.i'                                                    
      include 'datagasgal.i'
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
                                                                               
!                                                                              
!     N(R) - projected NFW; use eq. (41), (42) in Lokas & Mamon (2001)         
!                                                                              
      function sigmar1(tt)                                                     
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      parameter (pig=3.1415926535897932d0)                                     
      include 'paramcomm.i'                                                    
      include 'datagasgal.i'
      t=tt/rc                                                                  
                                                                               
!     no interlopers considered                                                
                                                                               
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
!                                                                              
!     The factor outside the integral for the                                  
!     sigma_r^2 formulae                                                       
!                                                                              
      function sr2out(t)                                                       
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      parameter (pig=3.1415926535897932d0)                                     
      include 'paramcomm.i'                                                    
      include 'datagasgal.i'
!      external gammarec                                                        
                                                                               
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then                       
         c=r200/rc                                                             
         gc=1./(dlog(c+1)-c/(c+1))                                             
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)         
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then                   
         xnu=rc/(2.*pig*t)/(t+rc)**3.                                          
      else                                                                     
         xnu=-1.d0/dsqrt(pig*1.d0)*gammarec(-al+0.5d0)/gammarec(-al+1.d0)*al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)                          
      endif                                                                    
                                                                               
      if (kani.eq.0) then                                                      
!                                                                              
!     Convert from cbe=beta' to bec=beta                                       
!     in the case of constant anisotropy                                       
!                                                                              
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         app=1.d-2                                                             
!                                                                              
!     Radial anisotropy                                                        
!                                                                              
         if(dabs(bec-1.d0).lt.app) then                                        
            sr2out=1./(t*t)                                                    
!                                                                              
!     Isotropic case                                                           
!                                                                              
         elseif (dabs(bec).lt.app) then                                        
            sr2out=1.                                                          
!                                                                              
!     General constant anisotropy                                              
!                                                                              
         else                                                                  
            sr2out=t**(-2.*bec)                                                
         endif                                                                 
!                                                                              
!     Osipkov Merritt                                                          
!                                                                              
      elseif (kani.eq.2) then                                                  
         sr2out=1./(t*t+cbe*cbe)                                               
!                                                                              
!     Mamon Lokas anisotropy profiles                                          
!                                                                              
      elseif (kani.eq.1) then                                                  
         sr2out=1./(t+cbe)                                                     
!                                                                              
!     Simplified Wojtak anisotropy profile                                     
!                                                                              
      elseif (kani.eq.3) then                                                  
         sr2out=(t**1.5+cbe**1.5)**(-4./3.)                                    
!                                                                              
!     Simplified Tiret (modified ML)                                           
!                                                                              
      elseif (kani.eq.4) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         sr2out=(t+rm2)**(-2.*bec)                                             
!                                                                              
!     modfied Tiret (non-zero central anisotropy)                              
!                                                                              
      elseif (kani.eq.5) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values               
         sr2out=(t+rm2)**(-4.*bec)*t**(2.*bec)                                 
      else                                                                     
!                                                                              
!     Hansen+Moore relation for NFW mass profile                               
!                                                                              
         rhm=rc ! use nu(r)                                                    
!c         rhm=rs ! use rho(r)                                                 
         sr2out=t**(2.*(bhm-ahm))*(rhm+t)**(4.*bhm) ! My solution              
      endif                                                                    
!                                                                              
!     Divide by the density profile                                            
!                                                                              
      sr2out=sr2out/xnu                                                        
      return                                                                   
      end                                                                      
                                                                               
!                                                                              
!     Integrand function for the determination of                              
!     nu(r)*sigma_r^2(r), given beta(r) and M(r)                               
!     from eqs. A3, A4, A6 in Mamon & Lokas (2005)                             
!     or from eq.(A3) in Mamon, Biviano & Boue' (2013)                         
!                                                                              
!       modified by GAM to handle integration of ln r                          
!                                                                              
      function sr2int(alr)                                                     
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      parameter (pig=3.1415926535897932d0,grav=4.302e-9)                       
      dimension rsvalues(28),r100values(28)                                    
      include 'paramcomm.i'                                                   
      include 'datagasgal.i'

      data rsvalues/0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.30,0.35,0.40,0.45,0.50,1.00/                                      
      data r100values/1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476/           
!      external gammarec                                                        
!                                                                              
                                                                               
      gm200=r200*v200*v200                                                     
!                                                                              
                                                                               
      t = dexp(alr)                                                            
      if (kmp.eq.1) then                                                       
!                                                                              
!     DM mass profile is genealized-NFW                                        
!                                                                              
!                                                                              
!     If not requested, the mass profile becomes                               
!     identical to NFW with the normalisation                                  
!     set to v200^2*r200/G at r/r200=1                                         
!                                                                              
                                                                               
            bhyp1 = 3.d0-gam                                                   
            bhyp2 = 3.d0-gam                                                   
            ahyp = 4.d0-gam                                                    
                                                                               
            chyp = -r200/rs                                                    
            thyp = -t/rs                                                       
            call HYGFX(bhyp1,bhyp2,ahyp,chyp,hypgeoc)                          
            call HYGFX(bhyp1,bhyp2,ahyp,thyp,hypgeot)                          
                                                                               
            fac200g=(r200/rs)**(3.d0-gam) * hypgeoc                            
            gnfw   =(t/rs)**(3.d0-gam) * hypgeot                               
                                                                               
            xm=gnfw/fac200g                                                    
                                                                               
!$$$            fac200=1./(dlog(1.d0+r200/rs)-(r200/rs)/(1.d0+r200/rs))        
!$$$            xm=(dlog(1.d0+t/rs)-t/rs/(1.d0+t/rs))*fac200                   

      elseif (kmp.eq.2) then                                                   
!                                                                              
!     M(r) is Hernquist; no gas contribution allowed for the time being        
!     Normalisation set to 1 at r/r200=1                                       
!                                                                              
         fac200=(r200+rs)*(r200+rs)/(r200*r200)                                
         xm=t*t/(t+rs)/(t+rs)*fac200                                           
      elseif (kmp.eq.3) then                                                   
!                                                                              
!     M(r) is PIEMD; no gas contribution allowed for the time begin            
!     Normalisation set to 1 at r/r200=1                                       
!                                                                              
         fac200=1./(rcut*datan(r200/rcut)-rs*datan(r200/rs))                   
         xm=fac200*(rcut*datan(t/rcut)-rs*datan(t/rs))                         
      elseif (kmp.eq.4) then                                                   
!                                                                              
!     M(r) is Burkert; no gas contribution allowed for the time begin          
!     Normalisation set to 1 at r/r200=1                                       
!                                                                              
         trs=t/rs                                                              
         rvrs=r200/rs                                                          
         fac200=1./(dlog(1+rvrs*rvrs)+2.*dlog(1.+rvrs)-2.*datan(rvrs))         
         xm=fac200*(dlog(1+trs*trs)+2.*dlog(1.+trs)-2.*datan(trs))             
      elseif (kmp.eq.5) then                                                   
!                                                                              
!     M(r) is Soft Isoth Sph; no gas contribution allowed for the time begin   
!     Normalisation set to 1 at r/r200=1                                       
!                                                                              
         fac200=1./(r200-rs*datan(r200/rs))                                    
         xm=fac200*(t-rs*datan(t/rs))                                          
      else                                                                     
!                                                                              
!     M(r) is Einasto m=5; no gas contribution allowed for the time begin      
!     Normalisation set to 1 at r/r200=1                                       
!                                                                              
         eim=5.                                                                
         agp=3.*eim                                                            
         xgp=2.*eim*(r200/rs)**(1./eim)                                        
         call dincog(agp,xgp,gin,gim,gip200)                                   
         fac200=1./gip200                                                      
         xgp=2.*eim*(t/rs)**(1./eim)                                           
         call dincog(agp,xgp,gin,gim,gip)                                      
         xm=fac200*gip                                                         
!c         write(*,999) r200,rs,t,3.*eim,2.*eim*(r200/rs)**(1./eim),           
!c     &        2.*eim*(t/rs)**(1./eim),gip200,gip,xm                          
!c 999     format(4(1x,f5.2),2(1x,f6.2),3(1x,d11.3))                           
!c         eim=5.                                                              
!c         fac200=1./dgammp(3.*eim,2.*eim*(r200/rs)**(1./eim))                 
!c         xm=fac200*dgammp(3.*eim,2.*eim*(t/rs)**(1./eim))                    
      endif                                                                    
                                                                               
!     add the mass of the BCG in units of G*M200                               
!     (set the free param xmstarl to zero if this is not required)             
                                                                               
      xmbcg=xmstarl*xlumbcg*t/rjaf/(1.+t/rjaf)*grav/gm200

!     GAS the gas mass can be ignored with labelgas = 0 
      call hunt(rgas, ngas,t,jlo)
      fmassgas = xlabelgas*xmgas(jlo)*grav/gm200

!     GAL the gal mass can be ignored with labelgal = 0 
      call hunt(rgal, ngal,t,jlo1)
      fmassgal = xlabelgal*xmgal(jlo1) *grav/gm200     
                                                                               
      xm=xm+xmbcg+fmassgas+fmassgal                                                        
!                                                                              
!     nu(r) is beta-model, NFW or Hernquist                                    
!                                                                              
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then                       
         c=r200/rc                                                             
         gc=1./(dlog(c+1)-c/(c+1))                                             
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)         
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then                   
         xnu=rc/(2.*pig*t)/(t+rc)**3.                                          
      else                                                                     
         xnu=-1.d0/dsqrt(pig*1.d0)*gammarec(-al+0.5d0)/gammarec(-al+1.d0)*al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)                          
      endif                                                                    
                                                                               
!                                                                              
!     beta(r)                                                                  
!                                                                              
!                                                                              
!     Convert from cbe=beta' to bec=beta                                       
!     in the case of constant anisotropy                                       
!                                                                              
      if (kani.eq.0) then                                                      
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         app=1.d-2                                                             
!                                                                              
!     Radial anisotropy                                                        
!                                                                              
         if(dabs(bec-1.d0).lt.app) then                                        
            sr2int=xnu*xm                                                      
!                                                                              
!     Isotropic case                                                           
!                                                                              
         elseif (dabs(bec).lt.app) then                                        
            sr2int=xnu*xm/(t*t)                                                
!                                                                              
!     General constant anisotropy                                              
!                                                                              
         else                                                                  
            sr2int=xnu*xm*t**(2.*bec-2.)                                       
         endif                                                                 
!                                                                              
!     Osipkov Merritt                                                          
!                                                                              
      elseif (kani.eq.2) then                                                  
         sr2int=xnu*xm*(t*t+cbe*cbe)/(t*t)                                     
!                                                                              
!     Mamon Lokas anisotropy profile                                           
!                                                                              
      elseif (kani.eq.1) then                                                  
         sr2int=xnu*xm*(t+cbe)/(t*t)                                           
!                                                                              
!     Simplified Wojtak anisotropy profile                                     
!                                                                              
      elseif (kani.eq.3) then                                                  
         sr2int=xnu*xm*(t**1.5+cbe**1.5)**(4./3.)/(t*t)                        
!                                                                              
!     Simplified Tiret (modified ML profile)                                   
!                                                                              
      elseif (kani.eq.4) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         sr2int=xnu*xm*(t+rm2)**(2.*bec)/(t*t)                                 
!                                                                              
!     modified Tiret (non-zero central anisotropy)                             
!                                                                              
      elseif (kani.eq.5) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values               
         sr2int=xnu*xm*(t+rm2)**(4.*bec)*t**(-2.-2.*bec)                       
!                                                                              
!     Hansen+Moore relation for NFW mass profile                               
!                                                                              
      else                                                                     
         rhm=rc ! use nu(r)                                                    
!c         rhm=rs ! use rho(r)                                                 
         sr2int=xnu*xm*t**(2.*(ahm-bhm-1.))*(rhm+t)**(-4.*bhm)                 
      endif                                                                    
                                                                               
      sr2int = t*sr2int*gm200    ! integral in dlog and G*M scaling            
!                                                                              
      return                                                                   
      end                                                                      
!                                                                              
! GWENU returns integrand of g(R,vz) integrated over u = cosh-1 (r/R)          
!                                                                              
!     Integrand to derive the f(R,vlos)                                        
!     distribution function for given velocity                                 
!     vlos, given anisotropy and mass profile                                  
!     (from Gwenael Boue')                                                     
!                                                                              
      double precision function gwenu(u)                                       
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      dimension rjl(1), yrjl(1)                                                
      parameter (pig=3.1415926535897932d0)                                     
      include 'paramcomm.i'                                                    
      include 'sr.i'   
      include 'datagasgal.i'                                                        
!      external gammarec                                                        
!                                                                              
      t = xmin*dcosh(u)                                                        
                                                                               
!c      write(*,*) ' vel and error is',vj,ej                                   
                                                                               
!                                                                              
!     tracer density nu                                                        
!                                                                              
      if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then                       
         c=r200/rc                                                             
         gc=1./(dlog(c+1)-c/(c+1))                                             
         xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)         
      elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then                   
         xnu=rc/(2.*pig*t)/(t+rc)**3.                                          
      else                                                                     
         xnu=-1.d0/dsqrt(pig*1.d0)*gammarec(-al+0.5d0)/gammarec(-al+1.d0)*al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)                          
      endif                                                                    
!                                                                              
!     Interpolate to get the sigma_r(t)                                        
!                                                                              
      rjl(1)=dlog(t)                                                           
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)                   
      if (ier .eq. 34) then                                                    
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',rjl(1)                                                               
      endif                                                                    
      sigr=dexp(yrjl(1))                                                       
                                                                               
!                                                                              
!     beta(r)                                                                  
!                                                                              
!     Convert from cbe=beta' to bec=beta                                       
!     in the case of constant anisotropy                                       
!                                                                              
      if (kani.eq.0) then                                                      
         bec=1.-1./(cbe*cbe)                                                   
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
!                                                                              
!     Osipkov Merritt                                                          
!                                                                              
      elseif (kani.eq.2) then                                                  
         bec=t*t/(t*t+cbe*cbe)                                                 
!                                                                              
!     Mamon Lokas                                                              
!                                                                              
      elseif (kani.eq.1) then                                                  
         bec=0.5*t/(t+cbe)                                                     
!                                                                              
!     Simplified Wojtak                                                        
!                                                                              
      elseif (kani.eq.3) then                                                  
         bec=t**1.5/(t**1.5+cbe**1.5)                                          
!                                                                              
!     Simplified Tiret or modified ML                                          
!                                                                              
      elseif (kani.eq.4) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=(1.-1./(cbe*cbe))*t/(t+rm2)                                       
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
!                                                                              
!     modified Tiret (non-zero central anisotropy)                             
!                                                                              
      elseif (kani.eq.5) then                                                  
         rm2=rs                     ! SIS or other unspecified mass profiles
         if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
         if (kmp.eq.2) rm2=0.5*rs   ! Hernquist                                
         if (kmp.eq.4) rm2=1.521*rs ! Burkert                                  
         bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)                                 
         if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values                 
         if (bec.le.-1.) bec=-0.999d0  ! avoid unphysical values               
!                                                                              
!     Hansen+Moore relation for NFW mass profile                               
!                                                                              
      else                                                                     
         rhm=rc ! use nu(r)                                                    
!c         rhm=rs ! use rho(r)                                                 
         bec=ahm-bhm*(rhm+3.*t)/(rhm+t)                                        
      endif                                                                    
!                                                                              
!     sigma_z(R,r)                                                             
!                                                                              
      sigmaz=sigr*dsqrt(1.-bec*xmin*xmin/(t*t))                                
                                                                               
!     add velocity error                                                       
                                                                               
      sigmaz=sqrt(sigmaz*sigmaz+ej*ej)                                         
                                                                               
!                                                                              
!     Integrand function (note that sigmaz is in units of v200, just           
!     as vj)                                                                   
                                                                               
         gwenu=2.d0*xnu/dsqrt(2.*pig)/sigmaz*dcosh(u)*dexp(-vj*vj/(2.*sigmaz*sigmaz))                                 
                                                                               
                                                                               
      return                                                                   
      end                                                                      
                                                                               
!                                                                              
!     temp                                                                     
!                                                                              
      subroutine temp(t1,sigrt1,sigrt2,sigzt1t2,bec)                           
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      parameter (pig=3.1415926535897932d0)                                     
      real*8 rjl(1), yrjl(1)                                                   
      include 'paramcomm.i'                                                    
      include 'sr.i'                                                           
      include 'datagasgal.i'
!                                                                              
      t2 = 2.*t1                                                               
!                                                                              
!     Interpolate to get the sigma_r(t)                                        
!                                                                              
      rjl(1)=dlog(t1)                                                          
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)                   
      if (ier .eq. 34) then                                                    
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=',rjl(1)                                                               
      endif                                                                    
      sigrt1=dexp(yrjl(1))                                                     
                                                                               
      rjl(1)=dlog(t2)                                                          
      call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)                   
      if (ier .eq. 34) then                                                    
         print *,' ICSEVU in GWENU: max xris=', xris(ninterp), ' r=', rjl(1)                                                              
      endif                                                                    
      sigrt2=dexp(yrjl(1))                                                     
                                                                               
!                                                                              
!     beta(r)                                                                  
!                                                                              
!     Convert from cbe=beta' to bec=beta                                       
!     in the case of constant anisotropy                                       
!                                                                              
      bec=1.-1./(cbe*cbe)                                                      
      if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values                     
!                                                                              
!     sigma_z(R,r)                                                             
!                                                                              
      sigzt1t2=sigrt1*dsqrt(1.-bec*t1*t1/(t2*t2))                              
                                                                               
      return                                                                   
      end                                                                      
                                                                               
!   imsl routine name   - dcadre                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - vax/double                                           
!                                                                              
!   latest revision     - june 1, 1982                                         
!                                                                              
!   purpose             - numerical integration of a function using            
!                           cautious adaptive romberg extrapolation            
!                                                                              
!   usage               - function dcadre (f,a,b,aerr,rerr,error,ier)          
!                                                                              
!   arguments    dcadre - estimate of the integral of f(x) from a to b.        
!                           (output).                                          
!                f      - a single-argument real function subprogram           
!                           supplied by the user. (input)                      
!                           f must be declared external in the                 
!                           calling program.                                   
!                a,b    - the two endpoints of the interval of                 
!                           integration. (input)                               
!                aerr   - desired absolute error in the answer. (input)        
!                rerr   - desired relative error in the answer. (input)        
!                error  - estimated bound on the absolute error of             
!                           the output number, dcadre. (output)                
!                ier    - error parameter. (output)                            
!                         warning error(with fix)                              
!                           ier = 65 implies that one or more                  
!                             singularities were successfully handled.         
!                           ier = 66 implies that, in some                     
!                             subinterval(s), the estimate of the              
!                             integral was accepted merely because the         
!                             estimated error was small, even though no        
!                             regular behavior was recognized.                 
!                         terminal error                                       
!                           ier = 131 indicates failure due to                 
!                             insufficient internal working storage.           
!                           ier = 132 indicates failure due to                 
!                             too much noise in the function (relative         
!                             to the given error requirements) or              
!                             due to an ill-behaved integrand.                 
!                           ier = 133 indicates that rerr is greater           
!                             than 0.1, or rerr is less than 0.0, or           
!                             rerr is too small for the precision of           
!                             the machine.                                     
!                                                                              
!   precision/hardware  - single and double/h32                                
!                       - single/h36,h48,h60                                   
!                                                                              
!   reqd. imsl routines - uertst,ugetio                                        
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   remarks  1.  dcadre can, in many cases, handle jump                        
!                discontinuities. see document reference for full              
!                details.                                                      
!            2.  the relative error parameter rerr must be in the              
!                interval (0.0,0.1) inclusively. for example,                  
!                rerr = 0.1 indicates that the estimate of the                 
!                integral is to be correct to one digit, whereas               
!                rerr = .0001 calls for four digits of accuracy.               
!                if dcadre determines that the relative accuracy               
!                requirements cannot be satisfied, ier is set to               
!                133 (rerr should be large enough that, when added             
!                to 100.0, the result is a number greater than                 
!                100.0).                                                       
!            3.  the absolute error parameter, aerr, should be non-            
!                negative. in order to give a reasonable value for             
!                aerr, the user must know the approximate magnitude            
!                of the integral being computed. in many cases it is           
!                satisfactory to use aerr = 0.0. in this case, only            
!                the relative error requirement is satisfied in the            
!                computation.                                                  
!            4.  even when ier is not equal to 0, dcadre returns the           
!                best estimate that has been computed.                         
!                quoting from the document reference- a very cautious          
!                man would accept dcadre only if ier is 0 or 65. the           
!                merely reasonable man would keep the faith even if            
!                ier is 66. the adventurous man is quite often right           
!                in accepting dcadre even if ier is 131 or 132.                
!            5.  dcadre may return wrong answers if f has a periodic           
!                factor with high frequency and the interval (a,b)             
!                contains an integral number of periods. in this case          
!                the easiest fix is to divide the interval into two            
!                subintervals (a,c) and (c,b) such that neither                
!                contains an integral number of periods (pick c at             
!                random), and call dcadre to integrate over each               
!                subinterval.                                                  
!                                                                              
!   copyright           - 1978 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      double precision function dcadre (f,a,b,aerr,rerr,error,ier)             
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
!                                  specifications for local variables          
      dimension            ibegs(30)                                           
      real*8 length,jumptl                                                     
      dimension   t(10,10),r(10),ait(10),dif(10),rn(4),ts(2049)                
      dimension   begin(30),finis(30),est(30)                                  
      logical            h2conv,aitken,right,reglar,reglsv(30)                 
      data               aitlow,h2tol,aittol,jumptl,maxts,maxtbl,mxstge/1.1d0,.15d0,.1d0,.01d0,2049,10,30/                   
      data               rn(1),rn(2),rn(3),rn(4)/.7142005d0,.3466282d0,.843751d0,.1263305d0/           
      data               zero,p1,half,one,two,four,fourp5,ten,hun/0.0d0,0.1d0,0.5d0,1.0d0,2.0d0,4.0d0,4.5d0,10.0d0,100.0d0/                                 
!                                  first executable statement                  
      ier = 0                                                                  
      cadre = zero                                                             
      error = zero                                                             
      curest = zero                                                            
      vint = zero                                                              
      length = dabs(b-a)                                                       
      if (length .eq. zero) go to 215                                          
      if (rerr .gt. p1 .or. rerr .lt. zero) go to 210                          
      hrerr = rerr+hun                                                         
      if (aerr .eq. zero .and. hrerr .le. hun) go to 210                       
      errr = rerr                                                              
      erra = dabs(aerr)                                                        
      stepmn = length/(two**mxstge)                                            
      stepnm = dmax1(length,dabs(a),dabs(b))*ten                               
      stage = half                                                             
      istage = 1                                                               
      fnsize = zero                                                            
      prever = zero                                                            
      reglar = .false.                                                         
!                                  the given interval of integration           
!                                    is the first interval considered.         
      beg = a                                                                  
      fbeg = f(beg)*half                                                       
      ts(1) = fbeg                                                             
      ibeg = 1                                                                 
      edn = b                                                                  
      fend = f(edn)*half                                                       
      ts(2) = fend                                                             
      iend = 2                                                                 
    5 right = .false.                                                          
!                                  investigation of a particular               
!                                    subinterval begins at this point.         
   10 step = edn - beg                                                         
      astep =  dabs(step)                                                      
      if (astep .lt. stepmn) go to 205                                         
      hrerr = stepnm+astep                                                     
      if (hrerr .eq. stepnm) go to 205                                         
      t(1,1) = fbeg + fend                                                     
      tabs = dabs(fbeg) + dabs(fend)                                           
      l = 1                                                                    
      n = 1                                                                    
      h2conv = .false.                                                         
      aitken = .false.                                                         
   15 lm1 = l                                                                  
      l = l + 1                                                                
!                                  calculate the next trapezoid sum,           
!                                    t(l,1), which is based on *n2* + 1        
!                                    equispaced points. here,                  
!                                    n2 = n*2 = 2**(l-1).                      
      n2 = n+n                                                                 
      fn = n2                                                                  
      istep = (iend - ibeg)/n                                                  
      if (istep .gt. 1) go to 25                                               
      ii = iend                                                                
      iend = iend + n                                                          
      if (iend .gt. maxts) go to 200                                           
      hovn = step/fn                                                           
      iii = iend                                                               
      fi = one                                                                 
      do 20 i=1,n2,2                                                           
         ts(iii) = ts(ii)                                                      
         ts(iii-1) = f(edn - fi * hovn)                                        
         fi = fi+two                                                           
         iii = iii-2                                                           
         ii = ii-1                                                             
   20 continue                                                                 
      istep = 2                                                                
   25 istep2 = ibeg + istep/2                                                  
      sum = zero                                                               
      sumabs = zero                                                            
      do 30 i=istep2,iend,istep                                                
         sum = sum + ts(i)                                                     
         sumabs = sumabs + dabs(ts(i))                                         
   30 continue                                                                 
      t(l,1) = t(l-1,1)*half+sum/fn                                            
      tabs = tabs*half+sumabs/fn                                               
      n = n2                                                                   
!                                  get preliminary value for *vint*            
!                                    from last trapezoid sum and update        
!                                    the error requirement *ergoal*            
!                                    for this subinterval.                     
      it = 1                                                                   
      vint = step*t(l,1)                                                       
      tabtlm = tabs*ten                                                        
      fnsize = dmax1(fnsize,dabs(t(l,1)))                                      
      ergl = astep*fnsize*ten                                                  
      ergoal = stage*dmax1(erra,errr*dabs(curest+vint))                        
!                                  complete row l and column l of *t*          
!                                    array.                                    
      fextrp = one                                                             
      do 35 i=1,lm1                                                            
         fextrp = fextrp*four                                                  
         t(i,l) = t(l,i) - t(l-1,i)                                            
         t(l,i+1) = t(l,i) + t(i,l)/(fextrp-one)                               
   35 continue                                                                 
      errer = astep*dabs(t(1,l))                                               
!                                  preliminary decision procedure              
!                                    if l = 2 and t(2,1) = t(1,1),             
!                                    go to 135 to follow up the                
!                                    impression that intergrand is             
!                                    straight line.                            
      if (l .gt. 2) go to 40                                                   
      hrerr = tabs+p1*dabs(t(1,2))                                             
      if (hrerr .eq. tabs) go to 135                                           
      go to 15                                                                 
!                                  caculate next ratios for                    
!                                    columns 1,...,l-2 of t-table              
!                                    ratio is set to zero if difference        
!                                    in last two entries of column is          
!                                    about zero                                
   40 do 45 i=2,lm1                                                            
         diff = zero                                                           
         hrerr = tabtlm+dabs(t(i-1,l))                                         
         if (hrerr .ne. tabtlm) diff = t(i-1,lm1)/t(i-1,l)                     
         t(i-1,lm1) = diff                                                     
   45 continue                                                                 
      if (dabs(four-t(1,lm1)) .le. h2tol) go to 60                             
      if (t(1,lm1) .eq. zero) go to 55                                         
      if (dabs(two-dabs(t(1,lm1))) .lt. jumptl) go to 130                      
      if (l .eq. 3) go to 15                                                   
      h2conv = .false.                                                         
      if (dabs((t(1,lm1)-t(1,l-2))/t(1,lm1)) .le. aittol) go to 75             
   50 if (reglar) go to 55                                                     
      if (l .eq. 4) go to 15                                                   
      hrerr = ergl+errer                                                       
   55 if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 175                   
      go to 145                                                                
!                                  cautious romberg extrapolation              
   60 if (h2conv) go to 65                                                     
      aitken = .false.                                                         
      h2conv = .true.                                                          
   65 fextrp = four                                                            
   70 it = it + 1                                                              
      vint = step*t(l,it)                                                      
      errer = dabs(step/(fextrp-one)*t(it-1,l))                                
      if (errer .le. ergoal) go to 160                                         
      hrerr = ergl+errer                                                       
      if (hrerr .eq. ergl) go to 160                                           
      if (it .eq. lm1) go to 125                                               
      if (t(it,lm1) .eq. zero) go to 70                                        
      if (t(it,lm1) .le. fextrp) go to 125                                     
      if (dabs(t(it,lm1)/four-fextrp)/fextrp .lt. aittol) fextrp = fextrp*four                                              
      go to 70                                                                 
!                                  integrand may have x**alpha type            
!                                    singularity                               
!                                    resulting in a ratio of *sing*  =         
!                                    2**(alpha + 1)                            
   75 if (t(1,lm1) .lt. aitlow) go to 175                                      
      if (aitken) go to 80                                                     
      h2conv = .false.                                                         
      aitken = .true.                                                          
   80 fextrp = t(l-2,lm1)                                                      
      if (fextrp .gt. fourp5) go to 65                                         
      if (fextrp .lt. aitlow) go to 175                                        
      if (dabs(fextrp-t(l-3,lm1))/t(1,lm1) .gt. h2tol) go to 175               
      sing = fextrp                                                            
      fextm1 = one/(fextrp - one)                                              
      ait(1) = zero                                                            
      do 85 i=2,l                                                              
         ait(i) = t(i,1) + (t(i,1)-t(i-1,1))*fextm1                            
         r(i) = t(1,i-1)                                                       
         dif(i) = ait(i) - ait(i-1)                                            
   85 continue                                                                 
      it = 2                                                                   
   90 vint = step*ait(l)                                                       
      errer = errer*fextm1                                                     
      hrerr = ergl+errer                                                       
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 95                    
      ier = max0(ier,65)                                                       
      go to 160                                                                
   95 it = it + 1                                                              
      if (it .eq. lm1) go to 125                                               
      if (it .gt. 3) go to 100                                                 
      h2next = four                                                            
      singnx = sing+sing                                                       
  100 if (h2next .lt. singnx) go to 105                                        
      fextrp = singnx                                                          
      singnx = singnx+singnx                                                   
      go to 110                                                                
  105 fextrp = h2next                                                          
      h2next = four*h2next                                                     
  110 do 115 i=it,lm1                                                          
         r(i+1) = zero                                                         
         hrerr = tabtlm+dabs(dif(i+1))                                         
         if (hrerr .ne. tabtlm) r(i+1) = dif(i)/dif(i+1)                       
  115 continue                                                                 
      h2tfex = -h2tol*fextrp                                                   
      if (r(l) - fextrp .lt. h2tfex) go to 125                                 
      if (r(l-1)-fextrp .lt. h2tfex) go to 125                                 
      errer = astep*dabs(dif(l))                                               
      fextm1 = one/(fextrp - one)                                              
      do 120 i=it,l                                                            
         ait(i) = ait(i) + dif(i)*fextm1                                       
         dif(i) = ait(i) - ait(i-1)                                            
  120 continue                                                                 
      go to 90                                                                 
!                                  current trapezoid sum and resulting         
!                                    extrapolated values did not give          
!                                    a small enough *errer*.                   
!                                    note -- having prever .lt. errer          
!                                    is an almost certain sign of              
!                                    beginning trouble with in the func-       
!                                    tion values. hence, a watch for,          
!                                    and control of, noise should              
!                                    begin here.                               
  125 fextrp = dmax1(prever/errer,aitlow)                                      
      prever = errer                                                           
      if (l .lt. 5) go to 15                                                   
      if (l-it .gt. 2 .and. istage .lt. mxstge) go to 170                      
      erret = errer/(fextrp**(maxtbl-l))                                       
      hrerr = ergl+erret                                                       
      if (erret .gt. ergoal .and. hrerr .ne. ergl) go to 170                   
      go to 15                                                                 
!                                  integrand has jump (see notes)              
  130 hrerr = ergl+errer                                                       
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 170                   
!                                  note that  2*fn = 2**l                      
      diff = dabs(t(1,l))*(fn+fn)                                              
      go to 160                                                                
!                                  integrand is straight line                  
!                                    test this assumption by comparing         
!                                    the value of the integrand at             
!                                    four *randomly chosen* points with        
!                                    the value of the straight line            
!                                    interpolating the integrand at the        
!                                    two end points of the sub-interval.       
!                                    if test is passed, accept *vint*          
  135 slope = (fend-fbeg)*two                                                  
      fbeg2 = fbeg+fbeg                                                        
      do 140 i=1,4                                                             
         diff = dabs(f(beg+rn(i)*step) - fbeg2-rn(i)*slope)                    
         hrerr = tabtlm+diff                                                   
         if(hrerr .ne. tabtlm) go to 155                                       
  140 continue                                                                 
      go to 160                                                                
!                                  noise may be dominant feature               
!                                    estimate noise level by comparing         
!                                    the value of the integrand at             
!                                    four *randomly chosen* points with        
!                                    the value of the straight line            
!                                    interpolating the integrand at the        
!                                    two endpoints. if small enough,           
!                                    accept *vint*                             
  145 slope = (fend-fbeg)*two                                                  
      fbeg2 = fbeg+fbeg                                                        
      i = 1                                                                    
  150 diff = dabs(f(beg+rn(i)*step) - fbeg2-rn(i)*slope)                       
  155 errer = dmax1(errer,astep*diff)                                          
      hrerr = ergl+errer                                                       
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 175                   
      i = i+1                                                                  
      if (i .le. 4) go to 150                                                  
      ier = 66                                                                 
!                                  intergration over current sub-              
!                                    interval successful                       
!                                    add *vint* to *dcadre* and *errer*        
!                                    to *error*, then set up next sub-         
!                                    interval, if any.                         
  160 cadre = cadre + vint                                                     
      error = error + errer                                                    
      if (right) go to 165                                                     
      istage = istage - 1                                                      
      if (istage .eq. 0) go to 220                                             
      reglar = reglsv(istage)                                                  
      beg = begin(istage)                                                      
      edn = finis(istage)                                                      
      curest = curest - est(istage+1) + vint                                   
      iend = ibeg - 1                                                          
      fend = ts(iend)                                                          
      ibeg = ibegs(istage)                                                     
      go to 180                                                                
  165 curest = curest + vint                                                   
      stage = stage+stage                                                      
      iend = ibeg                                                              
      ibeg = ibegs(istage)                                                     
      edn = beg                                                                
      beg = begin(istage)                                                      
      fend = fbeg                                                              
      fbeg = ts(ibeg)                                                          
      go to 5                                                                  
!                                  integration over current subinterval        
!                                    is unsuccessful. mark subinterval         
!                                    for further subdivision. set up           
!                                    next subinterval.                         
  170 reglar = .true.                                                          
  175 if (istage .eq. mxstge) go to 205                                        
      if (right) go to 185                                                     
      reglsv(istage+1) = reglar                                                
      begin(istage) = beg                                                      
      ibegs(istage) = ibeg                                                     
      stage = stage*half                                                       
  180 right = .true.                                                           
      beg = (beg+edn)*half                                                     
      ibeg = (ibeg+iend)/2                                                     
      ts(ibeg) = ts(ibeg)*half                                                 
      fbeg = ts(ibeg)                                                          
      go to 10                                                                 
  185 nnleft = ibeg - ibegs(istage)                                            
      if (iend+nnleft .ge. maxts) go to 200                                    
      iii = ibegs(istage)                                                      
      ii = iend                                                                
      do 190 i=iii,ibeg                                                        
         ii = ii + 1                                                           
         ts(ii) = ts(i)                                                        
  190 continue                                                                 
      do 195 i=ibeg,ii                                                         
         ts(iii) = ts(i)                                                       
         iii = iii + 1                                                         
  195 continue                                                                 
      iend = iend + 1                                                          
      ibeg = iend - nnleft                                                     
      fend = fbeg                                                              
      fbeg = ts(ibeg)                                                          
      finis(istage) = edn                                                      
      edn = beg                                                                
      beg = begin(istage)                                                      
      begin(istage) = edn                                                      
      reglsv(istage) = reglar                                                  
      istage = istage + 1                                                      
      reglar = reglsv(istage)                                                  
      est(istage) = vint                                                       
      curest = curest + est(istage)                                            
      go to 5                                                                  
!                                  failure to handle given integra-            
!                                    tion problem                              
  200 ier = 131                                                                
      go to 215                                                                
  205 ier = 132                                                                
      go to 215                                                                
  210 ier = 133                                                                
  215 cadre = curest + vint                                                    
  220 dcadre = cadre                                                           
 9000 continue                                                                 
!cc      if (ier .ne. 0) call uertst (ier,6hdcadre)                            
 9005 return                                                                   
      end                                                                      
                                                                               
!   imsl routine name   - icsccu                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - vax/double                                           
!                                                                              
!   latest revision     - june 1, 1980                                         
!                                                                              
!   purpose             - cubic spline interpolation                           
!                           (easy-to-use version)                              
!                                                                              
!   usage               - call icsccu (x,y,nx,c,ic,ier)                        
!                                                                              
!   arguments    x      - vector of length nx containing the abscissae         
!                           of the nx data points (x(i),y(i)) i=1,...,         
!                           nx. (input) x must be ordered so that              
!                           x(i) .lt. x(i+1).                                  
!                y      - vector of length nx containing the ordinates         
!                           (or function values) of the nx data points.        
!                           (input)                                            
!                nx     - number of elements in x and y. (input) nx            
!                           must be .ge. 2.                                    
!                c      - spline coefficients. (output) c is an nx-1 by        
!                           3 matrix. the value of the spline                  
!                           approximation at t is                              
!                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)         
!                           where x(i) .le. t .lt. x(i+1) and                  
!                           d = t-x(i).                                        
!                ic     - row dimension of matrix c exactly as                 
!                           specified in the dimension statement in            
!                           the calling program. (input)                       
!                ier    - error parameter. (output)                            
!                         terminal error                                       
!                           ier = 129, ic is less than nx-1.                   
!                           ier = 130, nx is less than 2.                      
!                           ier = 131, input abscissa are not ordered          
!                             so that x(1) .lt. x(2) ... .lt. x(nx).           
!                                                                              
!   precision/hardware  - single and double/h32                                
!                       - single/h36,h48,h60                                   
!                                                                              
!   reqd. imsl routines - uertst,ugetio                                        
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   copyright           - 1980 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      subroutine icsccu (x,y,nx,c,ic,ier)                                      
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      dimension x(nx),y(nx),c(ic,3)                                            
!                                  specifications for local variables          
      dimension cnx(3)                                                         
!                                  first executable statement                  
      nm1 = nx-1                                                               
      ier = 129                                                                
      if (ic .lt. nm1) go to 9000                                              
      ier = 130                                                                
      if (nx .lt. 2) go to 9000                                                
      ier = 131                                                                
      if (nx .eq. 2) go to 45                                                  
!                                  compute not-a-knot spline                   
      do 5 m = 2,nm1                                                           
         mm1=m-1                                                               
         c(m,2) = x(m)-x(mm1)                                                  
         if (c(m,2).le.0.0d0) go to 9000                                       
         c(m,3) = (y(m)-y(mm1))/c(m,2)                                         
    5 continue                                                                 
      cnx(2) = x(nx)-x(nm1)                                                    
      if (cnx(2).le.0.0d0) go to 9000                                          
      cnx(3) = (y(nx)-y(nm1))/cnx(2)                                           
      ier = 0                                                                  
      nm2 = nx-2                                                               
      if (nx .gt. 3) go to 10                                                  
      c(1,3) = cnx(2)                                                          
      c(1,2) = c(2,2)+cnx(2)                                                   
      c(1,1) = ((c(2,2)+2.d0*c(1,2))*c(2,3)*cnx(2)+c(2,2)**2*cnx(3))/c(1,2)                                                                  
      go to 20                                                                 
   10 c(1,3) = c(3,2)                                                          
      c(1,2) = c(2,2)+c(3,2)                                                   
      c(1,1) = ((c(2,2)+2.d0*c(1,2))*c(2,3)*c(3,2)+c(2,2)**2*c(3,3))/c(1,2)                                                                  
      do 15 m=2,nm2                                                            
         mp1=m+1                                                               
         mm1=m-1                                                               
         g = -c(mp1,2)/c(mm1,3)                                                
         c(m,1) = g*c(mm1,1)+3.d0*c(m,2)*c(mp1,3)+3.d0*c(mp1,2)*c(m,3)         
         c(m,3) = g*c(mm1,2)+2.d0*c(m,2)+2.d0*c(mp1,2)                         
   15 continue                                                                 
   20 g = -cnx(2)/c(nm2,3)                                                     
      c(nm1,1) = g*c(nm2,1)+3.d0*c(nm1,2)*cnx(3)+3.d0*cnx(2)*c(nm1,3)          
      c(nm1,3) = g*c(nm2,2)+2.d0*c(nm1,2)+2.d0*cnx(2)                          
      if (nx.gt.3) go to 25                                                    
      cnx(1)=2.d0*cnx(3)                                                       
      cnx(3)=1.d0                                                              
      g=-1.d0/c(nm1,3)                                                         
      go to 30                                                                 
   25 g = c(nm1,2)+cnx(2)                                                      
      cnx(1) = ((cnx(2)+2.d0*g)*cnx(3)*c(nm1,2)+cnx(2)**2*(y(nm1)-y(nx-2))/c(nm1,2))/g                                             
      g = -g/c(nm1,3)                                                          
      cnx(3) = c(nm1,2)                                                        
   30 cnx(3) = g*c(nm1,2)+cnx(3)                                               
      cnx(1) = (g*c(nm1,1)+cnx(1))/cnx(3)                                      
      c(nm1,1) = (c(nm1,1)-c(nm1,2)*cnx(1))/c(nm1,3)                           
      do 35 jj=1,nm2                                                           
         j = nm1-jj                                                            
         c(j,1) = (c(j,1)-c(j,2)*c(j+1,1))/c(j,3)                              
   35 continue                                                                 
      do 40 i=2,nm1                                                            
         im1 = i-1                                                             
         dtau = c(i,2)                                                         
         divdf1 = (y(i)-y(im1))/dtau                                           
         divdf3 = c(im1,1)+c(i,1)-2.d0*divdf1                                  
         c(im1,2) = (divdf1-c(im1,1)-divdf3)/dtau                              
         c(im1,3) = divdf3/dtau**2                                             
   40 continue                                                                 
      dtau = cnx(2)                                                            
      divdf1 = (y(nx)-y(nm1))/dtau                                             
      divdf3 = c(nm1,1)+cnx(1)-2.d0*divdf1                                     
      c(nm1,2) = (divdf1-c(nm1,1)-divdf3)/dtau                                 
      c(nm1,3) = divdf3/dtau**2                                                
      go to 9005                                                               
   45 if (x(1) .ge. x(2)) go to 9000                                           
      ier = 0                                                                  
      c(1,1) = (y(2)-y(1))/(x(2)-x(1))                                         
      c(1,2) = 0.0d0                                                           
      c(1,3) = 0.0d0                                                           
      go to 9005                                                               
 9000 continue                                                                 
      call uertst(ier,6hicsccu)                                                
 9005 return                                                                   
      end                                                                      
                                                                                                                                                           
!                                                                              
!   imsl routine name   - icsevu                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - vax/double                                           
!                                                                              
!   latest revision     - january 1, 1978                                      
!                                                                              
!   purpose             - evaluation of a cubic spline                         
!                                                                              
!   usage               - call icsevu(x,y,nx,c,ic,u,s,m,ier)                   
!                                                                              
!   arguments    x      - vector of length nx containing the abscissae         
!                           of the nx data points (x(i),y(i)) i=1,...,         
!                           nx (input). x must be ordered so that              
!                           x(i) .lt. x(i+1).                                  
!                y      - vector of length nx containing the ordinates         
!                           (or function values) of the nx data points         
!                           (input).                                           
!                nx     - number of elements in x and y (input).               
!                           nx must be .ge. 2.                                 
!                c      - spline coefficients (input). c is an nx-1 by         
!                           3 matrix.                                          
!                ic     - row dimension of matrix c exactly as                 
!                           specified in the dimension statement               
!                           in the calling program (input).                    
!                           ic must be .ge. nx-1.                              
!                u      - vector of length m containing the abscissae          
!                           of the m points at which the cubic spline          
!                           is to be evaluated (input).                        
!                s      - vector of length m (output).                         
!                           the value of the spline approximation at           
!                           u(i) is                                            
!                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)         
!                           where x(j) .le. u(i) .lt. x(j+1) and               
!                           d = u(i)-x(j).                                     
!                m      - number of elements in u and s (input).               
!                ier    - error parameter (output).                            
!                         warning error                                        
!                           ier = 33, u(i) is less than x(1).                  
!                           ier = 34, u(i) is greater than x(nx).              
!                                                                              
!   precision/hardware  - single and double/h32                                
!                       - single/h36,h48,h60                                   
!                                                                              
!   reqd. imsl routines - uertst,ugetio                                        
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   remarks  1.  the routine assumes that the abscissae of the nx              
!                data points are ordered such that x(i) is less than           
!                x(i+1) for i=1,...,nx-1. no check of this condition           
!                is made in the routine. unordered abscissae will cause        
!                the algorithm to produce incorrect results.                   
!            2.  the routine generates two warning errors. one error           
!                occurs if u(i) is less than x(1), for some i in the           
!                the interval (1,m) inclusively. the other error occurs        
!                if u(i) is greater than x(nx), for some i in the              
!                interval (1,m) inclusively.                                   
!            3.  the ordinate y(nx) is not used by the routine. for            
!                u(k) .gt. x(nx-1), the value of the spline, s(k), is          
!                given by                                                      
!                 s(k)=((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)         
!                where d=u(k)-x(nx-1).                                         
!                                                                              
!   copyright           - 1978 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      subroutine icsevu  (x,y,nx,c,ic,u,s,m,ier)                               
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      dimension x(nx),y(nx),c(ic,3),u(m),s(m)                                  
      data               i/1/,zero/0.0d0/                                      
!                                  first executable statement                  
      jer = 0                                                                  
      ker = 0                                                                  
      if (m .le. 0) go to 9005                                                 
      nxm1 = nx-1                                                              
      if (i .gt. nxm1) i = 1                                                   
!                                  evaluate spline at m points                 
      do 40 k=1,m                                                              
!                                  find the proper interval                    
         d = u(k)-x(i)                                                         
         if (d) 5,25,15                                                        
    5    if (i .eq. 1) go to 30                                                
         i = i-1                                                               
         d = u(k)-x(i)                                                         
         if (d) 5,25,20                                                        
   10    i = i+1                                                               
         d = dd                                                                
   15    if (i .ge. nx) go to 35                                               
         dd = u(k)-x(i+1)                                                      
         if (dd .ge. zero) go to 10                                            
         if (d .eq. zero) go to 25                                             
!                                  perform evaluation                          
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)                            
         go to 40                                                              
   25    s(k) = y(i)                                                           
         go to 40                                                              
!                                  warning - u(i) .lt. x(1)                    
   30    jer = 33                                                              
         go to 20                                                              
!                                  if u(i) .gt. x(nx) - warning                
   35    if (dd .gt. zero) ker = 34                                            
         d = u(k)-x(nxm1)                                                      
         i = nxm1                                                              
         go to 20                                                              
   40 continue                                                                 
      ier = max0(jer,ker)                                                      
 9000 continue                                                                 
      if (jer .gt. 0) call uertst(jer,6hicsevu)                                
      if (ker .gt. 0) call uertst(ker,6hicsevu)                                
 9005 return                                                                   
      end                                                                      
!   imsl routine name   - uerset                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - vax/single                                           
!                                                                              
!   latest revision     - january 1, 1978                                      
!                                                                              
!   purpose             - set message level for imsl routine uertst            
!                                                                              
!   usage               - call uerset (level,levold)                           
!                                                                              
!   arguments    level  - new value for message level. (input)                 
!                           output from imsl routine uertst is                 
!                           controlled selectively as follows,                 
!                             level = 4 causes all messages to be              
!                                       printed,                               
!                             level = 3 messages are printed if ier is         
!                                       greater than 32,                       
!                             level = 2 messages are printed if ier is         
!                                       greater than 64,                       
!                             level = 1 messages are printed if ier is         
!                                       greater than 128,                      
!                             level = 0 all message printing is                
!                                       suppressed.                            
!                levold - previous message level. (output)                     
!                                                                              
!   precision/hardware  - single/all                                           
!                                                                              
!   reqd. imsl routines - uertst,ugetio                                        
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   copyright           - 1978 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      subroutine uerset (level,levold)                                         
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
                                                                               
!                                  first executable statement                  
      levold = level                                                           
      call uertst (levold,6huerset)                                            
      return                                                                   
      end                                                                      
!   imsl routine name   - uertst                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - ibm/single                                           
!                                                                              
!   latest revision     - january 1, 1978                                      
! **                    - May 1982 by sjr                                      
!                                                                              
!                                                                              
!   purpose             - print a message reflecting an error condition        
!                                                                              
!   usage               - call uertst (ier,name)                               
!                                                                              
!   arguments    ier    - error parameter. (input)                             
!                           ier = i+j where                                    
!                             i = 128 implies terminal error,                  
!                             i =  64 implies warning with fix, and            
!                             i =  32 implies warning.                         
!                             j = error code relevant to calling               
!                                 routine.                                     
!                name   - a six character literal string giving the            
!                           name of the calling routine. (input)               
!                                                                              
!   precision/hardware  - single/all                                           
!                                                                              
!   reqd. imsl routines - ugetio                                               
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   remarks      the error message produced by uertst is written               
!                onto the standard output unit. the output unit                
!                number can be determined by calling ugetio as                 
!                follows..   call ugetio(1,nin,nout).                          
!                the output unit number can be changed by calling              
!                ugetio as follows..                                           
!                                nin = 0                                       
!                                nout = new output unit number                 
!                                call ugetio(3,nin,nout)                       
!                see the ugetio document for more details.                     
!                                                                              
!   copyright           - 1978 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      subroutine uertst (ier,name)                                             
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      dimension name(3),namset(3),nameq(3)                                     
      character*2        namsez(3),namez(3),iez                                
!                                                                              
      equivalence (namset,namsez),(nameq,namez),(ieq,iez)                      
!                                                                              
!      data               namsez/2hue,2hrs,2het/                               
      data               namsez/'ue','rs','et'/                                
      data               namez/'  ','  ','  '/                                 
!                                  first executable statement                  
      data               level/4/,ieqdf/0/,iez/'='/                            
      if (ier.gt.999) go to 25                                                 
      if (ier.lt.-32) go to 55                                                 
      if (ier.le.128) go to 5                                                  
      if (level.lt.1) go to 30                                                 
!                                  print terminal message                      
      call ugetio(1,nin,iounit)                                                
      if (ieqdf.eq.1) write(iounit,35) ier,nameq,ieq,name                      
      if (ieqdf.eq.0) write(iounit,35) ier,name                                
      go to 30                                                                 
    5 if (ier.le.64) go to 10                                                  
      if (level.lt.2) go to 30                                                 
!                                  print warning with fix message              
      call ugetio(1,nin,iounit)                                                
      if (ieqdf.eq.1) write(iounit,40) ier,nameq,ieq,name                      
      if (ieqdf.eq.0) write(iounit,40) ier,name                                
      go to 30                                                                 
   10 if (ier.le.32) go to 15                                                  
!                                  print warning message                       
      if (level.lt.3) go to 30                                                 
      call ugetio(1,nin,iounit)                                                
      if (ieqdf.eq.1) write(iounit,45) ier,nameq,ieq,name                      
      if (ieqdf.eq.0) write(iounit,45) ier,name                                
      go to 30                                                                 
   15 continue                                                                 
!                                  check for uerset call                       
      do 20 i=1,3                                                              
         if (name(i).ne.namset(i)) go to 25                                    
   20 continue                                                                 
      levold = level                                                           
      level = ier                                                              
      ier = levold                                                             
      if (level.lt.0) level = 4                                                
      if (level.gt.4) level = 4                                                
      go to 30                                                                 
   25 continue                                                                 
      if (level.lt.4) go to 30                                                 
!                                  print non-defined message                   
      call ugetio(1,nin,iounit)                                                
      if (ieqdf.eq.1) write(iounit,50) ier,nameq,ieq,name                      
!      if (ieqdf.eq.0) write(iounit,50) ier,name                                
   30 ieqdf = 0                                                                
      return                                                                   
   35 format(19h *** terminal error,10x,7h(ier = ,i3, 20h) from imsl routine ,3a2,a1,3a2)                               
   40 format(36h *** warning with fix error  (ier = ,i3, 20h) from imsl routine ,3a2,a1,3a2)                               
   45 format(18h *** warning error,11x,7h(ier = ,i3, 20h) from imsl routine ,3a2,a1,3a2)                               
   50 format(20h *** undefined error,9x,7h(ier = ,i5, 20h) from imsl routine ,3a2,a1,3a2)                               
!                                  save p for p = r case                       
!                                    p is the page name                        
!                                    r is the routine name                     
   55 ieqdf = 1                                                                
      do 60 i=1,3                                                              
   60 nameq(i) = name(i)                                                       
   65 return                                                                   
      end                                                                      
                                                                               
!   imsl routine name   - ugetio                                               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
!   computer            - ibm/single                                           
!                                                                              
!   latest revision     - january 1, 1978                                      
!                                                                              
!   purpose             - to retrieve current values and to set new            
!                           values for input and output unit                   
!                           identifiers.                                       
!                                                                              
!   usage               - call ugetio(iopt,nin,nout)                           
!                                                                              
!   arguments    iopt   - option parameter. (input)                            
!                           if iopt=1, the current input and output            
!                           unit identifier values are returned in nin         
!                           and nout, respectively.                            
!                           if iopt=2 (3) the internal value of                
!                           nin (nout) is reset for subsequent use.            
!                nin    - input unit identifier.                               
!                           output if iopt=1, input if iopt=2.                 
!                nout   - output unit identifier.                              
!                           output if iopt=1, input if iopt=3.                 
!                                                                              
!   precision/hardware  - single/all                                           
!                                                                              
!   reqd. imsl routines - none required                                        
!                                                                              
!   notation            - information on special notation and                  
!                           conventions is available in the manual             
!                           introduction or through imsl routine uhelp         
!                                                                              
!   remarks      each imsl routine that performs input and/or output           
!                operations calls ugetio to obtain the current unit            
!                identifier values. if ugetio is called with iopt=2 or 3       
!                new unit identifier values are established. subsequent        
!                input/output is performed on the new units.                   
!                                                                              
!   copyright           - 1978 by imsl, inc. all rights reserved.              
!                                                                              
!   warranty            - imsl warrants only that imsl testing has been        
!                           applied to this code. no other warranty,           
!                           expressed or implied, is applicable.               
!                                                                              
!-----------------------------------------------------------------------       
!                                                                              
      subroutine ugetio(iopt,nin,nout)                                         
!                                  specifications for arguments                
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      data               nind/5/,noutd/6/                                      
!                                  first executable statement                  
      if (iopt.eq.3) go to 10                                                  
      if (iopt.eq.2) go to 5                                                   
      if (iopt.ne.1) go to 9005                                                
      nin = nind                                                               
      nout = noutd                                                             
      go to 9005                                                               
    5 nind = nin                                                               
      go to 9005                                                               
   10 noutd = nout                                                             
 9005 return                                                                   
      end                                                                      
                                                                               
                                                                               
        SUBROUTINE DINCOG(A,X,GIN,GIM,GIP)                                     
!                                                                              
!       ===================================================                    
!       Purpose: Compute the incomplete gamma function                         
!                r(a,x), â(a,x) and P(a,x)                                     
!       Input :  a   --- Parameter ( a ó 170 )                                 
!                x   --- Argument                                              
!       Output:  GIN --- r(a,x)                                                
!                GIM --- â(a,x)                                                
!                GIP --- P(a,x)                                                
!       Routine called: GAMMA for computing â(x)                               
!       ===================================================                    
!                                                                              
        implicit real*8 (a-h,o-z)                                              
        implicit integer*4 (i-n)                                               
                                                                               
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   
        XAM=-X+A*DLOG(X)                                                       
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN                                   
           WRITE(*,*)'a and/or x too large'                                    
           STOP                                                                
        ENDIF                                                                  
        IF (X.EQ.0.0) THEN                                                     
           GIN=0.0                                                             
           CALL GAMMAJIN(A,GA)                                                 
           GIM=GA                                                              
           GIP=0.0                                                             
        ELSE IF (X.LE.1.0+A) THEN                                              
           S=1.0D0/A                                                           
           R=S                                                                 
           DO K=1,60                                                           
              R=R*X/(A+K)                                                      
              S=S+R                                                            
              IF (DABS(R/S).LT.1.0D-15) GO TO 15                               
           enddo                                                               
15         GIN=DEXP(XAM)*S                                                     
           CALL GAMMAJIN(A,GA)                                                 
           GIP=GIN/GA                                                          
           GIM=GA-GIN                                                          
        ELSE IF (X.GT.1.0+A) THEN                                              
           T0=0.0D0                                                            
           DO K=60,1,-1                                                        
              T0=(K-A)/(1.0D0+K/(X+T0))                                        
           enddo                                                               
           GIM=DEXP(XAM)/(X+T0)                                                
           CALL GAMMAJIN(A,GA)                                                 
           GIN=GA-GIM                                                          
           GIP=1.0D0-GIM/GA                                                    
        ENDIF                                                                  
        END                                                                    
                                                                               
                                                                               
        SUBROUTINE GAMMAJIN(X,GA)                                              
!                                                                              
!       ==================================================                     
!       Purpose: Compute gamma function â(x)                                   
!       Input :  x  --- Argument of â(x)                                       
!                       ( x is not equal to 0,-1,-2,úúú)                       
!       Output:  GA --- â(x)                                                   
!       ==================================================                     
!                                                                              
        implicit real*8 (a-h,o-z)                                              
        implicit integer*4 (i-n)                                               
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   
        DIMENSION G(26)                                                        
        PIg=3.141592653589793D0                                                 
        IF (X.EQ.INT(X)) THEN                                                  
           IF (X.GT.0.0D0) THEN                                                
              GA=1.0D0                                                         
              M1=X-1                                                           
              DO 10 K=2,M1                                                     
10               GA=GA*K                                                       
           ELSE                                                                
              GA=1.0D+300                                                      
           ENDIF                                                               
        ELSE                                                                   
           IF (DABS(X).GT.1.0D0) THEN                                          
              Z=DABS(X)                                                        
              M=INT(Z)                                                         
              R=1.0D0                                                          
              DO 15 K=1,M                                                      
15               R=R*(Z-K)                                                     
              Z=Z-M                                                            
           ELSE                                                                
              Z=X                                                              
           ENDIF                                                               
           DATA G/1.0D0,0.5772156649015329D0,-0.6558780715202538D0, -0.420026350340952D-1,0.1665386113822915D0,-.421977345555443D-1,-.96219715278770D-2, .72189432466630D-2,-.11651675918591D-2, -.2152416741149D-3,.1280502823882D-3, -.201348547807D-4,-.12504934821D-5, .11330272320D-5,-.2056338417D-6, .61160950D-8,.50020075D-8, -.11812746D-8,.1043427D-9, .77823D-11,-.36968D-11, .51D-12,-.206D-13, -.54D-14, .14D-14, .1D-15/                          
           GR=G(26)                                                            
           DO 20 K=25,1,-1                                                     
20            GR=GR*Z+G(K)                                                     
           GA=1.0D0/(GR*Z)                                                     
           IF (DABS(X).GT.1.0D0) THEN                                          
              GA=GA*R                                                          
              IF (X.LT.0.0D0) GA=-PIg/(X*GA*DSIN(PIg*X))                         
           ENDIF                                                               
        ENDIF                                                                  
        RETURN                                                                 
        END                                                                    
                                                                               
                                                                               
!                                                                              
!     determine infinity as the radius where sigma_z/v_z leq 0.1               
!                                                                              
      subroutine rinf(t)                                                       
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      dimension rjl(1), yrjl(1)                                                
      parameter (pig=3.1415926535897932d0)                                     
      include 'paramcomm.i'                                                    
      include 'sr.i'   
      include 'datagasgal.i'                                                        
!      external gammarec                                                        
                                                                               
      do j=1,1000                                                              
         t=0.01d0*dfloat(j)+10.d0                                              
!                                                                              
!     nu(t)                                                                    
!                                                                              
         if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then                    
            c=r200/rc                                                          
            gc=1./(dlog(c+1)-c/(c+1))                                          
            xnu=1.d0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./pig/(r200*r200*r200)                                              
         elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then                
            xnu=rc/(2.*pig*t)/(t+rc)**3.                                       
         else                                                                  
            xnu=-1.d0/dsqrt(pig*1.d0)*gammarec(-al+0.5d0)/gammarec(-al+1.d0)*al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)                       
         endif                                                                 
!                                                                              
!     Interpolate to get the sigma_r(t)                                        
!                                                                              
         rjl(1)=dlog(t)                                                        
         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)                
      if (ier .eq. 34) then                                                    
         print *,' ICSEVU in RINF: max xris=', xris(ninterp), ' r=',rjl        
      endif                                                                    
         sigr=dexp(yrjl(1))                                                    
!                                                                              
!     beta(r)                                                                  
!                                                                              
!     Convert from cbe=beta' to bec=beta                                       
!     in the case of constant anisotropy                                       
!                                                                              
         if (kani.eq.0) then                                                   
            bec=1.-1./(cbe*cbe)                                                
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values               
!                                                                              
!     Osipkov Merritt                                                          
!                                                                              
         elseif (kani.eq.2) then                                               
            bec=t*t/(t*t+cbe*cbe)                                              
!                                                                              
!     Mamon Lokas                                                              
!                                                                              
         elseif (kani.eq.1) then                                               
            bec=0.5*t/(t+cbe)                                                  
!                                                                              
!     Simplified Wojtak                                                        
!                                                                              
         elseif (kani.eq.3) then                                               
            bec=t**1.5/(t**1.5+cbe**1.5)                                       
!                                                                              
!     Simplified Tiret or modified ML                                          
!                                                                              
         elseif (kani.eq.4) then                                               
            rm2=rs                     ! SIS or other unspecified mass profiles
            if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist                               
            if (kmp.eq.4) rm2=1.521*rs ! Burkert                               
            bec=(1.-1./(cbe*cbe))*t/(t+rm2)                                    
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values               
!                                                                              
!     modified Tiret (non-zero central anisotropy)                             
!                                                                              
         elseif (kani.eq.5) then                                               
            rm2=rs                     ! SIS or other unspecified mass profiles
            if (kmp.eq.1) rm2=rs*(2.-gam) ! gNFW
            if (kmp.eq.2) rm2=0.5*rs ! Hernquist                               
            if (kmp.eq.4) rm2=1.521*rs ! Burkert                               
            bec=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)                              
            if (bec.ge.1.) bec=0.999d0 ! avoid unphysical values               
            if (bec.le.-1.) bec=-0.999d0 ! avoid unphysical values             
!                                                                              
!     Hansen+Moore relation for NFW mass profile                               
!                                                                              
         else                                                                  
            rhm=rc ! use nu(r)                                                 
!c         rhm=rs ! use rho(r)                                                 
            bec=ahm-bhm*(rhm+3.*t)/(rhm+t)                                     
         endif                                                                 
!                                                                              
!     sigma_z(R,r)                                                             
!                                                                              
         sigmaz=sigr*dsqrt(1.-bec*xmin*xmin/(t*t))                             
!                                                                              
!     Ratio sigma_z/v_z                                                        
!                                                                              
         rat=dabs(sigmaz/vj)                                                   
         if (rat.le.0.1) return                                                
      enddo                                                                    
                                                                               
      return                                                                   
      end                                                                      
!                                                                              
!     N(R) - beta-function                                                     
!                                                                              
      function sigmar3(tt)                                                     
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
      include 'paramcomm.i'                                                    
      include 'datagasgal.i'
      t=tt/rc                                                                  
!cc      sigmar3=(1.d0+t*t)**al/(r200*r200)                                    
      sigmar3=(1.d0+t*t)**al                                                   
      return                                                                   
      end                                                                      
                                                                               
!                                                                              
        SUBROUTINE HYGFX(A,B,C,X,HF)                                           
!                                                                              
!       ====================================================                   
!       Purpose: Compute hypergeometric function F(a,b,c,x)                    
!       Input :  a --- Parameter                                               
!                b --- Parameter                                               
!                c --- Parameter, c <> 0,-1,-2,...                             
!                x --- Argument   ( x < 1 )                                    
! Gary Mamon: a b and c must be specified as separate variables                
!   in the calling function!                                                   
!       Output:  HF --- F(a,b,c,x)                                             
!       Routines called:                                                       
!            (1) r8_GAMMA for computing gamma function                         
!            (2) r8_PSI for computing psi function                             
!       ====================================================                   
!                                                                              
        implicit real*8 (a-h,o-z)                                              
        implicit integer*4 (i-n)                                               
                                                                               
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   
        LOGICAL L0,L1,L2,L3,L4,L5                                              
        PIg=3.141592653589793D0                                                 
        EL=.5772156649015329D0                                                 
        L0=C.EQ.INT(C).AND.C.LT.0.0                                            
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0                                 
        L2=A.EQ.INT(A).AND.A.LT.0.0                                            
        L3=B.EQ.INT(B).AND.B.LT.0.0                                            
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0                                      
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0                                      
        IF (L0.OR.L1) THEN                                                     
           WRITE(*,*)'The hypergeometric series is divergent'                  
           RETURN                                                              
        ENDIF                                                                  
        EPS=1.0D-15                                                            
        IF (X.GT.0.95) EPS=1.0D-8                                              
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN                             
           HF=1.0D0                                                            
           RETURN                                                              
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN                         
           gc=r8_gamma(c)                                                      
           gcab=r8_gamma(c-a-b)                                                
           gca=r8_gamma(c-a)                                                   
           gcb=r8_gamma(c-b)                                                   
!$$$           CALL GAMMA(C,GC)                                                
!$$$           CALL GAMMA(C-A-B,GCAB)                                          
!$$$           CALL GAMMA(C-A,GCA)                                             
!$$$           CALL GAMMA(C-B,GCB)                                             
           HF=GC*GCAB/(GCA*GCB)                                                
           RETURN                                                              
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN               
           G0=DSQRT(PIg)*2.0D0**(-A)                                            
           g1=r8_gamma(c)                                                      
           g2=r8_gamma(1.0D0+A/2.0-B)                                          
           g3=r8_gamma(0.5D0+0.5*A)                                            
!$$$           CALL GAMMA(C,G1)                                                
!$$$           CALL GAMMA(1.0D0+A/2.0-B,G2)                                    
!$$$           CALL GAMMA(0.5D0+0.5*A,G3)                                      
           HF=G0*G1/(G2*G3)                                                    
           RETURN                                                              
        ELSE IF (L2.OR.L3) THEN                                                
           IF (L2) NM=INT(ABS(A))                                              
           IF (L3) NM=INT(ABS(B))                                              
           HF=1.0D0                                                            
           R=1.0D0                                                             
           DO 10 K=1,NM                                                        
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X                    
10            HF=HF+R                                                          
           RETURN                                                              
        ELSE IF (L4.OR.L5) THEN                                                
           IF (L4) NM=INT(ABS(C-A))                                            
           IF (L5) NM=INT(ABS(C-B))                                            
           HF=1.0D0                                                            
           R=1.0D0                                                             
           DO 15 K=1,NM                                                        
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X                
15            HF=HF+R                                                          
           HF=(1.0D0-X)**(C-A-B)*HF                                            
           RETURN                                                              
        ENDIF                                                                  
        AA=A                                                                   
        BB=B                                                                   
        X1=X                                                                   
        IF (X.LT.0.0D0) THEN                                                   
           X=X/(X-1.0D0)                                                       
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN                            
              A=BB                                                             
              B=AA                                                             
           ENDIF                                                               
           B=C-B                                                               
        ENDIF                                                                  
        IF (X.GE.0.75D0) THEN                                                  
           GM=0.0D0                                                            
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN                         
              M=INT(C-A-B)                                                     
              ga=r8_gamma(a)                                                   
              gb=r8_gamma(b)                                                   
              gc=r8_gamma(c)                                                   
              gam=r8_gamma(a+m)                                                
              gbm=r8_gamma(b+m)                                                
!$$$              CALL GAMMA(A,GA)                                             
!$$$              CALL GAMMA(B,GB)                                             
!$$$              CALL GAMMA(C,GC)                                             
!$$$              CALL GAMMA(A+M,GAM)                                          
!$$$              CALL GAMMA(B+M,GBM)                                          
              pa=r8_psi(a)                                                     
              pb=r8_psi(b)                                                     
!$$$              CALL PSI(A,PA)                                               
!$$$              CALL PSI(B,PB)                                               
              IF (M.NE.0) GM=1.0D0                                             
              DO 30 J=1,ABS(M)-1                                               
30               GM=GM*J                                                       
              RM=1.0D0                                                         
              DO 35 J=1,ABS(M)                                                 
35               RM=RM*J                                                       
              F0=1.0D0                                                         
              R0=1.0D0                                                         
              R1=1.0D0                                                         
              SP0=0.D0                                                         
              SP1=0.0D0                                                         
              IF (M.GE.0) THEN                                                 
                 C0=GM*GC/(GAM*GBM)                                            
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)                                
                 DO 40 K=1,M-1                                                 
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)              
40                  F0=F0+R0                                                   
                 DO 45 K=1,M                                                   
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K                
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)                           
                 DO 55 K=1,250                                                 
                    SP1=SP1+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))        
                    SM=0.0D0                                                   
                    DO 50 J=1,M                                                
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/(B+J+K-1.0)                                          
                    RP=PA+PB+2.0D0*EL+SP1+SM+DLOG(1.0D0-X)                      
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)          
                    F1=F1+R1*RP                                                
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60                  
55                  HW=F1                                                      
60               HF=F0*C0+F1*C1                                                
              ELSE IF (M.LT.0) THEN                                            
                 M=-M                                                          
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)                                 
                 C1=-(-1)**M*GC/(GAM*GBM*RM)                                   
                 DO 65 K=1,M-1                                                 
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)          
65                  F0=F0+R0                                                   
                 DO 70 K=1,M                                                   
70                  SP0=SP0+1.0D0/K                                            
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)                           
                 DO 80 K=1,250                                                 
                    SP1=SP1+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))        
                    SM=0.0D0                                                   
                    DO 75 J=1,M                                                
75                     SM=SM+1.0D0/(J+K)                                       
                    RP=PA+PB+2.0D0*EL+SP1-SM+DLOG(1.0D0-X)                      
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)              
                    F1=F1+R1*RP                                                
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85                  
80                  HW=F1                                                      
85               HF=F0*C0+F1*C1                                                
              ENDIF                                                            
           ELSE                                                                
              ga=r8_gamma(a)                                                   
              gb=r8_gamma(b)                                                   
              gc=r8_gamma(c)                                                   
              gca=r8_gamma(c-a)                                                
              gcb=r8_gamma(c-b)                                                
              gcab=r8_gamma(c-a-b)                                             
              gabc=r8_gamma(a+b-c)                                             
                                                                               
!$$$              CALL GAMMA(A,GA)                                             
!$$$              CALL GAMMA(B,GB)                                             
!$$$              CALL GAMMA(C,GC)                                             
!$$$              CALL GAMMA(C-A,GCA)                                          
!$$$              CALL GAMMA(C-B,GCB)                                          
!$$$              CALL GAMMA(C-A-B,GCAB)                                       
!$$$              CALL GAMMA(A+B-C,GABC)                                       
              C0=GC*GCAB/(GCA*GCB)                                             
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)                            
              HF=0.0D0                                                         
              R0=C0                                                            
              R1=C1                                                            
              DO 90 K=1,250                                                    
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)             
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))*(1.0-X)                                                   
                 HF=HF+R0+R1                                                   
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95                     
90               HW=HF                                                         
95            HF=HF+C0+C1                                                      
           ENDIF                                                               
        ELSE                                                                   
           A0=1.0D0                                                            
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.C.GT.B.AND.C.LT.2.0D0*B) THEN                                   
              A0=(1.0D0-X)**(C-A-B)                                            
              A=C-A                                                            
              B=C-B                                                            
           ENDIF                                                               
           HF=1.0D0                                                            
           R=1.0D0                                                             
           DO 100 K=1,250                                                      
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X                    
              HF=HF+R                                                          
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105                       
100           HW=HF                                                            
105        HF=A0*HF                                                            
        ENDIF                                                                  
        IF (X1.LT.0.0D0) THEN                                                  
           X=X1                                                                
           C0=1.0D0/(1.0D0-X)**AA                                              
           HF=C0*HF                                                            
        ENDIF                                                                  
        A=AA                                                                   
        B=BB                                                                   
!$$$        IF (K.GT.120) WRITE(*,115)                                         
!$$$115     FORMAT(1X,'Warning! You should check the accuracy')                
        Return                                                                 
        END                                                                    
                                                                               
                                                                               
      function r8_psi ( xx )                                                   
                                                                               
!*********************************************************************72       
!                                                                              
!c R8_PSI evaluates the function Psi(X).                                       
!                                                                              
!  Discussion:                                                                 
!                                                                              
!    This routine evaluates the logarithmic derivative of the                  
!    GAMMA function,                                                           
!                                                                              
!      PSI(X) = d/dX (GAMMA(X)) / GAMMA(X)                                     
!             = d/dX LN ( GAMMA(X) )                                           
!                                                                              
!    for real X, where either                                                  
!                                                                              
!      -XMAX1 < X < -XMIN  and X is not a negative integer),                   
!                                                                              
!    or                                                                        
!                                                                              
!      XMIN < X.                                                               
!                                                                              
!  Modified:                                                                   
!                                                                              
!    23 January 2008                                                           
!                                                                              
!  Author:                                                                     
!                                                                              
!    William Cody                                                              
!                                                                              
!  Reference:                                                                  
!                                                                              
!    William Cody, Anthony Strecok, Henry Thacher,                             
!    Chebyshev Approximations for the Psi Function,                            
!    Mathematics of Computation,                                               
!    Volume 27, Number 121, January 1973, pages 123-127.                       
!                                                                              
!  Parameters:                                                                 
!                                                                              
!    Input, double precision XX, the argument of the function.                 
!                                                                              
!    Output, double precision R8_PSI, the value of the function.               
!                                                                              
      implicit real*8 (a-h,o-z)                                                
      implicit integer*4 (i-n)                                                 
                                                                               
      dimension p1(9),p2(7),q1(8),q2(6)                                        
                                                                               
!                                                                              
!  Mathematical constants.  PIOV4 = pi / 4                                     
!                                                                              
      data zero /0.0d0/                                                        
      data fourth / 0.25d0/                                                    
      data half / 0.5d0 /                                                      
      data one / 1.0d0 /                                                       
      data three /3.0d0/                                                       
      data four /4.0d0/                                                        
      data piov4 /7.8539816339744830962d-01/                                   
!                                                                              
!  Machine-dependent constants                                                 
!                                                                              
      data xinf /1.70d+38/                                                     
      data xmin1 /5.89d-39/                                                    
      data xmax1 /3.60d+16/                                                    
      data xsmall /2.05d-09/                                                   
      data xlarge /2.04d+15/                                                   
!                                                                              
!  Zero of psi(x)                                                              
!                                                                              
      data x01 /187.0d0/                                                       
      data x01d /128.0d0/                                                      
      data x02 /6.9464496836234126266d-04/                                     
!                                                                              
!  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]           
!                                                                              
      data p1/4.5104681245762934160d-03,5.4932855833000385356d+00,3.7646693175929276856d+02,7.9525490849151998065d+03,7.1451595818951933210d+04,3.0655976301987365674d+05,6.3606997788964458797d+05,5.8041312783537569993d+05,1.6585695029761022321d+05/                                       
      data q1/9.6141654774222358525d+01,2.6287715790581193330d+03,2.9862497022250277920d+04,1.6206566091533671639d+05,4.3487880712768329037d+05,5.4256384537269993733d+05,2.4242185002017985252d+05,6.4155223783576225996d-08/             
!                                                                              
!  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)                  
!  for 3.0 < x.                                                                
!                                                                              
      data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01,-1.9784554148719218667d+01,-8.8100958828312219821d+00,-1.4479614616899842986d+00,-7.3689600332394549911d-02,-6.5135387732718171306d-21/                                      
      data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02,2.4736979003315290057d+02, 1.0742543875702278326d+02, 1.7463965060678569906d+01, 8.8427520398873480342d-01/           
                                                                               
      x = xx                                                                   
      w = abs ( x )                                                            
      aug = zero                                                               
!                                                                              
!  Check for valid arguments, then branch to appropriate algorithm.            
!                                                                              
      if ( -x .ge. xmax1 .or. w .lt. xmin1 ) then                              
        r8_psi = xinf                                                          
        if ( zero .lt. x ) then                                                
          r8_psi = -xinf                                                       
        end if                                                                 
        return                                                                 
      end if                                                                   
                                                                               
      if ( x .ge. half ) then                                                  
        go to 200                                                              
!                                                                              
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)         
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.                    
!                                                                              
      else if ( w .le. xsmall ) then                                           
        aug = -one / x                                                         
        go to 150                                                              
      end if                                                                   
!                                                                              
!  Argument reduction for cotangent.                                           
!                                                                              
  100 continue                                                                 
                                                                               
      if ( x .lt. zero ) then                                                  
        sgn = piov4                                                            
      else                                                                     
        sgn = - piov4                                                          
      end if                                                                   
                                                                               
      w = w - aint ( w )                                                       
      nq = int ( w * four )                                                    
      w = four * ( w - dble ( nq ) * fourth )                                  
!                                                                              
!  W is now related to the fractional part of 4.0 * X.                         
!  Adjust argument to correspond to values in the first                        
!  quadrant and determine the sign.                                            
!                                                                              
      n = nq / 2                                                               
                                                                               
      if ( n + n .ne. nq ) then                                                
        w = one - w                                                            
      end if                                                                   
                                                                               
      z = piov4 * w                                                            
                                                                               
      if ( mod ( n, 2 ) .ne. 0 ) then                                          
        sgn = - sgn                                                            
      end if                                                                   
!                                                                              
!  Determine the final value for  -pi * cotan(pi*x).                           
!                                                                              
      n = ( nq + 1 ) / 2                                                       
      if ( mod ( n, 2 ) .eq. 0 ) then                                          
!                                                                              
!  Check for singularity.                                                      
!                                                                              
        if ( z .eq. zero ) then                                                
          r8_psi = xinf                                                        
          if ( zero .lt. x ) then                                              
            r8_psi = -xinf                                                     
          end if                                                               
          return                                                               
        end if                                                                 
                                                                               
        aug = sgn * ( four / tan ( z ) )                                       
                                                                               
      else                                                                     
        aug = sgn * ( four * tan ( z ) )                                       
      end if                                                                   
                                                                               
  150 continue                                                                 
                                                                               
      x = one - x                                                              
                                                                               
  200 continue                                                                 
!                                                                              
!  0.5 <= X <= 3.0.                                                            
!                                                                              
      if ( x .le. three ) then                                                 
                                                                               
        den = x                                                                
        upper = p1(1) * x                                                      
        do i = 1, 7                                                            
          den = ( den + q1(i) ) * x                                            
          upper = ( upper + p1(i+1) ) * x                                      
        end do                                                                 
        den = ( upper + p1(9) ) / ( den + q1(8) )                              
        x = ( x - x01 / x01d ) - x02                                           
        r8_psi = den * x + aug                                                 
        return                                                                 
                                                                               
      end if                                                                   
!                                                                              
!  3.0 < X.                                                                    
!                                                                              
      if ( x .lt. xlarge ) then                                                
        w = one / ( x * x )                                                    
        den = w                                                                
        upper = p2(1) * w                                                      
        do i = 1, 5                                                            
          den = ( den + q2(i) ) * w                                            
          upper = ( upper + p2(i+1) ) * w                                      
        end do                                                                 
        aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug             
      end if                                                                   
                                                                               
      r8_psi = aug + log ( x )                                                 
                                                                               
      return                                                                   
      end                                                                      
                                                                               
                                                                               
         function r8_gamma ( x )                                               
                                                                               
!*********************************************************************72       
!                                                                              
!c R8_GAMMA evaluates Gamma(X) for a real argument.                            
!                                                                              
!  Discussion:                                                                 
!                                                                              
!    This function was originally named DGAMMA.                                
!                                                                              
!    However, a number of Fortran compilers now include a library              
!    function of this name.  To avoid conflicts, this function was             
!    renamed R8_GAMMA.                                                         
!                                                                              
!    This routine calculates the GAMMA function for a real argument X.         
!    Computation is based on an algorithm outlined in reference 1.             
!    The program uses rational functions that approximate the GAMMA            
!    function to at least 20 significant decimal digits.  Coefficients         
!    for the approximation over the interval (1,2) are unpublished.            
!    Those for the approximation for 12 <= X are from reference 2.             
!                                                                              
!  Modified:                                                                   
!                                                                              
!    18 January 2008                                                           
!                                                                              
!  Author:                                                                     
!                                                                              
!    William Cody, Laura Stoltz                                                
!                                                                              
!  Reference:                                                                  
!                                                                              
!    William Cody,                                                             
!    An Overview of Software Development for Special Functions,                
!    in Numerical Analysis Dundee, 1975,                                       
!    edited by GA Watson,                                                      
!    Lecture Notes in Mathematics 506,                                         
!    Springer, 1976.                                                           
!                                                                              
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,                      
!    Charles Mesztenyi, John Rice, Henry Thatcher,                             
!    Christoph Witzgall,                                                       
!    Computer Approximations,                                                  
!    Wiley, 1968,                                                              
!    LC: QA297.C64.                                                            
!                                                                              
!  Parameters:                                                                 
!                                                                              
!    Input, double precision X, the argument of the function.                  
!                                                                              
!    Output, double precision R8_GAMMA, the value of the function.             
!                                                                              
         implicit real*8 (a-h,o-z)                                             
         implicit integer*4 (i-n)                                              
         dimension c(7),p(8),q(8)                                              
         logical parity                                                        
                                                                               
!                                                                              
!  Mathematical constants                                                      
!                                                                              
      data one /1.0D+00 /                                                      
      data half /0.5D+00/                                                      
      data twelve /12.0D+00/                                                   
      data two /2.0D+00 /                                                      
      data zero /0.0D+00/                                                      
      data sqrtpi /0.9189385332046727417803297D+00/                            
      data pig /3.1415926535897932384626434D+00/                                
!                                                                              
!  Machine dependent parameters                                                
!                                                                              
      data xbig / 171.624D+00 /                                                
      data xminin / 2.23D-308 /                                                
      data eps /2.22D-16/                                                      
      data xinf /1.79D+308/                                                    
!                                                                              
!  Numerator and denominator coefficients for rational minimax                 
!  approximation over (1,2).                                                   
!                                                                              
      data p/-1.71618513886549492533811d+00,2.47656508055759199108314d+01,-3.79804256470945635097577d+02,6.29331155312818442661052d+02,8.66966202790413211295064d+02,-3.14512729688483675254357d+04,-3.61444134186911729807069d+04,6.64561438202405440627855d+04/                                         
                                                                               
      data q/-3.08402300119738975254353D+01,3.15350626979604161529144D+02,-1.01515636749021914166146D+03,-3.10777167157231109440444D+03,2.25381184209801510330112D+04,4.75584627752788110767815D+03,-1.34659959864969306392456D+05,-1.15132259675553483497211D+05/                                         
!                                                                              
!  Coefficients for minimax approximation over (12, INF).                      
!                                                                              
      data c/-1.910444077728D-03,8.4171387781295D-04,-5.952379913043012D-04,7.93650793500350248D-04,-2.777777777777681622553D-03,8.333333333333333331554247D-02,5.7083835261D-03/                                                      
                                                                               
      parity = .false.                                                         
      fact = one                                                               
      n = 0                                                                    
      y = x                                                                    
!                                                                              
!  Argument is negative.                                                       
!                                                                              
      if ( y .le. zero ) then                                                  
                                                                               
        y = - x                                                                
        y1 = aint ( y )                                                        
        res = y - y1                                                           
                                                                               
        if ( res .ne. zero ) then                                              
                                                                               
          if ( y1 .ne. aint ( y1 * half ) * two ) then                         
            parity = .true.                                                    
          end if                                                               
                                                                               
          fact = - pig / sin ( pig * res )                                       
          y = y + one                                                          
                                                                               
        else                                                                   
                                                                               
          res = xinf                                                           
          r8_gamma = res                                                       
          return                                                               
                                                                               
        end if                                                                 
                                                                               
      end if                                                                   
!                                                                              
!  Argument is positive.                                                       
!                                                                              
      if ( y .lt. eps ) then                                                   
!                                                                              
!  Argument < EPS.                                                             
!                                                                              
        if ( xminin .le. y ) then                                              
          res = one / y                                                        
        else                                                                   
          res = xinf                                                           
          r8_gamma = res                                                       
          return                                                               
        end if                                                                 
                                                                               
      else if ( y .lt. twelve ) then                                           
                                                                               
        y1 = y                                                                 
!                                                                              
!  0.0 < argument < 1.0.                                                       
!                                                                              
        if ( y .lt. one ) then                                                 
                                                                               
          z = y                                                                
          y = y + one                                                          
!                                                                              
!  1.0 < argument < 12.0.                                                      
!  Reduce argument if necessary.                                               
!                                                                              
        else                                                                   
                                                                               
          n = int ( y ) - 1                                                    
          y = y - dble ( n )                                                   
          z = y - one                                                          
                                                                               
        end if                                                                 
!                                                                              
!  Evaluate approximation for 1.0 < argument < 2.0.                            
!                                                                              
        xnum = zero                                                            
        xden = one                                                             
        do i = 1, 8                                                            
          xnum = ( xnum + p(i) ) * z                                           
          xden = xden * z + q(i)                                               
        end do                                                                 
                                                                               
        res = xnum / xden + one                                                
!                                                                              
!  Adjust result for case  0.0 < argument < 1.0.                               
!                                                                              
        if ( y1 .lt. y ) then                                                  
                                                                               
          res = res / y1                                                       
!                                                                              
!  Adjust result for case 2.0 < argument < 12.0.                               
!                                                                              
        else if ( y .lt. y1 ) then                                             
                                                                               
          do i = 1, n                                                          
            res = res * y                                                      
            y = y + one                                                        
          end do                                                               
                                                                               
        end if                                                                 
                                                                               
      else                                                                     
!                                                                              
!  Evaluate for 12.0 <= argument.                                              
!                                                                              
        if ( y .le. xbig ) then                                                
                                                                               
          ysq = y * y                                                          
          sum = c(7)                                                           
          do i = 1, 6                                                          
            sum = sum / ysq + c(i)                                             
          end do                                                               
          sum = sum / y - y + sqrtpi                                           
          sum = sum + ( y - half ) * log ( y )                                 
          res = exp ( sum )                                                    
                                                                               
        else                                                                   
                                                                               
          res = xinf                                                           
          r8_gamma = res                                                       
          return                                                               
                                                                               
        end if                                                                 
                                                                               
      end if                                                                   
!                                                                              
!  Final adjustments and return.                                               
!                                                                              
      if ( parity ) then                                                       
        res = - res                                                            
      end if                                                                   
                                                                               
      if ( fact .ne. one ) then                                                
        res = fact / res                                                       
      end if                                                                   
                                                                               
      r8_gamma = res                                                           
                                                                               
      return                                                                   
      end                                                                      





!####Mio:calcolo likelihood e metto prior
    function Cl_Cash(this,CMB,Theory,DataParams)
    ! likelihood computation
    ! SZ nuisance in dataparams
    use, intrinsic :: ieee_arithmetic
    Class(ClLikelihood) :: this    !####Mio:parametri di likelihood
    Class (CMBParams):: CMB !####Mio:parametri di CMB
    Class(TCosmoTheoryPredictions), target :: Theory  !####Mio:parametri di Teoria
    real(mcp) DataParams(:)   
    real(mcp)  Cl_Cash
    REAL(DL) :: sum,sumvdp
    integer*4  npar,nparvdp
    real*8  parlik(6),f,fvdp,parvdp(6)

    Cl_Cash=logzero

!####Mio: do il valore dei parametri di CMB oppure Theory negli 'array' che entra nel modulo cosmology 
    !Mapping of nuisance parameters
  
    cosmopar%r200=DataParams(1)
    cosmopar%rjaf=DataParams(2)
    cosmopar%rs=DataParams(3)
    cosmopar%cbe=DataParams(4)
    cosmopar%gam=DataParams(5)
    cosmopar%xmstarl=DataParams(6)
    cosmopar%rabcg=DataParams(7)

    npar=6

    parlik(1)=dlog10(cosmopar%r200)
    parlik(2)=dlog10(cosmopar%rjaf)
    parlik(3)=dlog10(cosmopar%rs)
    parlik(4)=dlog10(cosmopar%cbe)
    parlik(5)=dlog10(cosmopar%gam)
    parlik(6)=dlog10(cosmopar%xmstarl)

    call mamlik(npar,parlik,f)
    sum=f

    nparvdp=6

    parvdp(1)=cosmopar%r200
    parvdp(2)=cosmopar%rs
    parvdp(3)=cosmopar%gam
    parvdp(4)=cosmopar%xmstarl
    parvdp(5)=cosmopar%rabcg
    parvdp(6)=cosmopar%rjaf

    call vdplik(nparvdp,parvdp,fvdp)
    sumvdp=fvdp

 !####Mio: qui dovra' arrivare la likelihood
    
    Cl_Cash=sum + sumvdp/2.
;    Cl_Cash=sum     ! neglect BCG VDP 

    ! Print*,'Param mam',10**(parlik)
    ! Print*,'Param vdp',parvdp
    Print*,'CL lnlike = ',Cl_Cash,sum,sumvdp,parvdp

    end function Cl_Cash


    end module clcounts

