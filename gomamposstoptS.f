      program gomamposstoptS

c     Version of 23/12/2025 - main features:
c     - modified gravity/DE models: Chameleon screening, Vainsthein 
c       screening, Refracted Gravity, Clustering Dark Energy
c     - for CS, possibility to run with mass models: NFW, Einasto,
c       Burkert, bNFW (it includes Hernquist), Isothermal and multi-component
c     - BCG dataset added. Maximum of 11 free parameters 
c     - gT gOM, and BP anisotropy models included
c     - flat prior and gaussian prior included, with the possibility of 
c       keeping flat priors on beta0 and betaInf while the others are
c       Gaussian
c    -  for gT, gOM, BP rbeta can be kept as a free parameter or fixed as
c       the r-2 of the mass profile (or dark matter profile for multi-component
c
c    - Gas and galaxies mass profiles can be passed in a parametric way,
c      along with their (gaussian) uncertainties. The parameters are treated 
c      as nuisance parameters in the MCMC module and optimized. They are saved
c      in an additional file called GasGal_mass.txt

c    - The Options file is deprecated. The grid/mcmc controls, as well as the
c      number of steps and lensing swithces are all in the input_pars.txt file

c    - Lensing chains added as a new file with several subcases depending on the
c      mass model used

c     ADDED THE gNFW PROFILE: option kmp=10
c     THE Exponent gamma is given by the dA3 parameter 

c     derived from gomamposst.f ...'opt' for optimization

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
c                 Aosta,   September 2019
c                 Aosta,   September 2021
c                 Aosta,   October 2021: removed copyrighted imsl routines;
c                                        substituted by John Burkardt's
c                                        routines (GNU public license) and
c                                        SLATEC routines by  Jones, R. E., (SNLA) 

c                 Aosta,   March 2022: few adjustments have been made,
c                                      cmake for installation, some bugs
c                                       fixed 

c     if needed (cmake doesn't work), use the script installation:
c     (change the permission of the script to make them executable)
c     ./script/script_Lib.sh
c     Compile with ./script/script_compile.sh
c     Execute as ./script/script_runmam.sh

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535d0,iga=25000,
     &     grav=4.302e-9,q0=-0.64,clight=2.99792458d5,npar=9)
      dimension di(iga),ve(iga),eve(iga),dns(iga),ra(iga),de(iga),
     &     rso(iga),vso(iga),bins(400),pars(27),   
     &     syb(1000),sybs(1000),vkb(1000),vkbs(1000),
     &     xx(iga),yy(iga),yboo(iga),wboo(iga),ww(iga),
     &     iw(iga), sigma(npar)
      character*75 fgroup,frnvn,fbfn,fbn,ffn,fmlv,fprop,line,
     &     fsb,fsf,fsft,fkb,fkf,fkft,fparmam, fbcg, fUme
      real*8 rjl(1), yrjl(1)
c   input files      
      character(len = 13) :: filename
      character(len=4) :: dat
      character(len=300) :: buffer, label
      character (len=300) :: massmod, animod, numod, amodbcg
      integer :: pos, pos1, pos2, posl1, posl2, posl, poscom
      integer, parameter :: fh = 15
      integer :: ios = 0
      integer :: lin = 0
      integer :: line2 = 0   !for reading Options 
      logical*4 :: read2=.False. !for r200
      logical*4 :: addrup=.False.
      logical*4 :: taddl=.False.
      logical*4 :: taddu=.False.
      logical*4 :: saddl=.False.
      logical*4 :: saddu=.False.
      
c     all file to be read by MAMPOSSt *************************
      include 'datarv.i'             !definition of data arrays
      include 'paramsoptS.i'         !all useful parameters
      include 'units.i'
      include 'sr.i'
      include 'vkr.i'
      include 'vlos.i'
      include 'tsallis.i'
      include 'datavdp.i'
      
      external fcn1,fcn2,fcn3,fa,fkii
      external sr2int,sr2out,sigmar1,sigmar2,sigmar3,gwenu
      external vr4nuint
      external sigmarnorm
      external gammarec
      external frLin
      external r_normal_ab 
      external fmg
      external f_nfw
      external gmg
      external g_nfw
      

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

c***********************************************************************
            
      write(*,*) ' '
      write(*,*) ' '
      write(*,*)    ' ******************************************'
      write(*,*)    ' ***                                    ***'
      if (ktsa.gt.0.5) then
         write(*,*) ' *** Welcome into Tsallis MG-MAMPOSSt!  ***'
      else
         write(*,*) ' *** Welcome into Gaussian MG-MAMPOSSt! ***'
      endif
      write(*,*)    ' ***                                    ***'
      write(*,*)    ' ******************************************'
      write(*,*) ' '
      write(*,*) ' '

      write(*,100)
 100  format(' Input files: '/ 
     &     ' (1) data (R [kpc], Vrf [km/s], errVrf [km/s], weights',/
     &     ' (2) input parameters for MAMPOSSt ')
      read(*,9) fgroup
      read(*,9) fparmam
 9    format(a)
      open(10,file=fgroup,status='old')   !data file 
      open(29,file=fparmam,status='old')  !parameter file (name is in gomamposst*.inp)

      write(*,196) fgroup
 196  format(//' Data-set is ',a50,//)

c     default parameters:
      nr200=0
      nrs=0
      nrc=0
      nbs=0
      ntmass=0
      nhone=0
      nb2=0
      nrbeta=0
      njaf=0
      nxmas=0
      nA3=0
      
      !Initial guess default values
      r200g=1.0d0
      r200=r200g
      rsg=1.0d0
      rs=rsg
      rcg=1.0d0
      rc=rcg
      cbeg=1.0d0
      cbe=cbeg
      rbetag=1.0
      rbeta=rbetag
      rjafg=0.01
      xmstarg=1.0
      
      dA3g=0.4
      
      cbe0g=1.2d0
      cbe0=cbe0g
      cb0low=0.0001d0
      cb0up=100.0d0
      rbetalow=0.0001d0
      rbetaup=10.0d0
      r2low = 0.001
      r2up = 100
      
      
      nweigh=0
      eim=3. !exponent for Einasto model
            
      kbcg=0
      kanbcg=0
      kdbcg=0 !not used anymore
      
c     baryonic component of gas: assumed to be a beta 
c     profile with  the following parameters:
      rho0gas = 1.0e13
      rsgas = 0.2d0
      
      srhoga = 0.0d0   !uncertainties in the gas parameters
      srsga = 0.0d0
      
C     baryonic components of galaxies: assumed to be a NFW with the 
C     following parameters:
      rhostar = 0.0e13
      rstar = 0.2d0
      
      srhost = 0.0d0  !uncertainties in the galaxy parameters
      srstar = 0.0d0

c     for beta =1 or beta = 2/3 in the beta model of the gas
      kalpha = 0 !Default: beta = 2/3

c
c     limiting value for the screening radius to be considered zero       
      alimi=0.0001
c     read file with parameters for MAMPOSSt
c
c     Exponential cutoff for the fifth force in Chameleon gravity?
      ncutof=1

c     Option used for the analysis of the FS2 survey on rc 
c     this sets rc = 4rs for Euclid clusters
      kfs2 = 1

c     Options to have flat priors in beta and beta0 while the other parameters 
c     have Gaussian priors. Default is no (0)  
      ipriorBeta0 = 0
      ipriorBeta = 0
	                  
c     Free rbeta for gT or gOM models
      
      ifreeBeta = 0
      
      !This is to add the probabiltiy distribution from lensing
      naddlens = 0
   
      !Default: nmcmc run 
      nmcmc = 1
      
      !Default: nmcmc run 
      Nsample = 100000
      
      !Default: no lensing
      nlens = 0
     
      nLogP=0 !Set equal to 0 to explore the phi space using the logarithm
               ! of the parameter. 
     
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

     
      do while (ios == 0)
         read(29, '(A)', iostat=ios) buffer
         if (ios == 0) then
            lin = lin + 1

        ! Find the first instance of whitespace.  Split label and data.
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
                    
            select case (label)
c*********** number of steps in the input parameters *******************            
            case ('nr200')  ! number of steps for r200 fit
                read2=.True.
                read(buffer, *, iostat=ios) nr200
                if (ios.eq.-1) then
                    ios=0
                    nr200=0
                endif
                pars(1)=nr200  
            case ('nrc')    ! number of steps for rc fit, scale radius of N(R)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces LfM, c_nu=c]
                            !   [if = -2 fit N(R) outside MAMPOSSt]
                read(buffer, *, iostat=ios) nrc
                if (ios.eq.-1) then
                    ios=0
                    nrc=0
                endif 
                pars(2)=nrc
            case ('nrs')    ! number of steps for rs fit, scale radius of M(r)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces MfL, c=c_nu]
                            !   [if = -2 forces LCDM, c=c(M)]
                read(buffer, *, iostat=ios) nrs
                if (ios.eq.-1) then
                    ios=0
                    nrs=0
                endif
                pars(3)=nrs 
            case ('nbeta')  ! number of steps for anisotropy parameter
                            !   [if = -1 forces a_ML=r_s]
                            !   [if = -2 forces Hansen+Moore]
                read(buffer, *, iostat=ios) nbs 
                if (ios.eq.-1) then
                    ios=0
                    nbs=0
                endif
                pars(4)=nbs 
            case ('nA1')    ! number of steps in the first MG parameter
                read(buffer, *, iostat=ios) ntmass
                if (ios.eq.-1) then
                    ios=0
                    ntmass=0
                endif 
                pars(5)=ntmass
            case ('nA2')    ! number of steps in the additional free parameter
                            ! could be the screening radius for f(R)-Hu Sawicki 
                            ! or Coupling constant Q for chameleon gravity
                read(buffer, *, iostat=ios) nhone
                if (ios.eq.-1) then
                    ios=0
                    nhone=0
                endif
                pars(6)=nhone
            case ('nA3')    ! number of steps in the additional free parameter
                            ! could be the exponent Q of Refracted gravity 
                            ! model or the gamma exponent in the gNFW case.
                read(buffer, *, iostat=ios) nA3
                if (ios.eq.-1) then
                    ios=0
                    nA3=0
                endif
 !               pars(10)=nhone
            case ('nbeta2')    ! number of steps in the additional free parameter
                            ! could be the screening radius for f(R)-Hu Sawicki 
                            ! or Coupling constant Q for chameleon gravity
                read(buffer, *, iostat=ios) nb2
                if (ios.eq.-1) then
                    ios=0
                    nb2=0
                endif
            case ('nrbeta')    ! number of steps in scale radius of the velocity
                               ! anisotropy profile 
                read(buffer, *, iostat=ios) nrbeta
                if (ios.eq.-1) then
                    ios=0
                    nrbeta=0
                endif


            case ('njaf')    ! number of steps in the jaffe radius BCG
                read(buffer, *, iostat=ios) njaf
                if (ios.eq.-1) then
                    ios=0
                    njaf=0
                endif

            case ('nxmas')    ! number of steps in the M/L BCG
                read(buffer, *, iostat=ios) nxmas
                if (ios.eq.-1) then
                    ios=0
                    nxmas=0
                endif

c***********************************************************************
c********** parameters input guess values ******************************
            case ('r200g')
                read(buffer, *, iostat=ios) r200g
                if (ios.eq.-1) then
                    ios=0
                    r200g=1.0d0 
                    stop('Error: input guess value required') 
                endif
                r200=r200g ! set r200 initial value to guess value
            case ('rcg')
                read(buffer, *, iostat=ios) rcg
                if (ios.eq.-1) then
                    ios=0
                    rcg=1.0d0
                    stop('Error: input guess value required') 
                endif
                rc=rcg ! set rc initial value to guess value
            case ('rsg')
                read(buffer, *, iostat=ios) rsg
                if (ios.eq.-1) then
                    ios=0
                    rsg=1.0d0
                    stop('Error: input guess value required') 
                endif
                rs=rsg ! set rs initial value to guess value
            case ('betag')
                read(buffer, *, iostat=ios) cbeg
                if (ios.eq.-1) then
                    ios=0
                    cbeg=1.0d0
                    stop('Error: input guess value required') 
                endif
                cbe=cbeg ! set cbe initial value to guess value
            case ('A1g')
                read(buffer, *, iostat=ios) tmassg
                if (ios.eq.-1) then
                    ios=0
                    tmassg=0.0d0
                    stop('Error: input guess value required') 
                endif
                
            case ('A2g')
                read(buffer, *, iostat=ios) screeg
                if (ios.eq.-1) then
                    ios=0
                    screeg=0.0d0
                    stop('Error: input guess value required') 
                endif
                
            case ('A3g')
                read(buffer, *, iostat=ios) dA3g
                if (ios.eq.-1) then
                    ios=0
                    dA3g=0.0d0
                    stop('Error: input guess value required') 
                endif                
            case ('beta2g')
                read(buffer, *, iostat=ios) cbe0g
                if (ios.eq.-1) then
                    ios=0
                    cbe0g=1.0d0
                    
                endif
                cbe0=cbe0g
            !Anisotropy scale radius (for BP model)
            case ('rbetag')
                read(buffer, *, iostat=ios) rbetag
                if (ios.eq.-1) then
                    ios=0
                    cbe0g=1.0d0
                    
                endif
                rbeta=rbetag
                                
            case ('rjafg')
                read(buffer, *, iostat=ios) rjafg
                if (ios.eq.-1) then
                    ios=0
                    rjafg=0.01d0
                    
                endif
                rjaf=rjafg     
                
            case ('xmstarg')
                read(buffer, *, iostat=ios) xmstarg
                if (ios.eq.-1) then
                    ios=0
                    xmstarg=1.0d0
                    
                endif
                xmstar=xmstarg 
c***********************************************************************
c********* Model Options ***********************************************

            case ('za')  ! average redshift of the cluster (needed to evaluate Hz)
                         ! since velocities are given as rest-frame in input file
                read(buffer, *, iostat=ios) za
                if (ios.eq.-1) then
                    ios=0
                    za=0.0
                endif
            case ('H0')  ! Hubble constant at z=0
                read(buffer, *, iostat=ios) h0
                if (ios.eq.-1) then
                    ios=0
                    h0=70.0
                endif
            case ('Olam')  ! Omega Lambda density parameter 
                read(buffer, *, iostat=ios) Olam
                if (ios.eq.-1) then
                    ios=0
                    Olam=0.7
                endif
            case ('Omegam')  ! Omega matter density parameter

                read(buffer, *, iostat=ios) Omegam
                if (ios.eq.-1) then
                    ios=0
                    Omegam=0.3
                endif
            case ('Rlow')  ! Inner radius for sample selection (Mpc)

                read(buffer, *, iostat=ios) rlowin
                if (ios.eq.-1) then
                    ios=0
                    rlowin=0.05
                endif
            case ('Rup')  ! Outer radius for sample selection (Mpc)

                read(buffer, *, iostat=ios) rupin
                if (ios.eq.-1) then
                    ios=0
                    rupin=2.0 
                    addrup=.True.
                    
                endif
            case ('N(R)')  ! model for number density profile N(R)
                           ! projected NFW / projected Hernquist / beta-model (1/2/3)

                read(buffer, *, iostat=ios) numod
                if (ios.eq.-1) then
                    ios=0
                    numod='pNFW'
                endif
                select case(numod)
                    case('pNFW')
                      knfit=1
                    case('pHer')
                      knfit=2
                    case('beta')
                      knfit=3
                    case default
                      print *, 'Invalid N(R), switching to default'
                      knfit=1
                end select
            case ('al')  ! N(R) exponent for beta-model (if knfit=3)

                read(buffer, *, iostat=ios) al
                if (ios.eq.-1) then
                    ios=0
                    al=-0.00
                endif
             
            case ('betaBCG')  !beta model for the bcg
             read(buffer, *, iostat=ios) amodbcg
                if (ios.eq.-1) then
                    ios=0
                    amodbcg='C'
                endif
                
                select case(amodbcg)
                    case('C')
                      kanbcg=0
                    case('OM')
                      kanbcg=1
                    case default
                      print *, 'Invalid betaBCG, switching to default'
                      kanbcg=0
                end select
                
            case ('rabcg')  ! N(R) exponent for beta-model (if knfit=3)

                read(buffer, *, iostat=ios) rabcg
                if (ios.eq.-1) then
                    ios=0
                    rabcg=20.0
                endif    
                             
            case ('xlumbcg')  ! N(R) exponent for beta-model (if knfit=3)

                read(buffer, *, iostat=ios) xlumbcg
                if (ios.eq.-1) then
                    ios=0
                    xlumbcg=1
                endif 
                !unit 1e11
                xlumbcg=xlumbcg*1e11
c                write(*,*) 'psdfopsodfpsodfpsdop', xlumbcg, rabcg

                                
            case ('M(r)')  ! rho(r) model: NFW/Hernquist/PIEMD/Burkert/SoftIS/Einasto_m=5/
                           ! mod_NFW linear f(R)/mod_NFW beyond Horndeski/mod_NFW general chameleon (1/2/3/4/5/6/7/8/9)

                read(buffer, *, iostat=ios) massmod
                if (ios.eq.-1) then
                    ios=0
                    massmod='NFW'
                endif 
                select case(massmod)
                    case('NFW')
                      kmp=1
                    case('Her')
                      kmp=2
                    case('PIEMD')
                      kmp=3
                    case('Bur')
                      kmp=4
                    case('SoftIS')
                      kmp=5
                    case('Eis')
                      kmp=6
                    !case('mNFW_LH')
                    !  kmp=7
                    case('mNFW_BH')
                      kmp=8   
                    case('mNFW_GC')
                      kmp=9
                    case('gNFW')
                      kmp=10    
                    case('mBur_GC')
                      kmp=11    
                    case('mEin_GC')
                      kmp=12
                    case('mIso_GC')
                      kmp=16                    
                    case('RG') !refracted gravity
                      kmp=13
                      kdbcg=0
                      kbary=1
                    case('BS') !boson star
                      kmp = 14
                      kbary = 1
                    case('mgNFW_GC') !modified gNFW -chameleon
                      kmp = 15
                    case('mTotal_GC') !multicomponent model - chameleon
                      kmp = 17
                    case('gNFW_CDE') !gNFW clustering DE
                      kmp = 18
                    case default
                      print *, 'Invalid mass mod, switching to default'
                      kmp=1
                end select
                pars(7)=kmp
            case ('Beta(r)')  ! Anisotropy model, beta'=constant, MamLok, OsiMer, 
                           !   simplified Wojtak, simplified Tiret, modified Tiret (0,1,2,3,4,5)
                read(buffer, *, iostat=ios) animod
                if (ios.eq.-1) then
                    ios=0
                    animod='C'
                endif 
                select case(animod)
                    case('C')
                      kani=0
                    case('ML')
                      kani=1
                    case('OM')
                      kani=2
                    case('gOM')
                      kani=21
                    case('WJ')
                      kani=3
                    case('T')
                      kani=4
                    case('gT')
                      kani=41
                    case('BP')
                      kani=8  !The Biviano&Pizzuti model
                    case('O')
                      kani=5                                                                                
                    case default
                      print *, 'Invalid number, switching to default'
                      kani=0
                end select
                pars(8)=kani
            case ('rcut')  ! PIEMD model rcut in Mpc (only used if kmp=3) 
                read(buffer, *, iostat=ios) rcut
                if (ios.eq.-1) then
                    ios=0
                    rcut=1.0
                endif
            case ('eim')  ! Einasto model exponent
                read(buffer, *, iostat=ios) eim
                if (ios.eq.-1) then
                    ios=0
                    eim=3.0
                endif
            case ('BCG')  !Include the BCG mass modeling?
                read(buffer, *, iostat=ios) kbcg
                if (ios.eq.-1) then
                    ios=0
                    kbcg=0
                endif
            case ('BCG_Data')  !Include the BCG dataset?
                read(buffer, *, iostat=ios) kbdt
                if (ios.eq.-1) then
                    ios = 0
                    kbdt = 0
                endif                
                
            case ('Baryons') !Include baryons? Only for BCG=1 or RG
                read(buffer, *, iostat=ios) kbary
                if (ios.eq.-1) then
                    ios=0
                    kbary=0
                endif
                
            case ('rho0gas')  !baryons critical density (beta  model)
                read(buffer, *, iostat=ios) rho0gas
                if (ios.eq.-1) then
                    ios=0
                    rho0gas=10.0d0
                endif
                rho0gas=rho0gas*1.0e13               

            case ('rsgas') !baryons scale radius (beta  model)
                read(buffer, *, iostat=ios) rsgas
                if (ios.eq.-1) then
                    ios=0
                    rsgas=0.12d0
                endif                           

            case ('srhoga') !baryons critical density error (beta  model)
                read(buffer, *, iostat=ios) srhoga
                if (ios.eq.-1) then
                    ios=0
                    srhoga=0.0d0
                endif 
                srhoga = srhoga*1.0e13

            case ('srsga') !baryons scale radius error (beta  model)
                read(buffer, *, iostat=ios) srsga
                if (ios.eq.-1) then
                    ios=0
                    srsga=0.0d0
                endif 
!   upgrade!
            case ('rhostar')  !baryons critical density (beta  model)
                read(buffer, *, iostat=ios) rhostar
                if (ios.eq.-1) then
                    ios=0
                    rhostar=10.0d0
                endif
                rhostar=rhostar*1.0e13               

            case ('rstar') !baryons scale radius (beta  model)
                read(buffer, *, iostat=ios) rstar
                if (ios.eq.-1) then
                    ios=0
                    rstar=0.12d0
                endif                           

            case ('srhostar') !baryons critical density error (beta  model)
                read(buffer, *, iostat=ios) srhost
                if (ios.eq.-1) then
                    ios=0
                    srhost=0.0d0
                endif 
                srhost = srhost*1.0e13

            case ('srstar') !baryons scale radius error (beta  model)
                read(buffer, *, iostat=ios) srstar
                if (ios.eq.-1) then
                    ios=0
                    srstar=0.0d0
                endif 
            case ('kbeta') ! exponent for beta model of gas: alpha = 2/3 or 1
                read(buffer, *, iostat=ios) kalpha
                if (ios.eq.-1) then
                    ios=0
                    kalpha=0
                endif 
	         case ('weights')    ! weights in the mass distribution
                read(buffer, *, iostat=ios) nweigh
                if (ios.eq.-1) then
                    ios=0
                    nweigh=0
                endif 	
            case ('FASTMODE')   
                read(buffer, *, iostat=ios) kbsp
                if (ios.eq.-1) then
                    ios=0
                    kbsp=0
                endif   
            case ('OPT')  ! optimization algorithm: 0/1/2=bobyqa/newuao/powell
                                            ! -1 skip optimization
                read(buffer, *, iostat=ios) kopt
                if (ios.eq.-1) then
                    ios=0
                    kopt=-1
                endif
            
            case ('FreeBeta')  ! optimization algorithm: 0/1/2=bobyqa/newuao/powell
                                            ! -1 skip optimization
                read(buffer, *, iostat=ios) ifreeBeta
                if (ios.eq.-1) then
                    ios=0
                    ifreeBeta=0
                endif
                
c   ******** now we read the prior type (flat or Gaussian) ************  
            case ('nmcmc')
                read(buffer, *, iostat=ios) nmcmc
                if (ios.eq.-1) then
                    ios=0
                    nmcmc=1
                endif  
            case ('Nsample')
                read(buffer, *, iostat=ios) Nsample
                if (ios.eq.-1) then
                    ios=0
                    Nsample=400000
                endif 
            
!              lensing                 
            case ('nlens')
                read(buffer, *, iostat=ios) nlens
                if (ios.eq.-1) then
                    ios=0
                    nlens=0
                endif 
                write(*,"('  Read nlens:   ',i10)") nlens 
                if (nlens.gt.1.or.nlens.lt.0) nlens=0
                !if (kmp.lt.7.or.kmp.eq.10) then
             !if(kmp.ne.7.and.kmp.ne.9.and.kmp.ne.8.and.kmp.ne.11) then
                 !nlens=0
                !endif
            case ('r200t')
                read(buffer, *, iostat=ios) r200t
                if (ios.eq.-1) then
                    ios=0
                    if (nlens.eq.1.and.naddlens.eq.0) then 
                      write(*,*) 'error: lensing mode requires an'
                      write(*,*) 'input value for r200. '
                      write(*,*) 'Switching to nlens=0'
                      nlens=0
                      r200t=0.0d0
                    endif
                endif 
                !write(*,"('  Read r200t:   ',f10.2)") r200t
            case ('rst')
                read(buffer, *, iostat=ios) rst
                if (ios.eq.-1) then
                    ios=0
                    if (nlens.eq.1.and.naddlens.eq.0) then 
                      write(*,*) 'error: lensing mode requires an'
                      write(*,*) 'input value for rs. '
                      write(*,*) 'Switching to nlens=0'
                      nlens=0
                      rst=0.0d0
                    endif
                endif
  
            case ('GaussianPrior')
                read(buffer, *, iostat=ios) ipriorGauss
                if (ios.eq.-1) then
                    ios=0
                    ipriorGauss=0
                endif
            case ('FlatPriorBeta')
                read(buffer, *, iostat=ios) ipriorBeta
                if (ios.eq.-1) then
                    ios=0
                    ipriorBeta=0
                endif
            case ('FlatPriorBeta0')
                read(buffer, *, iostat=ios) ipriorBeta0
                if (ios.eq.-1) then
                    ios=0
                    ipriorBeta0=0
                endif
c***********************************************************************
c*********** parameter limits ******************************************
            case ('r2low')  !r200 lower bound 
                read(buffer, *, iostat=ios) r2low
                if (ios.eq.-1) then
                    ios=0
                    r2low=0.0001d0
                endif
                pars(9)=r2low
            case ('r2up')  !r200 upper bound 
                                            
                read(buffer, *, iostat=ios) r2up
                if (ios.eq.-1) then
                    ios=0
                    r2up=100.0d0
                endif
                pars(10)=r2up
            case ('rclow')  !rc lower bound
                read(buffer, *, iostat=ios) rclow
                if (ios.eq.-1) then
                    ios=0
                    rclow=0.0001d0
                endif
                pars(11)=rclow
            case ('rcup')  !rc upper bound
                read(buffer, *, iostat=ios) rcup
                if (ios.eq.-1) then
                    ios=0
                    rcup=100.0d0
                endif
                pars(12)=rcup
            case ('rslow')  !rs lower bound
                read(buffer, *, iostat=ios) rslow
                if (ios.eq.-1) then
                    ios=0
                    rslow=0.0001d0
                endif
                pars(13)=rslow 
            case ('rsup')  !rs upper bound
                read(buffer, *, iostat=ios) rsup
                if (ios.eq.-1) then
                    ios=0
                    rsup=100.0d0
                endif 
                pars(14)=rsup
            case ('blow')  !beta lower bound
                read(buffer, *, iostat=ios) blow
                if (ios.eq.-1) then
                    ios=0
                    blow=0.0001d0
                endif
                pars(15)=blow
            case ('bup')  !beta upper bound
                read(buffer, *, iostat=ios) bup
                if (ios.eq.-1) then
                    ios=0
                    bup=100.0d0
                endif 
                pars(16)=bup 
            case ('A1low')  !first MG parameter lower bound
                read(buffer, *, iostat=ios) tmlow
                if (ios.eq.-1) then
                    ios=0
                    taddl=.True.
                    tmlow=0.0001d0
                endif
                pars(17)=tmlow
            case ('A1up')  !first MG parameter upper bound
                read(buffer, *, iostat=ios) tmup
                if (ios.eq.-1) then
                    ios=0
                    taddu=.True.                   
                    tmup=100.0d0
                endif
                pars(18)=tmup 
            case ('A2low')  !second MG parameter lower bound
                read(buffer, *, iostat=ios) scrlow
                if (ios.eq.-1) then
                    ios=0
                    saddl=.True.
                    scrlow=0.0001d0
                endif
                pars(19)=scrlow
            case ('A2up')  !second MG parameter upper bound
                read(buffer, *, iostat=ios) scrup
                if (ios.eq.-1) then
                    ios=0
                    saddu=.True.
                    scrup=100.0d0
                endif
                pars(20)=scrup 
                
           case ('A3low')  !third MG parameter lower bound
                read(buffer, *, iostat=ios) dA3low
                if (ios.eq.-1) then
                    ios=0
                    dA3low=0.0001d0
                endif

            case ('A3up')  !third MG parameter upper bound
                read(buffer, *, iostat=ios) dA3up
                if (ios.eq.-1) then
                    ios=0
                    dA3up=100.0d0
                endif

            case ('b2low')  !second beta lower bound
                read(buffer, *, iostat=ios) cb0low
                if (ios.eq.-1) then
                    ios=0
                    cb0low=0.0001d0
                endif
                
            case ('b2up')  !second beta parameter upper bound
                read(buffer, *, iostat=ios) cb0up
                if (ios.eq.-1) then
                    ios=0
                    cb0up=100.0d0
                endif

            case ('rbetalow')  !rbeta parameter upper bound
                read(buffer, *, iostat=ios) rbetalow
                if (ios.eq.-1) then
                    ios=0
                    rbetalow=0.00001d0
                endif
                
            case ('rbetaup')  !rbeta parameter upper bound
                read(buffer, *, iostat=ios) rbetaup
                if (ios.eq.-1) then
                    ios=0
                    rbetaup=10.0d0
                endif                
                
            case ('rjalow')  !rjaffe lower bounds
                read(buffer, *, iostat=ios) rjalow
                if (ios.eq.-1) then
                    ios=0
                    rjalow=0.0d0
                endif
            case ('rjaup')  !rjaffe upper bounds
                read(buffer, *, iostat=ios) rjaup
                if (ios.eq.-1) then
                    ios=0
                    rjaup=1.0d0
                endif

            case ('xmalow')  !M/L lower bounds
                read(buffer, *, iostat=ios) xmalow
                if (ios.eq.-1) then
                    ios=0
                    xmalow=0.0d0
                endif
            case ('xmaup')  !M/L upper bounds
                read(buffer, *, iostat=ios) xmaup
                if (ios.eq.-1) then
                    ios=0
                    xmaup=10.0d0
                endif
            case ('lensing_file')    
                read(buffer, *, iostat=ios) naddlens
                if (ios.eq.-1) then
                    ios=0
                    naddlens=0
                endif
            case ('nLogP')    
                read(buffer, *, iostat=ios) nLogP
                if (ios.eq.-1) then
                    ios=0
                    nLogP=1
                endif
            case default
            end select
         end if
      end do

      if (addrup) then
         if (read2) rupin=r200g !if not given, set the maximum fitting radius
                                !equal to the guess value of r200
      endif
      if (nbs.eq.-1.) kani=1 !   forced to MamLok if requested
      if (nbs.eq.-2.) kani=-1 !   if Hansen&Moore, beta(r) depends on rho(r)
      
      if(kmp.eq.8) then
        if (taddl) tmlow=-0.67
        if (taddu) tmup=8.0
        if (saddl) scrlow=-12.0
        if (saddu) scrup=12.0
      endif
   
      close(29)
      
      if (nrs.gt.-2) kfs2 = 0 !this is to avoid problems when not running over
      !Euclid mock data
       
      kmody=0 !modified gravity models that do not support bcg  
      if(kmp.eq.9) kmody=1
      if(kmp.eq.11.or.kmp.eq.12.or.kmp.eq.15.or.kmp.eq.16) kmody=1
      
      if(kmody.eq.1.or.kmp.eq.8) km2 = 1 !case of MG gravity with no BCG
      
      if(km2.eq.1.and.kbcg.eq.1) then
       write(*,*) 'BCG profile does not work with this models'
       write(*,*) 'excluding BCG from the analysis'
       kbdt = 0
       
       kbcg = 0
       
       kbary = 0 
       
      endif
      
      if(kmp.eq.17) then !The general chameleon case  always consider the bcg
        !and the baryons
        kbcg = 1
        kbary = 1
      endif

      
c************************************************************************

      
      if(r200g.gt.r2up.or.r200g.lt.r2low) then
         Stop('ERROR: GUESS VALUE r200 EXCEEDES PARAMETER LIMITS')
      endif 
      if(rsg.gt.rsup.or.rsg.lt.rslow) then
         Stop('ERROR: GUESS VALUE rs EXCEEDES PARAMETER LIMITS')
      endif 
      if(rcg.gt.rcup.or.rcg.lt.rclow) then
         Stop('ERROR: GUESS VALUE rc EXCEEDES PARAMETER LIMITS')
      endif 
      if(cbeg.gt.bup.or.cbeg.lt.blow) then
         Stop('ERROR: GUESS VALUE beta EXCEEDES PARAMETER LIMITS')
      endif 
            if(tmassg.gt.tmup.or.tmassg.lt.tmlow) then
         Stop('ERROR: GUESS VALUE A1 EXCEEDES PARAMETER LIMITS')
      endif 
            if(screeg.gt.scrup.or.screeg.lt.scrlow) then
         Stop('ERROR: GUESS VALUE A2 EXCEEDES PARAMETER LIMITS')
      endif
      if(dA3g.gt.dA3up.or.dA3g.lt.dA3low) then
         Stop('ERROR: GUESS VALUE A3 EXCEEDES PARAMETER LIMITS')
      endif
      
      if(cbe0g.lt.cb0low.or.cbe0g.gt.cb0up) then
         Stop('ERROR: GUESS VALUE cbe0g EXCEEDES PARAMETER LIMITS')
      endif
      
      if(rjafg.lt.rjalow.or.rjafg.gt.rjaup) then
         Stop('ERROR: GUESS VALUE rjafg EXCEEDES PARAMETER LIMITS')
      endif
      if(xmstarg.lt.xmalow.or.xmstarg.gt.xmaup) then
         Stop('ERROR: GUESS VALUE xmstarg EXCEEDES PARAMETER LIMITS')
      endif
      if(rbetag.lt.rbetalow.or.rbetag.gt.rbetaup) then
         Stop('ERROR: GUESS VALUE xmstarg EXCEEDES PARAMETER LIMITS')
      endif    

      tmass=tmassg          !set the mass to the guess value
      
      
c General chameleon screening: sub-case f(R)
      if (kmp.eq.9.and.nhone.lt.0) screeg=1./dsqrt(6.0d0)
      
      screen=screeg         !set screening to guess value

      dA3 = dA3g            !set A3 to guess value

c********* BCG equal to 0 **********************************************

        !set all BCG parameters to zero if not used
        
        if (kbcg.eq.0.and.kmp.ne.13) then
            nxmas = 0
            njaf = 0
            !kbary = 0
        endif
c***********************************************************************        
       
       !Define gas and gal mass
       rhog1 = rho0gas
       rsga1 = rsgas
       
       rhost1 = rhostar
       rst1 = rstar
       
c*********************************************************************** 
      
      write(*,731) animod, kani, massmod, kmp

 731  format(/,' You selected the anisotropy model ',A5,/,
     &         ' corresponding to the parameter kani= ',i2,/
     &         '  ',/
     &         ' You selected the mass  model ',A8,/,
     &     ' corresponding to the parameter kmp= ',i2,/)
   
      if (kmody.eq.1.or.kmp.gt.10) then
          if (tmassg.le.0) then
            stop('ERROR: inconsitent negative MG parameter')
          endif
          if (ipriorGauss.eq.0) then
           if (tmlow.le.0) tmlow=1.0e-4
           if (tmup.le.0) tmup=1.0e-4  
          endif
      endif
c exponent for f(R) Hu&Sawicki model (to be set manually)
      nhs=1
c Coupling constant for f(R)
      aQ=1./6
      
   
      
c     Stop if H&M beta(r) required and M(r) difft from NFW
c
      if (kani.lt.0.and.knfit.ne.1) then
         write(*,*) ' Hansen & Moore beta(r) currently '
         write(*,*) ' implemented only for pNFW N(R) model '
         stop
      endif
      
      kindt=0 !old features: add a unversal distribution for interlopers
      kgas=0  !old features: gas considered in a scaling-relation way

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
      read(*,9) fbcg !BCG file
      read(*,9) fUme !File from Keiichi
      
      open(20,file=frnvn,status='unknown')
      open(30,file=fbfn,status='unknown')
      open(40,file=fbn,status='unknown')
      open(50,file=ffn,status='unknown')
      open(60,file=fmlv,status='unknown')
      open(70,file=fsb,status='unknown')
      open(80,file=fsf,status='unknown')
      
      if (kbcg.eq.1.and.kbdt.eq.1) then
       
       open(81,file=fbcg,status='old') 
      endif
    
      !TEST 
      
     
      if (naddlens.eq.1) then
       
       write(*,*) 'Reading lensing from ', fUme
       open(82,file=fUme,status='old') 
      endif

!Information about the output files ************************************
!(Metadata
      write(30,*) 'record of N(R) model used'
             
      write(40,*) 'Binned project N(R) profile'
      write(40,*)  
      
      write(50,*) 'Best fit project N(R) profile'
      write(50,*) 'R[Mpc], N(R)' 
       
      WRITE(60,*) 'Table containing parameter values and -log(P)'
      WRITE(60,*) ''  
        
      write(70,*) 'Binned velocity dispersion profile'
      write(70,*) 'R[Mpc], VDP [km/s], err_up[km/s] err_low [km/s]'
      
      write(80,*) 'Fitted velocity dispersion profile'
      write(80,*) 'R[Mpc], VDP [km/s]'  

!***********************************************************************
      
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
!if weights are not given, read three columns in the data file 
      if (nweigh.eq.0) then
       read(10,*,end=111) dkpc,vkms,evkms 
       
      else
!Read four columns (weights is the last) 
       read(10,*,end=111) dkpc,vkms,evkms,wei
      endif
      j=j+1
      di(j)=dkpc/1.e3
      ve(j)=vkms
      eve(j)=evkms
      w(j)=1.0
      if (nweigh.ne.0) w(j)=wei
      
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

c     cosmological parameters
      omegal=Olam
      omega0=Omegam   
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
      write(*,*) ' After MG-MAMPOSSt: '
      write(*,*) ' '
      write(*,*) '    build output files of'
      write(*,*) '    - binned N(R), '
      write(*,*) '    - VDP of best-fit N(R),'
      write(*,*) '    - VDP of input N(R),' 
      write(*,*) '    - VDP solutions (''true'' values)'
      write(*,*) ' '

c     output a binned number density profile N(R)

      ibwt=0
      nbins=int(dsqrt(dfloat(nga)))
      do j=1,nga
         dns(j)=rso(j)
      enddo

      call sortp(dns,nga)
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
cc         facnorm=gc/(2.*pig)       
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

c     Here computes sigma_los 

      call dsort(rso,vso,nga,2)
      nbins=int(sqrt(float(nga)))/2.
      npbin=nga/nbins
      do j=1,nbins
         bins(j)=rso((j-1)*npbin)
         bins(j+1)=rso(j*npbin)
      enddo

      write(*,*) ' Using ',nbins,' bins for the VDP(R) '

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
            call sortp(syb,nboo)
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
      rinfinity=70.

      ius=80

      write(*,*) ' Evaluating expected VDP for '
      write(*,*) ' Max Lik solution values: '

 797  continue

      v200=10.*hz*r200

      write(*,516) r200,rc,r200/rc,rs,r200/rs,cbe,cbe0,rjaf,
     & xmstar,tmass,screen
 516  format('   r_200    = ',f6.3,
     &       /'   r_tracer = ',f6.3,' (c_tracer= ',f7.2,')'
     &       /'   r_mass   = ',f6.3,' (c_mass=   ',f7.2,')'
     &       /'   Anisotropy parameter = ',f8.4,
     &       /'   Second anisotropy parameter = ',f8.4,
     &       /'   BCG Jaffe radius = ',f8.4,
     &       /'   BCG star M/L = ',f8.4,
     &       /'   First MG parameter = ',f8.4,
     &       /'   Second MG parameter = ',f8.4,/)
c
c**************************************************************

      errrel=0.0005d0

c     evaluate sigma_r at some points (it will then interpolate)

      rlow=0.005d0    ! we want the profile in the inner region (12 Dec. 2014)

      do i=1,ninterp
         xx2 = dlog(2.d0*rinfinity)
         xx1 = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(ninterp-1)

         call dgaus8 (sr2int,xx1,xx2, errrel, risl, IERR)

         xmin=dexp(xx1)
         risok=dsqrt(risl*sr2out(xmin))
         if (risl.gt.1.8d195) risok=rismax 
         if (risl.le.0.d0) risok=rismin
         xris(i)=xx1
         yris(i)=dlog(risok)
         
      enddo
      write(*,*) ' sigma_r evaluated'
c
c compute spline coeffs for later interpolation of sigma_r
c
      call SPLINE_CUBIC_SET(ninterp,xris,yris,2,0.d0,2,0.d0,ypp2)
c    

      errrel=0.005d0
      errabs=0.

      do i=1,50
         xx1 = dlog(rlow)+dfloat(i-1)*
     &    dlog(1.001*rinfinity/rlow)/dfloat(50-1)
         xx2=dlog(2.*rinfinity)
         xmin=dexp(xx1)
         
         call dgaus8 (fa,xx1,xx2, errrel, ris2n, IERR)

         if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then 
            pnr=sigmar1(xmin)
         elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
            pnr=sigmar2(xmin)
         else
            pnr=sigmar3(xmin)
         endif
         ris=dsqrt(ris2n/pnr)
!         write (*,*) xmin, ','
!         write (*,*) ris, ','
         if (ris.eq.ris) write(ius,*) xmin,ris
      enddo
      close(ius)

      
      end
