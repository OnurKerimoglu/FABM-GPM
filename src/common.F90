#include "fabm_driver.h"

module gpm_common

   use fabm_types
   use fabm_standard_variables

   implicit none

   public

   real(rk),parameter :: CMass   = 12.011_rk
   real(rk),parameter :: qnRF = 14._rk/106._rk
   real(rk),parameter :: qpRF = 1._rk/106._rk
   real(rk),parameter :: qsRF = 15._rk/106._rk
   real(rk), parameter :: pi=acos(-1._rk)
   real(rk), parameter :: s2d = 86400.
   real(rk), parameter :: tres_pXn= 0.0    ! treshold value for phytoplankton losses and mortality
   real(rk), parameter :: TINY=1.0e-3      !a small number, e.g., used as the treshold value for exudation
   real(rk), parameter :: eps     = 1.0e-6 ! number smaller than "TINY", e.g, used as the treshold value for mortality
   
   ! !PUBLIC DERIVED TYPES:
   type type_env
     real(rk) :: temp,salt,doy,depth,zmax,I0,par,parmean
   end type
   
   !elements resolved
   type,public :: type_elms 
     real(rk) :: C
     real(rk) :: P     
     real(rk) :: N
     real(rk) :: Si
     !real(rk) :: Ccal
   end type
   
   ! dissolved inorganic material species resolved
   type, public :: type_dim 
     real(rk) :: C
     real(rk) :: P
     real(rk) :: NO3
     real(rk) :: NH4
     real(rk) :: Si
     !real(rk) :: Ccal
   end type
   
   ! data and procedures relevant for the nutrition status of organisms 
   type, extends(type_elms),public :: org_basic
     real(rk) :: QP     
     real(rk) :: QN
     real(rk) :: QPr
     real(rk) :: QNr
     real(rk) :: lim_dip
     real(rk) :: lim_dop
     real(rk) :: lim_p
     real(rk) :: lim_no3
     real(rk) :: lim_nh4
     real(rk) :: lim_n
     real(rk) :: limNP
     real(rk) :: nutlim
     logical  :: dop_uptake
     contains
     procedure :: get_nut_Q
     procedure :: get_nut_QLU
   end type
   
   !variables relevant for autotrophic processes
   type, extends(org_basic),public :: org_autotrophic
     real(rk) :: Chl
     real(rk) :: QChl
     real(rk) :: fI
     contains
     procedure :: get_fI
     procedure :: get_chl
     procedure :: get_PP_upt4fs
     procedure :: get_losses_exud
     procedure :: get_losses_mort_aut
   end type
   
   !variables relevant for heterotrophic processes
   type, extends(org_basic),public :: org_heterotrophic
     real(rk) ::  KzFact
     real(rk) ::  ftotC !,ftotP,ftotN
     contains
     procedure :: get_GRonMultiPrey
     procedure :: get_fMon
     procedure :: get_adjasef
     procedure :: get_losses_excr
     procedure :: get_losses_mort_het
   end type
   
   !structure to collect prey parameters (variable identifiers, stoichiometric ratios, etc)
   type,public                           :: prey_pars
      type (type_state_variable_id)      :: id_C,id_P,id_N
      type (type_dependency_id)          :: id_QPr,id_QNr
      type (type_diagnostic_variable_id) :: id_realpref
      real(rk)                           :: pref,C2Si
   end type
   
   !structure to collect prey data
   type,public                          :: prey_data
      real(rk),dimension(:),allocatable :: corpref
      real(rk),dimension(:),allocatable :: rpref
      real(rk),dimension(:),allocatable :: C
      real(rk),dimension(:),allocatable :: P      
      real(rk),dimension(:),allocatable :: N
      !real(rk),dimension(:),allocatable :: Chl !when dynamic
      real(rk),dimension(:),allocatable :: Si
      real(rk),dimension(:),allocatable :: grC
      real(rk),dimension(:),allocatable :: grP
      real(rk),dimension(:),allocatable :: grN
      !real(rk),dimension(:),allocatable :: grChl
      real(rk),dimension(:),allocatable :: grSi
      real(rk),dimension(:),allocatable :: Qr
      real(rk),dimension(:),allocatable :: QPr
      real(rk),dimension(:),allocatable :: QNr
      real(rk),dimension(:),allocatable :: weight
   end type
   
   ! a base model type to be potentially used for all 'life' modules
   type, extends(type_base_model),public :: type_GPMbase
     !state vars
     type (type_state_variable_id)      :: id_boundC,id_boundP,id_boundN

     !diagnostics
     type (type_diagnostic_variable_id) :: id_QP,id_QPr,id_QN,id_QNr,id_N2P
     type (type_dependency_id)          :: id_QPr_dep,id_QNr_dep
     type (type_dependency_id)          :: id_depth,id_temp
     type (type_global_dependency_id)   :: id_doy

     !Model parameters
     integer  :: metvel,metext,metIntSt,metTresp
     real(rk) :: kc,w,Q10,Tref
     real(rk) :: C2N,C2P
     real(rk) :: rmd,rmdq,frac_d2x
     real(rk) :: QPmax,QPmin,QNmax,QNmin
     
   end type

   !an intermediate 'life' structure, e.g., for uptake and recyling?

   ! a fabm-type to be used for describing potentially autotrophic organisms
   type, extends(type_GPMbase),public :: type_GPMaut
     !potentially shared
     !  Variable identifiers
     type (type_state_variable_id)      :: id_DIC,id_DIP,id_DINO3,id_DINH4,id_DISi
     type (type_state_variable_id)      :: id_SOC,id_DOC,id_DOP,id_DON
     type (type_state_variable_id)      :: id_det1C,id_det1P,id_det1N
     type (type_state_variable_id)      :: id_det2C,id_det2P,id_det2N,id_det2Si !, id_det2k
     type (type_state_variable_id)      :: id_O2
     
     !autotrophs
     !diagnostics
     type (type_diagnostic_variable_id) :: id_Closs,id_Ploss,id_Nloss
     type (type_diagnostic_variable_id) :: id_dPAR,id_NPPR,id_exudsoc
     type (type_diagnostic_variable_id) :: id_Cgain_A,id_Pgain_A,id_NO3gain_A,id_NH4gain_A
     type (type_diagnostic_variable_id) :: id_MuClim_A,id_Plim,id_Nlim,id_Silim
     type (type_diagnostic_variable_id) :: id_Chl,id_QChl
     !type (type_diagnostic_variable_id) :: id_pc_o2o = pc_dic = CGain_A
     type (type_dependency_id)          :: id_par,id_Chl_dep
     type (type_horizontal_dependency_id)::id_I_0

!    Model parameters
     logical  :: resolve_Si,dop_allowed !,resolve_cal,resolve_carb
     integer  :: metchl,metIresp,metCexc
     real(rk) :: kchl,gam
     real(rk) :: C2Si,C2Chl !,C2Ccal
     real(rk) :: vCmax,excess,Kp,Kno3,Knh4,Ksi
     real(rk) :: islope,Iopt,Imin
     real(rk) :: VPmax,VNmax
     !real(rk) :: c_max,rccalc_min,xkc_i,xkk_i,detach_min
   end type
   
   ! a fabm-type to be used for describing potentially mixotrophic organisms. 
   ! it extends from the type_autopar, and adds heterotrophic parameters
   type, extends(type_GPMaut),public :: type_GPMmixo
      !diagnostics
      type (type_diagnostic_variable_id) :: id_IngasC,id_IngasP,id_IngasN
      type (type_diagnostic_variable_id) :: id_IngC,id_IngP,id_IngN
      type (type_diagnostic_variable_id) :: id_asefC,id_asefP,id_asefN
      type (type_diagnostic_variable_id) :: id_IngunasC,id_IngunasP,id_IngunasN,id_IngunasSi
      type (type_diagnostic_variable_id) :: id_respC!,id_excrP,id_excrN
      !type (type_diagnostic_variable_id) :: id_o2o_zc = zc_dic= respC
      
!     Model parameters
      !heteroptrophic
      logical  :: dynpref !,resolve_Si
      integer  :: num_prey
      type(prey_pars),dimension(:),allocatable :: prpar
      real(rk) :: rmn,gmax,Kz,asefC,asefP,asefN,unas_detfrac
   end type

   
   ! Aggregate diagnostics for e.g., carbon budgets.
   type (type_bulk_standard_variable),parameter :: total_NPPR = type_bulk_standard_variable(name='total_NPPR',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
   !?
   
   contains 
   

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Grazing on Multiple Prey Items
!
! !INTERFACE:
   subroutine get_GronMultiPrey(org,mpar,prpar,fT,prdat,Ing)
   !
! !DESCRIPTION:
! Here, grazing on multiple prey is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_heterotrophic)     :: org
   type (type_GPMmixo),intent(in):: mpar 
   real(rk),intent(in)          :: fT
   type(prey_pars),intent(in)   :: prpar(:)
   type(prey_data)              :: prdat
   type (type_elms)             :: Ing
!  !LOCAL VARIABLES
   real(rk)                     :: prefmin !,KzFact,ftotC,ftotP,ftotN
   integer                      :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve the values of all grazed food and calculate the total food (perceived, based on pref values)
   prdat%corpref=0.0
   org%ftotC=0.0 !;ftotP=0.0;ftotN=0.0
   DO i=1,mpar%num_prey
     !correct the preference if required (indicated by a negative value):
     if (prpar(i)%pref .lt. 0) then
       !At nutrient replete conditions, some plankton species build colonies that reduce their edibility. 
       !To emulate this, (when pref <0) we assume that the preference (edibility) sigmoidially increases from 
       !prefmin (=pref/5 @Q=Qmax) to the originally specified (absolute) value (@Q=Qmin)
       prefmin=-1*prpar(i)%pref/5
       prdat%corpref(i)=prefmin+(-1*prpar(i)%pref-prefmin) * (1.0-1.0/(1.0+exp(- prdat%Qr(i)*10.0+5.0)))
       !write(*,'(A, I1, A, 3F13.9)') '-prey#',i,' fQprey, oldpref, newpref: ',prdat%Qr(i),prpar(i)%pref,prdat%rpref(i)
     else 
       prdat%corpref(i)=prpar(i)%pref
     end if
     
     !calculate the total pref-weighted food (needed only if prefs are dynamically adjusted)
     if (mpar%dynpref) then 
        org%ftotC=org%ftotC+prdat%corpref(i)*prdat%C(i)
        !org%ftotP=ftotP+prdat%corpref(i)*prdat%P(i)
        !org%ftotN=ftotN+prdat%corpref(i)*prdat%N(i)
     end if
     
   END DO
   
   !calculate grazing rate for each prey unit and total grazing  
   prdat%grC=0.0; prdat%grN=0.0;prdat%grP=0.0
   prdat%rpref=0.0
   DO i=1,mpar%num_prey
      !if the preferences are specified to change dynamically (e.g., as in Fasham 1990)
      if (mpar%dynpref) then 
        prdat%rpref(i)=prdat%corpref(i)*prdat%C(i)/org%ftotC
        !write(*,'(A,I1,A,3F8.4)'),' '//trim(GPMpl%name)//' prey#',i,' corpref,C/totC,rpref:',prdat%corpref(i),prdat%C(i)/ftotC,prdat%rpref(i)
      else
        prdat%rpref(i)=prdat%corpref(i)
      end if
      !weight of food-i among all available food
      org%KzFact=1.0
      prdat%weight(i)       = org%get_fMon(prdat%C(i),prdat%rpref(i),org%ftotC,mpar%Kz*org%KzFact)
      !grazing rate of the item
      prdat%grC(i)       = mpar%gmax*fT* prdat%weight(i)                 !molCprey/molCpred/d
      prdat%grP(i)       = prdat%grC(i)         * prdat%P(i)/prdat%C(i)  !molPprey/molCpred/d
      prdat%grN(i)       = prdat%grC(i)         * prdat%N(i)/prdat%C(i)  !molNprey/molCpred/d
      !prdat%grChl(i)     = prdat%grC(i)         * prdat%Chl(i)/prdat%C(i)  !molChlprey/molCpred/d
      if (mpar%resolve_Si) then
        prdat%grSi(i)      = prdat%grC(i)       * prdat%Si(i)/prdat%C(i) !molChlprey/molChlpred/d
        !write(*,*),' '//trim(mpar%name)//' prey#',i,'prSi,prC,grSi',prdat%Si(i),prdat%C(i),prdat%grSi(i)
      end if
      !molXprey/molCpred/d = molCprey/molCpred/d  * molXprey/molCprey
      
      !write(*,*),' '//trim(mpar%name)//' prey#',i,'grC,prC,pref,ftotC',prdat%grC(i),prdat%C(i),prdat%rpref(i),ftotC
      
      !write(*,*),' '//trim(GPMpl%name)//' prey#',i,'grC,grP,grN',prdat%grC(i),prdat%grP(i),prdat%grN(i)
   END DO
   
   !write(*,'(A, 2F13.10)'),' '//trim(GPMpl%name)//' grSi:',prdat%grSi  
   
   !calculate the ingestion rates (prdat%X*(1-fA))
   Ing%C=sum(prdat%grC)  !molCprey/molCpred/d
   Ing%P=sum(prdat%grP)  !molPprey/molCpred/d 
   Ing%N=sum(prdat%grN)  !molNprey/molCpred/d 
   if (mpar%resolve_Si) then
     Ing%Si=sum(prdat%grSi)  !molSiprey/molCpred/d 
   end if
   
   end subroutine get_gronmultiprey
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Monod formulation for zoophytoplankton grazing on multiple prey (as in Fasham 1990). 
!
! !INTERFACE:
   pure real(rk) function get_fmon(org,food,pref,foodtot,Kzeff)
!                                    
! !DESCRIPTION:
! !The function calculates only the grazing for a single one of those multiple prey items 
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_heterotrophic),intent(in) :: org
   !type (type_GPMmixo),intent(in) :: mpar
   real(rk), intent(in)          :: food,pref,Kzeff,foodtot!,!KzFact
!
!EOP
!-----------------------------------------------------------------------
!BOC

   get_fmon = food*pref / (Kzeff+foodtot)
   
   end function get_fmon
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust assimilation efficiencies
!
! !INTERFACE:
   subroutine get_adjasef(org,mpar,Ing,asef,Ingas,Ingunas)
   !
! !DESCRIPTION:
!  Here, the asssimlation efficiencies for each element is adjusted to maintain Qmin<Q<Qmax
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_heterotrophic)     :: org
   type (type_GPMmixo),intent(in) :: mpar 
   type(prey_data)              :: prdat
   type (type_elms)             :: Ing,asef,Ingas,Ingunas
!  !LOCAL VARIABLES
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !C
   asef%C=mpar%asefC
   
   if ((mpar%metIntSt .eq. 1) .and. (mpar%QPmax .ne. mpar%QPmin)) then !non-homeostatis:
     asef%P=mpar%asefP*(_ONE_-org%QPr)  !asef%P: decreases with increasing quota
     asef%N=mpar%asefN*(_ONE_-org%QNr)  !asef%P: decreases with increasing quota
   else                                 !homeostasis:
     asef%P=mpar%asefP
     if (asef%P*Ing%P .gt. asef%C*Ing%C*org%QP) then !P-surplus
       asef%P=asef%C*Ing%C*org%QP/Ing%P          !down-regulate asef%P according to asef%C
       !if (debw) write(*,'(A,2F13.10)'),' '//trim(mpar%name)//' (adj eP) P surplus.  aeC, aeP:',asef%C,asef%P
     else if (asef%P*Ing%P .lt. asef%C*Ing%C*org%QP) then !C-surplus
       asef%C=asef%P*Ing%P/(Ing%C*org%QP)        ! down-regulate asef%C according to asef%P
       !if (debw) write(*,'(A,2F13.10)'),' '//trim(mpar%name)//' (adj eP) C surplus.  aeC, aeP:',asef%C,asef%P
     endif !if C-surplus
                              !homeostasis:
     asef%N=mpar%asefN
     if (asef%N*Ing%N .gt. asef%C*Ing%C*org%QN) then !N-surplus
       asef%N=asef%C*Ing%C*org%QN/Ing%N              !down-regulate asef%N according to asef%C
       !if (debw) write(*,'(A,2F13.10)'),' '//trim(mpar%name)//' (adj eN) N surplus. aeC, aeN:',asef%C,asef%N
     else if (asef%N*Ing%N .lt. asef%C*Ing%C*org%QN) then !C-surplus
       asef%C=asef%N*Ing%N/(Ing%C*org%QN)            !down-regulate asef%C according to asef%N
       !if (debw) write(*,'(A,2F13.10)'),' '//trim(mpar%name)//' (adj eN) C surplus. aeC, aeN:',asef%C,asef%N
       !as asef%C was down-regulated, asef%P may need to be down-regulated as well

       if (asef%P*Ing%P .gt. asef%C*Ing%C*org%QP) then !P-surplus
         asef%P=asef%C*Ing%C*org%QP/Ing%P          !adjust asef%P according to asef%C
         !if (debw) write(*,'(A,2F13.10)'),' '//trim(mpar%name)//' (adj eP2) P surplus. aeC, aeP:',asef%C,asef%P
       endif !if P-surplus
     endif ! C-surplus
   end if !if homeostatis
     
   !C
   Ingas%C=org%C*Ing%C*asef%C ![molCprey/m3//d]
   Ingunas%C=org%C*Ing%C*(1.0-asef%C)
   !P
   Ingas%P=org%C*Ing%P*asef%P
   Ingunas%P=org%C*Ing%P*(1.0-asef%P)
   !write(*,'(A,2F7.4)'),' '//trim(mpar%name)//' QP,PgainH/CgainH:',QP,Ingas%P/Ingas%C
   !N
   Ingas%N=org%C*Ing%N*asef%N
   Ingunas%N=org%C*Ing%N*(1.0-asef%N)
   if (mpar%resolve_Si) then
     Ingunas%Si=org%C*Ing%Si*(1.0-0.0) !i.e., asefSi=0.0
     !write(*,'(A,2F10.6)'),' '//trim(mpar%name)//' IngunasSi:',Ingunas%Si
   end if
   !if (debw) write(*,'(A,2F7.4)'),' '//trim(mpar%name)//' QN,NgainH/CgainH:',QN,Ingas%N/Ingas%C
   
   end subroutine get_adjasef
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get nutrient Quotas, Limitations and Uptake rates
!
! !INTERFACE:
   subroutine get_nut_Q(org,bpar)
!
! !DESCRIPTION:
! Here, nutrient limited phytoplankton growth is formulated
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (org_basic)        :: org
   type (type_GPMbase),intent(in) :: bpar
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !if variation of internal stores is not resolved
   if (bpar%metIntSt .eq. 0) then

     !QUOTAS
     if (.true.) then !this is the default method
       org%QP=1./bpar%C2P
       org%P=org%C*org%QP
       org%QN=1./bpar%C2N
       org%N=org%C*org%QN
     else
       !QUOTAS (e.g., for debugging purposes)
       org%QP=org%P/org%C
       org%QPr=min(_ONE_,max(_ZERO_,(org%QP-bpar%QPmin)/(bpar%QPmax-bpar%QPmin)))
       org%QN=org%N/org%C
       org%QNr=min(_ONE_,max(_ZERO_,(org%QN-bpar%QNmin)/(bpar%QNmax-bpar%QNmin)))
     end if
     
   else if (bpar%metIntSt .eq. 1) then
     
     !QUOTAS
     org%QP=org%P/org%C
     org%QPr=min(_ONE_,max(_ZERO_,(org%QP-bpar%QPmin)/(bpar%QPmax-bpar%QPmin)))
     org%QN=org%N/org%C
     org%QNr=min(_ONE_,max(_ZERO_,(org%QN-bpar%QNmin)/(bpar%QNmax-bpar%QNmin)))
     
   end if
   end subroutine
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get nutrient Quotas, Limitations and Uptake rates
!
! !INTERFACE:
   subroutine get_nut_QLU(org,apar,fT,di,dom,lim,upt)
!
! !DESCRIPTION:
! Here, nutrient limited phytoplankton growth is formulated
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (org_basic)        :: org
   type (type_GPMaut),intent(in) :: apar
   real(rk),intent(in)           :: fT
   type(type_dim),intent(in)     :: di
   type (type_elms)              :: dom,lim
   type(type_dim)                :: upt
!  !LOCAL VARIABLES
   real(rk)                      :: mno3,mnh4
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !intermediate quantities required by multiple approaches
   
   !N limitation: no3,nh4
   !when NO3 and NH3 are not seperately resolved
   !org%lim_N = di%N/(di%N+apar%Kn)
   !when NO3 and NH3 are resolved
   mno3=di%NO3/apar%kno3
   mnh4=di%NH4/apar%knh4
   org%lim_no3= mno3 / (1+mno3+mnh4) !in EH, Q1x
   org%lim_nh4= mnh4 / (1+mno3+mnh4) !in EH, Q2x
   org%lim_N = org%lim_no3+org%lim_nh4
   !Note that lim%N = (mno3+mnh4)/(1+mno3+mnh4) < 1.0  (which should be the case)
   !other way of expressing lim_ terms:
   !!lim_no3=di%NO3/(di%NO3+apar%kno3*(1+di%NH4/apar%knh4))  
   !!lim_nh4=di%NH4/(di%NH4+apar%knh4*(1+di%NO3/apar%kno3))
     
   !P limitation: choose between dip/dop
   org%lim_dip= di%P/(di%P+apar%Kp)
   org%lim_P=org%lim_dip
   org%dop_uptake=.false.
   if (apar%dop_allowed) then
     org%lim_dop= dom%P/(dom%P+apar%Kp)
     if ((org%lim_dip<=0.7) .and. (org%lim_dop>org%lim_dip)) then 
       org%lim_P=org%lim_dop
       org%dop_uptake=.true.
     end if
   end if
     
   !C
   lim%C=1.0 ! todo: introduce a DIC limitation term to avoid growth when DIC=0
   !org%nutlim=lim%C
   
   !N&P
   !if variation of internal stores is not resolved
   if (apar%metIntSt .eq. 0) then

     !QUOTAS
     if (.true.) then !this is the default method
       org%QP=1./apar%C2P
       org%P=org%C*org%QP
       org%QN=1./apar%C2N
       org%N=org%C*org%QN
     else
       !QUOTAS (e.g., for debugging purposes)
       org%QP=org%P/org%C
       org%QPr=min(_ONE_,max(_ZERO_,(org%QP-apar%QPmin)/(apar%QPmax-apar%QPmin)))
       org%QN=org%N/org%C
       org%QNr=min(_ONE_,max(_ZERO_,(org%QN-apar%QNmin)/(apar%QNmax-apar%QNmin)))
     end if
     
     !LIMITATIONS
     lim%P=org%lim_P
     lim%N=org%lim_N
     
     !QR approximations
     org%QPr=lim%P
     org%QNr=lim%N
     
     !UPTAKE
     !At this stage, nut. uptake rates are not yet known; they will be based on the C-uptake (growth)
     
   !if variation of internal stores is  resolved according to droop approach
   else if (apar%metIntSt .eq. 1) then
     
     !QUOTAS
     org%QP=org%P/org%C
     org%QPr=min(_ONE_,max(_ZERO_,(org%QP-apar%QPmin)/(apar%QPmax-apar%QPmin)))
     org%QN=org%N/org%C
     org%QNr=min(_ONE_,max(_ZERO_,(org%QN-apar%QNmin)/(apar%QNmax-apar%QNmin)))
     
     !LIMITATIONS
     !fN=0 @ Q=Qmin and fN->1 as QP->inf
     lim%P= _ONE_ - apar%QPmin/org%QP
     lim%N= _ONE_ - apar%QNmin/org%QN
     
     !UPTAKE
     !P
     upt%P = org%C*fT*apar%VPmax*(1.0-org%QPr)*org%lim_P ![mmol/m3/d]
     !N
     !when no3 and nh4 are not resolved:
     !upt%N = org%C*fT*apar%VNmax*(1.0-org%QNr)*org%lim_N ![mmol/m3/d]
     !NO3,NH4 seperate:
     upt%NO3 = org%C*fT*apar%VNmax*(1.0-org%QNr)* org%lim_no3
     upt%NH4 = org%C*fT*apar%VNmax*(1.0-org%QNr)* org%lim_nh4
   else
     call apar%fatal_error('common.F90:get_fI','for '//trim(apar%name)// ' specified metIntSt option is not available')
   end if
   
   org%limNP=min(lim%P,lim%N)
   org%nutlim=org%limNP
   
   if (apar%resolve_Si) then
     !QUOTAS
     org%Si=org%C / apar%C2Si !for consistency
     !LIMITATION
     lim%Si=di%Si/(di%Si+apar%Ksi)
     org%nutlim=min(org%nutlim,lim%Si) !lim_nps in EH
     !UPTAKE
     !At this stage, nut. uptake rates are not yet known; they will be based on the C-uptake (growth)
   end if
   
      !cal
   !if (self%resolve_carb .and. self%resolve_cal) then
     !amount of calcite-C bound to the organism
     !org%Ccal=org%C / self%C2Ccal
     
     !L838- if coupled to carbonate model
     !if (.false.) then
     !  omega=_ONE_
       !omega=ch*ca/aksp !EH 840 ??
     !else
     !  omega=_ZERO_
     !end if
      
     !if (omega .gt. 1.0) then
     !  phylim_calc=(omega-1.0)/((omega-1.0)+xkk) !L856
     !else
     !  phylimc_calc=_ZERO_
     !end if
     !org%nutlim=min(org%nutlim,phylimc_calc) !lim_npc in EH
     
     !loss of calcite: L893 -
     !if a maximum of coccoliths is reached, then all surplus coccospheres are dettached,
            !else, only 10 % are dettached(Tyrell&Taylor)
     !rccalc_neu = (biop3c +f_FromTo(i_dic_p3c) -f_FromTo(i_p3c_doc)         &
     !              -f_FromTo(i_p3c_soc) -f_FromTo(i_p3c_d1c)         &
     !              -f_FromTo(i_p3c_d2c))/(biop3k+f_FromTo(i_dic_p3k))

     !if (rccalc_neu<rccalc_min) then
     !   plossk = biop3k +f_FromTo(i_dic_p3k) -((biop3c+f_FromTo(i_dic_p3c)  &
     !                   -f_FromTo(i_p3c_doc) -f_FromTo(i_p3c_soc)           &
     !                   -f_FromTo(i_p3c_d1c) -f_FromTo(i_p3c_d2c))/rccalc_min)
               !rccalc = rccalc_min
     !else
        !plossk = detach_min*biop3k + (f_FromTo(i_p3c_d1c)+f_FromTo(i_p3c_d2c))/rccalc_neu
     !   plossk = detach_min*biop3k + (f_FromTo(i_p3c_d1c)+f_FromTo(i_p3c_d2c))/rccalc
        !rccalc = rcc
   !end if
   
   end subroutine
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get mortality loss rates for heterotrophs
!
! !INTERFACE:  
  subroutine get_losses_mort_het(org,mpar,fT,mort)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_heterotrophic)       :: org
   type(type_GPMmixo), intent(in)  :: mpar
   real(rk), intent(in)           :: fT
   type(type_elms)                :: mort
! !LOCAL VARIABLES
   real(rk)                       :: fy
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   fy=_ONE_
   if (org%C .lt. eps) fy=_ZERO_
   mort%C = fy * (mpar%rmd*fT + mpar%rmdq*org%C) * org%C
   mort%P = mort%C*org%QP !pic_d1c+pic_d2c
   mort%N = mort%C*org%QN
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get mortality loss rates for autotrophs
!
! !INTERFACE:  
  subroutine get_losses_mort_aut(org,apar,fT,mort)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_autotrophic)        :: org
   type(type_GPMaut), intent(in)  :: apar
   real(rk), intent(in)          :: fT
   type(type_elms)               :: mort
! !LOCAL VARIABLES
   real(rk)                      :: fy
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   fy=_ONE_
   if (org%C .lt. eps) fy=_ZERO_
   mort%C = fy * (apar%rmd*fT + apar%rmdq*org%C) * org%C
   mort%P = mort%C*org%QP !pic_d1c+pic_d2c
   mort%N = mort%C*org%QN
   if (apar%resolve_Si) then
     mort%Si=mort%C/apar%C2Si
   end if
   !if (apar%resolve_carb) then
   !  if (apar%resolve_cal) then
       !mort%Ccal=max(_ZERO_,plossk) !p3k_d2k !L910
   !  else 
   !    if (.not. apar%resolve_Si) then !non-diatoms
         !psk_d2k (L.1489)
         !mort%Ccal = (mort%C + exud%C)/apar%C2Ccal !psk_d2k
        !else
         !mort%Ccal=0.0 !diatoms accumulate no caCO3, so no loss thereof
   !    end if !if resolve_Si
   !  end if !if resolve_cal
   !end if !if resolve_carb
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get exudation loss rates
!
! !INTERFACE:  
  subroutine get_losses_exud(org,apar,upt,exud)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_autotrophic)       :: org
   type(type_GPMaut), intent(in)  :: apar
   type(type_dim), intent(in)   :: upt
   type(type_elms)               :: exud
! !LOCAL VARIABLES
   real(rk)                      :: fy
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   !L796 -
   fy = _ONE_
   if (org%C .lt. TINY) fy = _ZERO_ ! treshold for mortality/loss niche (tres_pXn)
   exud%C = fy * apar%gam * upt%C  
   exud%P = exud%C * org%QP
   exud%N = exud%C * org%QN
   if (apar%resolve_Si) then
     !in EH exud%Si is not included, but without it, Si is not conserved
     exud%Si = exud%C/apar%C2Si
   end if
   !exud%calC?
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get excretion loss rates
!
! !INTERFACE:  
  subroutine get_losses_excr(org,mpar,fT,excr)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_heterotrophic)      :: org
   type(type_GPMmixo), intent(in) :: mpar
   real(rk), intent(in)          :: fT
   type(type_elms)               :: excr
! !LOCAL VARIABLES
   real(rk)                      :: fy
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   fy = _ONE_
   if (org%C .lt. TINY) fy = _ZERO_ ! treshold for mortality/loss niche (tres_pXn)
   excr%C=fy*mpar%rmn*fT*org%C
   excr%P=excr%C*org%QP
   excr%N=excr%C*org%QN
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get chlorophyll bulk and quota value
!
! !INTERFACE:  
  subroutine get_chl(org,apar,env)
!
! !DESCRIPTION:
! Here, light limited growth is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_autotrophic)       :: org
   type(type_GPMaut), intent(in) :: apar
   type(type_env), intent(in)    :: env
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   select case (apar%metchl)
     case default
       call apar%fatal_error('common.F90:get_fT','for '//trim(apar%name)// ' specified metchl option is not available')
     case (0)
       org%QChl=1./apar%C2Chl * 12.0
       !mgChl/mmolC = gChl/gC  * 12gC/1molC
       org%Chl=org%C*org%QChl !mgChl
     case (1) 
       org%QChl=(0.003 + max(0.0,0.0154*exp(0.05*env%temp)*exp(-0.059*env%parmean)*org%limNP))*12.0
       !todo: parmean=?
       !aut%QChl=12*
       org%Chl=org%C*org%QChl
       
      !case(2) !Geider 1997
      !case(3) !IA
   end select
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light limitation 
!
! !INTERFACE:
   subroutine get_fI(org,apar,fT,env)
!
! !DESCRIPTION:
! Here, light limited growth is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class(org_autotrophic)       :: org
   type(type_GPMaut), intent(in) :: apar
   real(rk), intent(in)          :: fT
   type(type_env),intent(in)     :: env
   real(rk)                      :: fI,mumax
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !some formulations of light limitation depend on mumax [/d]
   !if the fixed stoich. is employed, just convert the units of vCmax[/s]
   !if the droop approach is used, vCmax is in fact V(@Q=Qinf). calculate mumax=vcmax*f(Q=Qmax)
   if (apar%metIntSt .eq. 0) then
     mumax=apar%vCmax*s2d
   else if (apar%metIntSt .eq. 1) then
     mumax=apar%vCmax*s2d*(_ONE_ - apar%QPmin/apar%QPmax)
   end if

   if (env%par .gt. 0.0) then
    !mumax=max(1e-9,mumax0) !to avoid 0/0 when par=0
     select case (apar%metIresp)
      case default
        call apar%fatal_error('common.F90:get_fI','for '//trim(apar%name)// ' specified metIresp option is not available')
      case(1) !monod type, e.g., schwaderer et al
        org%fI=env%par/(env%par+mumax/apar%islope)
      case(2) !with light inhibition type, e.g., schwaderer et al
        org%fI=env%par/( env%par**2 *mumax/apar%islope/apar%Iopt**2+ env%par*(1-2*mumax/apar%islope/apar%Iopt) +mumax/apar%islope ) 
      case(3) !tanh, eg. merico & oguz
        org%fI=TANH(env%par*apar%islope)
      case(4) ! Light acclimation formulation based on surface light intensity (Steele 1962)
        org%fI=env%par/max(0.25*env%I0,apar%Imin) *exp(_ONE_-env%par/max(0.25*env%I0,apar%Imin))
      !case (5) !todo:based on theta 
      
      ! EH light limitation
      ! Iopt=get_Iopt()
      ! fPAR(j,k,i) = par(j,k,i)/Iopt(j,k,i)*exp(one-par(j,k,i)/Iopt(j,k,i))
       
     end select
   else
     org%fI=0.0 !when mumax=0 and par=0, fI can become 0/0=NaN
   end if
   
   end subroutine
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get nutrient Quotas, Limitations and Uptake rates
!
! !INTERFACE:
   subroutine get_PP_upt4fs(org,apar,fT,upt,exud_soc)
!
! !DESCRIPTION:
! Here, nutrient limited phytoplankton growth is formulated
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (org_autotrophic)      :: org
   type (type_GPMaut),intent(in) :: apar
   real(rk),intent(in)           :: fT
   type (type_dim)              :: upt
   real(rk)                      :: exud_soc
! !LOCAL VARIABLES
   real(rk)                      :: vC_fTfI,vC_fTfIfN
   
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !light and temp. limited carbon uptake (growth) rate
   !some definitions to facilitate comparison with EH:
   ! vC_fTfI: [mmolC/m3/d], bulk C uptake rate as f(T,I). In EH, this is the recurrent 'p1C*fT1*fPAR(I,vp1,p1c)', and in Lorkowski et al, it's the 'Gross C fixation)
   ! vC_fTfIfN: [mmolC/m3/d], bulk C uptake rate as f(T,I,nutlim). In EH, this is 'Redfield C fixation' (dic_pxc_red)
   ! Aupt%C: [mmolC/m3/d], final (to appear in RHS) bulk C uptake rate. In EH, this is the 'Effective C fixation' (dic_pxc) 
   ! loss_soc = [mmolC/m3/d], bulk loss rate in SOC form. In EH, this is 'SOC exudation' (p1c_soc)
   
   vC_fTfI = org%C *apar%vCmax*fT*org%fI !in EH: dic_p1c (first calculation, i.e., growth without nut lim )
   !write(*,*)'apar%vCmax,fT,org%fI',apar%vCmax,fT,org%fI
   vC_fTfIfN = vC_fTfI * org%nutlim !in EH: dic_p1c_red
   select case(apar%metCexc)
     case default
       call apar%fatal_error('common.F90:get_PP_upt4s','for '//trim(apar%name)// ' specified metCexc option is not available')
     case(0)  !no excess uptake   
       upt%C = vC_fTfIfN !+ self%excess(vC_fTfI - vC_fTfIfN)
       exud_soc = 0.0 
     case(1)  !as in ECOHAM (Lorkowski et al 2012)
       upt%C = vC_fTfIfN + apar%excess*(vC_fTfI - vC_fTfIfN)
       !exud_soc = vC_fTfI - vC_fTfIfN !original code from EH leads to C leakage from phytoplankton, and instable Q
       exud_soc = apar%excess*(vC_fTfI - vC_fTfIfN) !suggested fix
       !write(*,'(A,3F14.10)')'vC_fTfI,vC_fTfIfN,exud_soc',vC_fTfI,vC_fTfIfN,exud_soc
     !case(2) then as in RECOM (Schartau et al 2007)
   end select
   
   !calculate the balanced-uptake terms (Si, and potentially N and P too)
   if (apar%metIntSt .eq. 0) then
    upt%P= vC_fTfIfN * org%QP !/apar%%C2P
    !when NO3 and NH3 are not seperately resolved
    !upt%N= vC_fTfIfN * org%QN !/apar%%C2N
    !when NO3 and NH3 are resolved:
    upt%NO3=upt%C*org%QN*org%lim_no3 !in EH, lim_no3=Q1x 
    upt%NH4=upt%C*org%QN*org%lim_nh4 !in EH, lim_nh4=Q21
   end if
   !uptake of Si and Ccal is based on C uptake anyway.
   if (apar%resolve_Si) then
     upt%Si= vC_fTfIfN/apar%C2Si
   end if
   !if (env%par .gt. 0.0) write(*,'(A,4F14.8)')'vC,fT,fI,fn',apar%%vCmax,fT,fI,org%nutlim 
   !if (apar%%resolve_carb) then
   !  if (apar%%resolve_cal) then
       !(dic_p3k)
       !upt%Ccal = org%C*fT*apar%%c_max*phylim_calc !*fp3k (fp3k=?)
       !modify excess uptake: assume it's additionally limited by calcification lim. (why?)
       !exud_soc =   exud_soc * phylim_calc
     !else
     !  if (.not. apar%%resolve_Si) then
         !if not resolving Silicate and not being cocco, assume a fixed caCO3 uptake 
         !upt%Ccal=Aupt%C/apar%%C2Ccal !dic_psk
     !  else 
         !upt%Ccal=0.0 !otherwise no CaCO3 was used for skeleton building
     !  end if !if not resolve_Si
     !end if !if resolve_cal
   !end if !if resolvie_carb
   
   end subroutine
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Temperature function
!
! !INTERFACE:
   real(rk) function get_fT(bpar,temp)
!
! !DESCRIPTION:
! Here, temperature response functions are formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type(type_GPMbase), intent(in) :: bpar
   real(rk), intent(in)         :: temp
   !real(rk)                     :: get_fT
!
!EOP
!-----------------------------------------------------------------------
!BOC

   select case (bpar%metTresp)
     case default
       call bpar%fatal_error('common.F90:get_fT','for '//trim(bpar%name)// ' specified metTresp option is not available')
     case(0)
       get_fT=1.0
     case (1)
       !fTeh=exp( alog(1.5)*((env%temp-10.)/10.) )
       get_fT=bpar%Q10**((temp-bpar%Tref)/bpar%Tref)
     !case (2)
     !  get_fT=exp(-((temp-bpar%Topt)**2)/bpar%Tint**2)
   end select
   
   end function 
!EOC
!-----------------------------------------------------------------------

! required if light limitation is calculated as in Ecoham
!function get_Iopt
!    !from eco_meteo.F90:
!    r_Iopt =  0.25  ! relaxation (-time) for adapting Iopt to current light
!    dep = dep + 0.5*dz(j,k,i) ! depth at center of current layer
!    xeps = extw_field(j,k,i)+exts_field(j,k,i)+extp_field(j,k,i)
!    turb_up     = turb_low
!    turb_center = turb_up + xeps*dz(j,k,i)*0.5
!    turb_low    = turb_up + xeps*dz(j,k,i)
!    par(j,k,i)  = pafr*radsol(j,i)*exp(-turb_center)
!    !CL ! calculate mean daily PAR using mean surface radiation and current turbidity
!    !par_mean(j,k,i) = pafr*radsol_mean(j,i)*exp(-turb_center)
!    par_mean(j,k,i) = 0.4*radsol_mean(j,i)*exp(-turb_center) ! *0.4=factor for wm-1 to mol quanta m-2 d-1
!    adepth  = min(adepth_max,dep)
!    ! ---- MK version ----
!    !if (adepth .le. adepth_max)then
!    !   turb_adepth = turb_up + xeps*min(dz(j,k,i)*0.5,adepthmax-adepth)
!    !endif
!    ! ---- jp revision ----
!    if (k.eq.1) xeps1 = xeps
!    turb_adepth = xeps1*adepth
!    ! --------------------
!    Iavg_adepth = pafr*radsol_mean(j,i)*exp(-turb_adepth)
!    dIopt_dt(j,k,i) = r_Iopt*(min(Iopt_max,max(Iopt_min,Iavg_adepth))-Iopt(j,k,i))
!    dep = dep + 0.5*dz(j,k,i)  ! depth at bottom of current layer
!     if (abs(opt_irr).le.1e-4) then
!      Iopt(:,:,:) = Iopt(:,:,:) + dIopt_dt(:,:,:)*dt
!     endif
!   end select


end module
