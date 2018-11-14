#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_zooplankton --- generic zooplankton model
!
! !INTERFACE:
   module gpm_zooplankton
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use gpm_common

   implicit none
   
!  default: all is private.
   private
   
! !PRIVATE DATA MEMBERS:
!
! !PUBLIC DERIVED TYPES:
   
   type, extends(type_EHmixo),public :: type_gpm_zooplankton
      contains
      
      procedure  :: initialize
      procedure :: do
      !procedure  :: get_vertical_movement !when implementing, e.g., active swimming
   end type
   
   !Local variables
   logical :: debug = .true.
   logical :: rightZ= .false. !writing values only one layer
   
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_zooplankton), intent(inout), target   :: self
   integer,                       intent(in)          :: configunit
!
! !LOCAL VARIABLES
!parameters of the self structure
   real(rk), parameter :: d_per_s = 1._rk/86400.
   integer                   :: i !prey counter
   character(len=2)          :: istr !i converted to character string
!EOP
!-----------------------------------------------------------------------
!BOC
   
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   
   !General Parameters
   call self%get_parameter(self%resolve_Si, 'resolve_Si',   '-',          'whether to resolve Si cycle',          default=.false.)
   
   call self%get_parameter(self%metext,  'metext',   '-',  'method for calculating light ext.',  default=0)
   call self%get_parameter(self%kc,      'kc',       'm^2/mmolC',   'specific light extinction',              default=0.03_rk )
   
   call self%get_parameter(self%metvel,   'metvel',   '-',  'method for calculating settling vel.',  default=0)
   call self%get_parameter(self%w,       'w',        'm/d',        'vertical velocity (<0 for sinking)',     default=0.0_rk, scale_factor=d_per_s)      
   
   call self%get_parameter(self%metTresp, 'metTresp',  '-',   'method to calcualte T resp.',   default=1)
   call self%get_parameter(self%Q10,     'Q10',     '-',    'Q10 for kinetic rates',       default=1.5_rk)
   call self%get_parameter(self%Tref,     'Tref',    'celcius',    'reference temperature for kinetic rates',default=10.0_rk)
   
   call self%get_parameter(self%rmd,     'rmd',     '/d',       'linear mortality rate -min',                  default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmdq,    'rmdq',    '/d',       'quadratic mortality rate',               default=0.001_rk, scale_factor=d_per_s)
   
   call self%get_parameter(self%metIntSt, 'metIntSt', '-', 'method for representing the internal states', default=0)
   call self%get_parameter(self%C2N,     'C2N',     'molC/molN',   'molar C:N ratio',  default=6.625_rk) !rf: 6.625
   call self%get_parameter(self%C2P,     'C2P',     'molC/molP',   'molar C:P ratio',  default=132.5_rk) !rf:106
   call self%get_parameter(self%QPmax,    'QPmax',    'molP/molC',    'Max. P Quota',                 default=0.04_rk) 
   call self%get_parameter(self%QPmin,    'QPmin',    'molP/molC',    'Subsistance P Quota',          default=0.01_rk) 
   call self%get_parameter(self%QNmax,    'QNmax',    'molN/molC',    'Max. N Quota',                 default=0.22_rk) 
   call self%get_parameter(self%QNmin,    'QNmin',    'molN/molC',    'Subsistance N Quota',          default=0.12_rk)
   
   ! HETEROTROPHIC PARAMETERS
   call self%get_parameter(self%rmn,     'rmn',     '/d',       'loss rate to nutrients',                 default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%gmax,        'gmax',        '/d',     'max. grazing rate',                    default=1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Kz,          'Kz',          'mmolC/m3','half saturation of ingestion rate',    default=0.00001_rk)
   call self%get_parameter(self%asefC,       'asefC',       '-',       'assimilation efficiency of the ingested C',default=0.75_rk)
   call self%get_parameter(self%asefP,       'asefP',       '-',       'assimilation efficiency of the ingested P',default=1.0_rk)   
   call self%get_parameter(self%asefN,       'asefN',       '-',       'assimilation efficiency of the ingested N',default=1.0_rk) 
   call self%get_parameter(self%unas_detfrac,'unas_detfrac','-',       'detritus fraction of the unassimmilated Ingestion',default=0.75_rk)
   call self%get_parameter(self%num_prey,    'num_prey',    '-',       'number of prey targets',               default=8)
   call self%get_parameter(self%dynpref,     'dynpref',     '-',       'dynamic preference switch',            default=.false.)
   
   ! Register state variables
   call self%register_state_variable(self%id_boundC,'C','mmolC/m^3','bound carbon', & 
                                    minimum=_ZERO_, specific_light_extinction=self%kc,vertical_movement=self%w*d_per_s)
   
   if (self%metIntSt .eq. 0) then
     !call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
     !                               minimum=_ZERO_)
     !call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
     !                               minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundC,scale_factor=1./self%C2P)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundC,scale_factor=1./self%C2N)
   else if (self%metIntSt .eq. 1) then  ! .or. (self%metIntSt .eq. 0) !(for debugging purposes)
     call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
                                    minimum=_ZERO_,vertical_movement=self%w*d_per_s)
     call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
                                    minimum=_ZERO_,vertical_movement=self%w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundP)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundN)
   end if
   

   ! Register links to external models, if variable names are provided in namelist. 
    allocate(self%prpar(self%num_prey))
    DO i=1,self%num_prey
     write (istr,'(i0)') i
     call self%get_parameter(self%prpar(i)%pref,'prey'//trim(istr)//'pref','-','preference for prey-'//trim(istr))
     !C
     call self%register_state_dependency(self%prpar(i)%id_C,'prey'//trim(istr)//'C')
     !P
     !todo: s.1) try to access self%prpar(i)%id_P, if available, setC2P=0. If not available, calc P from C2P
     !      s.2) try to access C2P, if not provided, set C2P=0, then try to access id_P 
     call self%register_state_dependency(self%prpar(i)%id_P,'prey'//trim(istr)//'P')
     if (self%prpar(i)%pref .lt. 0.0) then !if pref<0, real pref will be f(QPr(i)
       call self%register_dependency(self%prpar(i)%id_QPr,'prey'//trim(istr)//'QPr')
     end if
     !N
     call self%register_state_dependency(self%prpar(i)%id_N,'prey'//trim(istr)//'N')
     if (self%prpar(i)%pref .lt. 0.0) then !if pref<0, real pref will be f(QPr(i)
       call self%register_dependency(self%prpar(i)%id_QNr,'prey'//trim(istr)//'QNr')
     end if
     !Chl: todo
     if (self%resolve_Si) then
       !when C2Si is not explicitly provided, assume the prey is non-diatom, so Si2C=0 (C2Si=0 will mean the same)
       call self%get_parameter(self%prpar(i)%C2Si,'prey'//trim(istr)//'C2Si','-','C:Si ratio of prey-'//trim(istr),default=0.0_rk)
     end if
    END DO
   
   ! linking to DIM and POM pools are required for recycling
   !C
   call self%register_state_dependency(self%id_detC,'detC')
   !P
   call self%register_state_dependency(self%id_DIP,'DIP')
   call self%register_state_dependency(self%id_detP,'detP')
   !N
   call self%register_state_dependency(self%id_DIN,'DIN')
   call self%register_state_dependency(self%id_detN,'detN')
   !Si
   if (self%resolve_Si) then
     call self%register_state_dependency(self%id_DISi,'DISi')
     call self%register_state_dependency(self%id_detSi,'detSi')
   end if
   
   !Diagnostics
    
   !general
   call self%register_diagnostic_variable(self%id_QP,'QP','molP/molC', 'molar P:C ratio',         &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_QPr,'QPr','-', '(QP-QPmin)/(QPmax-QPmin)',         &
                                          output=output_time_step_averaged) 
   call self%register_diagnostic_variable(self%id_QN,'QN','molN/molC', 'molar N:C ratio',         &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_QNr,'QNr','-', '(QN-QNmin)/(QNmax-QNmin)',         &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_N2P,'N2P','molN/molP', 'molar N:P ratio',         &
                                          output=output_time_step_averaged)
                                          
   call self%register_diagnostic_variable(self%id_Closs,'Closs','mmolC/m^3/d', ' bulk C loss rate',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Ploss,'Ploss','mmolP/m^3/d', ' bulk P loss rate',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Nloss,'Nloss','mmolN/m^3/d', ' bulk N loss rate',   &
                                          output=output_time_step_averaged) 
   
   !heterotrophy                                       
   do i=1,self%num_prey
     write (istr,'(i0)') i
     call self%register_diagnostic_variable(self%prpar(i)%id_realpref,'real_pref_prey'//trim(istr),'-', 'realized pref for prey'//trim(istr),&
                                          output=output_time_step_averaged)
   end do
   
   call self%register_diagnostic_variable(self%id_IngC,'Ing_C','molC/molC/d', 'total sp. C grazing rate', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_IngunasC,'Ingunas_C','/d',  'unassimilated fraction of sp. ingestion of C', &
                                          output=output_time_step_averaged)                                          
   call self%register_diagnostic_variable(self%id_IngasC,'Ingas_C','/d', 'sp. C gain rate by heterotrophy',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_asefC,'asef_C','-',  'assimilation efficiency of C', &
                                          output=output_time_step_averaged)
   
   call self%register_diagnostic_variable(self%id_IngP,'Ing_P','molP/molC/d', 'total sp. P grazing rate', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_IngunasP,'Ingunas_P','/d',  'unassimilated fraction of sp. ingestion of P', &
                                          output=output_time_step_averaged)                                       
   call self%register_diagnostic_variable(self%id_IngasP,'Ingas_P','/d', 'sp. P gain rate by herbivory',   &
                                          output=output_time_step_averaged) 
   call self%register_diagnostic_variable(self%id_asefP,'asef_P','-',  'assimilation efficiency of P', &
                                          output=output_time_step_averaged)

   call self%register_diagnostic_variable(self%id_IngN,'Ing_N','molN/molC/d', 'total sp. N grazing rate', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_IngunasN,'Ingunas_N','/d',  'unassimilated fraction of sp. ingestion of N', &
                                          output=output_time_step_averaged)                                        
   call self%register_diagnostic_variable(self%id_IngasN,'Ingas_N','/d', 'sp. N gain rate by herbivory',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_asefN,'asef_N','-',  'assimilation efficiency of N', &
                                          output=output_time_step_averaged)
   
   call self%register_diagnostic_variable(self%id_IngunasSi,'Ingunas_Si','/d',  'unassimilated fraction of sp. ingestion of Si', &
                                          output=output_time_step_averaged) 
                                          
   ! Register environmental dependencies
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   
   !register the diagnostic var QPr also as a dependency, such that its value can be accessed
   call self%register_dependency(self%id_QPr_dep,'QPr', '-', 'relative QP')
   !register the diagnostic var QPr also as a dependency, such that its value can be accessed
   call self%register_dependency(self%id_QNr_dep,'QNr', '-', 'relative QN')
   
   return

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
   
! !INPUT PARAMETERS:
   class (type_gpm_zooplankton), intent(in)     :: self 
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: fT
   type (type_env)            :: env
   type (type_elms)           :: Ing,asef,Ingas,Ingunas,excr,mort
   type (org_heterotrophic)   :: org !quantities of the organismal system being modelled (state and derived values)
   integer                    :: i
   character(len=2)           :: istr !i converted to character string
   logical                    :: oneofPMZ,debw
   type(prey_data)            :: prdat
   allocate(prdat%corpref(self%num_prey))
   allocate(prdat%rpref(self%num_prey))
   allocate(prdat%weight(self%num_prey))
   allocate(prdat%C(self%num_prey))
   allocate(prdat%P(self%num_prey))
   allocate(prdat%N(self%num_prey))
   allocate(prdat%Si(self%num_prey))
   allocate(prdat%grC(self%num_prey))
   allocate(prdat%grP(self%num_prey))
   allocate(prdat%grN(self%num_prey))
   allocate(prdat%grSi(self%num_prey))
   allocate(prdat%Qr(self%num_prey))
   allocate(prdat%QPr(self%num_prey))
   allocate(prdat%QNr(self%num_prey))
   
   !TODO: resolve the case for metIntSt=0 of the prey?
   !TODO: query prey%id_N, idP,id_Si instead of resolve_N,P,Si
   !prdat: use dim property (?) instead of allocatable?
   
!EOP
!-----------------------------------------------------------------------
!BOC
   !if (debug) then
   !  _GET_GLOBAL_(self%id_doy,doy) !day of year
   !  write(*,'(A,2F7.3)')'doy:',doy
   !end if
   
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_
   
   !reset the rates
   Ing%C=0.0;Ingunas%C=0.0;Ingas%C=0.0;
   Ing%P=0.0;Ingunas%P=0.0;Ingas%P=0.0;
   Ing%N=0.0;Ingunas%N=0.0;Ingas%N=0.0;
   Ing%Si=0.0;Ingunas%Si=0.0;
   
   !-------------------------------------------------------------------------
   !PREPARE: retrieve variables
   !
   !_GET_GLOBAL_(self%id_doy,doy) !day of year
   _GET_(self%id_depth,env%depth)     ! depth
   _GET_(self%id_temp,env%temp)       ! temperature
   _GET_(self%id_par,env%par)         ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,env%I0)  ! surface short wave radiation

   !LOGICALS FOR DEBUGGING 
   if (debug .and. (env%depth .lt. 5.6 .and. env%depth .gt. 4.4)) then
     debw=.true.
   else 
     debw=.false.
   end if
   
   !if (debw) write(*,'(1A)',advance='no') ''
   
   ! RETRIEVE STATE VARIABLES and DEPENDENCIES
   ! retrieve bound elements
   _GET_(self%id_boundC,org%C) ! C bound to the organismal system
   if ((self%metIntSt .eq. 1)) then ! .or. (self%metIntSt .eq. 0) !(for debugging purposes)
     _GET_(self%id_boundP,org%P)
     _GET_(self%id_boundN,org%N) 
   end if
   
   !collect prey abundances
   !reset prdat%X, prdat%QX
   prdat%C=0.0
   prdat%P=0.0
   prdat%N=0.0
   prdat%Si=0.0
   prdat%QPr=0.0
   prdat%QNr=0.0
   prdat%Qr=0.0 !this is min(QPr,QNr)
   DO i=1,self%num_prey
     !C
     _GET_STATE_(self%prpar(i)%id_C,prdat%C(i))
     !P
     !TODO: if id_P present?
       _GET_STATE_(self%prpar(i)%id_P,prdat%P(i))
       !if pref(i)<0, then pref=f(Qr(i))
       if (self%prpar(i)%pref .lt. 0.0) then
         _GET_(self%prpar(i)%id_QPr,prdat%QPr(i))
       else 
         prdat%Qr(i)=1.0
       end if
     !end if
     !N
     !TODO: if id_N present?
       _GET_STATE_(self%prpar(i)%id_N,prdat%N(i))
       !if pref(i)<0, then pref=f(Qr(i))
       if (self%prpar(i)%pref .lt. 0.0) then !default=0.0, which means no Si in prey
         _GET_(self%prpar(i)%id_QNr,prdat%QNr(i))
         prdat%Qr(i)=min(prdat%QPr(i),prdat%QNr(i))
       else 
         prdat%Qr(i)=1.0
       end if
     !end if
     !Chl: Query id_Chl (if_dynamic?)
     !Si:
     if (self%resolve_Si) then
       if (self%prpar(i)%C2Si .gt. 0.0) then
         prdat%Si(i)=prdat%C(i)/self%prpar(i)%C2Si
       else
         prdat%Si(i)=0.0_rk
       end if
       !write(*,*),'prey#',i,' c2si,C,Si:',self%prpar(i)%C2Si,prdat%C(i),prdat%Si(i)
     end if
   END DO 
   !
   !END OF PREPARE
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !START OF CALCULATIONS: no other calculation outside this box
   !
   !Temperature function
   fT = get_fT(self%type_EHbase,env%temp)
   
   !Nutrient quotas
   call org%get_nut_Q(self%type_EHbase)
   
   !calculate grazing rate for each prey unit and total grazing 
   call org%get_GronMultiPrey(self%type_EHmixo,self%prpar,fT,prdat,Ing)
   !this gives grazrateX [molXprey/molXpred/d]
     
   !adjust the assimilation efficiencies to maintain  Qmin<Q<Qmax
   call org%get_adjasef(self%type_EHmixo,Ing,asef,Ingas,Ingunas)
   
   ! LOSSES
   call org%get_losses_excr(self%type_EHmixo,fT,excr)
   call org%get_losses_mort_het(self%type_EHmixo,fT,mort)   
   !
   !END OF CALCULATIONS: no other calculation outside this box
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !WRITE
   !
   ! SET RHS 
   !C
   _SET_ODE_(self%id_boundC, Ingas%C - excr%C - mort%C)
   !P
   if (self%metIntSt .eq. 1) then
     _SET_ODE_(self%id_boundP, Ingas%P - excr%P - mort%P) 
     _SET_ODE_(self%id_boundN, Ingas%N - excr%N - mort%N)
   end if
   
   !grazing targets
   DO i=1,self%num_prey
     !C
     _SET_ODE_(self%prpar(i)%id_C,-prdat%grC(i)*org%C) !molC/molC/d *molC/m3 =molC/m3/d !*(1.0-self%fracaut)
     !P
     !TODO: if id_P present?
       _SET_ODE_(self%prpar(i)%id_P,-prdat%grP(i)*org%C) !molP/molC/d *molC/m3 =molP/m3/d
     !end if
     !TODO: if id_N present?
       _SET_ODE_(self%prpar(i)%id_N,-prdat%grN(i)*org%C) !molN/molC/d *molC/m3 =molP/m3/d
     !end if
   END DO
   
   !recycling to nutrient pools
   !C
   _SET_ODE_(self%id_detC,mort%C)
   _SET_ODE_(self%id_detC,Ingunas%C*self%unas_detfrac)
   !P
   _SET_ODE_(self%id_DIP,excr%P)
   _SET_ODE_(self%id_detP,mort%P)
   _SET_ODE_(self%id_detP,Ingunas%P*self%unas_detfrac) 
   _SET_ODE_(self%id_DIP,Ingunas%P*(1.-self%unas_detfrac)) 
     !if (debw) write(*,'(A, 3F14.10)') 'rmd+rmdq,Ingunas%C,Ingunas%P', rmd+rmdq,Ingunas%C,Ingunas%P
   !N
   _SET_ODE_(self%id_DIN,excr%N)
   _SET_ODE_(self%id_detN,mort%P)
   _SET_ODE_(self%id_detN,Ingunas%N*self%unas_detfrac) 
   _SET_ODE_(self%id_DIN,Ingunas%N*(1.-self%unas_detfrac)) 
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'Ingunas%P', Ingunas%P*s2d
   !Si
   if (self%resolve_Si) then
     _SET_ODE_(self%id_detSi,Ingunas%Si*self%unas_detfrac)
     _SET_ODE_(self%id_DISi,Ingunas%Si*(1.-self%unas_detfrac))
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'Ingunas%Si', Ingunas%Si
   end if
   
   ! SET DIAGNOSTICS
   _SET_DIAGNOSTIC_(self%id_Closs,(mort%C+excr%C)*s2d)
   _SET_DIAGNOSTIC_(self%id_Ploss,(mort%P+excr%P)*s2d)
   _SET_DIAGNOSTIC_(self%id_Nloss,(mort%N+excr%N)*s2d)
   
   if ((self%metIntSt .eq. 0)) then
     !_SET_DIAGNOSTIC_(self%id_QP, org%QP)
     !_SET_DIAGNOSTIC_(self%id_QN, org%QN)
     _SET_DIAGNOSTIC_(self%id_QPr, org%QPr)
     _SET_DIAGNOSTIC_(self%id_QNr, org%QNr)
   else if ((self%metIntSt .eq. 1)) then
     _SET_DIAGNOSTIC_(self%id_QP, org%QP)
     _SET_DIAGNOSTIC_(self%id_QN, org%QN)
     _SET_DIAGNOSTIC_(self%id_QPr, org%QPr)
     _SET_DIAGNOSTIC_(self%id_QNr, org%QNr)
     _SET_DIAGNOSTIC_(self%id_N2P, org%N/org%P)
   end if
   !heterotrophy
   DO i=1,self%num_prey
    _SET_DIAGNOSTIC_(self%prpar(i)%id_realpref,prdat%rpref(i))	
   END DO
   _SET_DIAGNOSTIC_(self%id_IngC,Ing%C*s2d)
   _SET_DIAGNOSTIC_(self%id_IngasC, Ingas%C*s2d)
   _SET_DIAGNOSTIC_(self%id_IngunasC, Ingunas%C*s2d)
   _SET_DIAGNOSTIC_(self%id_asefC,asef%C)
   !P
   _SET_DIAGNOSTIC_(self%id_IngP,Ing%P*s2d)
   _SET_DIAGNOSTIC_(self%id_IngasP, Ingas%P*s2d)
   _SET_DIAGNOSTIC_(self%id_IngunasP, Ingunas%P*s2d)
   _SET_DIAGNOSTIC_(self%id_asefP,asef%P)
   !N
   _SET_DIAGNOSTIC_(self%id_IngN,Ing%N*s2d)
   _SET_DIAGNOSTIC_(self%id_IngasN, Ingas%N*s2d)
   _SET_DIAGNOSTIC_(self%id_IngunasN, Ingunas%N*s2d)
   _SET_DIAGNOSTIC_(self%id_asefN,asef%N)
   !Si
   if (self%resolve_Si) then
     _SET_DIAGNOSTIC_(self%id_IngunasSi, Ingunas%Si*s2d)
   end if
   !
   !END OF WRITE
   !-------------------------------------------------------------------------
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
   
   end subroutine do
!EOC
!-----------------------------------------------------------------------
end module gpm_zooplankton
