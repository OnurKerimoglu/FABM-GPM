#include "fabm_driver.h"

!avoid problems due to invalid diagnostic values:
!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_phytoplankton --- generic phytoplankton model
!
! !INTERFACE:
   module gpm_phytoplankton
!
!TODO:
! - diagnostic N and P even with constant stoichiometry (?)
! - self%metCexc -> general switch (aux)
! - self%resolve_Si -> general switch (aux)
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
   
   type, extends(type_GPMaut),public :: type_gpm_phytoplankton
      contains
      procedure  :: initialize
      procedure  :: do
      procedure  :: get_vertical_movement
      procedure  :: get_light_extinction
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
   class (type_gpm_phytoplankton), intent(inout), target   :: self
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
   call self%get_parameter(self%Idm_met, 'Idm_met',   '-',          'method to calculate daily average PAR',          default=0)
   call self%get_parameter(self%resolve_Si,  'resolve_Si',   '-',     'whether to resolve Si cycle',          default=.false.)
   call self%get_parameter(self%resolve_DIC,  'resolve_DIC',   '-',     'whether to resolve DIC',          default=.false.)
   call self%get_parameter(self%lim_Si,  'lim_Si',   '-',     'whether limited by Si',          default=.false.)
   
   call self%get_parameter(self%metext,  'metext',   '-',  'method for calculating light ext.',  default=0)
   call self%get_parameter(self%kc,      'kc',       'm^2/mmolC',  'C-specific light extinction',              default=0.0045_rk )
   call self%get_parameter(self%kc,      'kchl',     'm^2/gChl',   'chl-specific light extinction',              default=0.023_rk )
   
   call self%get_parameter(self%metvel,   'metvel',   '-',  'method for calculating settling vel.',  default=0)
   call self%get_parameter(self%w,       'w',        'm/d',        'vertical velocity (<0 for sinking)',     default=0.0_rk, scale_factor=d_per_s)  

   call self%get_parameter(self%metTresp, 'metTresp',  '-',   'method to calcualte T resp.',   default=1)
   call self%get_parameter(self%Q10,      'Q10',     '-',     'Q10 for kinetic rates fr',       default=1.5_rk)
   call self%get_parameter(self%Tref,     'Tref',    'celcius',    'reference temperature for aut. proc.',default=10.0_rk)

   call self%get_parameter(self%rmd,     'rmd',     '/d',      'linear mortality rate',                  default=0.035_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmdq,    'rmdq',    'm^3/mmolC/d', 'quadratic mortality rate',               default=0.01_rk, scale_factor=d_per_s) ! 0.01
   call self%get_parameter(self%frac_d2x, 'frac_d2x',     '-',        'fraction of fast sinking detritus in phyto mortality',     default=0.15_rk)

   call self%get_parameter(self%metIntSt, 'metIntSt', '-', 'method for representing the internal states', default=0)
   call self%get_parameter(self%C2N,     'C2N',     'molC/molN',   'molar C:N ratio',  default=6.625_rk) !rf: 6.625
   call self%get_parameter(self%C2P,     'C2P',     'molC/molP',   'molar C:P ratio',  default=132.5_rk) !rf:106
   call self%get_parameter(self%QPmax,    'QPmax',    'molP/molC',    'Max. P Quota',                 default=0.04_rk)
   call self%get_parameter(self%QPmin,    'QPmin',    'molP/molC',    'Subsistance P Quota',          default=0.01_rk)
   call self%get_parameter(self%QNmax,    'QNmax',    'molN/molC',    'Max. N Quota',                 default=0.22_rk)
   call self%get_parameter(self%QNmin,    'QNmin',    'molN/molC',    'Subsistance N Quota',          default=0.12_rk)

   ! Autotrophic Parameters
   call self%get_parameter(self%gam,     'gam',     '-',       'exudation fraction',       default=0.05_rk)
   call self%get_parameter(self%vCmax,    'vCmax',      '/d',       'maximum production rate',          default=1.1_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%VPmax,    'VPmax',     'mmolP/mmolC/d','Max. phosphorus uptake rate',          default=1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%VNmax,    'VNmax',     'mmolN/mmolC/d','Max. nitrogen uptake rate',          default=14.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%dop_allowed,   'dop_allowed',   '-',      'whether to DOP can be used instead of DIP',          default=.false.)
   call self%get_parameter(self%Kp,       'Kp',        'mmolP/m^3',   'half-saturation P concentration',  default=0.05_rk)
   call self%get_parameter(self%Kno3,     'Knh4',        'mmolN/m^3',   'half-saturation ammonium concentration',  default=4.2_rk)
   call self%get_parameter(self%Knh4,     'Kno3',        'mmolN/m^3',   'half-saturation nitrate concentration',  default=4.2_rk)
   call self%get_parameter(self%metIresp,'metIresp', '-',         'light response',                          default=1)
   call self%get_parameter(self%islope,   'islope',    'mmolC/m^2/W',  'slope of the P-I curve',                  default=0.05_rk) !when metIresp=1 & 3
   call self%get_parameter(self%Iopt,     'Iopt',      'W/m^2',  'half-saturation nutrient concentration',  default=100._rk) !when metIresp=2 & 4
   call self%get_parameter(self%islope_perchl,   'islope_perchl',    'gCm2/gChl/molQuanta',  'Chl specific slope of the P-I curve',                  default=10.0_rk) 
   call self%get_parameter(self%metCexc,   'metCexc',  '-',   'exc C uptake method',   default=0)
   call self%get_parameter(self%excess,   'excess',  '-',   'excess carbon assimilation fraction',   default=0.0_rk)
   call self%get_parameter(self%C2Si,     'C2Si',    'molC/molSi',  'molar C:Si ratio', default=5.76_rk)  !Brzezinski, 1985: 106/15=7.067
   call self%get_parameter(self%Ksi,      'Ksi',      'mmolSi/m^3',   'half-saturation Si concentration',  default=0.5_rk)
   call self%get_parameter(self%metchl,   'metchl', '-', 'method for representing chlorophyll', default=0)
   call self%get_parameter(self%Chl2C,    'Chl2C',   'gChl/gC',   'gChl:gC ratio', default=0.02_rk)
   call self%get_parameter(self%Chl2Cmax, 'Chl2Cmax','gChl/gC',   'gChl:gC ratio', default=0.05_rk)
   !call self%get_parameter(self%resolve_cal,  'resolve_cal',   '-',      'whether to resolve calcification',          default=.false.)
   !call self%get_parameter(self%C2Ccal,   'C2Ccal',   'molC/molC', 'molar ratio of organic to calcite C', default=5.76_rk)  !Brzezinski, 1985: 106/15=7.067
   
   ! Register state variables
   call self%register_state_variable(self%id_boundC,'C','mmolC/m^3','bound carbon', & 
                                    minimum=0.0_rk,no_river_dilution=.true.)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_boundC)
   if (self%metIntSt .eq. 0) then
     !call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
     !                               minimum=0.0_rk)
     !call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
     !                               minimum=0.0_rk)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundC,scale_factor=1./self%C2P)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundC,scale_factor=1./self%C2N)
   else if (self%metIntSt .eq. 1) then  ! .or. (self%metIntSt .eq. 0) !(for debugging purposes)
     call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
                                    minimum=0.0_rk,no_river_dilution=.true.)
     call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
                                    minimum=0.0_rk,no_river_dilution=.true.)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundP)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundN)
   end if
   !Si
   if (self%lim_Si) then
     !Si stored in plankton is not explicitly resolved
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_boundC,scale_factor=1./self%C2Si)
   end if
   
   !calC
   !if (self%resolve_calC) then
     !calC stored in plankton is not explicitly resolved
     !call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_phyC,scale_factor=1./self%C2Ccal)
   !end if
   
   ! linking to DIM and POM pools
   !C
   if (self%resolve_DIC) then
     call self%register_state_dependency(self%id_DIC,'DIC')
   end if
   call self%register_state_dependency(self%id_DOC,'DOC')
   if (self%metCexc .ne. 0) then
     call self%register_state_dependency(self%id_SOC,'SOC')
   end if
   call self%register_state_dependency(self%id_det1C,'det1C')
   call self%register_state_dependency(self%id_det2C,'det2C')
   !P
   call self%register_state_dependency(self%id_DIP,'DIP')
   call self%register_state_dependency(self%id_DOP,'DOP')
   call self%register_state_dependency(self%id_det1P,'det1P')
   call self%register_state_dependency(self%id_det2P,'det2P')
   !N
   call self%register_state_dependency(self%id_DINO3,'DINO3')
   call self%register_state_dependency(self%id_DINH4,'DINH4')
   call self%register_state_dependency(self%id_DON,'DON')
   call self%register_state_dependency(self%id_det1N,'det1N')
   call self%register_state_dependency(self%id_det2N,'det2N')
   !O2
   call self%register_state_variable(self%id_O2,'O2','mmol O2/m^3','O2 in water')
   !Si
   if (self%lim_Si) then
     call self%register_state_dependency(self%id_DISi,'DISi')
     call self%register_state_dependency(self%id_det2Si,'det2Si')
   end if
   !calC
   !if (self%resolve_calC) then
     ! call self%register_state_dependency(?)
   !end if
   
   !Diagnostics
    
   !general
   if ((self%metIntSt .eq. 0)) then 
     !call self%register_diagnostic_variable(self%id_QN,'QN','molN/molC', 'fixed molar N:C ratio',         &
     !                                     output=output_time_step_averaged)
     !call self%register_diagnostic_variable(self%id_QP,'QP','molP/molC', 'fixed molar P:C ratio',         &
     !                                     output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_QPr,'QPr','-', 'P limitation as estimated by Monod approach',         &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_QNr,'QNr','-', 'N limitation as estimated by Monod approach',         &
                                          output=output_instantaneous)                                          
   else if ((self%metIntSt .eq. 1)) then  
     call self%register_diagnostic_variable(self%id_QP,'QP','molP/molC', 'molar P:C ratio',         &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_QPr,'QPr','-', '(QP-QPmin)/(QPmax-QPmin)',         &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_QN,'QN','molN/molC', 'molar N:C ratio',         &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_QNr,'QNr','-', '(QN-QNmin)/(QNmax-QNmin)',         &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_N2P,'N2P','molN/molP', 'molar N:P ratio',         &
                                          output=output_instantaneous)
    else 
      call self%fatal_error('phytoplankton.F90/initialize:','for '//trim(self%name)// ' specified metIntSt  option is not available')
   end if
   
   call self%register_diagnostic_variable(self%id_Closs,'Closs','mmolC/m^3/d', ' bulk C loss rate',   &
                                           output=output_time_step_averaged) 
   call self%register_diagnostic_variable(self%id_Ploss,'Ploss','mmolP/m^3/d', ' bulk P loss rate',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Nloss,'Nloss','mmolN/m^3/d', ' bulk N loss rate',   &
                                          output=output_time_step_averaged) 
                                          
   ! Autototrophy
   call self%register_diagnostic_variable(self%id_Ilim,'limI','-', 'Light limitation',  &
                                          output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Plim,'limP','-', 'P limitation',  &
                                          output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_Nlim,'limN','-', 'N limitation',  &
                                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MuClim_A,'MuClimA','/d', 'contribution of autotrophy to C-lim growth',&
                                        output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NPPR,'NPPR','mmolC/m^3/d','-NPPR',              &
                                        output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Cgain_A,'CgainA','/d', 'sp. C gain rate by autotrophy',   &
                                        output=output_time_step_averaged)                                          
   call self%register_diagnostic_variable(self%id_Pgain_A,'PgainA','/d', 'sp. P gain rate by autotrophy',   &
                                        output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NO3gain_A,'NO3gainA','/d', 'sp. N gain rate by autotrophy',   &
                                        output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NH4gain_A,'NH4gainA','/d', 'sp. N gain rate by autotrophy',   &
                                        output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_exudsoc,'exudsoc','mmolC/m^3/d', ' bulk C exudation rate',   &
                                           output=output_time_step_averaged)  
                                           
   if (self%lim_Si) then
     call self%register_diagnostic_variable(self%id_Silim,'limSi','-', 'Si limitation',  &
                                          output=output_instantaneous)
   end if
  
   select case (self%metchl)
     case default 
       call self%fatal_error('phytoplankton.F90/initialize:','for '//trim(self%name)// ' specified metchl option is not available')
     case (0)
       call self%register_diagnostic_variable(self%id_diagChl,'Chl','mg/m^3', 'Chl estimated from a fixed Chl:C ratio',         &
                                          output=output_instantaneous)
       call self%add_to_aggregate_variable(total_chlorophyll,self%id_diagChl)
     case (1,2)
       call self%register_diagnostic_variable(self%id_diagChl,'Chl','mg/m^3', 'diagnostically calculated Chl concentration',         &
                                          output=output_instantaneous)  
       call self%add_to_aggregate_variable(total_chlorophyll,self%id_diagChl)
     case (20)
       call self%register_state_variable(self%id_boundChl,'Chl','mg/m^3','bound Chlorophyll', & 
                                    minimum=0.0_rk,no_river_dilution=.true.)
       call self%add_to_aggregate_variable(total_chlorophyll,self%id_boundChl)
       call self%register_diagnostic_variable(self%id_chlrho,'Rho_per_Tmax','-', 'regulatory factor for chl synthesis (rho/Tmax in G97 eq.4)', &
                                          output=output_instantaneous) 
   end select 
  
   if (self%metchl .ne. 0) then
    call self%register_diagnostic_variable(self%id_QChl,'QChl','gChl/gC', 'diagnostically calculated Chl:C ratio',         &
                                          output=output_instantaneous) 
   end if
   
   if (self%metvel .gt. 0) then
    call self%register_diagnostic_variable(self%id_sinkvel,'sinkvel','m s-1', 'sinking velocity',         &
                                          output=output_instantaneous) 
   end if
   
   !contribute to aggregate variables
   call self%add_to_aggregate_variable(total_NPPR,self%id_NPPR)
   
   !register the diagnostic var Qr and Chl also as dependencies, such that its value can be accessed
   call self%register_dependency(self%id_QPr_dep,'QPr', '-', 'relative QP')
   call self%register_dependency(self%id_QNr_dep,'QNr', '-', 'relative QN')
   call self%register_dependency(self%id_Chl_dep,'Chl', '-', 'Chl concentration')
   
   ! Register environmental dependencies
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   if (self%Idm_met .eq. 1) then
     call self%register_dependency(self%id_I0dm,'I0dm','W/m^2',       'daily mean I0')
     !I0dm is a diagnostic variable provided by another module specified in the yamlaveraged')
   else if (self%Idm_met .eq. 2) then
     call self%register_dependency(self%id_pardm,'PARdm','W/m^2',       'daily mean PAR')
     !PARdm is a diagnostic variable provided by another module specified in the yamlaveraged')
   end if

   return

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides
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
   class (type_gpm_phytoplankton), intent(in)     :: self 
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   logical                    :: debw
   real(rk)                   :: fT,exud_soc,vert_vel
   type (type_env)            :: env
   type (type_elms)           :: dom,Alim,exud,mort
   type (type_dim)            :: di,Aupt
   type (org_autotrophic)     :: org !quantities of the organismal system being modelled (state and derived values) 
   
!EOP
!-----------------------------------------------------------------------
!BOC
   
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   !Reset the rates
   !General
   mort%C=0.0; mort%P=0.0; mort%N=0.0; mort%Si=0.0
   !Autotrophy:
   Aupt%C=0.0; Aupt%P=0.0; Aupt%NH4=0.0; Aupt%NO3=0.0; Aupt%Si=0.0
   exud%C=0.0; exud%P=0.0; exud%N=0.0; exud%Si=0.0; exud_soc=0.0
   
   !-------------------------------------------------------------------------
   !PREPARE: retrieve variables
   !
   ! !Environmental Dependencies
   !_GET_GLOBAL_(self%id_doy,doy)     ! day of year
   _GET_(self%id_depth,env%depth)     ! depth
   _GET_(self%id_temp,env%temp)       ! temperature
   _GET_(self%id_par,env%par)         ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,env%I0)  ! surface short wave radiation
   if (self%Idm_met .eq. 1) then
     _GET_HORIZONTAL_(self%id_I0dm,env%Idm)  ! surface short wave radiation, daily output_time_step_averaged
   else if (self%Idm_met .eq. 2) then
     _GET_(self%id_pardm,env%Idm)  ! short wave radiation, daily output_time_step_averaged
   else
     env%Idm=-99.0_rk !should result in an error
   end if
   
   ! Debugging logic
   if (debug .and. (env%depth .lt. 1.0 .and. env%depth .gt. 0.0)) then
     debw=.true.
   else 
     debw=.false.
   end if
   
   !if (debw) write(*,'(1A)',advance='no') ''
   
   ! !Retrieve state variables
   
   ! retrieve bound elements
   _GET_(self%id_boundC,org%C) ! C bound to the organismal system
   if ((self%metIntSt .eq. 1)) then ! .or. (self%metIntSt .eq. 0) !(for debugging purposes)
     _GET_(self%id_boundP,org%P)
     _GET_(self%id_boundN,org%N) 
   end if 
   
   select case (self%metchl) 
     case (20) !for dynamical Chl options of metchl
       _GET_(self%id_boundChl,org%Chl)
       org%QChl =  org%Chl / (org%C * 12.0_rk)
       !gChl/gC = mgChl/m3 / (mmolC/m3 * 12mgC/mmolC)
       !write(*,*)'phyto L411. Chl,QChl:',org%Chl,org%QChl
   end select
   
   ! collect resources
   !if (self%resolve_DIC) then
   !_GET_(self%id_DIC,di%C)
   !end if
   _GET_(self%id_DIP,di%P)
   if (self%dop_allowed) then
     _GET_(self%id_DOP,dom%P)
   end if
   _GET_(self%id_DINO3,di%NO3) ! nitrate
   _GET_(self%id_DINH4,di%NH4) ! ammonium
   if (self%lim_Si) then
     _GET_(self%id_DISi,di%Si)
   end if
   !
   !END OF PREPARE
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !START OF CALCULATIONS: no other calculation outside this box
   !
   !Temperature function
   fT = get_fT(self%type_GPMbase,env%temp)
   
   !Nutrient quotas, limitation and uptake
   !todo: split the get_nutQ 
   call org%get_nut_QLU(self%type_GPMaut,fT,di,dom,Alim,Aupt)
   !write(*,*)'N: lim%P,lim%N,upt%P,upt%N:',Alim%P,Alim%N !,Aupt%P*s2d,Aupt%N*s2d
   
   !Light limitation
   !Chl is potentially needed to calculate light limitation
   select case (self%metchl) 
     case (0,1,2) !for diagnostic Chl options of metchl
       call org%get_chl(self%type_GPMaut,env)
   end select
   call org%get_fI(self%type_GPMaut,env)
   !if (env%par .gt. 0.0) write(*,*)'vc,is,io,im,fT,p',vC_atQmax,self%islope,self%Iopt,self%Imin,fT,env%par
   
   !Primary production (upt%C) and  uptake terms for elements resolved as fixed stoichiometry
   call org%get_PP_upt4fs(self%type_GPMaut,fT,Aupt,exud_soc)
   !obtain the chl. synth rates, of Chl is dynamic (if not, the subr will do nothing)
   call org%get_chldyn(self%type_GPMaut,env,Alim,Aupt)
   
   ! LOSSES
   call org%get_losses_exud(self%type_GPMaut,Aupt,exud)
   call org%get_losses_mort(self%type_GPMaut,1.0_rk,fT,mort)
   !
   !END OF CALCULATIONS: no other calculation outside this box
   !-------------------------------------------------------------------------
   
   !-------------------------------------------------------------------------
   !WRITE
   !
   ! RHS 
   !C
   !if (env%par .gt. 0.0) write(*,'(A,5F14.10)')'Cupt,Nupt,exsoc,exC,exN',Aupt%C*s2d,Aupt%N*s2d,exud_soc*s2d,exud%C*s2d,exud%N*s2d
   !if (env%par .gt. 0.0) write(*,'(A,3F14.8)')'Aupt%C:N,exud%C:N,mort%C:N:',(Aupt%C-exud_soc)/Aupt%N,(exud%C)/exud%N,mort%C/mort%N
   !write(*,'(A,2F14.10)')'C2N, RHS C/N',self%C2N,(Aupt%C - exud%C - mort%C - exud_soc)/(Aupt%N - exud%N - mort%N)
   !> **Phytoplankton Source Terms**
   !> dpC/dt = dic_pC - pC_DOC - pC_SOC - M^C
   _SET_ODE_(self%id_boundC, Aupt%C - exud%C - exud_soc- mort%C)
   if ((self%metIntSt .eq. 1)) then ! .or. (self%metIntSt .eq. 0)  !(for debugging purposes)
     !> dpP/dt = dip_pP - pC_DOP - M^P
     !> dpN/dt = no3_pN + nh4_pN - pC_DOP - M^N
     _SET_ODE_(self%id_boundP, Aupt%P - exud%P - mort%P)
     _SET_ODE_(self%id_boundN, Aupt%NO3+Aupt%NH4 - exud%N - mort%N)
   end if
   
   select case (self%metchl)
     case (20) !list here the dynamic Chlorophyll options
       !> dpChl/dt = \%rhochl*dic_pC -(pC_DOC + pC_SOC + M^C)*12*Theta
       _SET_DIAGNOSTIC_(self%id_chlrho,_REPLNAN_(org%chlrho))
       _SET_ODE_(self%id_boundChl, org%chlsynth - (exud%C + exud_soc + mort%C)*12.0_rk*org%QChl)
       !                             mgChl/m3/s - (mmolC/m3/s)*12.0mgC/mgChl *mgChl/mgC
       !write(*,*)'org%chlsynth', org%chlsynth,'loss:', - (exud%C + exud_soc + mort%C)*org%QChl,'RHS:',org%chlsynth - (exud%C + exud_soc + mort%C)*org%QChl
   end select
   
   !O2
   _SET_ODE_(self%id_O2,Aupt%C)
   
   ! Uptake Targets
   !C
   if (self%resolve_DIC) then
    _SET_ODE_(self%id_DIC,-Aupt%C )
   end if
   !P
   if (org%dop_uptake) then
     _SET_ODE_(self%id_DOP,-Aupt%P)
   else
     _SET_ODE_(self%id_DIP,-Aupt%P) 
   end if
   !N
   _SET_ODE_(self%id_DINO3,-Aupt%NO3) !*org%lim_no3/Alim%N)
   _SET_ODE_(self%id_DINH4,-Aupt%NH4) !*org%lim_nh4/Alim%N) 
   !write(*,*)'upt:',Aupt%NO3+Aupt%NH4,'build:',Aupt%C/self%C2N,'diff:',Aupt%NO3+Aupt%NH4-Aupt%C/self%C2N
   !Si
   if (self%lim_Si) then
     _SET_ODE_(self%id_DISi,-Aupt%Si)
   end if
   !if (self%resolve_carb) then
     !if (self%resolve_cal) then !coccos
     !  _SET_ODE_(self%id_DICcal, -Aupt%Ccal)
     !else 
     !  if (.not. self%resolve_Si) then !non-diatoms
     !    _SET_ODE_(self%id_DICcal, -Aupt%Ccal)
     !  end if !if resolve_Si
     !end if !if resolve_cal
   !end if !if resolve_carb
   
   ! Recycling to Nutrient Pools
   !C
   _SET_ODE_(self%id_DOC,exud%C)
   _SET_ODE_(self%id_det1C,mort%C*(1.0-self%frac_d2x))
   _SET_ODE_(self%id_det2C,mort%C*self%frac_d2x)
   if (self%metCexc .ne. 0) then
     _SET_ODE_(self%id_SOC,exud_soc)
   end if
   !P
   _SET_ODE_(self%id_DOP,exud%P)
   _SET_ODE_(self%id_det1P,mort%P*(1.0-self%frac_d2x))
   _SET_ODE_(self%id_det2P,mort%P*self%frac_d2x)
   !N
   _SET_ODE_(self%id_DON, exud%N)
   _SET_ODE_(self%id_det1N, mort%N*(1.0-self%frac_d2x))
   _SET_ODE_(self%id_det2N, mort%N*self%frac_d2x)
   !Si
   if (self%lim_Si) then
     _SET_ODE_(self%id_dISi,exud%Si)
     _SET_ODE_(self%id_det2Si,mort%Si)
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'IngSiunas', IngSiunas
     !in EH, L1473-1477: todo:ask
     !f_FromTo(i_p1s_d2s)=(f_FromTo(i_p1c_d1c)+f_FromTo(i_p1c_d2c))/rcs
     !dia_adds_loss=(f_FromTo(i_p1c_doc)+f_FromTo(i_p1c_z1c)+f_FromTo(i_p1c_z2c))/rcs
     !dia_ups=f_FromTo(i_n5s_p1s)
     !f_FromTo(i_n5s_p1s)=max(0.0,dia_ups-dia_adds_loss) !MK verstehe ich nicht !JP Vermindere n5s Aufnahme
     !f_FromTo(i_p1s_d2s)=f_FromTo(i_p1s_d2s)+max(0.0,dia_adds_loss-dia_ups)
   end if
   !if (self%resolve_carb) then
     !if (self%resolve_cal) then
      !_SET_ODE_(id_detCcal,mort%Ccal) (p3k_d2k)
     !else
       !if (.not. self%resolve_Si) then
         !p3k_d2k (L.1489)
         !=p2c_d1c + p2c_d2c + p2c_doc/q_c_cal
       !end if
     !end if
   !end if
   
   ! Diagnostics
   !General
   _SET_DIAGNOSTIC_(self%id_Closs,_REPLNAN_((mort%C+exud%C+exud_soc)*s2d))
   _SET_DIAGNOSTIC_(self%id_Ploss,_REPLNAN_((mort%P+exud%P)*s2d))
   _SET_DIAGNOSTIC_(self%id_Nloss,_REPLNAN_((mort%N+exud%N)*s2d))
   if ((self%metIntSt .eq. 0)) then
     !_SET_DIAGNOSTIC_(self%id_QP, org%QP)
     !_SET_DIAGNOSTIC_(self%id_QN, org%QN)
     _SET_DIAGNOSTIC_(self%id_QPr, _REPLNAN_(org%QPr))
     _SET_DIAGNOSTIC_(self%id_QNr, _REPLNAN_(org%QNr))
   else if ((self%metIntSt .eq. 1)) then
     _SET_DIAGNOSTIC_(self%id_QP, _REPLNAN_(org%QP))
     _SET_DIAGNOSTIC_(self%id_QN, _REPLNAN_(org%QN))
     _SET_DIAGNOSTIC_(self%id_QPr, _REPLNAN_(org%QPr))
     _SET_DIAGNOSTIC_(self%id_QNr, _REPLNAN_(org%QNr))
     _SET_DIAGNOSTIC_(self%id_N2P, _REPLNAN_(org%N/org%P))
   end if
   _SET_DIAGNOSTIC_(self%id_Plim,_REPLNAN_(Alim%P))
   _SET_DIAGNOSTIC_(self%id_Nlim,_REPLNAN_(Alim%N))
   !Autotrophy
   _SET_DIAGNOSTIC_(self%id_Ilim,_REPLNAN_(org%fI))
   _SET_DIAGNOSTIC_(self%id_MuClim_A,_REPLNAN_(Aupt%C/org%C*s2d))
   _SET_DIAGNOSTIC_(self%id_NPPR, _REPLNAN_((Aupt%C-exud%C-exud_soc)*s2d))
   _SET_DIAGNOSTIC_(self%id_Cgain_A, _REPLNAN_(Aupt%C*s2d))
   _SET_DIAGNOSTIC_(self%id_Pgain_A, _REPLNAN_(Aupt%P*s2d))
   _SET_DIAGNOSTIC_(self%id_NO3gain_A, _REPLNAN_(Aupt%NO3*s2d))
   _SET_DIAGNOSTIC_(self%id_NH4gain_A, _REPLNAN_(Aupt%NH4*s2d))
   _SET_DIAGNOSTIC_(self%id_exudsoc,_REPLNAN_(exud_soc*s2d))
   if (self%lim_Si) then
     _SET_DIAGNOSTIC_(self%id_Silim,_REPLNAN_(Alim%Si))
   end if
   _SET_DIAGNOSTIC_(self%id_diagChl, _REPLNAN_(org%Chl))
   if (self%metchl .ne. 0) then
     _SET_DIAGNOSTIC_(self%id_QChl, _REPLNAN_(org%QChl))
   end if
   select case (self%metvel)
     case default
       call self%fatal_error('phytoplankton.F90/get_vertical_movement','for '//trim(self%name)// ' specified metvel option is not available')
     case (0)
       vert_vel=self%w
     case (1)
      !vert_vel=self%w * (1-1/(1+exp(10*(.5-min(org%QPr,org%QNr)))))
      vert_vel = -self%w*(0.1_rk+0.9_rk*exp( -5._rk * min(org%QPr,org%QNr)))
     _SET_DIAGNOSTIC_(self%id_sinkvel,_REPLNAN_(vert_vel*s2d))
     !write(*,*)'org%limNP,exp( -5._rk * org%limNP),vert_vel:',org%limNP,exp( -5._rk * org%limNP),vert_vel*s2d
   end select
   !
   !END OF WRITE
   !-------------------------------------------------------------------------
   
   ! Leave spatial loops (if any)
   _LOOP_END_
   
   end subroutine do
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
!
! !INTERFACE:
   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
!
! !INPUT PARAMETERS:
   class (type_gpm_phytoplankton), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
!
! !LOCAL VARIABLES:
   real(rk)                   :: vert_vel,temp,par,I0
   real(rk)                   :: QPr,QNr !,yearday,par
   real(rk), parameter        :: s2d= 86400.
!
!EOP
!-----------------------------------------------------------------------
!BOC

   ! way to get the day of year
   !_GET_DEPENDENCY_SCALAR_(self%id_yearday,yearday)  ! day of year (diff(d/m/Y - 1/1/Y))

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   select case (self%metvel)
     case default
       call self%fatal_error('phytoplankton.F90/get_vertical_movement','for '//trim(self%name)// ' specified metvel option is not available')
     case (0)
       vert_vel=self%w
     case (1)  !sinking velocity decreases with increasing quota
       _GET_ (self%id_QPr_dep,QPr)
       _GET_ (self%id_QNr_dep,QNr)
      !vert_vel=self%w * (1-1/(1+exp(10*(.5-min(QPr,QNr)))))
      vert_vel = -self%w*(0.1_rk+0.9_rk*exp( -5._rk * min(QPr,QNr)))
   end select

   !Set these calculated vertical_movement values
   _SET_VERTICAL_MOVEMENT_(self%id_boundC,vert_vel)
   if (self%metIntSt .eq. 1) then
     _SET_VERTICAL_MOVEMENT_(self%id_boundP,vert_vel)
     _SET_VERTICAL_MOVEMENT_(self%id_boundN,vert_vel)
   end if
   
   !Chlorophyll
   select case (self%metchl) 
     case (20) !for dynamical Chl options of metchl
       _SET_VERTICAL_MOVEMENT_(self%id_boundChl,vert_vel)
   end select
   
   !_SET_DIAGNOSTIC_(self%id_vv , vert_vel*s2d)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_vertical_movement
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_gpm_phytoplankton), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                 :: orgC,orgChl
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Self-shading with explicit contribution from background plankton concentration.
   select case (self%metext)
    case default
     call self%fatal_error('phytoplankton.F90/get_light_extinction','for '//trim(self%name)// ' specified metext option is not available')
    case (0) ! extp_d093
     _GET_(self%id_boundC,orgC) ! carbon biomass
     _SET_EXTINCTION_(self%kc*orgC)
    case (1)
     _GET_(self%id_Chl_dep,orgChl) ! light extinction based on Chlorophyll
     _SET_EXTINCTION_(self%kchl*orgChl)
   end select
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------


end module gpm_phytoplankton
