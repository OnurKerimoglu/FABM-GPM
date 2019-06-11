#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_mixotroph --- generic mixotroph model
!
!TODO:
!
! !INTERFACE:
   module gpm_mixotroph
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
   
   type, extends(type_GPMmixo),public :: type_gpm_mixotroph
      contains
      
      procedure  :: initialize
      procedure :: do
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
   class (type_gpm_mixotroph), intent(inout), target   :: self
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
   call self%get_parameter(self%resolve_Si, 'resolve_Si',   '-',          'whether to resolve Si cycle',          default=.false.)
   call self%get_parameter(self%resolve_DIC,  'resolve_DIC',   '-',     'whether to resolve DIC',          default=.false.)
   call self%get_parameter(self%fracaut, 'fracaut',  '-',  'fraction of autotrophy', default=0.5_rk)
   call self%get_parameter(self%metext,  'metext',   '-',  'method for calculating light ext.',  default=0)
   call self%get_parameter(self%kc,      'kc',       'm^2/mmolC',   'specific light extinction',              default=0.03_rk )
   call self%get_parameter(self%kc,      'kchl',     'm^2/gChl',   'chl-specific light extinction',              default=0.023_rk )
   
   call self%get_parameter(self%metvel,   'metvel',   '-',  'method for calculating settling vel.',  default=0)
   call self%get_parameter(self%w,       'w',        'm/d',        'vertical velocity (<0 for sinking)',     default=0.0_rk, scale_factor=d_per_s)      
   
   call self%get_parameter(self%metTresp, 'metTresp',  '-',   'method to calcualte T resp.',   default=1)
   call self%get_parameter(self%Q10,     'Q10',     '-',    'Q10 for kinetic rates',       default=1.5_rk)
   call self%get_parameter(self%Tref,     'Tref',    'celcius',    'reference temperature for kinetic rates',default=10.0_rk)
   
   call self%get_parameter(self%rmd,     'rmd',     '/d',       'linear mortality rate -min',                  default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmdq,    'rmdq',    '/d',       'quadratic mortality rate',               default=0.001_rk, scale_factor=d_per_s)
   call self%get_parameter(self%frac_d2x, 'frac_d2x',     '-',        'fraction of fast sinking detritus in zoo mortality/excretion',     default=0.15_rk)
   
   call self%get_parameter(self%metIntSt, 'metIntSt', '-', 'method for representing the internal states', default=0)
   call self%get_parameter(self%C2N,     'C2N',     'molC/molN',   'molar C:N ratio',  default=6.625_rk) !rf: 6.625
   call self%get_parameter(self%C2P,     'C2P',     'molC/molP',   'molar C:P ratio',  default=132.5_rk) !rf:106
   call self%get_parameter(self%QPmax,    'QPmax',    'molP/molC',    'Max. P Quota',                 default=0.04_rk) 
   call self%get_parameter(self%QPmin,    'QPmin',    'molP/molC',    'Subsistance P Quota',          default=0.01_rk) 
   call self%get_parameter(self%QNmax,    'QNmax',    'molN/molC',    'Max. N Quota',                 default=0.22_rk) 
   call self%get_parameter(self%QNmin,    'QNmin',    'molN/molC',    'Subsistance N Quota',          default=0.12_rk)
   
   ! Heterotrophic Parameters
   call self%get_parameter(self%rmn,     'rmn',     '/d',       'loss rate to nutrients',                 default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%gmax,        'gmax',        '/d',     'max. grazing rate',                    default=1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Kz,          'Kz',          'mmolC/m3','half saturation of ingestion rate',    default=0.00001_rk)
   call self%get_parameter(self%asefC,       'asefC',       '-',       'assimilation efficiency of the ingested C',default=0.75_rk)
   call self%get_parameter(self%asefP,       'asefP',       '-',       'assimilation efficiency of the ingested P',default=1.0_rk)   
   call self%get_parameter(self%asefN,       'asefN',       '-',       'assimilation efficiency of the ingested N',default=1.0_rk) 
   call self%get_parameter(self%unas_detfrac,'unas_detfrac','-',       'detritus fraction of the unassimmilated Ingestion',default=0.75_rk)
   call self%get_parameter(self%num_prey,    'num_prey',    '-',       'number of prey targets',               default=8)
   call self%get_parameter(self%dynpref,     'dynpref',     '-',       'dynamic preference switch',            default=.false.)
   
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

   ! Register state variables
   call self%register_state_variable(self%id_boundC,'C','mmolC/m^3','bound carbon', & 
                                    minimum=_ZERO_, specific_light_extinction=self%kc,vertical_movement=self%w*d_per_s, &
                                    no_river_dilution=.true.)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_boundC)
   
   if (self%metIntSt .eq. 0) then
     !call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
     !                               minimum=_ZERO_)
     !call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
     !                               minimum=_ZERO_)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundC,scale_factor=1./self%C2P)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundC,scale_factor=1./self%C2N)
   else if (self%metIntSt .eq. 1) then  ! .or. (self%metIntSt .eq. 0) !(for debugging purposes)
     call self%register_state_variable(self%id_boundP,'P','mmolP/m^3','bound phosphorus', & 
                                    minimum=_ZERO_,vertical_movement=self%w*d_per_s,no_river_dilution=.true.)
     call self%register_state_variable(self%id_boundN,'N','mmolN/m^3','bound nitrogen', & 
                                    minimum=_ZERO_,vertical_movement=self%w*d_per_s,no_river_dilution=.true.)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_boundP)                             
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_boundN)
   end if
   
   ! Register links to external models, if variable names are provided in namelist. 
   if (self%fracaut .lt. 1.0) then
    allocate(self%prpar(self%num_prey))
    DO i=1,self%num_prey
     write (istr,'(i0)') i
     call self%get_parameter(self%prpar(i)%pref,'prey'//trim(istr)//'pref','-','preference for prey-'//trim(istr))
     !C
     call self%register_state_dependency(self%prpar(i)%id_C,'prey'//trim(istr)//'C')
     
     ! logic: access coupled the state variable 'prey(i)X', and parameter 'prey(i)C2X'. If prey(i)X is not available, prey(i)C2X will be used to calc X
     !P
     call self%get_parameter(self%prpar(i)%C2P, 'prey'//trim(istr)//'C2P',     '-',       'prey'//trim(istr)//' C:P ratio',  default=-1.0_rk)
     if (self%prpar(i)%C2P .lt. 0.0) then 
       call self%register_state_dependency(self%prpar(i)%id_P,'prey'//trim(istr)//'P','mmolP/m^3','bound phosphorus in prey'//trim(istr), required=.false.)
     end if 
     
     !N
     call self%get_parameter(self%prpar(i)%C2N, 'prey'//trim(istr)//'C2N',     '-',       'prey'//trim(istr)//' C:N ratio',  default=-1.0_rk)
     if (self%prpar(i)%C2N .lt. 0.0) then 
       call self%register_state_dependency(self%prpar(i)%id_N,'prey'//trim(istr)//'N','mmolN/m^3','bound nitrogen in prey'//trim(istr), required=.false.)
     end if
     
     !Chl
     call self%register_state_dependency(self%prpar(i)%id_Chl,'prey'//trim(istr)//'Chl','mg/m^3','bound Chlorophyll in prey'//trim(istr), required=.false.)
     
     if (self%resolve_Si) then
       !when C2Si is not explicitly provided, assume the prey is non-diatom, so Si2C=0 (C2Si=0 will mean the same)
       call self%get_parameter(self%prpar(i)%C2Si,'prey'//trim(istr)//'C2Si','-','C:Si ratio of prey-'//trim(istr),default=0.0_rk)
     end if
    END DO
   end if
   
   ! linking to DIM and POM pools are required for recycling
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
   if (self%resolve_Si) then
     call self%register_state_dependency(self%id_DISi,'DISi')
     call self%register_state_dependency(self%id_det2Si,'det2Si')
   end if
   
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
   call self%register_diagnostic_variable(self%id_NPPR,'NPPR','/d','-NPPR',              &
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
                                           
  select case (self%metchl)
     case default 
       call self%fatal_error('mixotroph.F90/initialize:','for '//trim(self%name)// ' specified metchl option is not available')
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
                                    minimum=_ZERO_,no_river_dilution=.true.)
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

   ! Heterotrophy
   if (self%fracaut .lt. 1.0) then
    do i=1,self%num_prey
     write (istr,'(i0)') i
     call self%register_diagnostic_variable(self%prpar(i)%id_realpref,'real_pref_prey'//trim(istr),'-', 'realized pref for prey'//trim(istr),&
                                          output=output_instantaneous)
    end do
   
    call self%register_diagnostic_variable(self%id_respC,'resp_C','mmolC/m^3/d', 'respiration rate', &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_IngC,'Ing_C','molC/molC/d', 'total sp. C grazing rate', &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_IngunasC,'Ingunas_C','/d',  'unassimilated fraction of sp. ingestion of C', &
                                          output=output_time_step_averaged)                                          
    call self%register_diagnostic_variable(self%id_IngasC,'Ingas_C','/d', 'sp. C gain rate by heterotrophy',   &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_asefC,'asef_C','-',  'assimilation efficiency of C', &
                                          output=output_instantaneous)
   
    call self%register_diagnostic_variable(self%id_IngP,'Ing_P','molP/molC/d', 'total sp. P grazing rate', &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_IngunasP,'Ingunas_P','/d',  'unassimilated fraction of sp. ingestion of P', &
                                          output=output_time_step_averaged)                                       
    call self%register_diagnostic_variable(self%id_IngasP,'Ingas_P','/d', 'sp. P gain rate by herbivory',   &
                                          output=output_time_step_averaged) 
    call self%register_diagnostic_variable(self%id_asefP,'asef_P','-',  'assimilation efficiency of P', &
                                          output=output_instantaneous)

    call self%register_diagnostic_variable(self%id_IngN,'Ing_N','molN/molC/d', 'total sp. N grazing rate', &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_IngunasN,'Ingunas_N','/d',  'unassimilated fraction of sp. ingestion of N', &
                                          output=output_time_step_averaged)                                        
    call self%register_diagnostic_variable(self%id_IngasN,'Ingas_N','/d', 'sp. N gain rate by herbivory',   &
                                          output=output_time_step_averaged)
    call self%register_diagnostic_variable(self%id_asefN,'asef_N','-',  'assimilation efficiency of N', &
                                          output=output_instantaneous)
    if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_IngunasSi,'Ingunas_Si','/d',  'unassimilated fraction of sp. ingestion of Si', &
                                          output=output_time_step_averaged) 
    end if
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
   class (type_gpm_mixotroph), intent(in)     :: self 
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   logical                    :: debw
   real(rk)                   :: fT,gmax,exud_soc,vert_vel
   type (type_env)            :: env
   type (type_elms)           :: Ing,asef,Ingas,Ingunas,excr,mort
   type (type_elms)           :: dom,Alim,exud
   type (type_dim)            :: di,Aupt
   type (org_mixotrophic)     :: org !quantities of the organismal system being modelled (state and derived values)
   integer                    :: i
   character(len=2)           :: istr !i converted to character string
   type(prey_data)            :: prdat
   if (self%fracaut .lt. 1.0) then
    allocate(prdat%corpref(self%num_prey))
    allocate(prdat%rpref(self%num_prey))
    allocate(prdat%weight(self%num_prey))
    allocate(prdat%C(self%num_prey))
    allocate(prdat%P(self%num_prey))
    allocate(prdat%N(self%num_prey))
    allocate(prdat%Chl(self%num_prey))
    allocate(prdat%Si(self%num_prey))
    allocate(prdat%grC(self%num_prey))
    allocate(prdat%grP(self%num_prey))
    allocate(prdat%grN(self%num_prey))
    allocate(prdat%grChl(self%num_prey))
    allocate(prdat%grSi(self%num_prey))
    allocate(prdat%Qr(self%num_prey))
    allocate(prdat%QPr(self%num_prey))
    allocate(prdat%QNr(self%num_prey))
   end if
   
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
   
   !Reset the rates
   !General
   mort%C=0.0;mort%P=0.0;mort%N=0.0;mort%Si=0.0
   !Heterotrophy:
   Ing%C=0.0;Ing%P=0.0;Ing%N=0.0;Ing%Si=0.0
   Ingas%C=0.0;Ingas%P=0.0;Ingas%N=0.0;Ingas%Si=0.0
   Ingunas%C=0.0;Ingunas%P=0.0;Ingunas%N=0.0;Ingunas%Si=0.0
   excr%C=0.0;excr%P=0.0;excr%N=0.0;excr%Si=0.0
   !Autotrophy:
   Aupt%C=0.0; Aupt%P=0.0; Aupt%NH4=0.0; Aupt%NO3=0.0; Aupt%Si=0.0
   exud%C=0.0; exud%P=0.0; exud%N=0.0; exud%Si=0.0; exud_soc=0.0
    
   
   !-------------------------------------------------------------------------
   !PREPARE: retrieve variables
   !
   !_GET_GLOBAL_(self%id_doy,doy) !day of year
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
   
   !LOGICALS FOR DEBUGGING 
   if (debug .and. (env%depth .lt. 1.0 .and. env%depth .gt. 0.0)) then
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
   
   select case (self%metchl) 
     case (20) !for dynamical Chl options of metchl
       _GET_(self%id_boundChl,org%Chl)
       org%QChl =  org%Chl / (org%C * 12.0_rk)
       !gChl/gC = mgChl/m3 / (mmolC/m3 * 12mgC/mmolC)
       !write(*,*)'phyto L411. Chl,QChl:',org%Chl,org%QChl
   end select
   
   ! collect resources
   ! if (self%resolve_DIC) then
   !_GET_(self%id_DIC,di%C) !
   ! end if
   _GET_(self%id_DIP,di%P)
   if (self%dop_allowed) then
     _GET_(self%id_DOP,dom%P)
   end if
   _GET_(self%id_DINO3,di%NO3) ! nitrate
   _GET_(self%id_DINH4,di%NH4) ! ammonium
   
   !collect prey abundances
   if (self%fracaut .lt. 1.0) then
    !reset prdat%X, prdat%QX
    prdat%C=0.0
    prdat%P=0.0
    prdat%N=0.0
    prdat%Chl=0.0
    prdat%Si=0.0
    prdat%QPr=0.0
    prdat%QNr=0.0
    prdat%Qr=0.0 !this is min(QPr,QNr)
    DO i=1,self%num_prey
     write (istr,'(i0)') i
     !C
     _GET_STATE_(self%prpar(i)%id_C,prdat%C(i))
     if (prdat%C(i) .lt. TINYPREYC) then
       prdat%C(i)=0.0 ! prdat%C(i)=0 will spare the prey
     end if
     !P
     if (_AVAILABLE_(self%prpar(i)%id_P)) then
       _GET_STATE_(self%prpar(i)%id_P,prdat%P(i))
       if (prdat%P(i) .lt. TINYPREYC/106._rk) then
           prdat%C(i)=0.0 ! prdat%C(i)=0 will spare the prey
         end if
     else
       if (self%prpar(i)%C2P .lt. 0.0_rk) then
         call self%fatal_error('gpm_mixotroph_do','for prey'//trim(istr)//', no explicit P variable was coupled, or no (>0) C2P par was provided')
       end if
       prdat%P(i)=prdat%C(i) / self%prpar(i)%C2P
     end if
     !N
     if (_AVAILABLE_(self%prpar(i)%id_N)) then
       _GET_STATE_(self%prpar(i)%id_N,prdat%N(i))
       if (prdat%N(i) .lt. TINYPREYC*16._rk/106._rk) then
         prdat%C(i)=0.0 ! prdat%C(i)=0 will spare the prey
       end if
     else
       if (self%prpar(i)%C2N .lt. 0.0_rk) then
         call self%fatal_error('gpm_mixotroph_do','for prey'//trim(istr)//', no explicit N variable was coupled, or no (>0) C2N par was provided')
       end if
       prdat%N(i)=prdat%C(i) / self%prpar(i)%C2N
     end if
     
     !if pref(i)<0, then pref=f(Qr(i))
     if (self%prpar(i)%pref .lt. 0.0) then !default=0.0, which means no Si in prey
       _GET_(self%prpar(i)%id_QPr,prdat%QPr(i))
       _GET_(self%prpar(i)%id_QNr,prdat%QNr(i))
       prdat%Qr(i)=min(prdat%QPr(i),prdat%QNr(i))
     end if

     !Chl:
     prdat%Chl(i)=0.0 !if not explicitly provided, is not needed anyway
     if (_AVAILABLE_(self%prpar(i)%id_Chl)) then
       _GET_STATE_(self%prpar(i)%id_Chl,prdat%Chl(i))
       if (prdat%Chl(i) .lt. TINYPREYC/20._rk) then
         prdat%C(i)=0.0 ! prdat%C(i)=0 will spare the prey
       end if
     end if
     
     !Si:
     prdat%Si(i)=0.0_rk
     if (self%resolve_Si) then
       if (self%prpar(i)%C2Si .gt. 0.0) then
         prdat%Si(i)=prdat%C(i)/self%prpar(i)%C2Si
       end if
       !write(*,*),'prey#',i,' c2si,C,Si:',self%prpar(i)%C2Si,prdat%C(i),prdat%Si(i)
     end if
    END DO
   end if
   !
   !END OF PREPARE
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !START OF CALCULATIONS: no other calculation outside this box
   !
   !Temperature function
   fT = get_fT(self%type_GPMbase,env%temp)
   
   ! Autotrophy
   if (self%fracaut .gt. 0.0) then
     !Nutrient quotas, limitation and uptake

     call org%get_nut_QLU(self%type_GPMaut,fT,di,dom,Alim,Aupt)
     !write(*,*)'N: lim%P,lim%N,upt%P,upt%N:',Alim%P,Alim%N !,Aupt%P*s2d,Aupt%N*s2d
   
     ! Light limitation
     !Chl is potentially needed to calculate light limitation
     select case (self%metchl) 
      case (0,1,2) !for diagnostic Chl options of metchl
       call org%get_chl(self%type_GPMaut,env)
      end select
     call org%get_fI(self%type_GPMaut,env)
     !if (env%par .gt. 0.0) write(*,*)'vc,is,io,im,fT,p',vC_atQmax,self%islope,self%Iopt,self%Imin,fT,env%par
     
     !Primary production (upt%C) and  uptake terms for elements resolved as fixed stoichiometry
     call org%get_PP_upt4fs(self%type_GPMaut,fT,Aupt,exud_soc)
     
     !modify the autotrophic gain/loss rates by fracaut
     Aupt%C=self%fracaut*Aupt%C
     Aupt%P=self%fracaut*Aupt%P
     Aupt%NO3=self%fracaut*Aupt%NO3
     Aupt%NH4=self%fracaut*Aupt%NH4
     
     !obtain the chl. synth rates, of Chl is dynamic (if not, the subr will do nothing)
     call org%get_chldyn(self%type_GPMaut,env,Alim,Aupt)
     
     !Aut. Losses
     call org%get_losses_exud(self%type_GPMaut,Aupt,exud)
   end if !else: all the default values (Aupt,exud) were set to 0 (above)

   !Heterotrophy
   if (self%fracaut .lt. 1.0) then
     if (self%fracaut .lt. 0.001) then !i.e., if purely heterotroph, get_nut_QLU was not already called
       !Nutrient quotas
       call org%get_nut_Q(self%type_GPMbase)
     end if
     !calculate grazing rate for each prey unit and total grazing
     
     !modify the heterotrophic gain rates (=prop to gmax) by (1-fracaut)
     gmax=self%gmax*(1-self%fracaut)
     
     call org%get_GronMultiPrey(self%type_GPMmixo,self%prpar,gmax,fT,prdat,Ing)
     !this gives grazrateX [molXprey/molXpred/d]
     
     !adjust the assimilation efficiencies to maintain  Qmin<Q<Qmax
     call org%get_adjasef(self%type_GPMmixo,Ing,asef,Ingas,Ingunas)
   
     ! Het. Losses
     call org%get_losses_excr(self%type_GPMmixo,fT,excr)
   end if !else: all the default values (Ing, Inas, Ingunas,excr) were set to 0 (above)
   
   ! General Losses
   call org%get_losses_mort(self%type_GPMaut,self%fracaut,fT,mort)  
   !
   !END OF CALCULATIONS: no other calculation outside this box
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
   !WRITE
   !
   ! SET RHS 
   !C
   _SET_ODE_(self%id_boundC,   Aupt%C - exud%C - exud_soc + Ingas%C - excr%C - mort%C)
   !P,N
   if (self%metIntSt .eq. 1) then
     _SET_ODE_(self%id_boundP, Aupt%P - exud%P + Ingas%P - excr%P - mort%P) 
     _SET_ODE_(self%id_boundN, Aupt%NO3+Aupt%NH4 - exud%N + Ingas%N - excr%N - mort%N)
   end if
   
   select case (self%metchl)
     case (20) !list here the dynamic Chlorophyll options
       _SET_DIAGNOSTIC_(self%id_chlrho,org%chlrho)
       _SET_ODE_(self%id_boundChl, org%chlsynth - (exud%C + exud_soc + mort%C)*12.0_rk*org%QChl)
       !                             mgChl/m3/s - (mmolC/m3/s)*12.0mgC/mgChl *mgChl/mgC
       !write(*,*)'org%chlsynth', org%chlsynth,'loss:', - (exud%C + exud_soc + mort%C)*org%QChl,'RHS:',org%chlsynth - (exud%C + exud_soc + mort%C)*org%QChl
   end select
   
   !grazing targets
   if (self%fracaut .lt. 1.0) then
    DO i=1,self%num_prey
     !C
     _SET_ODE_(self%prpar(i)%id_C,-prdat%grC(i)*org%C) !molC/molC/d *molC/m3 =molC/m3/d !*(1.0-self%fracaut)
     !P
     if (_AVAILABLE_(self%prpar(i)%id_P)) then
       _SET_ODE_(self%prpar(i)%id_P,-prdat%grP(i)*org%C) !molP/molC/d *molC/m3 =molP/m3/d
     end if
     !N
     if (_AVAILABLE_(self%prpar(i)%id_N)) then
       _SET_ODE_(self%prpar(i)%id_N,-prdat%grN(i)*org%C) !molN/molC/d *molC/m3 =molP/m3/d
     end if
     !Chl
     if (_AVAILABLE_(self%prpar(i)%id_Chl)) then
       _SET_ODE_(self%prpar(i)%id_Chl,-prdat%grChl(i)*org%C) !molChl/molC/d *molC/m3 =molChl/m3/d
     end if
    END DO
   end if
   
   !O2
   _SET_ODE_(self%id_O2,Aupt%C-excr%C)
   
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
   
   !Recycling to nutrient pools
   !C
   if (self%resolve_DIC) then
     _SET_ODE_(self%id_DIC,excr%C)
   end if
   _SET_ODE_(self%id_DOC,exud%C + Ingunas%C*(1.-self%unas_detfrac))
   _SET_ODE_(self%id_det1C,(Ingunas%C*self%unas_detfrac+mort%C)*(1.-self%frac_d2x))
   _SET_ODE_(self%id_det2C,(Ingunas%C*self%unas_detfrac+mort%C)*self%frac_d2x)
   if (self%metCexc .ne. 0) then
     _SET_ODE_(self%id_SOC,exud_soc)
   end if
   !P
   _SET_ODE_(self%id_DIP,excr%P)
   _SET_ODE_(self%id_DOP, exud%P + Ingunas%P*(1.-self%unas_detfrac)) 
   _SET_ODE_(self%id_det1P,(Ingunas%P*self%unas_detfrac+mort%P)*(1.-self%frac_d2x))
   _SET_ODE_(self%id_det2P,(Ingunas%P*self%unas_detfrac+mort%P)*self%frac_d2x)
   !if (debw) write(*,'(A, 3F14.10)') 'rmd+rmdq,Ingunas%C,Ingunas%P', rmd+rmdq,Ingunas%C,Ingunas%P
   !N
   _SET_ODE_(self%id_DINH4,excr%N)
   _SET_ODE_(self%id_DON,exud%N + Ingunas%N*(1.-self%unas_detfrac))
   _SET_ODE_(self%id_det1N,(Ingunas%N*self%unas_detfrac+mort%N)*(1.-self%frac_d2x))
   _SET_ODE_(self%id_det2N,(Ingunas%N*self%unas_detfrac+mort%N)*self%frac_d2x)
   !if (debw) write(*,'(2A, 2F14.10)') self%name, 'Ingunas%P', Ingunas%P*s2d
   !Si
   if (self%resolve_Si) then
     _SET_ODE_(self%id_DISi,Ingunas%Si*(1.-self%unas_detfrac))
     _SET_ODE_(self%id_det2Si,Ingunas%Si*self%unas_detfrac)
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'Ingunas%Si', Ingunas%Si
   end if
   
   ! SET DIAGNOSTICS
   !General
   _SET_DIAGNOSTIC_(self%id_Closs,(mort%C+excr%C+exud%C)*s2d)
   _SET_DIAGNOSTIC_(self%id_Ploss,(mort%P+excr%P+exud%P)*s2d)
   _SET_DIAGNOSTIC_(self%id_Nloss,(mort%N+excr%N+exud%N)*s2d)
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
   !Autotrophy
   if (self%fracaut .gt. 0.0) then
    _SET_DIAGNOSTIC_(self%id_Plim,Alim%P)
    _SET_DIAGNOSTIC_(self%id_Nlim,Alim%N)
    _SET_DIAGNOSTIC_(self%id_Ilim,org%fI)
    _SET_DIAGNOSTIC_(self%id_MuClim_A,Aupt%C/org%C*s2d)
    _SET_DIAGNOSTIC_(self%id_NPPR, (Aupt%C-exud%C-exud_soc)*s2d)
    _SET_DIAGNOSTIC_(self%id_Cgain_A, Aupt%C*s2d)
    _SET_DIAGNOSTIC_(self%id_Pgain_A, Aupt%P*s2d)
    _SET_DIAGNOSTIC_(self%id_NO3gain_A, Aupt%NO3*s2d)
    _SET_DIAGNOSTIC_(self%id_NH4gain_A, Aupt%NH4*s2d)
    _SET_DIAGNOSTIC_(self%id_exudsoc,exud_soc*s2d)
    _SET_DIAGNOSTIC_(self%id_diagChl, org%Chl)
    if (self%metchl .ne. 0) then
      _SET_DIAGNOSTIC_(self%id_QChl, org%QChl)
    end if
   end if
   select case (self%metvel)
     case default
       call self%fatal_error('phytoplankton.F90/get_vertical_movement','for '//trim(self%name)// ' specified metvel option is not available')
     case (0)
       vert_vel=self%w
     case (1)
      !vert_vel=self%w * (1-1/(1+exp(10*(.5-Qr))))
      vert_vel = -self%w*(0.1_rk+0.9_rk*exp( -5._rk * min(org%QPr,org%QNr)))
     _SET_DIAGNOSTIC_(self%id_sinkvel,vert_vel*s2d)
   end select
   
   !Heterotrophy
   if (self%fracaut .lt. 1.0) then
    DO i=1,self%num_prey
     _SET_DIAGNOSTIC_(self%prpar(i)%id_realpref,prdat%rpref(i))	
    END DO
    _SET_DIAGNOSTIC_(self%id_respC,excr%C*s2d)
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
   end if
   !
   !END OF WRITE
   !-------------------------------------------------------------------------
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
   
   end subroutine do
!EOC
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
!
! !INTERFACE:
   subroutine get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)
!
! !INPUT PARAMETERS:
   class (type_gpm_mixotroph), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
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
   _FABM_LOOP_BEGIN_

   select case (self%metvel)
     case default
       call self%fatal_error('mixotroph.F90/get_vertical_movement','for '//trim(self%name)// ' specified metvel option is not available')
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
   _FABM_LOOP_END_

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
   class (type_gpm_mixotroph), intent(in)     :: self
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
     call self%fatal_error('mixotroph.F90/get_light_extinction','for '//trim(self%name)// ' specified metext option is not available')
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


end module gpm_mixotroph
