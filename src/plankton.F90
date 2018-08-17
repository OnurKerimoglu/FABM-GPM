#include "fabm_driver.h"
#include "fabm.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_plankton --- HZG plankton (phyto,zoo,mixo) model
! The parameter 'fracaut' controls the fraction of autotorophy 
! (1.0 = pure autotrophy, 0.0 = pure heterotrophy) 
!
! !INTERFACE:
   module gpm_plankton
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none
   
!  default: all is private.
   private
   
! !PRIVATE DATA MEMBERS:
!
! !PUBLIC DERIVED TYPES:
   type prey_data
      real(rk),dimension(:),allocatable :: corpref
      real(rk),dimension(:),allocatable :: rpref
      real(rk),dimension(:),allocatable :: C
      real(rk),dimension(:),allocatable :: P      
      real(rk),dimension(:),allocatable :: N
      real(rk),dimension(:),allocatable :: grC
      real(rk),dimension(:),allocatable :: grP
      real(rk),dimension(:),allocatable :: grN
      real(rk),dimension(:),allocatable :: grSi
      real(rk),dimension(:),allocatable :: Qr
      real(rk),dimension(:),allocatable :: QPr
      real(rk),dimension(:),allocatable :: QNr
   end type
   
   type prey_pars
      type (type_state_variable_id)      :: id_C,id_P,id_N
      type (type_dependency_id)          :: id_QPr,id_QNr
      type (type_diagnostic_variable_id) :: id_realpref
      real(rk)                           :: pref,Si2C
   end type
   
   type, extends(type_base_model),public :: type_gpm_plankton
      !  Variable identifiers
      !state vars
      type (type_state_variable_id)      :: id_plaC,id_detC
      type (type_state_variable_id)      :: id_plaP,id_DIP,id_detP
      type (type_state_variable_id)      :: id_plaN,id_DIN,id_detN
      type (type_state_variable_id)      :: id_DISi,id_detSi
      !diagnostics
      type (type_diagnostic_variable_id) :: id_dPAR,id_NPPR,id_rmd
      !type (type_diagnostic_variable_id) :: id_diag1, id_diag2, id_diag3, id_diag4, id_diag5
      type (type_diagnostic_variable_id) :: id_Cgain,id_Closs,id_MuClim,id_Cgain_H,id_Cgain_A,id_IngC,id_asefC,id_IngCunas,id_MuClim_A
      type (type_diagnostic_variable_id) :: id_Pgain,id_Ploss,id_Plim,id_Pgain_H,id_Pgain_A,id_IngP,id_asefP,id_IngPunas,id_QP,id_QPr
      type (type_diagnostic_variable_id) :: id_Ngain,id_Nloss,id_Nlim,id_Ngain_H,id_Ngain_A,id_IngN,id_asefN,id_IngNunas,id_QN,id_QNr,id_N2P
      type (type_diagnostic_variable_id) :: id_Silim,id_IngSiunas
      !depndencies
      type (type_dependency_id)          :: id_QPr_dep,id_QNr_dep
      type (type_dependency_id)          :: id_par,id_depth,id_temp
      type (type_global_dependency_id)   :: id_doy
      type (type_horizontal_dependency_id)::id_I_0

!     Model parameters
      !general
      real(rk) :: fracaut,kc,w,Isc,Tsc,m0,rmn,rmd,rmdq,rmdH,rmdfZ_offset
      real(rk) :: QPmax,QPmin,QNmax,QNmin  
      logical  :: dynpref,resolve_P,resolve_N,resolve_Si
      integer  :: tempresp,lightresp,velmet
      !autotrophic
      real(rk) :: rmax,Q10,Tref,Topt,Tint,Ki,islope,Iopt,Imin,Kp,VPmax,Kn,VNmax,Ksi,Si2C
      integer  :: num_prey
      !heteroptrophic
      real(rk) :: Q10H,TrefH,gmax,Kz,asefC,asefP,asefN,unas_detfrac
      TYPE(prey_pars),dimension(:),allocatable :: prey      
      
      contains
      
      procedure  :: initialize
      procedure  :: do
      procedure  :: get_vertical_movement
      procedure  :: get_light_extinction
      
   end type
   
   !Local variables
   logical :: debug = .false.
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
   class (type_gpm_plankton), intent(inout), target   :: self
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
   call self%get_parameter(self%resolve_P,  'resolve_P',   '-',           'whether to resolve P cycle',          default=.true.) 
   call self%get_parameter(self%resolve_N,  'resolve_N',   '-',           'whether to resolve N cycle',          default=.false.) 
   call self%get_parameter(self%resolve_Si, 'resolve_Si',   '-',          'whether to resolve Si cycle',          default=.false.)
   call self%get_parameter(self%fracaut, 'fracaut',  '-',           'Fraction of Autotrophs',                default=0.5_rk)   
   call self%get_parameter(self%m0,      'm0',       'mmol/m^3',    'background concentration ',              default=0.0225_rk)
   call self%get_parameter(self%kc,      'kc',       'm^2/mmolC',   'specific light extinction',              default=0.03_rk )
   call self%get_parameter(self%velmet,  'velmet',   '-',           'velocity method',                        default=1)
   call self%get_parameter(self%w,       'w',        'm/d',       'vertical velocity (<0 for sinking)',     default=0.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Tsc,     'Tsc',      '-',           'scaling factor for w=f(T)',             default=8.0_rk ) !when velmet=3
   call self%get_parameter(self%Isc,     'Isc',      '-',           'scaling factor for w=f(O)',             default=10.0_rk) !when velmet=4
   call self%get_parameter(self%rmn,     'rmn',     '/d',       'loss rate to nutrients',                 default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmd,     'rmd',     '/d',       'linear mortality rate -min',                  default=0.05_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmdH,     'rmdH',     '/d',       'linear mortality rate -max',                default=0.2_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%rmdfZ_offset,     'rmdfZ_offset',     'm',       'off-set factor for rmd=fz',    default=70.0_rk)
   call self%get_parameter(self%rmdq,    'rmdq',    '/d',       'quadratic mortality rate',               default=0.001_rk, scale_factor=d_per_s)
   call self%get_parameter(self%tempresp, 'tempresp', '-',          'temperature response',                default=1)
   call self%get_parameter(self%Q10,      'Q10',      '-',          'Q10 for autotrophic processes',       default=2.0_rk)    
   call self%get_parameter(self%Tref,     'Tref',     'celcius',    'reference temperature for aut. proc.',default=20.0_rk)
   call self%get_parameter(self%Topt,     'Topt',     'celcius',    'optimal temperature',                 default=25._rk) !when tempresp=3|4 or velmet=3
   call self%get_parameter(self%Tint,     'Tint',     'celcius',    'tempresp=3 parameter',                default=15._rk) !when tempresp=3
   call self%get_parameter(self%QPmax,    'QPmax',    'molP/molC',    'Max. P Quota',                 default=0.04_rk) 
   call self%get_parameter(self%QPmin,    'QPmin',    'molP/molC',    'Subsistance P Quota',          default=0.01_rk) 
   call self%get_parameter(self%QNmax,    'QNmax',    'molN/molC',    'Max. N Quota',                 default=0.22_rk) 
   call self%get_parameter(self%QNmin,    'QNmin',    'molN/molC',    'Subsistance N Quota',          default=0.12_rk) 
   ! AUTOTROPHIC PARAMATERS (!when fracaut>0)
   call self%get_parameter(self%rmax,     'rmax',      '/d',       'maximum specific growth rate',          default=1.0_rk,  scale_factor=d_per_s) !when fracaut>0
   call self%get_parameter(self%lightresp,'lightresp', '-',         'light response',                          default=1) 
   call self%get_parameter(self%islope,   'islope',    'mmolC/m^2/W',  'slope of the P-I curve',                  default=0.05_rk) !when lightresp=1 & 3
   call self%get_parameter(self%Iopt,     'Iopt',      'W/m^2',  'half-saturation nutrient concentration',  default=100._rk) !when lightresp=2 & 4
   call self%get_parameter(self%Imin ,    'Imin',      'W/m^2',     'min. light intensity to be adjusted',     default=25._rk) !when lightresp=5
   call self%get_parameter(self%Kp,       'Kp',        'mmolP/m^3',   'half-saturation P concentration',  default=0.3_rk)
   call self%get_parameter(self%VPmax,    'VPmax',     'mmolP/mmolC/d','Max. phosphorus uptake rate',          default=1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Kn,       'Kn',        'mmolN/m^3',   'half-saturation N concentration',  default=4.2_rk)
   call self%get_parameter(self%VNmax,    'VNmax',     'mmolN/mmolC/d','Max. nitrogen uptake rate',          default=14.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Ksi,      'Ksi',      'mmolSi/m^3',   'half-saturation Si concentration',  default=0.5_rk)
   call self%get_parameter(self%Si2C,     'Si2C',     'molSi/molC',   'molar Si:C ratio',  default=0.0_rk) !15/106=0.1415 (Brzezinski, 1985)
   ! HETEROTROPHIC PARAMETERS
   call self%get_parameter(self%Q10H,        'Q10H',        '-',       'Q10 for heterotrophic processes',      default=2.0_rk)    
   call self%get_parameter(self%TrefH,       'TrefH',       'celcius', 'reference temperature for het. proc.', default=20.0_rk)
   call self%get_parameter(self%gmax,        'gmax',        '/d',     'max. grazing rate',                    default=1.0_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Kz,          'Kz',          'mmolC/m3','half saturation of ingestion rate',    default=0.00001_rk)
   call self%get_parameter(self%asefC,       'asefC',       '-',       'assimilation efficiency of the ingested C',default=0.75_rk)
   call self%get_parameter(self%asefP,       'asefP',       '-',       'assimilation efficiency of the ingested P',default=1.0_rk)   
   call self%get_parameter(self%asefN,       'asefN',       '-',       'assimilation efficiency of the ingested N',default=1.0_rk) 
   call self%get_parameter(self%unas_detfrac,'unas_detfrac','-',       'detritus fraction of the unassimmilated Ingestion',default=0.75_rk)
   call self%get_parameter(self%num_prey,    'num_prey',    '-',       'number of prey targets',               default=8)
   call self%get_parameter(self%dynpref,     'dynpref',     '-',       'dynamic preference switch',            default=.false.)
   
   ! Register state variables
   call self%register_state_variable(self%id_plaC,'C','mmolC/m^3','bound carbon', & 
                                    minimum=_ZERO_, specific_light_extinction=self%kc,vertical_movement=self%w*d_per_s)
   !P
   if (self%resolve_P) then                                    
     call self%register_state_variable(self%id_plaP,'P','mmolP/m^3','bound phosphorus', & 
                                    minimum=_ZERO_, vertical_movement=self%w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_plaP)
   end if
   !N
   if (self%resolve_N) then                                    
     call self%register_state_variable(self%id_plaN,'N','mmolN/m^3','bound nitrogen', & 
                                    minimum=_ZERO_, vertical_movement=self%w*d_per_s)
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_plaN)
   end if
   !Si
   if (self%resolve_Si) then
     !Si stored in phytoplankton is not explicitly resolved
     if (self%fracaut .eq. 1.0) then
       call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_plaC,scale_factor=self%Si2C)
     end if
   end if
   
   
!    !phyto-P, zoo-P and mixo-P
!    if (self%fracaut .gt. 0.9) then
!      !C
!      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_phyto_carbon',units='mmolC/m**3', conserved=.true.),self%id_plaC)
!      !P
!      if (self%resolve_P) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_phyto_phosphorus',units='mmolP/m**3', conserved=.true.),self%id_plaP)
!      end if
!      !N
!      if (self%resolve_N) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_phyto_nitrogen',units='mmolN/m**3', conserved=.true.),self%id_plaN)
!      end if 
!    else if (self%fracaut .lt. 0.1) then
!      !C
!      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_zoo_carbon',units='mmolC/m**3', conserved=.true.),self%id_plaC)
!      !P
!      if (self%resolve_P) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_zoo_phosphorus',units='mmolP/m**3', conserved=.true.),self%id_plaP)
!      end if
!      !N
!      if (self%resolve_N) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_zoo_nitrogen',units='mmolN/m**3', conserved=.true.),self%id_plaN)
!      end if
!    else 
!      !C
!      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_mixo_carbon',units='mmolC/m**3', conserved=.true.),self%id_plaC)
!      !P
!      if (self%resolve_P) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_mixo_phosphorus',units='mmolP/m**3', conserved=.true.),self%id_plaP)
!      end if
!      !P
!      if (self%resolve_N) then
!        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_mixo_nitrogen',units='mmolN/m**3', conserved=.true.),self%id_plaN)
!      end if
!    end if
   
   ! Register links to external models, if variable names are provided in namelist. 
   if (self%fracaut .lt. 0.999) then
    allocate(self%prey(self%num_prey))
    DO i=1,self%num_prey
     write (istr,'(i0)') i
     call self%get_parameter(self%prey(i)%pref,'prey'//trim(istr)//'pref','-','preference for prey-'//trim(istr))
     !C
     call self%register_state_dependency(self%prey(i)%id_C,'prey'//trim(istr)//'C')
     !P
     if (self%resolve_P) then
       call self%register_state_dependency(self%prey(i)%id_P,'prey'//trim(istr)//'P')
       if (self%prey(i)%pref .lt. 0.0) then !if pref<0, real pref will be f(QPr(i)
         call self%register_dependency(self%prey(i)%id_QPr,'prey'//trim(istr)//'QPr')
       end if
    end if
    !N
    if (self%resolve_N) then
       call self%register_state_dependency(self%prey(i)%id_N,'prey'//trim(istr)//'N')
       if (self%prey(i)%pref .lt. 0.0) then !if pref<0, real pref will be f(QPr(i)
         call self%register_dependency(self%prey(i)%id_QNr,'prey'//trim(istr)//'QNr')
       end if
    end if
    !Si
    if (self%resolve_Si) then
       !when Si2C is not explicitly provided, assume the prey is non-diatom, so Si2C=0
       call self%get_parameter(self%prey(i)%Si2C,'prey'//trim(istr)//'Si2C','-','Si:C ratio of prey-'//trim(istr),default=0.0_rk)
    end if
    END DO
   end if 
   
   ! linking to DIM and POM pools are required in any case
   !C
   call self%register_state_dependency(self%id_detC,'detC')
   !P
   if (self%resolve_P) then
     call self%register_state_dependency(self%id_DIP,'DIP')
     call self%register_state_dependency(self%id_detP,'detP')
   end if
   !P
   if (self%resolve_N) then
     call self%register_state_dependency(self%id_DIN,'DIN')
     call self%register_state_dependency(self%id_detN,'detN')
   end if
   !P
   if (self%resolve_Si) then
     call self%register_state_dependency(self%id_DISi,'DISi')
     call self%register_state_dependency(self%id_detSi,'detSi')
   end if 
   
   !Diagnostics
    
   !general
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m^2','-PAR',    &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rmd,'rmd','/d', ' specific linear mort. rate',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Cgain,'Cgain','mmolC/m^3/d', ' bulk C gain rate',   &
                                          output=output_time_step_averaged)  
   call self%register_diagnostic_variable(self%id_Closs,'Closs','mmolC/m^3/d', ' bulk C loss rate',   &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_MuClim,'MuClim','/d', 'C limited sp. growth rate',  &
                                          output=output_time_step_averaged)
   if (self%resolve_P) then                                       
   call self%register_diagnostic_variable(self%id_QP,'QP','molP/molC', 'molar P:C ratio',         &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_QPr,'QPr','-', '(QP-QPmin)/(QPmax-QPmin)',         &
                                          output=output_time_step_averaged) 
     call self%register_diagnostic_variable(self%id_Pgain,'Pgain','mmolP/m^3/d', ' bulk P gain rate',   &
                                          output=output_time_step_averaged)  
     call self%register_diagnostic_variable(self%id_Ploss,'Ploss','mmolP/m^3/d', ' bulk P loss rate',   &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_Plim,'Plim','-', 'P limitation',  &
                                          output=output_time_step_averaged)
   end if
   if (self%resolve_N) then
     call self%register_diagnostic_variable(self%id_QN,'QN','molN/molC', 'molar N:C ratio',         &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_QNr,'QNr','-', '(QN-QNmin)/(QNmax-QNmin)',         &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_Ngain,'Ngain','mmolN/m^3/d', ' bulk N gain rate',   &
                                          output=output_time_step_averaged)  
     call self%register_diagnostic_variable(self%id_Nloss,'Nloss','mmolN/m^3/d', ' bulk N loss rate',   &
                                          output=output_time_step_averaged) 
     call self%register_diagnostic_variable(self%id_Nlim,'Nlim','-', 'N limitation',  &
                                          output=output_time_step_averaged)
   end if
   if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_Silim,'Silim','-', 'Si limitation',  &
                                          output=output_time_step_averaged)
   end if
   if (self%resolve_N .and. self%resolve_P) then 
     call self%register_diagnostic_variable(self%id_N2P,'N2P','molN/molP', 'molar N:P ratio',         &
                                          output=output_time_step_averaged)
   end if
     
   !heterotrophy
   if (self%fracaut .lt. .999) then                                          
     do i=1,self%num_prey
       write (istr,'(i0)') i
       call self%register_diagnostic_variable(self%prey(i)%id_realpref,'real_pref_prey'//trim(istr),'-', 'realized pref for prey'//trim(istr),&
                                          output=output_time_step_averaged)
     end do
     call self%register_diagnostic_variable(self%id_IngC,'IngC','molC/molC/d', 'total sp. C grazing rate', &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_IngCunas,'IngCunas','/d',  'unassimilated fraction of sp. ingestion of C', &
                                          output=output_time_step_averaged)                                          
     call self%register_diagnostic_variable(self%id_Cgain_H,'CgainH','/d', 'sp. C gain rate by heterotrophy',   &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_asefC,'asefC','-',  'assimilation efficiency of C', &
                                          output=output_time_step_averaged)
     if (self%resolve_P) then
       call self%register_diagnostic_variable(self%id_IngP,'IngP','molP/molC/d', 'total sp. P grazing rate', &
                                          output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_IngPunas,'IngPunas','/d',  'unassimilated fraction of sp. ingestion of P', &
                                          output=output_time_step_averaged)                                       
       call self%register_diagnostic_variable(self%id_Pgain_H,'PgainH','/d', 'sp. P gain rate by herbivory',   &
                                          output=output_time_step_averaged) 
       call self%register_diagnostic_variable(self%id_asefP,'asefP','-',  'assimilation efficiency of P', &
                                          output=output_time_step_averaged)
     end if
     if (self%resolve_N) then
       call self%register_diagnostic_variable(self%id_IngN,'IngN','molN/molC/d', 'total sp. N grazing rate', &
                                          output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_IngNunas,'IngNunas','/d',  'unassimilated fraction of sp. ingestion of N', &
                                          output=output_time_step_averaged)                                        
       call self%register_diagnostic_variable(self%id_Ngain_H,'NgainH','/d', 'sp. N gain rate by herbivory',   &
                                          output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_asefN,'asefN','-',  'assimilation efficiency of N', &
                                          output=output_time_step_averaged)
     end if
   end if
   
   !autotrophy
   if (self%fracaut .gt. .001) then
     call self%register_diagnostic_variable(self%id_MuClim_A,'MuClimA','/d', 'contribution of autotrophy to C limited sp. growth',&
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_NPPR,'NPPR','/d','-NPPR',              &
                                          output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_Cgain_A,'CgainA','/d', 'sp. C gain rate by autotrophy',   &
                                          output=output_time_step_averaged)                                          
     if (self%resolve_P) then
       call self%register_diagnostic_variable(self%id_Pgain_A,'PgainA','/d', 'sp. P gain rate by autotrophy',   &
                                          output=output_time_step_averaged)
     end if
     if (self%resolve_N) then
       call self%register_diagnostic_variable(self%id_Ngain_A,'NgainA','/d', 'sp. N gain rate by autotrophy',   &
                                          output=output_time_step_averaged)                                          
     end if
   end if

   ! Register environmental dependencies
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   !P
   if (self%resolve_P) then 
     !register the diagnostic var QPr also as a dependency, such that its value can be accessed
     call self%register_dependency(self%id_QPr_dep,'QPr', '-', 'relative QP')
   end if
   !N
   if (self%resolve_N) then 
     !register the diagnostic var QPr also as a dependency, such that its value can be accessed
     call self%register_dependency(self%id_QNr_dep,'QNr', '-', 'relative QN')
   end if
   
   return

   end subroutine initialize
!EOC
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
   class (type_gpm_plankton), intent(in)     :: self 
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: par,I_0,temp,depth,doy
   real(rk)                   :: iopt,rmd,rmdfz,rmdq,rmn,fT,fTA,ftH
   real(rk)                   :: plaC,asefC,MuClim,nutlim
   real(rk)                   :: plaP,asefP,Plim,QPreal,QP,fQP,QPrel,DIP
   real(rk)                   :: plaN,asefN,Nlim,QNreal,QN,fQN,QNrel,DIN
   real(rk)                   :: Silim,DISi
   real(rk)                   :: CGain,CGainA,CGainH,IngC,IngCunas,MuClimA
   real(rk)                   :: PGain,PGainA,PGainH,IngP,IngPunas,uptP
   real(rk)                   :: NGain,NGainA,NGainH,IngN,IngNunas,uptN
   real(rk)                   :: IngSi,IngSiunas
   integer                    :: i
   character(len=2)           :: istr !i converted to character string
   real(rk), parameter        :: secs_pr_day = 86400.
   logical                    :: oneofPMZ,debw
   type(prey_data)            :: PR
   allocate(PR%corpref(self%num_prey))
   allocate(PR%rpref(self%num_prey))
   allocate(PR%C(self%num_prey))
   allocate(PR%N(self%num_prey))
   allocate(PR%P(self%num_prey))
   allocate(PR%grC(self%num_prey))
   allocate(PR%grN(self%num_prey))
   allocate(PR%grP(self%num_prey))
   allocate(PR%grSi(self%num_prey))
   allocate(PR%Qr(self%num_prey))
   allocate(PR%QPr(self%num_prey))
   allocate(PR%QNr(self%num_prey))
   
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
   IngC=0.0;IngCunas=0.0;CGainH=0.0;CGainA=0.0;CGain=0.0;MuClimA=0.0;
   IngP=0.0;IngPunas=0.0;PGainH=0.0;PGainA=0.0;PGain=0.0;uptP=0.0
   IngN=0.0;IngNunas=0.0;NGainH=0.0;NGainA=0.0;NGain=0.0;uptN=0.0
   IngSi=0.0;IngSiunas=0.0;
   
   !_GET_GLOBAL_(self%id_doy,doy) !day of year
   _GET_(self%id_depth,depth)  ! depth
   _GET_(self%id_temp,temp)  ! temperature
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   !LOGICALS FOR DEBUGGING 
   if (debug .and. (depth .lt. 5.6 .and. depth .gt. 4.4)) then
     debw=.true.
   else 
     debw=.false.
   end if
   
   !if (debw) write(*,'(1A)',advance='no') ''
   
   ! RETRIEVE STATE VARIABLES and DEPENDENCIES
   !C
   _GET_(self%id_plaC,plaC) ! plankton carbon
   !P
   if (self%resolve_P) then 
     _GET_(self%id_plaP,plaP) ! plankton phosphours
   end if
   !P
   if (self%resolve_N) then 
     _GET_(self%id_plaN,plaN) ! plankton nitrogen
   end if
   
   IF (self%fracaut .lt. 1.0) then !i.e., if there will be any heterotrophy
    !collect prey abundances
    !reset PR%X, PR%QX
    PR%C=0.0
    PR%P=0.0
    PR%N=0.0
    PR%QPr=0.0
    PR%QNr=0.0
    PR%Qr=0.0 !this is min(QPr,QNr)
    DO i=1,self%num_prey
      !C
      _GET_STATE_(self%prey(i)%id_C,PR%C(i))
      !P
      if (self%resolve_P) then 
        _GET_STATE_(self%prey(i)%id_P,PR%P(i))
        !if pref(i)<0, then pref=f(Qr(i))
        if (self%prey(i)%pref .lt. 0.0) then
          _GET_(self%prey(i)%id_QPr,PR%QPr(i))
          if (.not. self%resolve_N) then
            PR%Qr(i)=PR%QPr(i)
          end if
        else 
          PR%Qr(i)=1.0
        end if
      end if
      !N
      if (self%resolve_N) then 
        _GET_STATE_(self%prey(i)%id_N,PR%N(i))
        !if pref(i)<0, then pref=f(Qr(i))
        if (self%prey(i)%pref .lt. 0.0) then
          _GET_(self%prey(i)%id_QNr,PR%QNr(i))
          if (.not. self%resolve_P) then
            PR%Qr(i)=PR%QNr(i)
          else
            PR%Qr(i)=min(PR%QPr(i),PR%QNr(i))
          end if
        else 
          PR%Qr(i)=1.0
        end if
      end if
    END DO 
   END IF
   
   !CALCULATE RHS terms
   
   !temperature functions
   if (self%fracaut .eq. 1.0) then
    fT=fT_fun(self,temp,'a') !tempmeth-2 = f(Q10), Q10 for autotrophs
   else if (self%fracaut .eq. 0.0) then
    fT=fT_fun(self,temp,'h') !tempmethod-3=f(Q10H): Q10H for heterotrophs
   else
    fTA=fT_fun(self,temp,'a') !tempmeth-2 = f(Q10), Q10 for autotrophs
    fTH=fT_fun(self,temp,'h') !tempmethod-3=f(Q10H): Q10H for heterotrophs
    fT=(self%fracaut*fTA+(1-self%fracaut)*fTH)/2 ! weighted average of fTA and fTH: to be used for combined processes 
   end if
   !if (debw) write(*,*) 'fT:',fT
   
   !internal states
   !P
   if (self%resolve_P) then
     QP=plaP/plaC  
     QPreal=QP
     !min qp correction
!      if (QP .lt. self%QPmin) then
!        if (debw) write(*,'(1A, 2F13.9)'),' '//trim(self%name)//'QP pushed to QPmin :', QP, self%QPmin 
!        QP=self%QPmin 
!       end if
!       !max qp correction
!       if (QP .gt. self%QPmax) then
!         if (debw) write(*,'(1A, 2F13.9)'),' '//trim(self%name)//'QP pushed to QPmax :', QP, self%QPmax 
!         QP=self%QPmax 
!       end if
     !down-regulation of nutrient uptake rate with increasing quota is described as a linear function
     if (self%QPmax .ne. self%QPmin) then
       fQP=min(1.0,max(0.0,(self%QPmax - QP)/(self%QPmax - self%QPmin)))
       !relative quota: a proxy for 'health'
       QPrel=1-fQP !(QP-self%QPmin)/(self%QPmax - self%QPmin) 
       !write(*,'(A,4F13.10)'),' '//trim(self%name)//'QP,fQP,QPrel',QP,fQP,QPrel
     else
       fQP=1.0
       QPrel=1.0
     end if
   end if
   !N
   if (self%resolve_N) then
     QN=plaN/plaC
     QNreal=QN
     !down-regulation of nutrient uptake rate with increasing quota is described as a linear function
     !min mass correction
!      if (QN .lt. self%QNmin) then
!        if (debw) write(*,'(1A, 2F13.9)'),' '//trim(self%name)//'QN pushed to QNnim :', QN, self%QNmin 
!        QN=self%QNmin 
!      end if
!      !max qn correction
!      if (QN .gt. self%QNmax) then
!        if (debw) write(*,'(1A, 2F13.9)'),' '//trim(self%name)//'QN pushed to QNmax :', QN, self%QNmax 
!        QN=self%QNmax 
!      end if
     if (self%QNmax .ne. self%QNmin) then
       fQN=min(1.0,max(0.0,(self%QNmax - QN)/(self%QNmax - self%QNmin)))
       !relative quota: a proxy for 'health'
       QNrel=1-fQN
       !write(*,'(A,4F13.10)'),' '//trim(self%name)//'QN,fQN,QNrel',QN,fQN,QNrel
     else
       fQN=1.0
       QNrel=1.0
     end if
   end if
   
   !autotrophy: nutrient uptake & light limitation   
   if (self%fracaut .gt. 0.001) then
     !C
     !Light-limited growth rate
     MuClimA=self%fracaut*MufI_fun(self,par,I_0,fT)
     !P
     if (self%resolve_P) then
       _GET_(self%id_DIP,DIP) ! nutrients
       uptP = uptake_fun(self,fQP,DIP,self%Kp,self%VPmax,fT)
       PGainA=self%fracaut*uptP !molP/molC/d 
       !write(*,*)'(P) v, fQP, DIP',uptP,fQP,DIP
     end if
     !N
     if (self%resolve_N) then
       _GET_(self%id_DIN,DIN) ! nutrients
       uptN = uptake_fun(self,fQN,DIN,self%Kn,self%VNmax,fT)
       NGainA=self%fracaut*uptN !molN/molC/d 
       !write(*,*)'(N) v, fQN, DIN',uptN,fQN,DIN
     end if
     !Si
     if (self%resolve_Si) then
       _GET_(self%id_DISi,DISi) ! nutrients
       !Si is assumed to be not stored, so no need for calculating uptSi and SiGainA
     end if
   else
     MuClimA=0.0
     PGainA=0.0
     NGainA=0.0
   end if 
   
   !heterotrophy
   if (self%fracaut .lt. 0.999) then !heterotrophy
     !calculate grazing rate for each prey unit and total grazing 
     !call Graze_Multi_Prey(self,fT,foodC,foodP,foodN,Qrel_prey,grazrateC,grazrateP,grazrateN,realpref)
     call Graze_Multi_Prey(self,fT,PR)
     !this gives grazrateX [molXprey/molXpred/d]
     
     !calculate the ingestion rates (PR%X*(1-fA))
     !C
     IngC=(1.0-self%fracaut)*sum(PR%grC) !molCprey/molCpred/d
     !P
     if (self%resolve_P) then
       IngP=(1.0-self%fracaut)*sum(PR%grP)  !molPprey/molCpred/d 
     end if
     !N
     if (self%resolve_N) then
       IngN=(1.0-self%fracaut)*sum(PR%grN)  !molPprey/molCpred/d 
     end if
     !Si
     if (self%resolve_Si) then
       IngSi=(1.0-self%fracaut)*sum(PR%grSi)  !molSiprey/molCpred/d 
     end if
     
     !adjust the assimilation efficiencies to maintain  Qmin<Q<Qmax
     !C
     asefC=self%asefC
     !P
     if (self%resolve_P) then
       if (self%QPmax .ne. self%QPmin) then !non-homeostatis:
         asefP=self%asefP*fQP               !asefP: decreases with increasing quota
       else                                 !homeostasis:
         asefP=self%asefP
         if (asefP*IngP .gt. asefC*IngC*self%QPmax) then !P-surplus
           asefP=asefC*IngC*self%QPmax/IngP          !down-regulate asefP according to asefC
           !if (debw) write(*,'(A,2F13.10)'),' '//trim(self%name)//' (adj eP) P surplus.  aeC, aeP:',asefC,asefP
         else if (asefP*IngP .lt. asefC*IngC*self%QPmax) then !C-surplus
           asefC=asefP*IngP/(IngC*self%QPmax)        ! down-regulate asefC according to asefP
           !if (debw) write(*,'(A,2F13.10)'),' '//trim(self%name)//' (adj eP) C surplus.  aeC, aeP:',asefC,asefP
         endif !if C-surplus
       end if !if homeostatis
     end if !if resolve_P
     !N
     if (self%resolve_N) then
       if (self%QNmax .ne. self%QNmin) then      !non-homeostatis: 
         asefN=self%asefN*fQN                    !asefP: decreases with increasing quota
       else                                      !homeostasis:
         asefN=self%asefN
         if (asefN*IngN .gt. asefC*IngC*self%QNmax) then !N-surplus
           asefN=asefC*IngC*self%QNmax/IngN              !down-regulate asefN according to asefC
           !if (debw) write(*,'(A,2F13.10)'),' '//trim(self%name)//' (adj eN) N surplus. aeC, aeN:',asefC,asefN
         else if (asefN*IngN .lt. asefC*IngC*self%QNmax) then !C-surplus
           asefC=asefN*IngN/(IngC*self%QNmax)            !down-regulate asefC according to asefN
           !if (debw) write(*,'(A,2F13.10)'),' '//trim(self%name)//' (adj eN) C surplus. aeC, aeN:',asefC,asefN
           !as asefC was down-regulated, asefP may need to be down-regulated as well
           if (self%resolve_P .and. self%QPmax .eq. self%QPmin) then
             if (asefP*IngP .gt. asefC*IngC*self%QPmax) then !P-surplus
               asefP=asefC*IngC*self%QPmax/IngP          !adjust asefP according to asefC
               !if (debw) write(*,'(A,2F13.10)'),' '//trim(self%name)//' (adj eP2) P surplus. aeC, aeP:',asefC,asefP
             endif !if P-surplus
           end if !if resolve_P
         endif ! C-surplus
       end if !if homeostatis
     end if !resolve_N
     
     !C
     CGainH=IngC*asefC ![molCprey/molCpred /d]
     IngCunas=IngC*(1.0-asefC)
     !P
     if (self%resolve_P) then
       PGainH=IngP*asefP
       IngPunas=IngP*(1.0-asefP)
     end if
     !write(*,'(A,2F7.4)'),' '//trim(self%name)//' QP,PgainH/CgainH:',QP,PGainH/CGainH
     !N
     if (self%resolve_N) then
       NGainH=IngN*asefN
       IngNunas=IngN*(1.0-asefN)
     end if
     !Si
     if (self%resolve_Si) then
       IngSiunas=IngSi*(1.0-0.0) !i.e., asefSi=0.0
     end if
     !if (debw) write(*,'(A,2F7.4)'),' '//trim(self%name)//' QN,NgainH/CgainH:',QN,NGainH/CGainH
   else !noheterotrophy
     !C
     CGainH=0.0
     !P
     PGainH=0.0
     IngPunas=0.0
     !N
     NGainH=0.0
     IngNunas=0.0
     !Si
     IngSiunas=0.0
   end if
   
   !combined (autotrophy + heterotrophy): XGain
   !C-limited potential growth
   MuClim=MuClimA+CGainH
   !P
   nutlim=1.0
   if (self%resolve_P) then
     Plim=fQ_fun(QP,self%QPmax,self%QPmin) !Nutrient-limited growth rate
     nutlim=min(nutlim,Plim)
     PGain=PGainA+PGainH
   end if
   !N
   if (self%resolve_N) then
     Nlim=fQ_fun(QN,self%QNmax,self%QNmin) !Nutrient-limited growth rate
     nutlim=min(nutlim,Nlim)
     NGain=NGainA+NGainH
   end if
   !Si
   if (self%resolve_Si) then
     !Si is assumed to be not stored, therefore Si limitation is directly f(DISi) 
     Silim=DISi/(DISi+self%Ksi)
     nutlim=min(nutlim,Silim)
   end if
   
   if (self%fracaut .gt. 0.001) then !if any autotrophy involved
     CGain=min(self%rmax*fT*nutlim,MuClim)
     ! Autotrophic Cgain is total-CGainH (to potentially calculate DIC and DISi uptake)
     CGainA=max(0.0,CGain-CGainH)
   else
     CGain=CGainH
     CGainA=0.0
   end if
   
   !if (debw) write(*,'(A,6F16.12)'),' '//trim(self%name)//' QN,Ng/Cg,NgH,CgH,NgA,CgA:',QN,NGain/CGain,NGainH,CGainH,NGainA,CGainA
   
   ! mortality
   rmn=self%rmn*fT
   rmdq=self%rmdq*fT*plaC ! d-1 (mmolC/m3)-1 * mmolC/m3 = d-1
   !depth-dependent mortality rate
   rmdfz=self%rmd +(self%rmdH-self%rmd) * (1.0/(1.0+exp(0.1*(-depth+self%rmdfZ_offset))))
   rmd=rmdfz*fT
   !if (self%rmdH .NE. self%rmd) write(*, '(2A, 4F12.8)') self%name,'; depth,rmd,rmdH,rmdfz:', depth, self%rmd, self%rmdH, rmdfz
   
   ! SET RHS 
   !C
   _SET_ODE_(self%id_plaC, CGain*(plaC+self%m0) - plaC*(rmn + rmd + rmdq))
   !P
   if (self%resolve_P) then
     _SET_ODE_(self%id_plaP,  PGain*plaC + CGain*self%m0*QP  - plaP*(rmn + rmd + rmdq)) 
   end if
   !N
   if (self%resolve_N) then
     _SET_ODE_(self%id_plaN,  NGain*plaC + CGain*self%m0*QN  - plaN*(rmn + rmd + rmdq))
   end if

   !uptake target
   if (self%fracaut .gt. 0.001) then
     !C
     !_SET_ODE_(self%id_DIC,-CGainA*(plaC+self%m0) )
     !P
     if (self%resolve_P) then
       _SET_ODE_(self%id_DIP,-PGainA*(plaC+self%m0) ) !uptP*self%fracaut
     end if
     !N
     if (self%resolve_N) then
       _SET_ODE_(self%id_DIN,-NGainA*(plaC+self%m0) ) !uptN*self%fracaut
     end if
     !Si
     if (self%resolve_Si) then
       _SET_ODE_(self%id_DISi,-CGainA*self%Si2C*(plaC+self%m0) )
     end if
   end if
   
   !grazing targets
   if (self%fracaut .lt. 0.999) then
     if (debw) write(*,*),' '//trim(self%name)//' : PR%grC',PR%grC
     DO i=1,self%num_prey
       !C
       _SET_ODE_(self%prey(i)%id_C,-PR%grC(i)*plaC*(1.0-self%fracaut)) !molC/molC/d *molC/m3 =molC/m3/d
       !P
       if (self%resolve_P) then
         _SET_ODE_(self%prey(i)%id_P,-PR%grP(i)*plaC*(1.0-self%fracaut)) !molP/molC/d *molC/m3 =molP/m3/d
       end if
       !N
       if (self%resolve_N) then
         _SET_ODE_(self%prey(i)%id_N,-PR%grN(i)*plaC*(1.0-self%fracaut)) !molN/molC/d *molC/m3 =molP/m3/d
       end if
     END DO
   END IF
   
   !recycling to nutrient pools
   !C
   _SET_ODE_(self%id_detC,(rmd+rmdq)*plaC)
   _SET_ODE_(self%id_detC,IngCunas*self%unas_detfrac*plaC)
   !P
   if (self%resolve_P) then
     _SET_ODE_(self%id_DIP,rmn*plaP)
     _SET_ODE_(self%id_detP,(rmd+rmdq)*plaP)
     _SET_ODE_(self%id_detP,IngPunas*self%unas_detfrac*plaC) 
     _SET_ODE_(self%id_DIP,IngPunas*(1.-self%unas_detfrac)*plaC) 
     !if (debw) write(*,'(A, 3F14.10)') 'rmd+rmdq,IngCunas,IngPunas', rmd+rmdq,IngCunas,IngPunas
   end if
   !N
   if (self%resolve_N) then
     _SET_ODE_(self%id_DIN,rmn*plaN)
     _SET_ODE_(self%id_detN,(rmd+rmdq)*plaN)
     _SET_ODE_(self%id_detN,IngNunas*self%unas_detfrac*plaC) 
     _SET_ODE_(self%id_DIN,IngNunas*(1.-self%unas_detfrac)*plaC) 
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'IngPunas', IngPunas*secs_pr_day
   end if
   !Si
   if (self%resolve_Si) then
     if (self%fracaut .eq. 1.0) then
       _SET_ODE_(self%id_DISi,rmn*plaC*self%Si2C)
       _SET_ODE_(self%id_detSi,(rmd+rmdq)*plaC*self%Si2C)
     end if
     _SET_ODE_(self%id_detSi,IngSiunas*self%unas_detfrac*plaC)
     _SET_ODE_(self%id_DISi,IngSiunas*(1.-self%unas_detfrac)*plaC)
     !if (debw) write(*,'(2A, 2F14.10)') self%name, 'IngSiunas', IngSiunas
   end if
   
   ! SET DIAGNOSTICS
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_rmd,(rmd)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_Cgain,CGain*plaC*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_Closs,(plaC*(rmn+rmd+rmdq))*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_MuClim,MuClim*secs_pr_day)
   if (self%resolve_P) then
     _SET_DIAGNOSTIC_(self%id_QP, QPreal)
     _SET_DIAGNOSTIC_(self%id_QPr , QPrel)
     _SET_DIAGNOSTIC_(self%id_Pgain,PGain*plaC*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Ploss,(plaP*(rmn+rmd+rmdq))*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Plim,Plim)
   end if
   if (self%resolve_N) then
     _SET_DIAGNOSTIC_(self%id_QN, QNreal)
     _SET_DIAGNOSTIC_(self%id_QNr , QNrel)
     _SET_DIAGNOSTIC_(self%id_Ngain,NGain*plaC*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Nloss,(plaN*(rmn+rmd+rmdq))*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Nlim,Nlim)
   end if
   if (self%resolve_Si) then
     _SET_DIAGNOSTIC_(self%id_Silim,Silim)
   end if
   if (self%resolve_N .and. self%resolve_P) then
     _SET_DIAGNOSTIC_(self%id_N2P, plaN/plaP)
   end if
   
   !heterotrophy
   if (self%fracaut .lt. 0.999) then
     DO i=1,self%num_prey
      _SET_DIAGNOSTIC_(self%prey(i)%id_realpref,PR%rpref(i))	
     END DO
     _SET_DIAGNOSTIC_(self%id_Cgain_H, CGainH*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_IngC,IngC*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_IngCunas, IngCunas*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_asefC,asefC)
     if (self%resolve_P) then
       _SET_DIAGNOSTIC_(self%id_Pgain_H, PGainH*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_IngP,IngP*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_IngPunas, IngPunas*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_asefP,asefP)
     end if
     if (self%resolve_N) then
       _SET_DIAGNOSTIC_(self%id_Ngain_H, NGainH*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_IngN,IngN*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_IngNunas, IngNunas*secs_pr_day)
       _SET_DIAGNOSTIC_(self%id_asefN,asefN)
     end if
     !Si
     if (self%resolve_Si) then
       _SET_DIAGNOSTIC_(self%id_IngSiunas, IngSiunas*secs_pr_day)
     end if
   end if
   if (self%fracaut .gt. 0.001) then 
     _SET_DIAGNOSTIC_(self%id_MuClim_A,MuClimA*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_NPPR,(CGainA-rmn-rmd)*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Cgain_A, CGainA*secs_pr_day)
     if (self%resolve_P) then
       _SET_DIAGNOSTIC_(self%id_Pgain_A, PGainA*secs_pr_day)
     end if
     if (self%resolve_N) then
       _SET_DIAGNOSTIC_(self%id_Ngain_A, NGainA*secs_pr_day)
     end if
   end if
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
   
   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Grazing on Multiple Prey Items
!
! !INTERFACE:
   subroutine Graze_Multi_Prey(self,fT,PR)
   !
! !DESCRIPTION:
! Here, grazing on multiple prey is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_plankton), intent(in)     :: self
   real(rk),intent(in)          :: fT
   type(prey_data)              :: PR
   real(rk)                     :: prefmin,KzFact,ftotC,ftotP,ftotN
   integer                      :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   ! Retrieve the values of all grazed food and calculate the total food (perceived, based on pref values)
   PR%corpref=0.0
   ftotC=0.0 !;ftotP=0.0;ftotN=0.0
   DO i=1,self%num_prey
     !correct the preference if required (indicated by a negative value):
     if (self%prey(i)%pref .lt. 0) then
       !At nutrient replete conditions, some phytoplankton species build colonies that reduce their edibility. 
       !To emulate this, (when pref <0) we assume that the preference (edibility) sigmoidially increases from 
       !prefmin (=pref/5 @Q=Qmax) to the originally specified (absolute) value (@Q=Qmin)
       prefmin=-1*self%prey(i)%pref/5
       PR%corpref(i)=prefmin+(-1*self%prey(i)%pref-prefmin) * (1.0-1.0/(1.0+exp(- PR%Qr(i)*10.0+5.0)))
       !write(*,'(A, I1, A, 3F13.9)') '-prey#',i,' fQprey, oldpref, newpref: ',PR%Qr(i),self%prey(i)%pref,PR%rpref(i)
     else 
       PR%corpref(i)=self%prey(i)%pref
     end if
     
     !calculate the total pref-weighted food (needed only if prefs are dynamically adjusted)
     if (self%dynpref) then 
        ftotC=ftotC+PR%corpref(i)*PR%C(i)
        !ftotP=ftotP+PR%corpref(i)*PR%P(i)
        !ftotN=ftotN+PR%corpref(i)*PR%N(i)
     end if
     
   END DO
   
   !calculate grazing rate for each prey unit and total grazing  
   PR%grC=0.0; PR%grN=0.0;PR%grP=0.0
   PR%rpref=0.0
   DO i=1,self%num_prey
      KzFact=1.0
      !if the preferences are specified to change dynamically (e.g., as in Fasham 1990)
      if (self%dynpref) then 
        PR%rpref(i)=PR%corpref(i)*PR%C(i)/ftotC
        !write(*,'(A,I1,A,3F8.4)'),' '//trim(self%name)//' prey#',i,' corpref,C/totC,rpref:',PR%corpref(i),PR%C(i)/ftotC,PR%rpref(i)
      else
        PR%rpref(i)=PR%corpref(i)
      end if    
      PR%grC(i)=fpzMonMulti(self,PR%C(i),PR%rpref(i),ftotC,fT,KzFact) !molCprey/molCpred/d
      PR%grP(i)       = PR%grC(i)         * PR%P(i)/PR%C(i) 
      PR%grN(i)       = PR%grC(i)         * PR%N(i)/PR%C(i)
      PR%grSi(i)      = PR%grC(i)         * self%prey(i)%Si2C
      !molXprey/molCpred/d = molCprey/molCpred/d  * molXprey/molCprey
      
      !write(*,*),' '//trim(self%name)//' prey#',i,'grC,grP,grN',PR%grC(i),PR%grP(i),PR%grN(i)
   END DO
   
   !write(*,'(A, 2F13.10)'),' '//trim(self%name)//' grSi:',PR%grSi  
   end subroutine Graze_Multi_Prey
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Monod formulation for zooplankton grazing on multiple prey (as in Fasham 1990).
!
! !INTERFACE:
   pure real(rk) function fpzMonMulti(self,food,pref,foodtot,fT,KzFact)
!                                    
! !DESCRIPTION:
! !The function calculates only the grazing for a single one of those multiple prey items 
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_plankton), intent(in)     :: self
   real(rk), intent(in)          :: food,pref,foodtot,fT,KzFact
   real(rk)                      :: gmax
!
!EOP
!-----------------------------------------------------------------------
!BOC

   gmax=self%gmax*fT
   fpzMonMulti = gmax*food*pref / (self%Kz*KzFact+foodtot)
   
   end function fpzMonMulti
!EOC
!-----------------------------------------------------------------------
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function MufI_fun(self,par,I_0,fT)
!
! !DESCRIPTION:
! Here, light limited growth is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_plankton), intent(in)     :: self
   real(rk), intent(in)         :: par,I_0,fT
   real(rk)                     :: mumax,fI
   real(rk), parameter :: secs_pr_day = 86400.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !some formulations of light limitation below depend on mumax.
   !mu_max is unknown, but only mu_inf. calculate rmax=rinf*f(Q=Qmax) 
   mumax=self%rmax*secs_pr_day*(_ONE_ - self%QPmin/self%QPmax)
   
   select case (self%lightresp)
     case(1) !monod type, e.g., schwaderer et al
       fI=par/(par+mumax/self%islope)
     case(2) !with light inhibition type, e.g., schwaderer et al
       fI=par/( par**2 *mumax/self%islope/self%Iopt**2+ par*(1-2*mumax/self%islope/self%Iopt) +mumax/self%islope ) 
     case(3) !tanh, eg. merico & oguz
       fI=TANH(par*self%islope)
     case(4) ! Light acclimation formulation based on surface light intensity.
       fI=par/max(0.25*I_0,self%Imin) *exp(_ONE_-par/max(0.25*I_0,self%Imin))
   end select
   
   MufI_fun=self%rmax*fT*fI
  
   end function MufI_fun
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten or Droop formulations for nutrient limitation
!
! !INTERFACE:
   pure real(rk) function fQ_fun(Q,Qmax,Qmin)
!
! !DESCRIPTION:
! Here, nutrient limited plankton growth is formulated
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   real(rk), intent(in)         :: Q,Qmax,Qmin
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !fN=0 @ Q=Qmin and fN->1 as QP->inf
   if (Qmax .ne. Qmin) then
     fQ_fun= _ONE_ - Qmin/Q
   else !if Qmax=Qmin, no nutrient limitation
     fQ_fun= _ONE_
   end if
  
   end function fQ_fun
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:                            
   pure real(rk) function uptake_fun(self,fQ,nut,K,Vmax,fT)
!
! !DESCRIPTION:
! Here, nutrient uptake is formulated
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_plankton), intent(in)     :: self
   real(rk), intent(in)         :: fQ,nut,Vmax,K,fT
   real(rk)                     :: fN
!
!EOP
!-----------------------------------------------------------------------
!BOC
 
   !External nutrient limitation is described as a Monod function
   fN = nut/(K+nut)

   uptake_fun = fT*Vmax*fQ*fN

   end function uptake_fun
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function fT_fun(self,temp,trlev)
!
! !DESCRIPTION:
! Here, temperature response functions are formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_plankton), intent(in)     :: self
   real(rk), intent(in)         :: temp
   character*1, intent(in)      :: trlev
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !temp effect on rmax
   select case (self%tempresp)
     case (1) !no temp effect
       fT_fun=1
     case (2) !Q10 
       if (trlev == 'a') then !autotrophy
         fT_fun=self%Q10**((temp-self%Tref)/10)
       else if (trlev == 'h') then
         fT_fun=self%Q10H**((temp-self%Tref)/10)
       end if
     case (3) ! Lancelot et al 2002
       fT_fun=exp(-((temp-self%Topt)**2)/self%Tint**2)

    end select
 end function fT_fun
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
   class (type_gpm_plankton), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
!
! !LOCAL VARIABLES:
   real(rk)                   :: vert_vel,temp,par,I0,QPr !,yearday,par
   real(rk), parameter        :: secs_pr_day = 86400.
!
!EOP
!-----------------------------------------------------------------------
!BOC

! way to get the day of year
!_GET_DEPENDENCY_SCALAR_(self%id_yearday,yearday)  ! day of year (diff(d/m/Y - 1/1/Y))

! Enter spatial loops (if any)
_FABM_LOOP_BEGIN_

select case (self%velmet)
  case (0)
   vert_vel=self%w
  case (1)  !sinking velocity decreases with increasing quota
   _GET_ (self%id_QPr_dep,QPr)
   vert_vel=self%w * (1-1/(1+exp(10*(.5-QPr))))
  case (2) !velocity is corrected by the viscosity, calculated as a function of temperature
   _GET_  (self%id_temp,temp)  ! temperature
   !in the expression below, denominator (10**..) gives the viscosity ratio at the ambient and reference (T=20 oC) temperatures (Kestin et al. 1978). So the overall expression is: mu(20)/mu(T) * vel
   vert_vel=self%w*1/(10**(((20-temp)/(temp+96))*(1.2378 - 1.303e-3*(20-temp) + 3.06e-6*(20-temp)**2 + 2.55e-8*(20-temp)**3)))  
  case (3) !active swimming towards the prescribed optimal temperature
   _GET_  (self%id_temp,temp)  ! temperature
   vert_vel=self%w * TANH((self%Topt-temp)/self%Tsc) ! when the expression is negative (temp>Topt), velocity is negative (sinking) self%Topt
  case (4) !active swimming towards the prescribed optimal light intensity
   _GET_   (self%id_par,par)  ! light
   _GET_HORIZONTAL_ (self%id_I_0,I0)  ! light
   if (I0 .gt. 0) then !when it's day light, swim to the optimal light intensity
    vert_vel=self%w * TANH((self%Iopt-par)/self%Isc) ! when the expression is negative (par>Iopt), velocity is negative (sinking) 
   !else  !otherwise sink down with max. speed ?
   ! vert_vel=-1.0*self%w
   end if
end select

!Set these calculated vertical_movement values 
_SET_VERTICAL_MOVEMENT_(self%id_plaC,vert_vel) 
_SET_VERTICAL_MOVEMENT_(self%id_plaP,vert_vel) 

!_SET_DIAGNOSTIC_(self%id_vv , vert_vel*secs_pr_day)

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
   class (type_gpm_plankton), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: plaC
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_plaC,plaC) ! plankton carbon biomass

   ! Self-shading with explicit contribution from background plankton concentration.
   _SET_EXTINCTION_(self%kc*plaC)
   !_SET_EXTINCTION_(0.01)
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
! vertical integrals can be calculated here, eg: integral primary production
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
     class (type_gpm_plankton), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_SURFACE_
!
!EOP
!-----------------------------------------------------------------------
!BOC
     _HORIZONTAL_LOOP_BEGIN_

        !_GET_HORIZONTAL_(self%id_totN_vertint,tot_vi_N)
        !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_vertint_diag,_REPLNAN_(tot_vi_N))
        !_GET_HORIZONTAL_(self%id_totC_vertint,tot_vi_C)
        !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_vertint_diag,_REPLNAN_(tot_vi_C))

      _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC
!-----------------------------------------------------------------------
   
end module gpm_plankton
