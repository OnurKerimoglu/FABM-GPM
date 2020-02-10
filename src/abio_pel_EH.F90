#include "fabm_driver.h"
#include "fabm.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gpm_abio_pel_EH --- GPM pelagic abiotic model (based on ECOHAMv5): DI(C,NH4,NO3,P,Si),DO(C,N,P),bac(C),det1(C,N,P),det2(C,N,P,Si)
!
!TODO:
! - self%metCexc -> general switch (aux)
! - self%resolve_Si -> general switch (aux)
!
! !INTERFACE:
   module gpm_abio_pel_EH
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_expressions

!  default: all is private.
   private
   
   real(rk), parameter :: pi=acos(-1._rk)
   real(rk), parameter :: s2d = 86400.
   
   ! !PUBLIC DERIVED TYPES:
   type,public :: type_env
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
   
   ! gpm_abio_pel_EH 
   type, extends (type_base_model),public :: type_gpm_abio_pel_EH
!     Variable identifiers
      type (type_state_variable_id)     :: id_DIC,id_DIP,id_DINO3,id_DINH4,id_DISi
      type (type_state_variable_id)     :: id_DOC,id_DOP,id_DON,id_SOC
      type (type_state_variable_id)    :: id_det1C,id_det1P,id_det1N
      type (type_state_variable_id)    :: id_det2C,id_det2P,id_det2N,id_det2Si !, id_det2k
      type (type_state_variable_id)    :: id_bacC,id_bacP,id_bacN
      type (type_state_variable_id)    :: id_O2
      type (type_dependency_id)         :: id_temp,id_salt
      type (type_horizontal_dependency_id)         :: id_wind,id_zmax
      type (type_global_dependency_id)  :: id_doy
      type (type_diagnostic_variable_id) :: id_det1QP,id_det1QN
      type (type_diagnostic_variable_id) :: id_det2QP,id_det2QN,id_det2QSi
      !flux diagnostics
      type (type_diagnostic_variable_id) :: id_det1_dom_C,id_det1_dom_P,id_det1_dom_N
      type (type_diagnostic_variable_id) :: id_det2_dom_C,id_det2_dom_P,id_det2_dom_N
      type (type_diagnostic_variable_id) :: id_det2_dim_Si,id_soc_doc
      type (type_diagnostic_variable_id) :: id_dom_dim_C,id_dom_dim_P,id_dom_dim_N
      type (type_diagnostic_variable_id) :: id_O2_Sat,id_O2_percSat,id_o2o_bac,id_o2o_n4n,id_n4n_n3n,id_n3n_nn2
      !bacteria
      type (type_diagnostic_variable_id) :: id_bacQP,id_bacQN
      type (type_diagnostic_variable_id) :: id_dom_bac_C,id_dom_bac_P,id_dom_bac_N
      type (type_diagnostic_variable_id) :: id_dim_bac_C,id_dim_bac_P,id_dim_bac_N
      type (type_diagnostic_variable_id) :: id_bac_dim_C,id_bac_dim_P,id_bac_dim_N
      type (type_horizontal_diagnostic_variable_id)  :: id_air_o2o !,id_air_o2c
      type (type_horizontal_diagnostic_variable_id)  :: id_DIP_Sflux,id_DINO3_Sflux,id_DINH4_Sflux,id_DISi_Sflux
      !light
      type (type_dependency_id)         :: id_par, id_par_dmean,id_env_k0
      type (type_diagnostic_variable_id)  :: id_dpar, id_dpar_dmean,id_k0,id_kw
      type (type_horizontal_dependency_id)  :: id_I0,id_I0_dmean
      type (type_horizontal_diagnostic_variable_id)  :: id_dI0, id_dI0_dmean

!     Model parameters
      !common
      real(rk) :: Q10,Tref,cwind
      real(rk) :: a_water,a_minfr,a_fz
      logical  :: flux_diags,resolve_Si,resolve_DIC,resolve_bac,resolve_bacnp
      !Surface Fluxes
      integer  :: metCexc,metairsea,metairseaEH,Idm_met,kwFzmaxMeth
      !Det
      real(rk) :: det1_w,det2_w,det1_kc,det2_kc,r_det1dec,r_det2dec,r_detdec_fC
      real(rk) :: dom_kc,soc_rate,rhoC,rhoN,rhoP,xknit 
      real(rk) :: b_exc,b_xk4,b_xkp,b_v,b_rcn,b_rcp
      integer  :: det_velmet
      !DIM
      real(rk) :: DIP_Sflux,DINO3_Sflux,DINH4_Sflux,DISi_Sflux
   
      contains
      
      procedure  :: initialize
      procedure  :: do
      procedure  :: do_surface
      procedure  :: get_light_extinction
      procedure  :: get_vertical_movement
      
   end type
   
   !Local variables
   LOGICAL :: debug = .false.
   
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the mat-model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the parameters and variables are registered det_with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel_EH), intent(inout), target  :: self
   integer,                  intent(in)              :: configunit
!
! !LOCAL VARIABLES:
!real(rk), parameter :: s2d = 1._rk/86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our odet_wn derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   !general
   call self%get_parameter(self%Idm_met, 'Idm_met',   '-',          'method to calculate daily average PAR',          default=0)
   call self%get_parameter(self%flux_diags,  'flux_diags',  '-',         'whether to output flux diagnostics',   default=.false.)
   call self%get_parameter(self%resolve_Si,  'resolve_Si',  '-',         'whether to resolve Si cycle',          default=.false.)
   call self%get_parameter(self%resolve_DIC,  'resolve_DIC',  '-',         'whether to resolve DIC',          default=.false.)
   call self%get_parameter(self%metCexc,   'metCexc',  '-',              'method to calculate excess C-uptake, 0:no excess utpake; 1: as in EH',    default=0)
   call self%get_parameter(self%metairsea, 'metairsea',  '-',            'method to resolve air-sea fluxes', default=1)
   call self%get_parameter(self%metairseaEH, 'metairseaEH',  '-',        'EcoHam method to resolve air-sea fluxes', default=1)
   call self%get_parameter(self%cwind, 'cwind',  'm/s',        'Constant wind speed', default=-9.0_rk)
   call self%get_parameter(self%kwFzmaxMeth, 'kwFzmaxMeth',  '-',          'background extinction method.',         default=0)
   call self%get_parameter(self%a_water      ,'a_water',  '/m', 'max. background attenuation coefficient',      default=2.0_rk)
   call self%get_parameter(self%a_minfr      ,'a_minfr',  '-', 'heuristic depth-dep attenuation',     default=0.1_rk)
  call self%get_parameter(self%a_fz         ,'a_fz',  '-', 'depth dependent turbidity gradient',        default=10.0_rk)
   call self%get_parameter(self%dom_kc,     'dom_kc',      'm^2/mmolC',  'C-specific light extinction of dom',           default=0.01_rk)
   call self%get_parameter(self%det1_kc,     'det1_kc',      'm^2/mmolC',  'C-specific light extinction of slow sinking (small) detritus',           default=0.03_rk)
   call self%get_parameter(self%det2_kc,     'det2_kc',      'm^2/mmolC',  'C-specific light extinction of fast sinking (large) detritus',           default=0.03_rk)
   call self%get_parameter(self%det_velmet, 'det_velmet',  '-',          'velocity method of detritus',         default=1)
   call self%get_parameter(self%det1_w,      'det1_w',       'm/d',        'vertical velocity for slow sinking detritus (<0 for sinking)',  default=-0.4_rk, scale_factor=1./s2d)
   call self%get_parameter(self%det2_w,      'det2_w',       'm/d',        'vertical velocity for fast sinking detritus (<0 for sinking)',  default=-10.0_rk, scale_factor=1./s2d)
   call self%get_parameter(self%Q10,        'Q10',         '-',          'Q10 for bacterial processes',         default=2.0_rk)
   call self%get_parameter(self%Tref,       'Tref',        'celcius',    'reference temperature for bac. proc.',default=10.0_rk)
   call self%get_parameter(self%r_det1dec,   'r_det1dec',    'd-1',       'detritus-1 (slow sinking) decay rate',  default=0.12_rk,  scale_factor=1./s2d) !in EH, this is xmu4n
   call self%get_parameter(self%r_det2dec,   'r_det2dec',    'd-1',       'detritus-2 (fast sinking) decay rate',  default=0.1_rk,  scale_factor=1./s2d) !in EH, this is xmu5n
   call self%get_parameter(self%r_detdec_fC,'r_detdec_fC', '0',          'detritus decay rate multiplier for C',  default=0.85_rk) !in EH, this is rxmu4c
   call self%get_parameter(self%soc_rate,   'soc_rate',    'd-1',        'decay rate of SOC to DOC',  default=0.00274_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%rhoC,   'rhoC',    'd-1',        'doc rem. rate',  default=0.1_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%rhoP,   'rhoP',    'd-1',        'dop rem. rate',  default=0.1_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%rhoN,   'rhoN',    'd-1',        'don rem. rate',  default=0.1_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%xknit,      'xknit',   'd-1',       'nitrificaiton rate', default=0.02_rk,scale_factor=1./s2d)
   ! bacterial rates
   call self%get_parameter(self%resolve_bac, 'resolve_bac', '-','whether to resolve bacterial remineralization', default=.false.)
   call self%get_parameter(self%resolve_bacnp, 'resolve_bacnp', '-','resolve bacterial N&P explicitly', default=.false.)
   call self%get_parameter(self%b_exc,   'b_exc',    'd-1',     'excretion rate of bac',  default=0.2_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%b_v,     'b_v',    'd-1',       'max upt rate of bac',  default=0.5_rk,  scale_factor=1./s2d)
   call self%get_parameter(self%b_xk4,   'b_xk4',  'mmol N/m^3', 'half sat conc for DIN upt of bac',  default=0.2_rk)
   call self%get_parameter(self%b_xkp,   'b_xkp',  'mmol P/m^3', 'half sat conc for DOP upt of bac',  default=0.02_rk)
   call self%get_parameter(self%b_rcn,   'b_rcn',    'molC/molN','molar C:N  ratio of bac',  default=5.0_rk)
   call self%get_parameter(self%b_rcp,   'b_rcp',    'molC/molP','molar C:P  ratio of bac',  default=50.0_rk)
   call self%get_parameter(self%DIP_Sflux,  'DIP_Sflux',   'mmol/m^2/d', 'DIP flux at surface',                 default=0.0_rk, scale_factor=1./s2d)
   call self%get_parameter(self%DINO3_Sflux,  'DINO3_Sflux',   'mmol/m^2/d', 'NO3 flux at surface',                 default=0.0_rk, scale_factor=1./s2d)
   call self%get_parameter(self%DINH4_Sflux,  'DINH4_Sflux',   'mmol/m^2/d', 'NH4 flux at surface',                 default=0.0_rk, scale_factor=1./s2d)
   call self%get_parameter(self%DISi_Sflux, 'DISi_Sflux',   'mmol/m^2/d','DISi flux at surface',                 default=0.0_rk, scale_factor=1./s2d)
   
   ! Register state variables
   !C
   if (self%resolve_DIC) then
    call self%register_state_variable(self%id_DIC,'DIC','mmol C/m^3','dissolved inorganic C', & 
                                    minimum=0.0_rk,no_river_dilution=.true.)
    call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_DIC)
   end if
   call self%register_state_variable(self%id_DOC,'DOC','mmol C/m^3','dissolved organic C', & 
                                    minimum=0.0_rk, specific_light_extinction=self%dom_kc,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_DOC)
   call self%register_state_variable(self%id_det1C,'det1C','mmol C/m^3','C in slow sinking detritus', & 
                                    minimum=0.0_rk, specific_light_extinction=self%det1_kc,vertical_movement=self%det1_w*s2d,&
                                    no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_det1C)
   call self%register_state_variable(self%id_det2C,'det2C','mmol C/m^3','C in fast sinking detritus', & 
                                    minimum=0.0_rk, specific_light_extinction=self%det2_kc,vertical_movement=self%det2_w*s2d,&
                                    no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_det2C)

   if (self%metCexc .ne. 0) then
     call self%register_state_variable(self%id_SOC,'SOC','mmol C/m^3','semi-dissolved organic C', & 
                                    minimum=0.0_rk, specific_light_extinction=self%dom_kc,no_river_dilution=.true.)
     call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_SOC)
   end if
   
   !P
   call self%register_state_variable(self%id_DIP,'DIP','mmolP/m^3','dissolved inorganic P', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DIP)
   call self%register_state_variable(self%id_DOP,'DOP','mmolP/m^3','dissolved organic P', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DOP)
   call self%register_state_variable(self%id_det1P,'det1P','mmolP/m^3','P in slow sinking detritus', & 
                                    minimum=0.0_rk,vertical_movement=self%det1_w*s2d,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_det1P)
   call self%register_state_variable(self%id_det2P,'det2P','mmolP/m^3','P in fast sinking detritus', & 
                                    minimum=0.0_rk,vertical_movement=self%det2_w*s2d,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_det2P)

   !N
   call self%register_state_variable(self%id_DINO3,'DINO3','mmolN/m^3','NO3', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DINO3)
   call self%register_state_variable(self%id_DINH4,'DINH4','mmolN/m^3','NH4', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DINH4)
   call self%register_state_variable(self%id_DON,'DON','mmol N/m^3','dissolved organic N', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DON)
   call self%register_state_variable(self%id_det1N,'det1N','mmolN/m^3','N in slow sinking detritus', & 
                                    minimum=0.0_rk,vertical_movement=self%det1_w*s2d,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_det1N)
   call self%register_state_variable(self%id_det2N,'det2N','mmolN/m^3','N in fast sinking detritus', & 
                                    minimum=0.0_rk,vertical_movement=self%det2_w*s2d,no_river_dilution=.false.)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_det2N)
   
   !O2
   call self%register_state_variable(self%id_O2,'O2','mmol O2/m^3','O2',no_river_dilution=.true.)!, & minimum=0.0_rk)
   
   !Si
   if (self%resolve_Si) then     
     call self%register_state_variable(self%id_DISi,'DISi','mmol Si/m^3','dissolved inorganic silica', & 
                                    minimum=0.0_rk,no_river_dilution=.false.)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_DISi)
     call self%register_state_variable(self%id_det2Si,'det2Si','mmol Si/m^3','(fast sinking) detrital silica', & 
                                    minimum=0.0_rk,vertical_movement=self%det2_w*s2d,no_river_dilution=.false.)
     call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_det2Si)
   end if 
   
   if (self%resolve_bac) then 
     call self%register_state_variable(self%id_bacC,'bacC','mmol C/m^3','carbon bound to bacteria', &
                                    minimum=0.0_rk,no_river_dilution=.true.)
     call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_bacC)
     if (self%resolve_bacnp) then
       call self%register_state_variable(self%id_bacP,'bacP','mmol P/m^3','phosphorus bound to bacteria', &
                                    minimum=0.0_rk,no_river_dilution=.true.)
       call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_bacP)
       call self%register_state_variable(self%id_bacN,'bacN','mmol N/m^3','nitrogen bound to bacteria', &
                                    minimum=0.0_rk,no_river_dilution=.true.)
       call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_bacN)
     else
       call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_bacC,scale_factor=1./self%b_rcp)
       call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_bacC,scale_factor=1./self%b_rcn)
     end if
   end if
   
   ! Register link to external pools
   ! Register diagnostic variables
   !C
   !call self%register_horizontal_diagnostic_variable(self%id_air_o2c,'air_o2c','mmolC/m^2/d', 'DIC (CO2) air flux',   &
   !                                       output=output_time_step_averaged)

   call self%register_diagnostic_variable(self%id_det1QP,'det1_QP','molP/molC', 'molar P:C ratio of slow sinking detritus',   &
                                          output=output_instantaneous)                                  
   call self%register_diagnostic_variable(self%id_det2QP,'det2_QP','molP/molC', 'molar P:C ratio of fast sinking detritus',   &
                                          output=output_instantaneous)
   call self%register_horizontal_diagnostic_variable(self%id_DIP_Sflux,'dip_f_surf','mmolP/m^2/d', 'DIP air flux',   &
                                          output=output_time_step_averaged)                                      
   !N
   call self%register_diagnostic_variable(self%id_det1QN,'det1_QN','molN/molC', 'molar N:C ratio of slow sinking detritus',   &
                                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_det2QN,'det2_QN','molN/molC', 'molar N:C ratio of fast sinking detritus',   &
                                          output=output_instantaneous)
   call self%register_horizontal_diagnostic_variable(self%id_DINO3_Sflux,'dino3_f_surf','mmolN/m^2/d', 'NO3 aire flux',   &
                                          output=output_time_step_averaged)
   call self%register_horizontal_diagnostic_variable(self%id_DINH4_Sflux,'dinh4_f_surf','mmolN/m^2/d', 'NH4 air flux',   &
                                          output=output_time_step_averaged)
   
   !O2
   !call self%register_diagnostic_variable(self%id_O2_Sat,   'O2_Sat','mmolO2/m^3','saturation-O2', &
   !                                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_O2_percSat,   'O2_percSat','%','O2 percent saturation', &
                                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_o2o_bac,   'o2o_bac','mmol/m^3/d','O2 consumption by bacteria', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_o2o_n4n,   'o2o_n4n','mmol/m^3/d','O2 consumption by nitrification', &
                                          output=output_time_step_averaged)                                       
   call self%register_diagnostic_variable(self%id_n4n_n3n,   'n4n_n3n','mmolN/m^2/d','nitrification rate', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_n3n_nn2,   'n3n_nn2','mmolN/m^2/d','denitrification rate', &
                                          output=output_time_step_averaged)  
   call self%register_horizontal_diagnostic_variable(self%id_air_o2o,'air_o2o','mmol/m^2/d', ' O2 air flux',   &
                                          output=output_time_step_averaged)                                          
   
   !Si
   if (self%resolve_Si) then
     call self%register_diagnostic_variable(self%id_det2QSi,'det2_QSi','molSi/mol C', 'molar Si:C ratio of fast sinking detritus',   &
                                          output=output_instantaneous)
     call self%register_horizontal_diagnostic_variable(self%id_DISi_Sflux,'DISi_f_surf','mmolSi/m^2/d', 'DISi air flux',   &
                                          output=output_time_step_averaged)
   end if
   if (self%resolve_bac .and. self%resolve_bacnp) then
     call self%register_diagnostic_variable(self%id_bacQP,'bac_QP','molP/molC', 'molar P:C ratio of bac',   &
                                          output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_bacQN,'bac_QN','molN/molC', 'molar N:C ratio of bac',   &
                                          output=output_instantaneous)                                           
   end if
   
   if (self%flux_diags) then
     call self%register_diagnostic_variable(self%id_det1_dom_C,'det1_dom_C','mmolC/m^3/d', 'molar C flux from det1 to DOM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_det1_dom_P,'det1_dom_P','mmolP/m^3/d', 'molar P flux from det1 to DOM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_det1_dom_N,'det1_dom_N','mmolN/m^3/d', 'molar N flux from det1 to DOM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_det2_dom_C,'det2_dom_C','mmolC/m^3/d', 'molar C flux from det2 to DOM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_det2_dom_P,'det2_dom_P','mmolP/m^3/d', 'molar P flux from det2 to DOM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_det2_dom_N,'det2_dom_N','mmolN/m^3/d', 'molar N flux from det2 to DOM', &
                                            output=output_time_step_averaged)
     if (self%resolve_Si) then
       call self%register_diagnostic_variable(self%id_det2_dim_Si,'det2_dim_Si','mmolSi/m^3/d', 'molar Si flux from det2 to DIM', &
                                            output=output_time_step_averaged)
     end if
     call self%register_diagnostic_variable(self%id_soc_doc,'soc_doc','mmolC/m^3/d', 'molar C flux from SOC to DOC', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dom_dim_C,'dom_dim_C','mmolC/m^3/d', 'molar C flux from DOM to DIM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dom_dim_P,'dom_dim_P','mmolP/m^3/d', 'molar P flux from DOM to DIM', &
                                            output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_dom_dim_N,'dom_dim_N','mmolN/m^3/d', 'molar N flux from DOM to DIM', &
                                            output=output_time_step_averaged)
     if (self%resolve_bac) then                                            
       call self%register_diagnostic_variable(self%id_dom_bac_C,'dom_bac_C','mmolC/m^3/d', 'molar C flux from DOM to bacteria', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_dom_bac_P,'dom_bac_P','mmolP/m^3/d', 'molar P flux from DOM to bacteria', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_dom_bac_N,'dom_bac_N','mmolN/m^3/d', 'molar N flux from DOM to bacteria', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_dim_bac_C,'dim_bac_C','mmolC/m^3/d', 'molar C flux from DIM to bacteria', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_dim_bac_P,'dim_bac_P','mmolP/m^3/d', 'molar P flux from DIM to bacteria', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_dim_bac_N,'dim_bac_N','mmolN/m^3/d', 'molar N flux from DIM to bacteria', &
                                            output=output_time_step_averaged)    
       call self%register_diagnostic_variable(self%id_bac_dim_C,'bac_dim_C','mmolC/m^3/d', 'molar C flux from bacteria to DIM', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_bac_dim_P,'bac_dim_P','mmolP/m^3/d', 'molar P flux from bacteria to DIM', &
                                            output=output_time_step_averaged)
       call self%register_diagnostic_variable(self%id_bac_dim_N,'bac_dim_N','mmolN/m^3/d', 'molar N flux from bacteria to DIM', &
                                            output=output_time_step_averaged)                                             
     end if                                       
   end if
   
   ! Register conserved quantities

   ! Register environmental dependencies
   
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity) !psu
   !if cwind is not set to a positive value, it is expected from the physical driver
   if (self%cwind .lt. 0.0) then
     write(*,*)'EH/abio_pel_EH: Constant wind (cwind) not provided (or assigned a negative value). Expected from the physical driver'
     call self%register_horizontal_dependency(self%id_wind,standard_variables%wind_speed) !m/s
   else
     write(*,*)'EH/abio_pel_EH: Prescribed constant wind speed (cwind=',self%cwind,') will be used for calculating O2 exchange'
   end if
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_zmax,standard_variables%bottom_depth)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_diagnostic_variable(self%id_dpar,'PAR','W/m^2',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_kw,'kw','/m',       'background attenuation_coefficient')
   call self%register_dependency(self%id_env_k0,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
   call self%register_diagnostic_variable(self%id_k0,'k0','/m',       'attenuation coefficient before setting attenuation within EH')
   call self%register_dependency(self%id_I0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_horizontal_diagnostic_variable(self%id_dI0,'I0','W/m^2',       'PAR at the surface')
   if (self%Idm_met .eq. 1) then
     call self%register_dependency(self%id_I0_dmean,temporal_mean(self%id_I0,period=1._rk*86400._rk,resolution=1._rk))
     call self%register_horizontal_diagnostic_variable(self%id_dI0_dmean,'I0dm','W/m^2',       'daily mean PAR at the surface')
   else if (self%Idm_met .eq. 2) then
     call self%register_dependency(self%id_par_dmean,temporal_mean(self%id_par,period=1._rk*86400._rk,resolution=1._rk))
     call self%register_diagnostic_variable(self%id_dpar_dmean,'PARdm','W/m^2',       'daily mean PAR')
   end if

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model
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
   class (type_gpm_abio_pel_EH), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   !Deps   
   type(type_env)             :: env
   !State variables
   type(type_elms)            :: det1,det2,dom,bac 
   type(type_dim)             :: di
   !Rates
   type(type_elms)            :: det1_dom,det2_dom !decay of detritus to dom
   type(type_dim)             :: dom_dim !remineralization of dom to di
   type(type_elms)            :: dom_bac !takeup of dom and dim by bacteria
   type(type_dim)             :: bac_dim,dim_bac !excretion and takeup of dim by bac by bacteria
   real(rk)                   :: fT,O2,soc,par,par_dmean
   real(rk)                   :: o2_sat,n4n_n3n,n3n_nn2,o2o_bac,o2o_n4n,d2s_dim_Si,soc_doc
!EOP
!-----------------------------------------------------------------------
!BOC
   if (debug) then
     _GET_GLOBAL_(self%id_doy,env%doy) !day of year
     write(*,'(A,2F9.5)')'doy:',env%doy
   end if
   
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   !-------------------------------------------------------------------------
   !PREPARE: retrieve variables
   !
   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,env%temp)  ! temperature
   _GET_(self%id_salt, env%salt)  ! salinity (psu)
   
   ! Retrieve current (local) state variable values.
   !C
   _GET_(self%id_det1C,det1%C) ! slow sinking detritus C
   _GET_(self%id_det2C,det2%C) ! fast sinking detritus C
   _GET_(self%id_DOC,dom%C) ! DOC
   !if (self%resolve_DIC) then
   !_GET_(self%id_DIC,di%C)   ! di C
   !end if
   if (self%metCexc .ne. 0) then
     _GET_(self%id_SOC,soc) ! soc
   end if
   !P
   _GET_(self%id_det1P,det1%P) ! slow sinking detritus C
   _GET_(self%id_det2P,det2%P) ! fast sinking detritus C
   _GET_(self%id_DOP,dom%P) ! DOP
   _GET_(self%id_DIP,di%P)   ! di P
   !N
   _GET_(self%id_det1N,det1%N) ! slow sinking detritus C
   _GET_(self%id_det2N,det2%N) ! fast sinking detritus C
   _GET_(self%id_DON,dom%N) ! DON
   _GET_(self%id_DINO3,di%NO3)   ! NO3
   _GET_(self%id_DINH4,di%NH4)   ! NH4
   !O2
   _GET_(self%id_O2,O2)   ! O2
   !Si
   if (self%resolve_Si) then
     det1%Si=0.0 !i.e., no such thing
     _GET_(self%id_det2Si,det2%Si) ! fast sinking detritus Si
   end if
   if (self%resolve_bac) then
     _GET_(self%id_bacC,bac%C) ! C bound to bac
     if (self%resolve_bacnp) then
       _GET_(self%id_bacP,bac%P) ! P bound to bac
       _GET_(self%id_bacN,bac%N) ! N bound to bac
     else
       bac%P=bac%C/self%b_rcp
       bac%N=bac%C/self%b_rcn
     end if
   end if
   !
   !END OF PREPARE
   !-------------------------------------------------------------------------
   
   !-------------------------------------------------------------------------
   !START OF CALCULATIONS: no other calculation outside this box
   !
   !Temperature function
   fT = get_fQ10(self,env%temp)
   
   !Detritus decay: det_dom: mmol/m3/d
   call get_detdec(self,det1,self%r_det1dec,fT,det1_dom) !r_det1dec=xmu4n in EH
   call get_detdec(self,det2,self%r_det2dec,fT,det2_dom) !r_det2dec=xmu5n in EH
   if (self%resolve_Si) then
     d2s_dim_Si= self%r_det2dec * fT * det2%Si * self%r_detdec_fC/10. !in EH, this is d2s_n5s (L1479)
   end if
   
   !Decay of soc to doc
   soc_doc=get_soc_doc(self,soc,fT)
   
   !Direct remineralization of dom to dim
   call get_rem(self,fT,dom,dom_dim)
   
   !Bacterial uptake and release: dom_bac,dim_bac,bac_dim: mmol/m3/d
   if (self%resolve_bac) then
     call get_bac_rates(self,bac,dom,di,fT,dom_bac,dim_bac,bac_dim)
   else
     bac_dim%C=0.0 !needed to calculate o2 consumption
   end if
   
   !Nitfication, Denitrification, Oxygen consumption
   call get_nit_denit_o2cons(self,di,O2,fT,bac_dim%C,dom_dim%C,n4n_n3n,n3n_nn2,o2o_bac,o2o_n4n)
   
   !O2 sat (required for diagnostics)
   o2_sat=osat_weiss(env%temp,env%salt)  ! saturation oxygen in mmolO2/m**3
   !
   !END OF CALCULATIONS: no other calculation outside this box
   !-------------------------------------------------------------------------
   
   !-------------------------------------------------------------------------
   !WRITE
   !
   ! Set temporal derivatives
   !detC
   _SET_ODE_(self%id_det1C,-det1_dom%C)
   _SET_ODE_(self%id_det1P,-det1_dom%P)
   _SET_ODE_(self%id_det1N,-det1_dom%N)
   _SET_ODE_(self%id_det2C,-det2_dom%C)
   _SET_ODE_(self%id_det2P,-det2_dom%P)
   _SET_ODE_(self%id_det2N,-det2_dom%N)
   if (self%resolve_Si) then
     _SET_ODE_(self%id_det2Si,-d2s_dim_Si)
   end if
   
   !soc
   _SET_ODE_(self%id_SOC,-soc_doc)
   
   !dom
   _SET_ODE_(self%id_DOC, det1_dom%C+det2_dom%C-dom_dim%C+soc_doc)
   _SET_ODE_(self%id_DOP, det1_dom%P+det2_dom%P-dom_dim%P)
   _SET_ODE_(self%id_DON, det1_dom%N+det2_dom%N-dom_dim%NH4)
   
   !dim
   if (self%resolve_DIC) then
     _SET_ODE_(self%id_DIC, dom_dim%C)
   end if
   _SET_ODE_(self%id_DIP, dom_dim%P)
   _SET_ODE_(self%id_DINO3, n4n_n3n-n3n_nn2) !dom_dim%NO3=0.0
   _SET_ODE_(self%id_DINH4, dom_dim%NH4-n4n_n3n)
   if (self%resolve_Si) then
     _SET_ODE_(self%id_DISi, d2s_dim_Si) !this is s2d_n5s in EH
   end if
   
   !O2
   _SET_ODE_(self%id_O2,-o2o_bac-o2o_n4n)
   
   if (self%resolve_bac) then
     !dom
     _SET_ODE_(self%id_DOC, -dom_bac%C)
     _SET_ODE_(self%id_DOP, -dom_bac%P)
     _SET_ODE_(self%id_DON, -dom_bac%N)
     
     !bac
     _SET_ODE_(self%id_bacC, dom_bac%C+dim_bac%C-bac_dim%C) 
     if (self%resolve_bacnp) then
       _SET_ODE_(self%id_bacP, dom_bac%P+dim_bac%P-bac_dim%P)
       _SET_ODE_(self%id_bacN, dom_bac%N+dim_bac%NH4-bac_dim%NH4)
     end if
     
     !dim
     if (self%resolve_DIC) then
       _SET_ODE_(self%id_DIC, bac_dim%C) ! no dim_bac%C
     end if
     _SET_ODE_(self%id_DIP, bac_dim%P-dim_bac%P) 
     _SET_ODE_(self%id_DINH4, bac_dim%NH4-dim_bac%NH4)
   end if
   
   ! Export diagnostic variables
   !P
   _SET_DIAGNOSTIC_(self%id_det1QP,det1%P/det1%C)
   _SET_DIAGNOSTIC_(self%id_det2QP,det2%P/det2%C)
   !N
   _SET_DIAGNOSTIC_(self%id_det1QN,det1%N/det1%C)
   _SET_DIAGNOSTIC_(self%id_det2QN,det2%N/det2%C)
   !O
   !_SET_DIAGNOSTIC_(self%id_O2_Sat,o2_sat)
   _SET_DIAGNOSTIC_(self%id_O2_percSat,100*O2/o2_sat)
   _SET_DIAGNOSTIC_(self%id_o2o_bac,o2o_bac*s2d)
   _SET_DIAGNOSTIC_(self%id_o2o_n4n,o2o_n4n*s2d)
   _SET_DIAGNOSTIC_(self%id_n4n_n3n,n4n_n3n*s2d)
   _SET_DIAGNOSTIC_(self%id_n3n_nn2,n3n_nn2*s2d)
   !Si
   if (self%resolve_Si) then
     _SET_DIAGNOSTIC_(self%id_det2QSi,det2%Si/det2%C)
   end if
   if (self%resolve_bac .and. self%resolve_bacnp) then
     _SET_DIAGNOSTIC_(self%id_bacQP,bac%P/bac%C)
     _SET_DIAGNOSTIC_(self%id_bacQN,bac%N/bac%C)
   end if
   
   if (self%flux_diags) then
   !id_dim_bac_C,id_bac_dim_C
     _SET_DIAGNOSTIC_(self%id_det1_dom_C,det1_dom%C*s2d)
     _SET_DIAGNOSTIC_(self%id_det1_dom_P,det1_dom%P*s2d)
     _SET_DIAGNOSTIC_(self%id_det1_dom_N,det1_dom%N*s2d)
     _SET_DIAGNOSTIC_(self%id_det2_dom_C,det2_dom%C*s2d)
     _SET_DIAGNOSTIC_(self%id_det2_dom_P,det2_dom%P*s2d)
     _SET_DIAGNOSTIC_(self%id_det2_dom_N,det2_dom%N*s2d)
     _SET_DIAGNOSTIC_(self%id_det2_dim_Si,d2s_dim_Si*s2d)
     _SET_DIAGNOSTIC_(self%id_soc_doc,soc_doc*s2d)
     _SET_DIAGNOSTIC_(self%id_dom_dim_C,dom_dim%C*s2d)
     _SET_DIAGNOSTIC_(self%id_dom_dim_P,dom_dim%P*s2d)
     _SET_DIAGNOSTIC_(self%id_dom_dim_N,dom_dim%NH4*s2d)
     if (self%resolve_bac) then
       _SET_DIAGNOSTIC_(self%id_dom_bac_C,dom_bac%C*s2d)
       _SET_DIAGNOSTIC_(self%id_dom_bac_P,dom_bac%P*s2d)
       _SET_DIAGNOSTIC_(self%id_dom_bac_N,dom_bac%N*s2d)
       _SET_DIAGNOSTIC_(self%id_dim_bac_C,dim_bac%C*s2d)
       _SET_DIAGNOSTIC_(self%id_dim_bac_P,dim_bac%P*s2d)
       _SET_DIAGNOSTIC_(self%id_dim_bac_N,dim_bac%NH4*s2d)
       _SET_DIAGNOSTIC_(self%id_bac_dim_C,bac_dim%C*s2d)
       _SET_DIAGNOSTIC_(self%id_bac_dim_P,bac_dim%P*s2d)
       _SET_DIAGNOSTIC_(self%id_bac_dim_N,bac_dim%NH4*s2d)
     end if
   end if
   
   ! Light
   _GET_(self%id_par,par) ! local photosynthetically active radiation
   _SET_DIAGNOSTIC_(self%id_dpar, par) ! mol/m2/d-1
   if (self%Idm_met .eq. 2) then
     _GET_(self%id_par_dmean,par_dmean)  ! surface short wave radiation, daily average
     _SET_DIAGNOSTIC_(self%id_dpar_dmean, max(0.0, par_dmean)) ! W/m2
   end if 
      
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
! !IROUTINE: Nitfication, Denitrification, O2 Consumption
!
! !INTERFACE:  
  real(rk) function get_soc_doc(self,soc,fT)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   real(rk),intent(in)                     :: soc,fT
!EOP
!-----------------------------------------------------------------------
!BOC
! 
   
  select case (self%metCexc)
    case default
      call self%fatal_error('abio_pel_EH.F90/get_soc_doc','for '//trim(self%name)// ' specified metCexc option is not available')
    case (0)
      get_soc_doc=0.0
    case (1) !EH
      get_soc_doc = soc*self%soc_rate !*iexcess (this is obsolete)
    !todo: *fT?
    !case (2) !Schartau et al 2007
  end select
   
  end function
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nitfication, Denitrification, O2 Consumption
!
! !INTERFACE:  
   subroutine get_nit_denit_o2cons(self,di,O2,fT,bac_dic,doc_dic,n4n_n3n,n3n_nn2,o2o_bac,o2o_n4n)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   type(type_dim), intent(in)   :: di
   real(rk),intent(in)           :: O2,fT,bac_dic,doc_dic
   real(rk)                      :: n4n_n3n,n3n_nn2,o2o_bac,o2o_n4n
! !LOCAL VARIABLES
   real(rk)                      :: osw,nsw,anitdep
   !
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   !oxic-anoxic switch (EH L.589-592)
   osw=1.0
   if (O2 <= 0.) osw=0.0
   nsw=1.0
   if (di%NO3 <= 0.1) nsw=0.0
   
   !nitrification
   anitdep=1.0 !todo: understand what this is (L1632)
   n4n_n3n = osw*self%xknit*di%NH4*fT*anitdep
   
   !denitrification
   n3n_nn2 = 0.5 * (1.0-osw)*nsw*(bac_dic/self%b_rcn+doc_dic/self%b_rcn)
   
   !oxic nitrification
   o2o_n4n = 2.0* n4n_n3n 
   
   !remineralization (bacterial or not)
   o2o_bac = osw*(bac_dic+doc_dic) &
           + (1.-osw)*(1-nsw)*(bac_dic+doc_dic)
   
   
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bacterial uptake
!
! !INTERFACE:  
  subroutine get_bac_rates (self,bac,dom,di,fT,dom_bac,dim_bac,bac_dim)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   type(type_elms),intent(in)    :: bac
   type(type_elms),intent(in)    :: dom
   type(type_dim),intent(in)     :: di
   real(rk), intent(in)          :: fT
   type(type_dim)               :: dim_bac,bac_dim
   type(type_elms)               :: dom_bac
! !LOCAL VARIABLES
   real(rk)          :: fdon,fdop,fbc,fbn,fdc_diff
   real(rk)          :: fbpmax,fbpreq,fbpdiff
   real(rk)          :: tiny_p,tiny_n
   real(rk),parameter :: tres_n1p=0.1_rk
   real(rk),parameter :: tres_n4n=0.1_rk
   real(rk),parameter :: tiny_c=1e-1
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   !reset all
   dom_bac%C=0.0; dom_bac%N=0.0; dom_bac%P=0.0
   dim_bac%C=0.0; dim_bac%NH4=0.0; dim_bac%P=0.0 !dim_bac%NO3=0.0; 
   bac_dim%C=0.0; bac_dim%NH4=0.0; bac_dim%P=0.0 !bac_dim%NO3=0.0;
   
   tiny_p=tiny_c/self%b_rcp
   tiny_n=tiny_c/self%b_rcn
   
   !uptake of DOM
   fdon = dom%N/(self%b_xk4+dom%N) !xk4
   fdop = dom%P/(self%b_xkp+dom%P)  !xkpb
   
     dom_bac%N = fT*self%b_v*fdon*bac%N !mmolN/m3/d
     dom_bac%C = dom_bac%N * dom%C/dom%N !mmolC/m3/d
   
   if (dom%P .gt. tiny_p) then
     dom_bac%P=fT*self%b_v*fdop*bac%P !mmolP/m3/d
   else
     dom_bac%P = 0.0_rk
   end if   
   
   !uptake/excretion of DIM (new_bac)
   !N&C
   !excretion
   if (bac%N .gt. tiny_n .and. bac%C .gt. tiny_c) then
     bac_dim%C=fT*self%b_exc*bac%C !b_exc=xmu3
     bac_dim%NH4=fT*self%b_exc*bac%N
   end if
   !uptake
   !dim_bac%NH4=0.0_rk !EH: L1528 
   fbc=dom_bac%C - bac_dim%C !-bac_z%C 
   fbn=dom_bac%N  - bac_dim%NH4 !+ dim_bac%NH4-bac_z%N
   fdc_diff=fbc - fbn*self%b_rcn
   if (fdc_diff .gt. 0.0_rk) then !more C than N
     if (di%NH4 <= tres_n4n) then 
       dim_bac%NH4=0.0_rk
       bac_dim%C=bac_dim%C+fdc_diff
       !write(*,*)'abioP,L758: fbc-fbn*bC/bN',dom_bac%C-bac_dim%C - (dom_bac%N+ dim_bac%NH4-bac_dim%NH4)*self%b_rcn
     else !balance by taking up more NH4:
       !original: I don't see why it should balance
       !dim_bac%NH4= dom_bac%C/self%b_rcn - dom_bac%N
       !new:
       dim_bac%NH4 = fdc_diff/self%b_rcn
       !write(*,*)'abioP,L764: fbc-fbn*bC/bN',dom_bac%C-bac_dim%C - (dom_bac%N+ dim_bac%NH4-bac_dim%NH4)*self%b_rcn
     end if
   else !more N than C: excrete NH4
     bac_dim%NH4 = bac_dim%NH4 + fdc_diff/self%b_rcn 
     !write(*,*)'abioP,L768: fbc-fbn*bC/bN',dom_bac%C-bac_dim%C - (dom_bac%N+ dim_bac%NH4-bac_dim%NH4)*self%b_rcn
   end if
   
   !P
   fbpreq = (dom_bac%N + dim_bac%NH4 - bac_dim%NH4) * self%b_rcn/self%b_rcp - dom_bac%P
   if (fbpreq .le. 0.0_rk) then !excrete excess P
     bac_dim%P = abs(fbpreq)
   else !excess C, try balancing by taking up from the dip pool
     fbpmax=self%b_v * di%P
     if (di%P <= tres_n1p) fbpmax=0.0_rk
     dim_bac%P=min(fbpreq,fbpmax)
     if (fbpreq .gt. fbpmax) then !excrete excess N and C
       fbpdiff=fbpreq - fbpmax
       bac_dim%NH4 = bac_dim%NH4 + fbpdiff * self%b_rcp/self%b_rcn
       bac_dim%C = bac_dim%C + fbpdiff * self%b_rcp
     end if
   end if
   
   !if ( self%resolve_bacnp ) then
   !  write(*,*)'abioP,L827: fbc-fbn*bC/bN',dom_bac%C-bac_dim%C - (dom_bac%N+ dim_bac%NH4-bac_dim%NH4)*self%b_rcn
   !  write(*,*)'abioP,L828: fbc-fbp*bC/bP',dom_bac%C-bac_dim%C - (dom_bac%P+dim_bac%P-bac_dim%P)*self%b_rcp
   !else 
   !  write(*,'(A,3F7.3)')'abiP.l830: bac_domN_net, dom_bacN_net, sum:',&
   !         s2d*(dom_bac%C)/self%b_rcn,& !+dim_bac%C-bac_dim%C)/self%b_rcn,&
   !         s2d*(-dom_bac%N),& !+bac_dim%NH4-dim_bac%NH4),&
   !         s2d*(dom_bac%C)/self%b_rcn+s2d*(-dom_bac%N)
            
   ! write(*,'(A,3F7.3)')'abiP.l836: bac_dimN_net, dim_bacN_net, sum:',&
   !         s2d*(dim_bac%C-bac_dim%C)/self%b_rcn,&
   !         s2d*(bac_dim%NH4-dim_bac%NH4),&
   !         s2d*((dim_bac%C-bac_dim%C)/self%b_rcn + (bac_dim%NH4-dim_bac%NH4))
   ! end if        
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: direct remineralization of dom to dim
!
! !INTERFACE:  
   subroutine get_rem(self,fT,dom,dom_dim)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   real(rk), intent(in)          :: fT
   type(type_elms), intent(in)   :: dom
   type(type_dim)               :: dom_dim
   
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   !assume that remineralization (dom->dim) is resolved EITHER through bacteria OR directly
   if (self%resolve_bac) then
     dom_dim%C=0.0
     dom_dim%P=0.0
     dom_dim%NO3=0.0
     dom_dim%NH4=0.0
   else
     dom_dim%C=self%rhoC*fT*dom%C 
     dom_dim%P=self%rhoP*fT*dom%P
     dom_dim%NH4=self%rhoN*fT*dom%N
     dom_dim%NO3=0.0 !remineralization produces only ammounium
   end if

  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get detritus decay rate
! !INTERFACE:  
   subroutine get_detdec(self,det,sr,fT,br)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:   
   type (type_gpm_abio_pel_EH), intent(in) :: self
   type(type_elms),intent(in)    :: det
   real(rk), intent(in)          :: sr !specific rate
   real(rk), intent(in)          :: fT
   type(type_elms)               :: br !bulk rates
! !LOCAL VARIABLES
   !real(rk)                      :: fy
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   
   br%P  = sr * fT * det%P !sr = r_det1dec, r_det2dec, which, in EH are xmu4n, xmu5n
   br%N  = sr * fT * det%N
   br%C  = sr * fT * det%C  * self%r_detdec_fC !in EH, r_detdec_fC = rxmu4c
  end subroutine 
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Temperature function
!
! !INTERFACE:
   pure real(rk) function get_fQ10(self,temp)
!
! !DESCRIPTION:
! Here, temperature response functions are formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   real(rk), intent(in)         :: temp!,Q10,Tref
!
!EOP
!-----------------------------------------------------------------------
!BOC

   get_fQ10=self%Q10**((temp-self%Tref)/self%Tref)
 
   end function 
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Weiss formula for the saturation oxygen (osat)
!
! !INTERFACE:
   real(rk) function osat_weiss(t,s)
!
! !DESCRIPTION:
! Weiss formula for the saturation oxygen (osat), copied from gotm/ergom
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  type (type_gotm_ergom), intent(in) :: self
  real(rk), intent(in)                 :: t,s
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  real(rk)                 :: tk
  real(rk)                 :: aa1=-173.4292_rk
  real(rk)                 :: aa2=249.6339_rk
  real(rk)                 :: a3=143.3483_rk
  real(rk)                 :: a4=-21.8492_rk
  real(rk)                 :: b1=-0.033096_rk
  real(rk)                 :: b2=0.014259_rk
  real(rk)                 :: b3=-0.001700_rk
  real(rk)                 :: kelvin=273.16_rk
  real(rk)                 :: mol_per_liter=44.661_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   tk=(t+kelvin)*0.01_rk
   osat_weiss=exp(aa1+aa2/tk+a3*log(tk)+a4*tk    &
              +s*(b1+(b2+b3*tk)*tk))*mol_per_liter
   return
   end function osat_weiss
!EOC
!-----------------------------------------------------------------------   

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Right hand sides of surface exchange model
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
! !IROUTINE: Process interaction betdet_ween benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically imples variable units * m/s
! [bottom fluxes] for the pelagic, and variable units/s [temporal derivatives] for the benthos.
! Positive values denote state variable increases, negative values state variable decreases.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel_EH), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: temp,salt,wnd,sc,p_vel,O2,O2_sat,air_o2o !,!air_c2o
   real(rk) :: I0,I0_dmean
!EOP
!-----------------------------------------------------------------------
!BOC
   
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_
   
   ! nut flux due to rivers, surface runoff, deposition
   !P
   _SET_SURFACE_EXCHANGE_(self%id_DIP, self%DIP_Sflux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DIP_Sflux, self%DIP_Sflux)
   !N
   _SET_SURFACE_EXCHANGE_(self%id_DINO3, self%DINO3_Sflux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DINO3_Sflux, self%DINO3_Sflux)
   _SET_SURFACE_EXCHANGE_(self%id_DINH4, self%DINH4_Sflux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DINH4_Sflux, self%DINH4_Sflux)
   !Si
   if (self%resolve_Si) then
     _SET_SURFACE_EXCHANGE_(self%id_DISi, self%DISi_Sflux)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DISi_Sflux, self%DISi_Sflux)
   end if
   
   !CO2:
   !call EH_air_sea_co2(ierr,pCO2a) !requires resolution of CO2
   !_SET_SURFACE_EXCHANGE_(self%id_DIC,air_o2c) 
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_air_o2c,air_o2c)
   
   !O2
   _GET_(self%id_O2, O2)   ! sea water dissolved oxygen in mmolO2/m**3
   _GET_(self%id_temp, temp)  ! water temperature (C)
   _GET_(self%id_salt, salt)  ! salinity (psu)
   !if not set a constant value, get wind speed
   if (self%cwind .lt. 0.0) then
     _GET_HORIZONTAL_(self%id_wind,wnd)
   else
     wnd=self%cwind
   end if
   !write(*,*)'wind:',wnd
   
   select case (self%metairsea)
     case default
       call self%fatal_error('abio_pel_EH.F90/do_surface','for '//trim(self%name)// ' specified metairsea option is not available')
     case (0)
       air_o2o=0.0_rk
       
     case (1) !as in EH
       call EH_air_sea_o2(self,O2,temp,salt,wnd,air_o2o)
       
     case (2) !ERGOM
       call ERGOM_air_sea_o2(self,O2,temp,salt,wnd,air_o2o)
       
   end select
   _SET_SURFACE_EXCHANGE_(self%id_O2, air_o2o)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_air_o2o, air_o2o*s2d) ! converts mmol/m2.s to mmol/m2.d
  
   !Light
   _GET_HORIZONTAL_(self%id_I0,I0)  ! surface short wave radiation
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dI0, I0) ! W/m2
   if (self%Idm_met .eq. 1) then
     _GET_HORIZONTAL_(self%id_I0_dmean,I0_dmean)  ! surface short wave radiation, daily average
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dI0_dmean, max(0.0, I0_dmean)) ! W/m2
   end if
   
   !Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: calculation of air-sea exchange of O2
! (copied from gotm/ergom)
!
!-----------------------------------------------------------------------
!
! !INTERFACE:

   subroutine ERGOM_air_sea_o2(self,O2,temp,salt,wnd,air_o2o)
   
   implicit none
   
   ! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   real(rk),     intent(in)     :: O2,temp,salt,wnd
! !OUTPUT PARAMETERS:
   real(rk)                     :: O2_Sat,air_o2o
! !LOCAL VARIABLES
   real(rk)                     :: sc,p_vel
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !o2 saturation
   O2_Sat = osat_weiss(temp,salt)  ! saturation oxygen in mmolO2/m**3

   !exchange coefficient (piston velocity)
   sc=1450.+(1.1*temp-71.0_rk)*temp
   if (wnd .gt. 13.0_rk) then
     p_vel = 5.9_rk*(5.9_rk*wnd-49.3_rk)/sqrt(sc)
   else
     if (wnd .lt. 3.6_rk) then
       p_vel = 1.003_rk*wnd/(sc)**(0.66_rk)
     else
       p_vel = 5.9_rk*(2.85_rk*wnd-9.65_rk)/sqrt(sc)
     end if
   end if
   if (p_vel .lt. 0.05_rk) then
     p_vel = 0.05_rk
   end if
   p_vel = p_vel/s2d
      
   air_o2o =p_vel*(O2_Sat-O2)
   !write (*,'(A,5(F10.5))') 'pvel*s2d,O2_sat,O2,air_o2o=', p_vel,O2_sat,O2,air_o2o*s2d
   
   end subroutine ERGOM_air_sea_o2
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: calculation of air-sea exchange of O2
! (COPIED FROM ECOHAMv5: eco_meteo.f90: 1221-)
!
!-----------------------------------------------------------------------
!
! !INTERFACE:

   subroutine EH_air_sea_o2(self,O2O,temp,salt,u10,air_o2o)
!
! !USES:
   !use mod_var
   !use mod_flux
   implicit none
!
! !INPUT PARAMETERS:
   type (type_gpm_abio_pel_EH), intent(in) :: self
   real(rk),     intent(in)     :: O2O,temp,salt,u10
! !OUTPUT PARAMETERS:
   real(rk)                     :: air_o2o
! !LOCAL VARIABLES:
   real          :: transfer_velocity, sc
   real          :: x1, x2, beta_dm, tabs !,fas
   integer       :: i, j, k, k0, kdz
   !from eco_par.f90:
   real, parameter    :: temzer=273.16

   ! Wanninkhof, R., 1992 JGR 97:7373-7382 Schmidt Number for oxygen
   real, parameter    :: scao=1953.4
   real, parameter    :: scbo=128.0
   real, parameter    :: scco=3.9918
   real, parameter    :: scdo=0.050091

   !Unesco 1986 Solubility of oxygen [mumol dm-1] (improvement of Weiss 1970)
   real, parameter    :: c1o=-135.90205
   real, parameter    :: c2o= 1.575701e5
   real, parameter    :: c3o=-6.642308e7
   real, parameter    :: c4o= 1.243800e10
   real, parameter    :: c5o=-8.621949e11
   real, parameter    :: d1o= 0.017674
   real, parameter    :: d2o=-10.754
   real, parameter    :: d3o= 2140.7
!EOP
!-----------------------------------------------------------------------
!BOC

   !air_sea_mode=1  !Wanninkhof 92
   !air_sea_mode=2  !Wanninkhof & McGillis 99
   !air_sea_mode=3  !Nightingale et al 2000

   transfer_velocity = 0.0

   ! Schmidt Number
   sc = ( (-scdo*temp +scco) *temp -scbo )*temp +scao
   ! transfer velocity
   select case(self%metairseaEH)
     case default
       call self%fatal_error('abio_pel_EH.F90/EH_air_sea_o2','for '//trim(self%name)// ' specified metairseaEH option is not available')
     case (1) ! Wanninkhof 92
       transfer_velocity = 0.3*u10*u10/sqrt(sc/660.)
     case (2) ! Wanninkhof & McGillis 99
       ! correction 0.0280 to 0.0283 (7.7.04)
       transfer_velocity = 0.0283*u10*u10*u10/sqrt(sc/660.)
     case (3) ! Nightingale et al 2000
       transfer_velocity = (2.5*((4.9946e-4*temp+1.6256e-2)      &
                            *temp+.5246)+0.222*u10*u10 &
                            +0.333*u10)/sqrt(sc/660.)
    end select 
    transfer_velocity = 0.24*transfer_velocity      ! unit is cm/h must be m/d
    !solubility mumol dm-3 (tabs in kelvin)
    tabs = temp + temzer
    beta_dm = exp( ((((c5o/tabs+c4o)/tabs+c3o)/tabs+c2o)/tabs+c1o)  &
                       -salt*((d3o/tabs+d2o)/tabs+d1o) )
    x1 = beta_dm
    x2 = max(0.0,O2O)                    ! [mmol m-3]
    air_o2o = transfer_velocity*(x1-x2)/s2d         ! [mmol m-2 s-1]
    !commented out in the original code:
    !sst(j,1,i,io2o) = sst(j,1,i,io2o) + f_from_to(j,1,i,i_air_o2o) ! [mmol m-2 d-1]
    ! fas = transfer_velocity*(x1-x2)/s2d                 ! [mmol m-2 1-1]
    ! do k = 1, kdz
    !   f_from_to(j,k,i,i_air_o2o) = fas/dz_air_sea  ! [mmol m-3 d-1]
    !   sst(j,k,i,io2o) = sst(j,k,i,io2o) + f_from_to(j,k,i,i_air_o2o)
    ! enddo

   end subroutine EH_air_sea_o2
!EOC
!-----------------------------------------------------------------------

! !    This requires the chemie module:
! !-----------------------------------------------------------------------
! !BOP
! !
! ! !IROUTINE: calculation of air-sea exchange of CO2
! ! (COPIED FROM ECOHAM5: eco_meteo.f90, L.1313-)
! !
! !-----------------------------------------------------------------------
! !
! ! !INTERFACE:
!    subroutine EH_air_sea_co2(pco2,airo2c)
! !
! ! !USES:
!    implicit none
! !
! ! !USES:
!    use mod_hydro,  only : temp, salt, rho
!    use mod_chemie, only : ch, itdic, io2c
!    use mod_var
!    use mod_flux
! 
!    implicit none
! !
! ! !INPUT PARAMETERS:
! !
! ! !OUTPUT PARAMETERS:
!    real(rk), intent(in)     :: pco2
!    real(rk), intent(inout)  :: airo2c
! !
! ! !LOCAL VARIABLES:
!    real          :: dz_air_sea, transfer_velocity, sc
!    real          :: x1, x2, tabs, alphas !,fas
!    integer       :: i, j, k, k0, kdz
! !
! !-----------------------------------------------------------------------

!    !air_sea_mode=1  !Wanninkhof 92
!    !air_sea_mode=2  !Wanninkhof & McGillis 99
!    !air_sea_mode=3  !Nightingale et al 2000

!    transfer_velocity = 0.0
! 
!    ! Schmidt Number
!    sc = ( (-scdo*temp(j,1,i) +scco) *temp(j,1,i) -scbo )*temp(j,1,i) +scao
!    ! transfer velocity
!    if (air_sea_mode==1) then  ! Wanninkhof 92
!      transfer_velocity=(2.5*((4.9946e-4*temp(j,1,i)+1.6256e-2)     &
!                                        *temp(j,1,i)+.5246)+0.3*u10(j,i)      &
!                                        *u10(j,i))/sqrt(sc/660.)
!    endif
!    if (air_sea_mode==2) then  ! Wanninkhof & McGillis 99
!      transfer_velocity=(2.5*((4.9946e-4*temp(j,1,i)+1.6256e-2)     &
!                                        *temp(j,1,i)+.5246)+0.0280*u10(j,i)   &
!                                        *u10(j,i)*u10(j,i))/sqrt(sc/660.)
!    endif
!    if (air_sea_mode==3) then  ! Nightingale et al 2000
!      transfer_velocity=(2.5*((4.9946e-4*temp(j,1,i)+1.6256e-2)     &
!                                        *temp(j,1,i)+.5246)+0.222*u10(j,i)    &
!                                        *u10(j,i)+0.333*u10(j,i))/sqrt(sc/660.)
!    endif
!    transfer_velocity = 0.24*transfer_velocity      ! unit is cm/h must be m/d
!    !solubility (tabs in kelvin)
!    tabs = temp(j,1,i) + temzer
!    alphas = exp( -60.2409 + 9345.17/tabs + 23.3585*alog(tabs/100.)  &
!                  + salt(j,1,i)*(.023517 - .023656*tabs/100.        &
!                  + 1.e-8*47.036*tabs*tabs) )      ! alphas: [mumol/kg/muatm]
!    !x1 = alphas *pco2a(j,i) *rho(j,1,i) ! NOTE: in ECOHAM "rho" is kg/l = (kg/m3)/1000 !!!
!    x1 = alphas *pco2a *rho(j,1,i) ! NOTE: in ECOHAM "rho" is kg/l = (kg/m3)/1000 !!!
!    x2 = ch(j,1,i,io2c)            ! [mmol m-3]
!    f_from_to(j,1,i,i_air_o2c) = transfer_velocity*(x1-x2)         ! [mmol m-2 d-1]
!    sst(j,1,i,idic) = sst(j,1,i,idic) + f_from_to(j,1,i,i_air_o2c) ! [mmol m-2 d-1]
   
!    airo2c=0.0
!    return
!    end subroutine air_sea_co2
! !EOC
! !-----------------------------------------------------------------------

! !-----------------------------------------------------------------------
! !BOP
! !
! ! !IROUTINE: determine time and lat dependend pco2(air) following Drange's PhD
! ! (COPIED FROM ECOHAM5)
! !
! !-----------------------------------------------------------------------
! !
! ! !INTERFACE:
!    subroutine EH_airpco2(td_,pco2a)
! !
! ! !USES:
!    implicit none
! !
! ! !INPUT PARAMETERS:
!    real(rk), intent(in)     :: td_ !day of year
!    real(rk)                 :: pco2a
! ! !LOCAL VARIABLES:
!    !integer       :: i, j
!    integer       :: itd,site
!    real          :: pco2a_ns(12), off
! !
! !-----------------------------------------------------------------------
! #include "call-trace.inc"
!    
!    site=1 !0: NorthSea, 1: Indo Pacific
!    select case (site)
!      case (1)
!        !North Sea
!        pco2a_ns(1)  = 359.0
!        pco2a_ns(2)  = 360.0
!        pco2a_ns(3)  = 361.0
!        pco2a_ns(4)  = 362.0
!        pco2a_ns(5)  = 360.0
!        pco2a_ns(6)  = 357.0
!        pco2a_ns(7)  = 351.0
!        pco2a_ns(8)  = 347.0
!        pco2a_ns(9)  = 347.0
!        pco2a_ns(10) = 352.0
!        pco2a_ns(11) = 357.0
!        pco2a_ns(12) = 360.0
!      case (2)
!        !values for Indo-Pacific low lat
!        pco2a_ns(1)  = 360.4 ! 359.
!        pco2a_ns(2)  = 360.6 ! 360.
!        pco2a_ns(3)  = 360.8 ! 361.
!        pco2a_ns(4)  = 361.0 ! 362.
!        pco2a_ns(5)  = 360.6 ! 360.
!        pco2a_ns(6)  = 360.0 ! 357.
!        pco2a_ns(7)  = 358.8 ! 351.
!        pco2a_ns(8)  = 358.0 ! 347.
!        pco2a_ns(9)  = 358.0 ! 347.
!        pco2a_ns(10) = 359.0 ! 352.
!        pco2a_ns(11) = 360.0 ! 357.
!        pco2a_ns(12) = 360.6 ! 360.
!    end select
! 
!    ! traditional formula ECOHAM4
!    !if (year<=1995) then
!    !   off = real(year-1995)*1.32
!    !else
!    !   off = real(year-1995)*2.00
!    !endif
! 
!    ! improved formula by Bernhard Mayer
!    if (year<=1987) then
!       off = real(year-1995)*1.32
!    else if (year<=1992) then
!       off = real(year-1995)*1.05
!    else if (year<=1995) then
!       off = real(year-1995)*1.32
!    else
!       off = real(year-1995)*2.00
!    endif
! 
!    itd = min(12,nint(td_/30.)+1)
!    ! commented out as long as no spatial variation is included
!    !do i = iStartCalc, iEndCalc
!    !   do j = jStartCalc, jEndCalc
!    !      pco2a(j,i) = pco2a_ns(itd)+off  !pCO2 in muatm
!    !   enddo
!    !enddo
!    pco2a = pco2a_ns(itd)+off  !pCO2 in muatm
! 
!    return
! 
!    end subroutine airpco2
! !EOC
! !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
!
! !INTERFACE:
   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_gpm_abio_pel_EH), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
!
! !LOCAL VARIABLES:
   real(rk)                   :: vert_vel1,vert_vel2,temp,fv !,yearday,par
!
!EOP
!-----------------------------------------------------------------------
!BOC

! det_way to get the day of year
!_GET_DEPENDENCY_SCALAR_(self%id_yearday,yearday)  ! day of year (diff(d/m/Y - 1/1/Y))

  ! Enter spatial loops (if any)
  _LOOP_BEGIN_

  select case (self%det_velmet)
    case default
      call self%fatal_error('abio_pel_EH.F90/get_vertical_movement','for '//trim(self%name)// ' specified det_velmet option is not available')
    case (1) 
      vert_vel1=self%det1_w
      vert_vel2=self%det2_w
    case (2) !velocity is corrected by the viscosity, calculated as a function of temperature
      _GET_  (self%id_temp,temp)  ! temperature
      !in the expression belodet_w, denominator (10**..) gives the viscosity ratio at the ambient and reference (T=20 oC) temperatures (Kestin et al. 1978). So the overall expression is: vel*mu(20)/mu(T)
      fv=1/(10**(((20-temp)/(temp+96))*(1.2378 - 1.303e-3*(20-temp) + 3.06e-6*(20-temp)**2 + 2.55e-8*(20-temp)**3))) 
      vert_vel1=self%det1_w*fv
      vert_vel2=self%det2_w*fv
  end select

  !Set these calculated vertical_movement values for the current p 
  _SET_VERTICAL_MOVEMENT_(self%id_det1C,vert_vel1)
  _SET_VERTICAL_MOVEMENT_(self%id_det1P,vert_vel1) 
  _SET_VERTICAL_MOVEMENT_(self%id_det1N,vert_vel1)
  _SET_VERTICAL_MOVEMENT_(self%id_det2C,vert_vel2)
  _SET_VERTICAL_MOVEMENT_(self%id_det2P,vert_vel2) 
  _SET_VERTICAL_MOVEMENT_(self%id_det2N,vert_vel2)
  if (self%resolve_Si) then
    _SET_VERTICAL_MOVEMENT_(self%id_det2Si,vert_vel2)
  end if
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
! !USES:
   implicit none
!   
! !INPUT PARAMETERS:
   class(type_gpm_abio_pel_EH), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: det1C,det2C
   real(rk)                     :: k0
   real(rk)                     :: kw,fz,ft,zmax,doy,F,A,B,L
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   ! Background walues
   _GET_(self%id_env_k0,k0)
   _SET_DIAGNOSTIC_(self%id_k0,k0)
   
   ! Retrieve current (local) state variable values.
   _GET_(self%id_det1C,det1C) ! detritus1-carbon
   _GET_(self%id_det2C,det2C) ! detritus2-carbon

   ! Self-shading det_with explicit contribution from detritus-C concentration.
   _SET_EXTINCTION_(self%det1_kc*det1C)
   _SET_EXTINCTION_(self%det2_kc*det2C)
   !_SET_EXTINCTION_(0.01)
   
   !background attenuation
   !default values
   fz=1.0_rk 
   ft=1.0_rk 
   if (self%kwFzmaxMeth .eq. 0) then
     !no explicit spm attenuation. Do nothing
     fz=0.0_rk; ft=0.0_rk
   else if (self%kwFzmaxMeth .eq. 1) then
    _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
    !f(z) exponential convergence to the 10% of 'self%a_water' with depth
    fz=self%a_minfr + (1.0_rk-self%a_minfr)*exp(-zmax/11.0)
   else if (self%kwFzmaxMeth .eq. 2) then
    _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
    !f(z)=sigmoidal function of depth with an upper plateau (100%) at 0-10 m and a lower (10%) for 30+
    fz=self%a_minfr+(1.0_rk-self%a_minfr)*(1.0-1.0/(1.0_rk+exp(self%a_fz-zmax*0.5)))
   else if (self%kwFzmaxMeth .eq. 3) then
    
    _GET_GLOBAL_ (self%id_doy,doy) !day of year
    _GET_HORIZONTAL_(self%id_zmax, zmax)  ! max depth
     
    !f(t)=sinusoidal function of day of year with minimum occuring during summer
    !original parameters below were fitted by J.Maerz using the scanfish data from the German Bight
    !A=8.036; B=9.78; L=102.42
    !adjusted parameters: results in a damped, lagged response:
    A=6.0
    B=12.0
    L=85.0
    F=0.05
    ft= F*(A*sin(2.0*doy*Pi/365.0 +2.0*L*Pi/365.0)+B)
    !write (*,'(A, F7.6)') 'ft term: ',0.05*(A*sin(2.0*doy*Pi/365.0 +2.0*L*Pi/365.0)+B)
    
    !f(z)=sigmoidal function of depth with an upper plateau (100%) at 0-10 m and a lower (10%) for 30+
    fz=self%a_minfr+(1.0-self%a_minfr)*(1.0-1.0/(1.0+exp(-zmax*0.5_rk+self%a_fz)))
    
    !write (*,'(A, F7.6)') 'fz term: ',(1.0-self%a_minfr)*(1.0-1.0/(1.0+exp(-zmax*0.5+10))) 
   end if
   
   kw=self%a_water*fz*ft
   
   _SET_EXTINCTION_(kw)
   _SET_DIAGNOSTIC_(self%id_kw,kw)
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------

   end module gpm_abio_pel_EH
