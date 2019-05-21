#! /bin/bash

#runs, plots and stores the results of various model versions
#example call: './run_EH-GPM_tests.sh 1

rootdir='2019_05_21_odemet2_n20_2014'
run=1

#if no argument provided, cycle through all model versions
if [ "$#" -lt 1 ]; then
  vset=(1 2 3)
else
  vset=($1)
fi

for v in ${vset[@]}; do 
if [ $v == 1 ]; then
 echo -e "\nfixstochl:"
 if [ $run == 1 ]; then
  ln -sf fabm-EH_AP_AS-GPM_P_Z_fixedSt.yaml fabm.yaml 
  ln -sf output_fixstochl.yaml output.yaml
  ln -sf restart_fixtochl_20140101.nc restart.nc
  gotm_fabm_onur >1D-40m_GPM-EH_fixstochl_err.txt 2>1D-40m_GPM-EH_fixstochl_log.txt
  mkdir -p ${rootdir}/fixstochl
  mv 1D-40m_GPM-EH_fixstochl* ${rootdir}/fixstochl
  #mv restart.nc ${rootdir}/fixstochl/restart_fixtochl_20140101.nc
  cp fabm-EH_AP_AS-GPM_P_Z_fixedSt.yaml ${rootdir}/fixstochl/
  echo "fixstochl: stored the results in: ${rootdir}/fixstochl"
 fi
 python plot_1D-EH-GPM.py ${rootdir}/fixstochl/1D-40m_GPM-EH_fixstochl_dm.nc #1 1 3
 tail -n17 ${rootdir}/fixstochl/1D-40m_GPM-EH_fixstochl_log.txt | head -n9
elif [ $v == 2 ]; then
 echo -e "\nfixsto_empchl:"
 if [ $run == 1 ]; then
  ln -sf fabm-EH_AP_AS-GPM_P_Z_fixedSt_empChl.yaml fabm.yaml 
  ln -sf output_fixsto_empchl.yaml output.yaml
  ln -sf restart_fixsto_empchl_20140101.nc restart.nc
  gotm_fabm_onur >1D-40m_GPM-EH_fixsto_empchl_err.txt 2>1D-40m_GPM-EH_fixsto_empchl_log.txt
  mkdir -p ${rootdir}/fixsto_empchl
  mv 1D-40m_GPM-EH_fixsto_empchl* ${rootdir}/fixsto_empchl
  #mv restart.nc ${rootdir}/fixsto_empchl/restart_fixto_empchl_20140101.nc
  cp fabm-EH_AP_AS-GPM_P_Z_fixedSt_empChl.yaml ${rootdir}/fixsto_empchl/
  echo "fixsto_empchl: stored the results in: ${rootdir}/fixsto_empchl"
 fi
 python plot_1D-EH-GPM.py ${rootdir}/fixsto_empchl/1D-40m_GPM-EH_fixsto_empchl_dm.nc #1 1 3
 tail -n17 ${rootdir}/fixsto_empchl/1D-40m_GPM-EH_fixsto_empchl_log.txt | head -n9
elif [ $v == 3 ]; then
 echo -e "\nvarsto:"
 if [ $run == 1 ]; then
  ln -sf fabm-EH_AP_AS-GPM_P_Z.yaml fabm.yaml 
  ln -sf output_varsto.yaml output.yaml
  ln -sf restart_varsto_20140101.nc restart.nc
  gotm_fabm_onur >1D-40m_GPM-EH_varsto_err.txt 2>1D-40m_GPM-EH_varsto_log.txt
  mkdir -p ${rootdir}/varsto
  mv 1D-40m_GPM-EH_varsto* ${rootdir}/varsto
  #mv restart.nc ${rootdir}/varsto/restart_varsto_20140101.nc
  cp fabm-EH_AP_AS-GPM_P_Z.yaml ${rootdir}/varsto/
  echo "varsto: stored the results in: ${rootdir}/varsto"
 fi
 python plot_1D-EH-GPM.py ${rootdir}/varsto/1D-40m_GPM-EH_varsto_dm.nc # 1 1 1
 tail -n17 ${rootdir}/varsto/1D-40m_GPM-EH_varsto_log.txt | head -n9
fi
done
