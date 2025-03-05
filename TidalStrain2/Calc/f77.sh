#!/bin/sh

#f77 -o conv_JST2UTC conv_JST2UTC.f

#f77 -o make_inputdata_X make_inputdata_X.f day_1901.f trans_magnitude.f csepform.f day2year.f year2day.f

#-f77 -o earthtide_mod earthtide_mod.f

f77 -o calcu_idep_ninf calcu_idep_ninf.f

#make -f make.loadgreenf3_mod
#-f77 -o loadgreenf3_mod loadgreenf3_mod.f sphlgd.f

f77 -o sum_both_list sum_both_list.f

f77 -o select_Lame select_Lame.f

f77 -o STRN2dVdTdS STRN2dVdTdS.f strn2stress.f

f77 -o dTdS2dCFF dTdS2dCFF.f

#--- 2020/06/29
f77 -o estimate_tidal_phase_level estimate_tidal_phase_level.f
f77 -o STRN2STRESS6 STRN2STRESS6.f
f77 -o calc_principal_stress calc_principal_stress.f
