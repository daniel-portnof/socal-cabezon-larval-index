#V3.30.12.00-safe;_2018_08_01;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_11.6
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#C growth parameters are estimated
#C spawner-recruitment bias adjustment Not tuned For optimality
#_data_and_control_files: simple_dat.ss // simple_ctl.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
2 #_Nblock_Patterns
1 1 #_Blocks_per_pattern
2000  2018    #Accounts for change in commerical regs (size limits)
2004  2018    #Accounts for change in recreational regs (size & bag limits)
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#  autogen
0 0 0 0 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
# 
#
#_Available timevary codes
#_Block types: 0: Pblock=Pbase*exp(TVP); 1: Pblock=Pbase+TVP; 2: Pblock=TVP; 3: Pblock=Pblock(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=Pbase*exp(TVP*env(y));  2: P(y)=Pbase+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=env(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=not implemented
0 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.01 0.5 0.293 -1.238726 0.26 3 2 0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 0 40 0.0001 14 99 0 -3 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 20 90 55.8049 58 99 0 -2 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.05 0.3 0.21318 0.185 99 0 -3 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.01 0.5 0.157348 0.2 99 0 -3 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.01 0.5 0.162416 0.1 99 0 -3 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 9.204e-06 9.204e-06 0.8 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 -3 4 3.187 3.187 3 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 -3 40 34.6 35 10 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -3 3 -0.7 -0.7 0.8 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 1 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 0 0 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
 0.01 0.6 0.274182 -0.9245811 0.275 3 2 0 0 0 0 0 0 0 # NatM_p_1_Mal_GP_1
# Sex: 2  BioPattern: 1  Growth
 0 40 0.0001 13 99 0 -3 0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
 20 90 42.3549 44.1 10 6 -2 0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
 0.05 0.6 0.329223 0.35 0.045 6 -4 0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1
 0.01 0.5 0.428774 0.3 99 0 -3 0 0 0 0 0 0 0 # CV_young_Mal_GP_1
 0.01 0.5 0.0709591 0.2 99 0 -4 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -3 3 1.163e-05 1.163e-05 0.8 0 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Mal_GP_1
 -3 4 3.118 3.19085 3 0 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution  
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_month_1
#  Cohort growth dev base
0.1 10  1 1 1 1 -1  0 0 0 0 0 0 0 # CohortGrowDev
# Movement                            
# Age Error from  parameters                      
# catch multiplier                          
# fraction  female, by  GP                      
1.00E-06  0.999999  0.5 0.5 0 0.5 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR     PR_SD   PR_type    PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             3            31          5.5          12            10    0           1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1          0.7           0.7          0.05    1          -4          0          0          0          0          0          0          0 # SR_BH_steep
             0             2          0.7           0.8           0.8    0          -4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5            0             0             1    0          -4          0          0          0          0          0          0          0 # SR_regime
             0             0            0             0             1    0         -99          0          0          0          0          0          0          0 # SR_autocorr
1 #do_recdev:  0=none; 1=devvector; 2=simple deviations
1970 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 5 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1970 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1977 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2011 #_last_yr_fullbias_adj_in_MPD
 2017 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.45 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1971R 1972R 1973R 1974R 1975R 1976R 1977R 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002F 2003F 2004F 2005F 2006F 2007F 2008F 2009F 2010F 2011F
#  0.1268 -0.0629442 0.0998014 -0.174095 0.0306484 0.714818 -0.0228752 0.00379775 0.261267 0.173626 0.0900049 -0.226622 -0.439888 -0.312088 0.393112 0.551725 0.21891 0.154932 -0.384786 0.596744 -0.682432 -0.273424 -0.829665 0.365024 -0.605267 0.455103 1.11072 -0.546499 -0.656469 0.171606 -0.301581 0 0 0 0 0 0 0 0 0 0
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.2 # F ballpark
-2018 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
5 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
4  # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 # 0    1    0     0.0001   0    99     -1   #1  InitF_FISHERY1_Comm_Non-Live
 # 0    1    0     0.0001   0    99     -1   #2  InitF_FISHERY2_Comm_Live
 # 0    1    0     0.0001   0    99     -1   #3  InitF_FISHERY3_Rec_MM
 # 0    1    0     0.0001   0    99     -1   #4  InitF_FISHERY4_Rec_BB
 # 0    1    0     0.0001   0    99     -1   #5  InitF_FISHERY5_Rec_PBR
 # 0    1    0     0.0001   0    99     -1   #6  InitF_FISHERY6_Rec_CPFV
 # 0    1    0     0.0001   0    99     -1   #7  InitF_FISHERY7_Ghost
#2011 2022
# F rates by fleet
# Yr:  1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# FISHERY 0 0.00211081 0.010609 0.0107037 0.0217063 0.0333334 0.0459509 0.0599453 0.0757167 0.107737 0.146876 0.162531 0.180868 0.202893 0.230365 0.266192 0.314644 0.338215 0.354481 0.356016 0.338877 0.238035 0.242891 0.250688 0.26355 0.283377 0.227156 0.238194 0.247552 0.252337 0.253174 0.0129829 0.0279253 0.038022 0.0447387 0.0493313 0.0527091 0.0554663 0.0579281 0.0602317 0.0624094
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
 5  1  0   1   0   1 #8  InitF_SURVEY1_CPFV 1960-1999
-9999 0 0 0 0 0
# #
# #_Q_parms(if_any);Qunits_are_ln(q)
# #_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
-15 15  -3.00892  0 1 0 -1  0 0 0 0 0 0 0 # LnQ_base_Depl
0 5 0  0.28  99  0 5  0 0 0 0 0 0 0 # Q_extraSD_CA_REC(3)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline 
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 0 0 #1  Comm_Non-Live
 24 0 0 0 #2  Comm_Live
 24 0 0 0 #3  Rec_Boat
 24 0 0 0 #5  Rec_Shore
 5 0 0 4  #8  CPFV 1960-1999
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+nages+1; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 10 0 0 0 #1  Comm_Non-Live
 10 0 0 0 #2  Comm_Live
 10 0 0 0 #3  Rec_Shore
 10 0 0 0 #5  Rec_Boat
 10 0 0 0 #8  CPFV_CPUE

#
#_LO    HI  INIT    PRIOR   PR_SD  PR_Type PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
#FLEET  1 NOT THE LIVE  FISHERY FISH                          
             10            91       34.3909       36.0888            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_Comm.Dead(1)
           -15             3      -6.02952     -0.969098            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_top_logit_Comm.Dead(1)
            -4            12       3.20995        3.5648            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_Comm.Dead(1)
           -10             6             4             4            99             0         -3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_Comm.Dead(1)
          -999            15          -999           -10            99             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_Comm.Dead(1)
            -5            15            10            10            99             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_Comm.Dead(1)
                                     
#FLEET  2 LIVE  FISH  FISHERY                             
  10  91  38.282  38  0.05  0 2  0 0 0 0 0.5 1 2 # PEAK      
  -5  3 -0.572  -0.3  0.05  0 3  0 0 0 0 0.5 1 2 # TOP:_width  of  plateau 
  -4  12  3.367 2.3 0.05  0 3  0 0 0 0 0.5 1 2 # Asc_width     
  -10 6 -1.997  3 0.05  0 3  0 0 0 0 0.5 1 2 # Desc_width      
  -999 5 -999  -5  0.05  0 -2  0 0 0 0 0.5 0 0 # INIT:_selectivity_at_fist_bin     
  -15  15  6.243 -3.5  0.05  0 2  0 0 0 0 0.5 1 2 # FINAL:_selectivity_at_last_bin      
                                      
#FLEET  3 Shore                                
  10  91  10.8  10.8  0.05  0 2 0 0 0 0 0.5 0 0 # PEAK      
  -15  3 -0.5  -0.5  0.05  0 3 0 0 0 0 0.5 0 0 # TOP:_width  of  plateau 
  -4  12  -0.6  -0.6  0.05  0 3 0 0 0 0 0.5 0 0 # Asc_width     
  -10 6 4 -1  0.05  0 -3  0 0 0 0 0.5 0 0 # Desc_width      
  -999 15  -999 -5  0.05  0 -2  0 0 0 0 0.5 0 0 # INIT:_selectivity_at_fist_bin     
  -15  15  10  -0.3  0.05  0 -2  0 0 0 0 0.5 0 0 # FINAL:_selectivity_at_last_bin      
                                                                            
#FLEET  4 Boat                                  
  10  91  45  45  0.05  0 2 0 0 0 0 0.5 2 2 # PEAK      
  -5  3 -2.8  -2.8  0.05  0 3 0 0 0 0 0.5 2 2 # TOP:_width  of  plateau 
  -4  12  4.1 4.1 0.05  0 3 0 0 0 0 0.5 2 2 # Asc_width     
  -10 6 4 -1  0.05  0 -3  0 0 0 0 0.5 0 0 # Desc_width      
  -999 15  -999 -5  0.05  0 -2  0 0 0 0 0.5 0 0 # INIT:_selectivity_at_fist_bin     
  -15  15  10  -0.3  0.05  0 -2  0 0 0 0 0.5 0 0 # FINAL:_selectivity_at_last_bin      
                                    
#Survey 1:00  CPFV  CPUE  1960-1999                             
  1 44  1 1 10  0 -3  0 0 0 0 0.5 0 0 #min  Len bin - fixed
  1 100 44  50  10  0 -4  0 0 0 0 0.5 0 0 #max  Len bin fixed 
                                      

#_no timevary selex parameters
#_timevary selex parameters
#FLEET  2 LIVE  FISH  FISHERY             
  10  91  39.0307 38  0.05  1 3  # PEAK    
  -5  3 -0.601767 -0.3  0.05  1 4  # TOP:_width  of  plateau
  -4  12  2.02299 2.3 0.05  1 4  # Asc_width   
  -10 6 2.56787 3 0.05  1 4  # Desc_width    
  -5  15  -0.816397 -3.5  0.05  1 3  # FINAL:_selectivity_at_last_bin    
                      
#FLEET  6 Boat                  
  10     91     45     45   99    0    5   # PEAK
 -5    3   -1.02507   -1.02507    99    0    5
 -4.0   12.0    4.1    4.1   99    0    5   # Asc_width
# -10 6 4 -1  0.05  1 -3  # Desc_width    
# -5  15  10  -0.3  0.05  1 -2  # FINAL:_selectivity_at_last_bin    
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# no timevary parameters
#
#
# Input variance adjustments factors: 
#_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
#1   8   0.28
#1   9   0.37
#1   10  0.29
#1   11  0.75
#1   12  0.34
#1   13  0.47
#1   14  0.79
4   2   0.067
4   3   1
4   4   0.553
# 4   5   1.3
# 4   6   1.2
# 5   7   1.3
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 3 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 # 1 5 1 0 1 #_CPUE/survey:_6 PISCO SMURF
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 0 0 0 #_CPUE/survey:_1
#  1 1 1 1 #_CPUE/survey:_2
#  1 1 1 1 #_CPUE/survey:_3
#  1 1 1 1 #_lencomp:_1
#  1 1 1 1 #_lencomp:_2
#  0 0 0 0 #_lencomp:_3
#  1 1 1 1 #_agecomp:_1
#  1 1 1 1 #_agecomp:_2
#  0 0 0 0 #_agecomp:_3
#  1 1 1 1 #_size-age:_1
#  1 1 1 1 #_size-age:_2
#  0 0 0 0 #_size-age:_3
#  1 1 1 1 #_init_equ_catch
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
#  0 0 0 0 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting 

999

