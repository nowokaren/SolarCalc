# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 12:30:24 2022

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
"""
import datetime
import gmat_set_orbit as gso
import orbital_mec as om
import os
import numpy as np

# --------------- Variables de entrada --------------- 

R_Tierra = 6371 # km
name =   "Prueba"
Epoch =  '21 Dec 2023 00:00:00.000'
SMA =  R_Tierra + 702 #669 #7035
ECC =  0.001 # 0.00001839 #0.0012
INC =    98.2 #98
RAAN =  64.3 # 27.5638 #90
AOP =    0 #90
TA =   0
ARB =  "Earth"    # AttitudeReferenceBody
BAV = [1,0,0]     # BodyAlignmentVectorX, BodyAlignmentVectorY, BodyAlignmentVectorZ
ACT = "Velocity"  # AttitudeConstraintType
BCV = [0,0,1]     # BodyConstraintVectorX, BodyConstraintVectorY, BodyConstraintVectorZ
step = 0.0001     # 0.5    # Days Calcular
duration_run = 0.25  #30      # Days Calcular
script_name = "orbit"
paneles = [[-1 , 1, 1], [-0.66795061,  0.74419228,  0.00445376]]

# -------------------------------------------------------- 

# --------------- Opciones Avanzadas ---------------------
force_overwrite = True
Attitude = "NadirPointing"
corte_subrun = 1 # días
scale = 1000
model_file = "\'"+os.getcwd()+"\\CubeSat3U-salientes.3ds"+ "\'"
values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, *BAV, ACT, *BCV, model_file]
title = gso.set_title(*values, Attitude)
file = "Runs/" + name + "_"+ title
scale_actitud = 700
scale_sol = 2000
scale_panel = 2000
lines = True  #Agregar lineas al gráfico de la simulación larga.
# ---------------------------------------------------------


if duration_run > corte_subrun:
    print(f"---- Periodo de tiempo de órbita simulada mayor a {corte_subrun} día(s) → Se realizaran dos simulaciones.")
    # ---- Simulación 1 - Larga ----
    gso.set_script(script_name, file, values, step, duration_run, force_overwrite = force_overwrite)
    gso.run_gmat(script_name, file)
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orb, position, sun_vect, euler_angs, t_orb_relativo, period, num_orb, raan, aop, ta, inc, ecc, sma, X_orb, Y_orb, Z_orb = om.data_orbit(file, paneles)
    index, cos, dist_sun, ax1, axes = om.cosenos(euler_angs, position, t_orb, X_orb, Y_orb, Z_orb, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, sun_vect=sun_vect, graph = True)
    
    
    # Reporte de settings
    gso.print_report(file, name, values, period, num_orb[:,1][-1], paneles, duration_run, step)
    
    # Reportes de angulos y cosenos
    cos, paneles, n_max1, Z_max1, plano_max_s1 = om.add_plano_max(file, cos, paneles, case = "simulación1")
    om.report_cos(file, t_orb, dist_sun, cos)
    om.report_ang(file, t_orb, dist_sun, cos)
    cos = om.add_shadow(cos)
    om.graph_cos(t_orb, cos, file, dist_sun, line = lines)
    
    
    # Eclipses
    start, stop, duration = om.data_eclip(file)                           # Toma los datos del file eclipse de GMAT   
    om.eclip_duration(start, duration, t_orb,  file)
    
    # -------- Simulación 2: Contiene al peor eclipse -------
    startWC, endWC, start_run, Epoch, i, duration_run2, step = om.subrunWC_params(start, stop, duration, t_orb, period)
    values = [Epoch, sma[i], ecc[i], inc[i], raan[i], aop[i], ta[i], ARB, *BAV, ACT, *BCV]
    gso.set_script(script_name, file, values, step, duration_run2, case = "WC")
    gso.run_gmat(script_name, file, case ="WC")
    
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orbWC, positionWC, sun_vectWC, euler_angsWC, t_orb_relativoWC, period, num_orbWC, raanWC, aopWC, taWC, incWC, eccWC, smaWC, X_orbWC, Y_orbWC, Z_orbWC = om.data_orbit(file, paneles, case = "WC")
    indexWC, cosWC, dist_sunWC, ax1, axes = om.cosenos(euler_angsWC, positionWC, t_orbWC, X_orbWC, Y_orbWC, Z_orbWC, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, sun_vect=sun_vectWC, graph = True, case = "WC")
    eclipI, eclip = om.getpos_eclip(t_orbWC, positionWC, start, stop, duration, num_orbWC[:,1], file, ax1, axes = axes, graph = True, case = "WC")  # Calcula las posiciones del satétlite en las que está eclipsado
    cosWC = om.add_eclip(cosWC, eclipI, eclip_stamp = np.nan)
    cos, paneles, n_max2, Z_max2, plano_max_s2 = om.add_plano_max(file, cosWC, paneles, graph = True, case = "simulacion2")
    
    # Reportes de angulos y cosenos
    om.report_cos(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    om.report_ang(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    
    # Reporte de settings
    gso.print_report(file, name, values, period, num_orbWC[:,1][-1], paneles, duration_run2, step, case = "WC")
    
else: 
    print(f"---- Periodo de tiempo de órbita simulada menor a {corte_subrun} día(s) → Se realiza una simulación.")
    # ---- Simulación única ----
    gso.set_script(script_name, file, values, step, duration_run, force_overwrite = force_overwrite)
    gso.run_gmat(script_name, file)
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orbWC, positionWC, sun_vectWC, euler_angsWC, t_orb_relativoWC, period, num_orbWC, raanWC, aopWC, taWC, incWC, eccWC, smaWC, X_orbWC, Y_orbWC, Z_orbWC = om.data_orbit(file, paneles)
    indexWC, cosWC, dist_sunWC, ax1, axes = om.cosenos(euler_angsWC, positionWC, t_orbWC, X_orbWC, Y_orbWC, Z_orbWC, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, graph = True, sun_vect=sun_vectWC)
    start, stop, duration = om.data_eclip(file, raise_exit = False)                           # Toma los datos del file eclipse de GMAT
    if duration != []:
        endWC = stop[duration.index(max(duration))]
        eclipI, eclip = om.getpos_eclip(t_orbWC, positionWC, start, stop, duration, num_orbWC[:,1], file, ax1, axes, graph = True, case = "WC")  # Calcula las posiciones del satétlite en las que está eclipsado
        cosWC = om.add_eclip(cosWC, eclipI, eclip_stamp = np.nan)
        
        cosWC, paneles, n_max1, Z_max1, plano_max_s1 = om.add_plano_max(file, cosWC, paneles, case = "simulación1")
        
        # Reportes de angulos y cosenos
        om.report_cos(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
        om.report_ang(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    else:     
        cosWC, paneles, n_max1, Z_max1, plano_max_s1 = om.add_plano_max(file, cosWC, paneles, case = "simulación1")
    om.graph_cos(t_orbWC, om.add_shadow(cosWC), file, dist_sunWC, line = True)
    gso.print_report(file, name, values, period, num_orbWC[:,1][-1], paneles, duration_run, step,)
    


if duration != []:
    # Órbita del peor caso (grafico y reporte)
    startWC = endWC - datetime.timedelta(seconds = period)
    t_orbWC, cosWC, dist_sunWC = om.select_timedelta(file, startWC, endWC, case ="WC")
    cosWC, paneles, n_max3, Z_max3, plano_max_s3 = om.add_plano_max(file, cosWC, paneles, case = "periodoWC")
    om.report_cos(file, t_orbWC, dist_sunWC, cosWC, case = "periodoWC")
    om.report_ang(file, t_orbWC, dist_sunWC, cosWC, case = "periodoWC")
    cosWC = om.add_shadow(cosWC)
    om.graph_cos(t_orbWC, cosWC, file, dist_sunWC, case = "periodoWC", line = True, period = period, eclip_dur = max(duration))
    

