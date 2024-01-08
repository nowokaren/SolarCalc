# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 12:30:24 2022

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 1.21
Fecha: 21/3/2023
"""
import datetime
import gmat_set_orbit as gso
import orbital_mec as om
import numpy as np
import configparser
import os
import copy as cp

ini_path = [i for i in os.listdir("../") if i[-4:] == ".ini"]
settings = configparser.ConfigParser()
settings.read("../"+ini_path[0])

name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, step, duration_run, stepWC, duration_runWC, script_name, paneles, start_periodWC, end_periodWC, startWC, endWC = gso.set_params_mec_orbit(settings["mecanica orbital"])
force_overwrite, corte_subrun, scale_actitud, scale_sol, scale_panel, lines, model_file = gso.set_params_op_avanz(settings["opciones avanzadas"])

if script_name == "NadirPointing":
    ARB, BAV, ACT, BCV = gso.att_NadirPointing(settings["Attitude NadirPointing"])
    values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, *BAV, ACT, *BCV, model_file]
    title = gso.set_title_nadir(*values)
elif script_name == "ThreeAxisKinematic":
    DCM1, DCM2, DCM3, AVX, AVY, AVZ = gso.att_ThreeAxisKinematic(settings["Attitude ThreeAxisKinematic"])
    values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, *DCM1, *DCM2, *DCM3, AVX, AVY, AVZ, model_file]
    title = gso.set_title_tak(*values)
    
file = "../Runs/" + name + "_"+ title
name_plano_max = []
    
if duration_run > corte_subrun:
    print(f"---- Periodo de tiempo de órbita simulada mayor a {corte_subrun} día(s) → Se realizaran dos simulaciones.")
    # ---- Simulación 1 - Larga ----
    gso.set_script(script_name, file, values, step, duration_run, ini_path, force_overwrite = force_overwrite)
    gso.run_gmat(script_name, file)
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orb, position, sun_vect, euler_angs, t_orb_relativo, period, num_orb, raan, aop, ta, inc, ecc, sma, X_orb, Y_orb, Z_orb = om.data_orbit(file)
    index, cos, dist_sun, ax1, axes = om.cosenos(euler_angs, position, t_orb, X_orb, Y_orb, Z_orb, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, sun_vect=sun_vect, graph = True)
    
    
    # Reporte de settings
    gso.print_report(file, name, values, period, num_orb[:,1][-1], axes, duration_run, step, name_plano_max)
    
    # Reportes de angulos y cosenos
    cos, axes, paneles, n_max1, Z_max1, plano_max_s1, name_plano_max = om.add_plano_max(file, cos, axes, paneles, name_plano_max, case = "simulación1")
    om.report_cos(file, t_orb, dist_sun, cos)
    om.report_ang(file, t_orb, dist_sun, cos)
    cos = om.add_shadow(cos)
    om.graph_cos(t_orb, cos, file, dist_sun, line = lines)
    
    
    # Eclipses
    start, stop, duration = om.data_eclip(file)                           # Toma los datos del file eclipse de GMAT   
    om.eclip_duration(start, duration, t_orb,  file)
    
    # -------- Simulación 2: Contiene al peor eclipse -------
    startWC, endWC, start_run, Epoch, i, duration_run2, step = om.subrunWC_params(start, stop, duration, startWC, endWC, t_orb, stepWC, duration_runWC, period)
    
    if script_name == "NadirPointing":
        values = [Epoch, sma[i], ecc[i], inc[i], raan[i], aop[i], ta[i], ARB, *BAV, ACT, *BCV, model_file]
    elif script_name == "ThreeAxisKinematic":
        values = [Epoch, sma[i], ecc[i], inc[i], raan[i], aop[i], ta[i], *DCM1, *DCM2, *DCM3, AVX, AVY, AVZ, model_file]
    gso.set_script(script_name, file, values, step, duration_run2, ini_path, case = "WC")
    gso.run_gmat(script_name, file, case ="WC")
    
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orbWC, positionWC, sun_vectWC, euler_angsWC, t_orb_relativoWC, period, num_orbWC, raanWC, aopWC, taWC, incWC, eccWC, smaWC, X_orbWC, Y_orbWC, Z_orbWC = om.data_orbit(file, paneles, case = "WC")
    indexWC, cosWC, dist_sunWC, ax1 = om.cosenos(euler_angsWC, positionWC, t_orbWC, X_orbWC, Y_orbWC, Z_orbWC, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, sun_vect=sun_vectWC, graph = True, case = "WC")
    eclipI, eclip = om.getpos_eclip(t_orbWC, positionWC, start, stop, duration, num_orbWC[:,1], file, ax1, axes = axes, graph = True, case = "WC")  # Calcula las posiciones del satétlite en las que está eclipsado
    cosWC = om.add_eclip(cosWC, eclipI, eclip_stamp = np.nan)
    cosWC, axesWC, paneles, n_max2, Z_max2, plano_max_s2, name_plano_max = om.add_plano_max(file, cosWC, axes, paneles, name_plano_max, graph = True, case = "simulacion2")
    
    # Reportes de angulos y cosenos
    om.report_cos(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    om.report_ang(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    
    # Reporte de settings
    gso.print_report(file, name, values, period, num_orbWC[:,1][-1], axesWC, duration_run2, step, name_plano_max, max(duration), case = "WC")
    
else: 
    print(f"---- Periodo de tiempo de órbita simulada menor a {corte_subrun} día(s) → Se realiza una simulación.")
    # ---- Simulación única ----
    gso.set_script(script_name, file, values, step, duration_run, ini_path, force_overwrite = force_overwrite)
    gso.run_gmat(script_name, file)
    
    # Datos de la órbita y calculo de cosenos y ángulos
    t_orbWC, positionWC, sun_vectWC, euler_angsWC, t_orb_relativoWC, period, num_orbWC, raanWC, aopWC, taWC, incWC, eccWC, smaWC, X_orbWC, Y_orbWC, Z_orbWC = om.data_orbit(file, paneles)
    indexWC, cosWC, dist_sunWC, ax1, axes= om.cosenos(euler_angsWC, positionWC, t_orbWC, X_orbWC, Y_orbWC, Z_orbWC, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, graph = True, sun_vect=sun_vectWC)
    start, stop, duration = om.data_eclip(file, raise_exit = False)                           # Toma los datos del file eclipse de GMAT
    if duration != []:
        endWC = stop[duration.index(max(duration))]
        eclipI, eclip = om.getpos_eclip(t_orbWC, positionWC, start, stop, duration, num_orbWC[:,1], file, ax1, axes, graph = True, case = "WC")  # Calcula las posiciones del satétlite en las que está eclipsado
        cosWC = om.add_eclip(cosWC, eclipI, eclip_stamp = np.nan)
        
        cosWC, axesWC, paneles, n_max1, Z_max1, plano_max_s1, name_plano_max = om.add_plano_max(file, cosWC, axesWC, paneles, name_plano_max, case = "simulación1")
        
        # Reportes de angulos y cosenos
        om.report_cos(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
        om.report_ang(file, t_orbWC, dist_sunWC, cosWC, case = "WC")
    else:     
        cosWC, axesWC, paneles, n_max1, Z_max1, plano_max_s1, name_plano_max = om.add_plano_max(file, cosWC, axesWC, paneles, name_plano_max, case = "simulación1")
    om.graph_cos(t_orbWC, om.add_shadow(cosWC), file, dist_sunWC, line = True)
    gso.print_report(file, name, values, period, num_orbWC[:,1][-1], axesWC, duration_run, step, name_plano_max,  max(duration))
    stepWC = cp.copy(step)


if duration != []:
    # Órbita del peor caso (grafico y reporte)
    if abs(start_periodWC) < 0.01:
        startWC = endWC - datetime.timedelta(seconds = period)
    else:
        startWC = start_periodWC
        endWC = end_periodWC

    t_orbWC3, cosWC3, dist_sunWC3 = om.select_timedelta(file, startWC, endWC, case ="WC")
    cosWC3, axes3, paneles3, n_max3, Z_max3, plano_max_s3, name_plano_max3 = om.add_plano_max(file, cosWC3, axesWC, paneles, name_plano_max, case = "periodoWC")
    om.report_cos(file, t_orbWC3, dist_sunWC3, cosWC3, case = "periodoWC")
    om.report_ang(file, t_orbWC3, dist_sunWC3, cosWC3, case = "periodoWC")
    cosWC3 = om.add_shadow(cosWC3)
    om.graph_cos(t_orbWC3, cosWC3, file, dist_sunWC3, case = "periodoWC", line = True, period = period, eclip_dur = max(duration))
    gso.print_report(file, name, values, period, num_orbWC[:,1][-1], axes3, duration_run2, stepWC, name_plano_max, max(duration), case = "periodoWC")
    
# Si se desea tomar un "pedazo" de los datos de cosenos del resultado de la simulacion 1 o 2 (WC), además del caso del peor período, se puede usar la siguiente función en la consola:
# t_orb4, cos4, dist_sun4, axes4, paneles4, n_max4, Z_max4, plano_max_s4, name_plano_max4 = gso.time_interval(file, name, values, num_orbWC, duration_run2, stepWC, startWC, endWC, axes3, paneles3, period = period, name_plano_max = name_plano_max3, eclip_dur = max(duration), case_in = "WC", case_out = "pruebita" )