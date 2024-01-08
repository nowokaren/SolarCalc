# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:52:57 2023

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 0.9
Fecha: 24/2/23
"""

import potencia_v0_9 as pot
import numpy as np
import configparser
import os

ini_path = [i for i in os.listdir("../") if i[-4:] == ".ini"]
dat = configparser.ConfigParser()
dat.read("../"+ini_path[0])
mec_orb = dat["mecanica orbital"]
op_av = dat ["opciones avanzadas"]
cell = dat ["celda solar"]
mision = dat ["mision"]
deg_E = dat ["degradaciones EOL"]
deg_B = dat ["degradaciones BOL"]

name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV, ACT, BCV, step, duration_run, script_name, paneles = pot.set_params_mec_orbit(dat["mecanica orbital"])
values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, *BAV, ACT, *BCV]
title = pot.set_title(*values, "NadirPointing")
file = "../Runs/" + name + "_"+ title
report = open(file+f"/{name}-resumen.txt", "r")
texto = report.readlines()
if "------------- PowerCalc -------------\n" in texto:
    report.close()
    print("---- Se están por pisar resultados de PowerCalc obtenidos anteriormente.")
    with open(file+f"/{name}-resumen.txt", "w") as report:
        for line in texto[:texto.index("------------- PowerCalc -------------\n")+1]:
            report.write(line)
    report = open(file+f"/{name}-resumen.txt", "a")
else:
    report.close()
    report = open(file+f"/{name}-resumen.txt", "a")
    report.write("\n")
    report.write("\n")
    report.write("------------- PowerCalc -------------")
    report.write("\n")


IVcurve_file, area_cell, Voc, Isc, Vmp, Imp, T_Vmp, T_Imp = pot.set_params_celda(cell)
T, porcentaje_celdas, consumo_file, case, caras = pot.set_params_condiciones(mision, report)

Imp, Pmp = pot.curvaIV(Vmp, IVcurve_file, graph = False)
Imp_B, Vmp_B, Pmp_B = pot.add_deg(deg_B, Imp, Vmp)
Pmp_B = pot.add_degT(Imp_B, Vmp_B, T, T_Imp, T_Vmp)
Imp_E, Vmp_E, Pmp_E = pot.add_deg(deg_E, Imp_B, Vmp_B)
Pmp_E = pot.add_degT(Imp_E, Vmp_E, T, T_Imp, T_Vmp)  

fechas_gen, duracion_gen, dist_sol, cos, eclip_gen, head, datos_gen = pot.read_cos(file, report, case = case)
datos = pot.read_resumen(file+f"/{name}-resumen.txt", ["Periodo", "Duración del eclipse más largo"])
fechas_cons, duracion_cons, P_cons, eclip_cons, datos_cons = pot.perfil_consumo(consumo_file, 1, fechas_gen[0], report, graph = "Consumo")

cos_sum = [sum(cos[:,i]) for i in range(len(cos[0,:]))]
posiciones = pot.orden(head, cos_sum, report)

end_eclip_cons = np.where(eclip_cons != 0)[0][0]
P_cons = pot.desfasaje_eclipse(P_cons, end_eclip_cons)

n = duracion_gen/duracion_cons
if n > 1 :     # Simulacion más larga que el consumo
    fechas_cons, P_cons = pot.to_periodic(n, fechas_cons, P_cons, fechas_gen[0], duracion_gen, report, graph = "Potencia")
    graph_P_tot = "Potencia"
else:          # Consumo más largo que la simulación
    graph_P_tot = "Consumo"
    print("---- La duración del perfil de consumo ingresado es mayor a la duración del período que se desea analizar.")
    report.write("\n")
report.write("---- La duración del perfil de consumo ingresado es mayor a la duración del período que se desea analizar.")
E_cons = pot.energy(fechas_cons, P_cons, report, case = "Consumo")

# Caso LABORATORIO 28°C ----------------------------------------
print(f"\n----------------------Laboratorio-------(P = {np.round(Pmp,2)})--------------------")
report.write("\n")
report.write(f"\n----------------------Laboratorio-------(P = {np.round(Pmp,2)})--------------------")

P_faces_L = pot.pot_faces(fechas_gen, cos, Pmp, head, graph = "Potencia de las caras", case = "Laboratorio")
E_faces_L, E_faces_L_tot = pot.energy_faces(fechas_gen, P_faces_L, head, area_cell, report, case = "Laboratorio")

gen_tot_L, cell_faces_L, sup_cells_L = pot.distribution_cells(E_cons[-1], E_faces_L_tot, head, posiciones, area_cell, report, caras = caras, distribucion = porcentaje_celdas, case = "Laboratorio")

P_tot_L = pot.pot_total(fechas_gen, P_faces_L,  fechas_cons, P_cons, cell_faces_L, area_cell, report, graph = graph_P_tot, case = "Laboratorio")
E_tot_L = pot.energy_total(fechas_gen, E_faces_L, cell_faces_L, report, graph = "Energía", case = "Laboratorio")

tiempo_bis, balance_L = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_L, report, graph = "Energía", case = "Laboratorio")

# Caso EN ORBITA / BOL: Temperatura de la órbita.-------------
print(f"\n----------------------BOL (T = {T}) (P = {np.round(Pmp_B,2)})---------------------")
report.write("\n")
report.write(f"\n----------------------BOL (T = {T}) (P = {np.round(Pmp_B,2)})---------------------")
P_faces_B = pot.pot_faces(fechas_gen, cos, Pmp_B, head, graph = "Potencia de las caras", case = "BOL")
E_faces_B, E_faces_B_tot = pot.energy_faces(fechas_gen, P_faces_B, head,area_cell, report, case = "BOL")

gen_tot_B, cell_faces_B, sup_cells_B = pot.distribution_cells(E_cons[-1], E_faces_B_tot, head, posiciones, area_cell, report, caras = caras, distribucion = porcentaje_celdas, case = "BOL")

P_tot_B = pot.pot_total(fechas_gen, P_faces_B,  fechas_cons, P_cons, cell_faces_B, area_cell, report, graph = graph_P_tot, case = "BOL")
E_tot_B = pot.energy_total(fechas_gen, E_faces_B, cell_faces_B, report, graph = "Energía", case = "BOL")

tiempo_bis, balance_B = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_B, report, graph = "Energía", case = "BOL")


# Caso EOL: ---------------------------------
print(f"\n----------------------EOL (T = {T}) (P = {np.round(Pmp_E,2)})---------------------")
report.write("\n")
report.write(f"\n----------------------EOL (T = {T}) (P = {np.round(Pmp_E,2)})---------------------")
P_faces_E = pot.pot_faces(fechas_gen, cos, Pmp_E, head, graph = "Potencia de las caras", case = "EOL")
E_faces_E, E_faces_E_tot = pot.energy_faces(fechas_gen, P_faces_E, head,area_cell, report, case = "EOL")

gen_tot_E, cell_faces_E, sup_cells_E = pot.distribution_cells(E_cons[-1], E_faces_E_tot, head, posiciones, area_cell, report, caras = caras, distribucion = porcentaje_celdas, case = "EOL")

P_tot_E = pot.pot_total(fechas_gen, P_faces_E,  fechas_cons, P_cons, cell_faces_E, area_cell, report, graph = graph_P_tot, case = "EOL")
E_tot_E = pot.energy_total(fechas_gen, E_faces_E, cell_faces_E, report, graph = "Energía", case = "EOL")

tiempo_bis, balance_E = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_E, report, graph = "Energía", case = "EOL")
    
report.close()