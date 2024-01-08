# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:52:57 2023

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 110
Fecha: 7/3/23
"""

import potencia as pot
import numpy as np
import configparser
import os
import copy as cp

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
# report = pot.start_report(file, name)
dic = pot.read_resumen(file+"/"+name+"-resumen.txt", ["- Dirección óptima -  caso periodoWC", "- Dirección óptima -  caso simulacion2", "- Dirección óptima -  caso simulación1"])
n_cases = len(dic) # periodo, simu1, simu2

IVcurve_file, area_cell, Voc, Isc, Vmp, Imp, T_Vmp, T_Imp = pot.set_params_celda(cell)
calc_cons, T, porcentaje_celdas, consumo_file, case, caras, cells, distrib_type, sup_disp, power_case, report, distribucion= pot.set_params_condiciones(mision, file)

Imp, Pmp = pot.curvaIV(Vmp, IVcurve_file, graph = False)
Imp_B, Vmp_B, Pmp_B = pot.add_deg(deg_B, Imp, Vmp)
Pmp_B = pot.add_degT(Imp_B, Vmp_B, T, T_Imp, T_Vmp)
Imp_E, Vmp_E, Pmp_E = pot.add_deg(deg_E, Imp_B, Vmp_B)
Pmp_E = pot.add_degT(Imp_E, Vmp_E, T, T_Imp, T_Vmp)  

fechas_gen, duracion_gen, dist_sol, cos, eclip_gen, head, datos_gen = pot.read_cos(file, report, case = case)
datos = pot.read_resumen(file+f"/{name}-resumen.txt", ["Periodo", "Duración del eclipse más largo"])
file+= "/Power_Calc/"+power_case+"/" 
fechas_cons, duracion_cons, P_cons, eclip_cons, datos_cons = pot.perfil_consumo(consumo_file, 1, fechas_gen[0], report, file, graph = "Consumo")
 

cos_sum = [sum(cos[:,i]) for i in range(len(cos[0,:]))]
posiciones = pot.orden(head, cos_sum, report)

end_eclip_cons = np.where(eclip_cons != 0)[0][0]
P_cons = pot.desfasaje_eclipse(P_cons, end_eclip_cons)

n = duracion_gen/duracion_cons
if n > 1 :     # Simulacion más larga que el consumo
    fechas_cons, P_cons = pot.to_periodic("el consumo", n, fechas_cons, P_cons, fechas_gen[0], duracion_gen, report, graph = "Potencia")
    graph_P_tot = "Potencia"
else:          # Consumo más largo que la simulación
    graph_P_tot = "Consumo"
    for i in range(len(cos[0,:])):
        fechas_gen, cos_aux = pot.to_periodic("el coseno", n**(-1), fechas_gen, cos[:,i], fechas_cons[0], duracion_cons, report, graph = "Potencia")
        if i == 0:
            cos_new = cos_aux
        else:
            cos_new = np.c_[cos_new, cos_aux]
    cos = cp.copy(cos_new)
    print("---- La duración del perfil de consumo ingresado es mayor a la duración del período que se desea analizar.")
    report.write("\n")
report.write("---- La duración del perfil de consumo ingresado es mayor a la duración del período que se desea analizar.")
E_cons = pot.energy(fechas_cons, P_cons, report, case = "Consumo")


# Caso LABORATORIO 28°C ----------------------------------------
pot.balance_case("Laboratorio", Pmp, file, report)
P_faces_L = pot.pot_faces(fechas_gen, cos, Pmp, head, file, graph = "Potencia de las caras", case = "Laboratorio")
E_faces_L, E_faces_L_tot = pot.energy_faces(fechas_gen, P_faces_L, head, area_cell, report, file, graph = "Energía caras" ,case = "Laboratorio")

if calc_cons:
    cell_faces_L, sup_cells_L = pot.calc_distribution_cell(E_faces_L_tot, distrib_type, head, cells, area_cell, report, sup_disp)
else:
    cell_faces_L, sup_faces_L = pot.distribution_cells(E_cons[-1], E_faces_L_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")

P_tot_L = pot.pot_total(fechas_gen, P_faces_L,  fechas_cons, P_cons, cell_faces_L, area_cell, float(datos["Periodo"]), report, file, graph = graph_P_tot, case = "Laboratorio")
E_tot_L = pot.energy_total(fechas_gen, E_faces_L, cell_faces_L, report, file, graph = "Energía", case = "Laboratorio")

tiempo_bis, balance_L = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_L, report, file, graph = "Energía", case = "Laboratorio")

# Caso EN ORBITA / BOL: Temperatura de la órbita.-------------
pot.balance_case("BOL", Pmp_B, file, report)
P_faces_B = pot.pot_faces(fechas_gen, cos, Pmp_B, head, file, graph = "Potencia de las caras", case = "BOL")
E_faces_B, E_faces_B_tot = pot.energy_faces(fechas_gen, P_faces_B, head,area_cell, report, file, graph = "Energía caras", case = "BOL")

if calc_cons:
    cell_faces_B, sup_cells_B = pot.calc_distribution_cell(E_faces_B_tot, distrib_type, head, cells, area_cell, report, sup_disp, posiciones)
else:
    cell_faces_B, sup_faces_B = pot.distribution_cells(E_cons[-1], E_faces_B_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")

P_tot_B = pot.pot_total(fechas_gen, P_faces_B,  fechas_cons, P_cons, cell_faces_B, area_cell, float(datos["Periodo"]), report, file, graph = graph_P_tot, case = "BOL")
E_tot_B = pot.energy_total(fechas_gen, E_faces_B, cell_faces_B, report, file, graph = "Energía", case = "BOL")

tiempo_bis, balance_B = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_B, report, file, graph = "Energía", case = "BOL")


# Caso EOL: ---------------------------------
pot.balance_case("EOL", Pmp_E, file, report)
P_faces_E = pot.pot_faces(fechas_gen, cos, Pmp_E, head, file, graph = "Potencia de las caras", case = "EOL")
E_faces_E, E_faces_E_tot = pot.energy_faces(fechas_gen, P_faces_E, head,area_cell, report, file, graph = "Energía caras", case = "EOL")

if calc_cons:
    cell_faces_E, sup_cells_E = pot.calc_distribution_cell(E_faces_E_tot, distrib_type, head, cells, area_cell, report, sup_disp, posiciones)
else:
    cell_faces_E, sup_faces_E = pot.distribution_cells(E_cons[-1], E_faces_E_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")


P_tot_E = pot.pot_total(fechas_gen, P_faces_E,  fechas_cons, P_cons, cell_faces_E, area_cell, float(datos["Periodo"]), report, file, graph = graph_P_tot, case = "EOL")
E_tot_E = pot.energy_total(fechas_gen, E_faces_E, cell_faces_E, report, file, graph = "Energía", case = "EOL")

tiempo_bis, balance_E = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_E, report, file, graph = "Energía", case = "EOL")
    
report.close()