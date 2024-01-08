# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:52:57 2023

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 1.14
Fecha: 23/3/23
"""

import potencia_v111 as pot
import numpy as np
import configparser
import os
import copy as cp
import pandas as pd

ini_path = [i for i in os.listdir("../") if i[-4:] == ".ini"]
dat = configparser.ConfigParser()
dat.read("../"+ini_path[0])
mec_orb = dat["mecanica orbital"]
op_av = dat ["opciones avanzadas"]
cell = dat ["celda solar"]
mision = dat ["mision"]
deg_E = dat ["degradaciones EOL"]
deg_B = dat ["degradaciones BOL"]
Potencia = pd.DataFrame()
Energia = pd.DataFrame()


name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, step, duration_run, script_name, paneles = pot.set_params_mec_orbit(dat["mecanica orbital"])

if script_name == "NadirPointing":
    ARB, BAV, ACT, BCV = pot.att_NadirPointing(dat["Attitude NadirPointing"])
    values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, *BAV, ACT, *BCV]
    title = pot.set_title_nadir(*values)
elif script_name == "ThreeAxisKinematic":
    DCM1, DCM2, DCM3, AVX, AVY, AVZ = pot.att_ThreeAxisKinematic(dat["Attitude ThreeAxisKinematic"])
    values = [Epoch, SMA, ECC, INC, RAAN, AOP, TA, *DCM1, *DCM2, *DCM3, AVX, AVY, AVZ]
    title = pot.set_title_tak(*values)
    
file = "../Runs/" + name + "_"+ title
# report = pot.start_report(file, name)
dic_paneles = pot.read_resumen(file+"/"+name+"-resumen.txt", ["Panel "])
len_paneles = len(dic_paneles.keys())
dic_cases = pot.read_resumen(file+"/"+name+"-resumen.txt", ["- Dirección óptima -  caso periodoWC", "- Dirección óptima -  caso simulacion2", "- Dirección óptima -  caso simulación1",])
n_cases = len(dic_cases.keys()) # periodo, simu1, simu2

IVcurve_file, area_cell, Voc, Jsc, Vmp, Jmp, T_Vmp, T_Jmp = pot.set_params_celda(cell)
calc_cons, T, porcentaje_celdas, consumo_file, case, caras, cells, distrib_type, sup_disp, power_case, report, distribucion= pot.set_params_condiciones(mision, file, n_cases, len_paneles)

Jmp, Pmp = pot.curvaIV(Vmp, IVcurve_file, graph = False)
Jmp_B, Vmp_B, Pmp_B = pot.add_deg(deg_B, Jmp, Vmp)
Pmp_B = pot.add_degT(Jmp_B, Vmp_B, T, T_Jmp, T_Vmp)
Jmp_E, Vmp_E, Pmp_E = pot.add_deg(deg_E, Jmp_B, Vmp_B)
Pmp_E = pot.add_degT(Jmp_E, Vmp_E, T, T_Jmp, T_Vmp)  

fechas_gen, duracion_gen, dist_sol, cos, eclip_gen, head, datos_gen = pot.read_cos(file, report, case = case)
datos = pot.read_resumen(file+f"/{name}-resumen.txt", ["Periodo", "Duración del eclipse más largo"])
file+= "/Power_Calc/"+power_case+"/" 

if calc_cons:
    fechas_cons = cp.copy(fechas_gen)
    duracion_cons = pot.to_tiempo(fechas_cons)
    P_cons = np.zeros(len(fechas_cons))
    graph_P_tot = "Potencia"
else:
    fechas_cons, duracion_cons, P_cons, eclip_cons, datos_cons= pot.perfil_consumo(consumo_file, 1, fechas_gen[0], report, file, graph = "Consumo")
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
    E_cons, Consumo = pot.energy(fechas_cons, P_cons, report = report, case = "Consumo")
    Consumo.to_csv(file+'data_Consumo.csv', index=False)

cos_sum = [sum(cos[:,i]) for i in range(len(cos[0,:]))]
posiciones = pot.orden(head, cos_sum, report)




# Caso LABORATORIO 28°C ----------------------------------------
pot.balance_case("Laboratorio", Pmp, file, report)
P_faces_L, Potencia = pot.pot_faces(fechas_gen, cos, Pmp, head, area_cell, file, Potencia, graph = "Potencia de las caras", case = "Laboratorio")
E_faces_L, E_faces_L_tot, Energia = pot.energy_faces(fechas_gen, P_faces_L, head, area_cell, report, file, Energia, graph = "Energía caras" ,case = "Laboratorio")

if calc_cons:
    cell_faces_L, sup_cells_L = pot.calc_distribution_cell(E_faces_L_tot, distrib_type, head, cells, area_cell, report, sup_disp)
else:
    cell_faces_L, sup_faces_L = pot.distribution_cells(E_cons[-1], E_faces_L_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")

P_tot_L, Potencia  = pot.pot_total(fechas_gen, P_faces_L,  fechas_cons, P_cons, cell_faces_L, area_cell, float(datos["Periodo"]), report, file, Potencia, graph = graph_P_tot, case = "Laboratorio")
E_tot_L, Energia = pot.energy_total(fechas_gen, E_faces_L, cell_faces_L, report, file, Energia, graph = "Energía", case = "Laboratorio")

if calc_cons == False:
    tiempo_bis, balance_L, Energia  = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_L, report, file, Energia, graph = "Energía", case = "Laboratorio")

# Caso EN ORBITA / BOL: Temperatura de la órbita.-------------
pot.balance_case("BOL", Pmp_B, file, report)
P_faces_B, Potencia  = pot.pot_faces(fechas_gen, cos, Pmp_B, head, area_cell, file, Potencia, graph = "Potencia de las caras", case = "BOL")
E_faces_B, E_faces_B_tot, Energia  = pot.energy_faces(fechas_gen, P_faces_B, head,area_cell, report, file, Energia, graph = "Energía caras", case = "BOL")

if calc_cons:
    cell_faces_B, sup_cells_B = pot.calc_distribution_cell(E_faces_B_tot, distrib_type, head, cells, area_cell, report, sup_disp, posiciones)
else:
    cell_faces_B, sup_faces_B = pot.distribution_cells(E_cons[-1], E_faces_B_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")

P_tot_B, Potencia  = pot.pot_total(fechas_gen, P_faces_B,  fechas_cons, P_cons, cell_faces_B, area_cell, float(datos["Periodo"]), report, file, Potencia, graph = graph_P_tot, case = "BOL")
E_tot_B, Energia  = pot.energy_total(fechas_gen, E_faces_B, cell_faces_B, report, file, Energia, graph = "Energía", case = "BOL")

if calc_cons == False:
    tiempo_bis, balance_B, Energia  = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_B, report, file, Energia, graph = "Energía", case = "BOL")


# Caso EOL: ---------------------------------
pot.balance_case("EOL", Pmp_E, file, report)
P_faces_E, Potencia  = pot.pot_faces(fechas_gen, cos, Pmp_E, head, area_cell, file, Potencia, graph = "Potencia de las caras", case = "EOL")
E_faces_E, E_faces_E_tot, Energia  = pot.energy_faces(fechas_gen, P_faces_E, head, area_cell, report, file, Energia, graph = "Energía caras", case = "EOL")

if calc_cons:
    cell_faces_E, sup_cells_E = pot.calc_distribution_cell(E_faces_E_tot, distrib_type, head, cells, area_cell, report, sup_disp, posiciones)
else:
    cell_faces_E, sup_faces_E = pot.distribution_cells(E_cons[-1], E_faces_E_tot, head, area_cell, report, distrib_type, caras, distribucion, case = "Laboratorio")


P_tot_E, Potencia  = pot.pot_total(fechas_gen, P_faces_E,  fechas_cons, P_cons, cell_faces_E, area_cell, float(datos["Periodo"]), report, file, Potencia, graph = graph_P_tot, case = "EOL")
E_tot_E, Energia  = pot.energy_total(fechas_gen, E_faces_E, cell_faces_E, report, file, Energia, graph = "Energía", case = "EOL")

if calc_cons == False:
    tiempo_bis, balance_E, Energia  = pot.balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot_E, report, file, Energia, graph = "Energía", case = "EOL")
    
report.close()
Potencia.to_csv(file+'data_Potencia.csv', sep ="\t" )
Energia.to_csv(file+'data_Energia.csv', sep = "\t")