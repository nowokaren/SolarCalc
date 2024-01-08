# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:38:22 2023

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 101
Fecha: 27/2/23
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import time
import copy as cp
import ast
import os

def to_tiempo(fechas):
    tiempo = [(x-fechas[0]).total_seconds() for x in fechas]
    return tiempo

def to_fechas(tiempo, start_date):
    fechas = [start_date+datetime.timedelta(seconds = t) for t in tiempo]
    return fechas

def set_params_mec_orbit(section_ini):
    mo = section_ini
    R_Tierra = 6371 # km
    name =  mo.get("Nombre") # "Prueba"
    Epoch = mo.get("Epoch") # '21 Dec 2023 00:00:00.000'
    SMA =  R_Tierra + mo.getfloat( "Altura") #702 #669 #7035
    ECC = mo.getfloat( "ECC") # 0.001 # 0.00001839 #0.0012
    INC =  mo.getfloat( "INC") #  98.2 #98
    RAAN = mo.getfloat( "RAAN") # 64.3 # 27.5638 #90
    AOP = mo.getfloat( "AOP") #   0 #90
    TA = mo.getfloat("TA") # 0
    ARB =  mo.get("AttitudeReferenceBody") # "Earth"    # AttitudeReferenceBody
    BAV =  ast.literal_eval(mo.get( "BodyAlignmentVector")) #[1,0,0]     # BodyAlignmentVectorX, BodyAlignmentVectorY, BodyAlignmentVectorZ
    ACT =  mo.get("AttitudeConstraintType") #"Velocity"  # AttitudeConstraintType
    BCV =  ast.literal_eval(mo.get("BodyConstraintVector")) #[0,0,1]     # BodyConstraintVectorX, BodyConstraintVectorY, BodyConstraintVectorZ
    step = mo.getfloat("mecanica orbital", "step") # 0.0001     # 0.5    # Days Calcular
    duration_run = mo.get("duration_run") # 0.25  #30      # Days Calcular
    script_name = mo.get("script_name") #"orbit"
    paneles = ast.literal_eval(mo.get("paneles")) # [[-1 , 1, 1], [-0.66795061,  0.74419228,  0.00445376]]
    return name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV, ACT, BCV, step, duration_run, script_name, paneles

def set_params_celda(cell):
    IVcurve_file 	=	cell.get("IVcurve_file")
    area_cell 	=	cell.getfloat("area_cell")	# cm2
    Voc		=	cell.getfloat("Voc")			# V
    Isc		=	cell.getfloat("Isc")*area_cell	# A
    Vmp		=	cell.getfloat("Vmp")			# V
    Imp		=	cell.getfloat("Imp")*area_cell		# A
    T_Vmp		=	cell.getfloat("Tc_Vmp")		# V/°C
    T_Imp		=	cell.getfloat("Tc_Imp")*area_cell	# A/°C
    return IVcurve_file, area_cell, Voc, Isc, Vmp, Imp, T_Vmp, T_Imp

def set_params_condiciones(cond, report):
    calc_cons = cond.getboolean("calc_cons")
    if calc_cons:
        print("---- Cálculo del consumo posible según la distribución de celdas por caras definida por la usuaria.")
        report.write("---- Cálculo del consumo posible según la distribución de celdas por caras definida por la usuaria.")
    else:
       print("---- Cálculo del la distribución de celdas por caras según el perfil de consumo definido por la usuaria.")
       report.write("---- Cálculo del la distribución de celdas por caras según el perfil de consumo definido por la usuaria.")
    T		=	cond.getfloat("T")          # °C
    porcentaje_celdas =	ast.literal_eval(cond.get("porcentaje_celdas"))   # %
    consumo_file	=	cond.get("consumo_file")
    case		=	cond.get("case")
    caras       = ast.literal_eval(cond.get("caras"))
    input_cells = cond.get("input_cells")
    cells = ast.literal_eval(cond.get("cells"))
    sup_disp = ast.literal_eval(cond.get("sup_disp"))
    if abs(sum(porcentaje_celdas)-1)> 1e-4:
        print("--!-- La suma de los porcentajes de las celdas es distinto de 1!")
        report.write("\n")
        report.write("--!-- La suma de los porcentajes de las celdas es distinto de 1!")
    if len(caras) > 1:
        if len(porcentaje_celdas) != len(caras):
            print("--!-- Indicar correctamente qué porcentaje de celdas del total de en cada cara")
            report.write("\n")
            report.write("--!-- Indicar correctamente qué porcentaje de celdas del total de en cada cara")
    return calc_cons, T, porcentaje_celdas, consumo_file, case, caras, cells, input_cells, sup_disp

def set_params_degradacion(deg):
    I_mismatch	=	deg.getfloat("I_mismatch")
    I_assembly	=	deg.getfloat("I_assembly")
    I_radiation	=	deg.getfloat("I_radiation")
    V_assembly	=	deg.getfloat("V_assembly")
    V_radiation	=	deg.getfloat("V_radiation")
    I_UV		=	deg.getfloat("I_UV")
    I_debris	=	deg.getfloat("I_debris")
    I_contamin	=	deg.getfloat("I_contamin")
    return I_mismatch, I_assembly, I_radiation, I_radiation, V_assembly, V_radiation, V_radiation, V_radiation, I_UV, I_debris, I_contamin

def curvaIV(Vx, file_name, graph = None):
    curvaIV = pd.read_excel(file_name, sheet_name='Hoja1' )
    
    I = curvaIV["I"]
    V = curvaIV["V"]
    P = curvaIV["P"]
    Ix = np.interp(Vx, V, I)
    Px = np.interp(Vx, V, P)
    if graph != None:
        plt.figure(graph)
        plt.plot(V, I, "-o", ms = 2, label = "Corriente")
        plt.plot(V, P, "-o", ms = 2, label = "Potencia")
        plt.legend()
        plt.grid()
        plt.plot(Vx,Ix, "or")
        plt.plot(Vx,Px, "or")
    return Ix, Px

def set_title(Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV0, BAV1, BAV2, ACT, BCV0, BCV1, BCV2, Attitude = "NadirPointing"):
    for i  in [BAV0, BAV1, BAV2, BCV0, BCV1, BCV2,]:
       if i<1 and i>0:
           i=str(i)[2:]
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA)}-ECC{str(ECC)[2:]}-INC{np.round(INC)}-RAAN{np.round(RAAN)}-AOP{np.round(AOP)}-TA{TA}_{Attitude[:3]}{ARB[0]}{BAV0}{BAV1}{BAV2}-{ACT[0]}{BCV0}{BCV1}{BCV2}"
    return title

def read_cos(file, report, case = 1):
    name_cos = "cosenos"
    if case!=1:
        name_cos += case
    datos = pd.read_csv(file+"/angulos/"+name_cos+".txt", sep = "\t")
    head = datos.columns[2:]
    fechas = [datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in datos["Fecha"].values]
    tiempo = to_tiempo(fechas)
    print(f"El período simulado dura {tiempo[-1]} segundos.")
    report.write("\n")
    report.write(f"El período simulado dura {tiempo[-1]} segundos.")
    dist_sol = datos["Distancia"].values
    cos = datos[head].values
    eclip_gen = [j[1]!=j[1] for j in cos]
    cos[np.isnan(cos)] = 0
    cos[cos<0]=0
    duracion = tiempo[-1]
    return fechas, duracion, dist_sol, cos, eclip_gen, head, datos

def read_resumen(file, datos):
    f = open(file, "r")
    dic = {}
    for i in f:
        d_i = i.split("\t")[0]
        for d in datos:
            if d in d_i:
                dic[d] = np.round(float(i.split("\t")[1]))
    return dic
    
def perfil_consumo(file_name, hoja,  start_date, report, file, graph = None):
    consumo = pd.read_excel(file_name, sheet_name=f'Hoja {hoja}' )
    tiempo = consumo["Time [hrs]"].values*60*60 # segundos
    fechas = to_fechas(tiempo, start_date)
    duracion = tiempo[-1]
    print(f"El perfil de consumo dura {np.round(tiempo[-1],2)} segundos.")
    report.write("\n")
    report.write(f"El perfil de consumo dura {np.round(tiempo[-1],2)} segundos.")
    potencia = consumo["W"].values
    eclip = consumo["Luz"].values
    eclip [eclip == "E"] = 0
    eclip [eclip == "G"] = 1
    if graph != None:
        plt.figure(graph, figsize = (10,5))
        plt.title( "Perfil de consumo")
        plt.xlabel("Tiempo [seg]")
        plt.ylabel("Potencia [W]")
        plt.plot(fechas, potencia, "-o", ms = 2, label = "Consumo")
        plt.legend()
        plt.grid()
        plt.savefig(file+"/consumo.png")
    return fechas, duracion, potencia, eclip, consumo

def pot_faces(fechas, cos, Px, head, file, graph = None, case = ""): 
    P_faces = np.array([Px*cosi for cosi in cos])
    if graph != None:
        for j, i in enumerate(head):
            plt.figure(graph + f"{case}", figsize = (15,5))
            plt.title(f"{case} Potencia de las caras")
            plt.plot(fechas, P_faces[:,j], "-o", ms = 2, label = i[3:])
            plt.xlabel("Fecha")
            plt.ylabel("Potencia [W/celda]  (Area celda = {area_cell} cm2)")
            plt.legend()
            plt.grid()
            plt.savefig(file+f"/{case}/pot_faces.png")
    return P_faces

def pot_total(fechas, P_faces,  fechas_cons, P_cons, cell_faces, area_cell, report, file, graph = None, case = ""):
    tiempo = to_tiempo(fechas)
    P_faces = np.array(P_faces)
    P_total = [np.sum(np.array(P_faces[i,:])*cell_faces) for i in range(len(P_faces[:,0]))]
    P_gen_media = np.trapz(P_total, tiempo)/tiempo[-1]
    print (f"Potencia media generada = {np.round(P_gen_media,2)} W/órbita")
    report.write("\n")
    report.write(f"Potencia media generada = {np.round(P_gen_media,2)} W/órbita")
    if graph != None:    
        plt.figure(graph+f" {case}", figsize = (15,5))
        plt.plot(fechas, P_total, "-o", ms = 2, label = "Suma paneles")
        plt.title(f"{case} Potencia - Cantidad de paneles en cada cara = "+str(cell_faces))
        plt.xlabel("Fecha")
        plt.ylabel("Potencia [W]")
        if case != None:
            plt.plot(fechas_cons, P_cons, "-o", ms = 2, label = "Consumo")
            plt.grid()
        plt.legend()
        plt.savefig(file+f"/{case}/pot_total.png")
    return P_total

def energy(fechas, potencia, report = None, graph = None, case = None,):
    '''case: (str) Caso particular que se está calculando. Se imprimirá en la consola.'''
    E = [0]
    integral=0
    if type(fechas[0]) == datetime.datetime:
        tiempo = to_tiempo(fechas)
    else:
        tiempo = cp.copy(fechas)
    for k in range(len(potencia)-1):    
        integral += np.trapz(potencia[k:k+2],tiempo[k:k+2])
        E.append(integral)
    if graph != None:
        plt.figure(graph)
        plt.plot(tiempo, E, "-o", ms = 2, label = case)
        plt.grid()
        plt.legend()
    if case != None:
        print(f'\n{case} ({np.round(tiempo[-1],2)}) = {np.round(E[-1]/3600000,2)} kWh')
        report.write("\n")
        report.write(f'\n{case} ({np.round(tiempo[-1],2)}) = {np.round(E[-1]/3600000,2)} kWh')
        P_cons_media = np.trapz(potencia, tiempo)/tiempo[-1]
        print (f"Potencia media consumida = {np.round(P_cons_media,2)} W")
        report.write("\n")
        report.write(f"Potencia media consumida = {np.round(P_cons_media,2)} W")
    return np.array(E)
    
def energy_faces(fechas, P_faces, head, area_cell, report, file, graph = None, case = ""):
    E_faces = []
    E_faces_tot = []
    tiempo = to_tiempo(fechas)
    print(f"\nEnergía total generada en la simulación de {np.round(tiempo[-1],2)} seg: ({case})")
    report.write("\n")
    report.write(f"\nEnergía total por celda generada en la simulación de {np.round(tiempo[-1],2)} seg: ({case})")
    for j, i in enumerate(head):
        E = np.array(energy(fechas, P_faces[:,j])) # Jolues
        E_faces.append(E)
        print("{:>10} {:>10} J {:>10} kWh".format(i[3:], np.round(E[-1],2), np.round(E[-1]/3600000,6)))
        report.write("\n")
        report.write("{:>10} {:>10} J {:>10} kWh".format(i[3:], np.round(E[-1],2), np.round(E[-1]/3600000,6)))
        E_faces_tot.append(E[-1])
        if graph != None:
            plt.figure(graph+f" {case}", figsize = (15,5))
            plt.title(f"{case} Energía generada en el periodo por cada cara")
            plt.plot(fechas, E, "-o", ms = 2, label = i[3:])
            plt.legend()
            plt.grid()
            plt.xlabel("Fecha")
            plt.ylabel(f"Energía [J/celda]  (Area celda = {area_cell} cm2)")
            plt.savefig(file+f"/{case}/energy_faces.png")
    print("\n")
    report.write("\n")
    return np.array(E_faces), np.array(E_faces_tot)

def energy_total(fechas, E_faces, cell_faces, report, file, graph = None, case = ""):
    E_faces = np.array(E_faces)
    E_tot = [np.sum(np.array(E_faces[:,i])*cell_faces) for i in range(len(E_faces[0]))]
    print(f"\nEnergía generada ({case}) = {np.round(E_tot[-1]/3600000, 2)} kWh")
    report.write("\n")
    report.write(f"\nEnergía generada ({case}) = {np.round(E_tot[-1]/3600000, 2)} kWh")
    if graph != None:    
        plt.figure(graph+ f" {case}", figsize = (15,5))
        plt.plot(fechas, E_tot, "-o", ms = 2, label = "Energía total")
        plt.title(f"{case} Energía. Celdas en cada cara = "+str(cell_faces))
        plt.xlabel("Fecha")
        plt.ylabel("Energía [J]")
        plt.grid()
        plt.savefig(file+f"/{case}/energy_total.png")
    return np.array(E_tot)

def to_periodic(new, n, fechas0, potencia0, start_date,  tiempo_max, report, graph = None):
    '''tiempo_max: tiempo máximo que se desea que dure la nueva función periodica.'''
    tiempo0 = np.array(to_tiempo(fechas0))
    tiempo_aux = cp.copy(tiempo0)
    t0 = time.process_time()
    n = int(np.ceil(n))
    for i in range(1,n):
        last = tiempo0 + tiempo_aux[-1]
        tiempo_aux = np.concatenate((tiempo_aux, last))
    tiempo1 = [t for t in tiempo_aux if t < tiempo_max]
    potencia1 = np.tile(potencia0, n)[:len(tiempo1)]
    fechas = to_fechas(tiempo1, start_date)
    t1 = time.process_time()
    print(f"El proceso de agrandar periodicamente {new} duró: {t1-t0}")
    report.write("\n")
    report.write(f"El proceso de agrandar periodicamente {new} duró: {t1-t0}")
    if graph != None:
        plt.figure(graph)
        plt.plot(fechas, potencia1, "-o", ms = 2, label = "Consumo periódico")
        plt.grid()
    return fechas, potencia1      

def desfasaje_eclipse(P_cons, index):
    potencia_aux = np.zeros(len(P_cons))
    potencia_aux[:len(P_cons[index:])] = P_cons[index:]
    potencia_aux[len(P_cons[index:]):] = P_cons[:index]
    P_cons = np.copy(potencia_aux)
    return P_cons

def orden(head, E_faces_tot, report): 
    E_faces_tot_ordenado = sorted(E_faces_tot, reverse=True) 
    posiciones = [] 
    for i in E_faces_tot: 
        posiciones.append(E_faces_tot_ordenado.index(i)) 
    print("\nOrden de generación de las caras: (Mayor a menor)\n")
    report.write("\n")
    report.write("\nOrden de generación de las caras: (Mayor a menor)\n")
    line = ""
    for i in head[posiciones][:-1]:
        line += i[-7:] + " > "
    line += head[posiciones][-1][-7:]+"\n"
    print(line)
    report.write("\n")
    report.write(line)
    return posiciones

def distribution_cells(consumo_total, gen_faces, head, posiciones, area_cell, report, caras = 1, distribucion = 1, case = ""):
    '''Si len(caras) != 1: hay que elegir valor para distribución distinto de 1'''
    gen_tot = []
    cell_faces = np.zeros(len(posiciones))
    sup_cells = np.zeros(len(posiciones))
    for i in caras:
        if abs(gen_faces[i])< 1e-4:
            print(f"--!-- La cara {i} elegidas para poner paneles tiene radiación nula.")
            report.write("\n")
            report.write(f"--!-- La cara {i} elegidas para poner paneles tiene radiación nula.")
            raise SystemExit
    report.write("\n")
    if len(caras) > 1:
        print('{:^10} {:^5} {:^15} {:^15} {:^15}'.format('Panel', "(%)", "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        report.write("\n")
        report.write('{:^10} {:^5} {:^15} {:^15} {:^15}'.format('Panel', "(%)", "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        for j, i in enumerate(distribucion):
            gen_panel = consumo_total*i
            if caras == "Max":
                i_cara = posiciones.index(j)
            else:
                i_cara = caras[j]
            cell_per_panel = gen_panel/gen_faces[i_cara]
            cell_faces[i_cara] = np.ceil(cell_per_panel)
            sup_cells [i_cara] = cell_faces[i_cara] * area_cell *0.0001 # m^2
            print('{:^10} {:^5} {:^15} {:^15} {:^15}'.format(f'{head[i_cara][-7:]}', i*100, np.round(gen_panel/3600000,2), int(np.ceil(cell_per_panel)), np.round(sup_cells[i_cara],2)))
            report.write("\n")
            report.write('{:^10} {:^5} {:^15} {:^15} {:^15}'.format(f'{head[i_cara][-7:]}', i*100, np.round(gen_panel/3600000,2), int(np.ceil(cell_per_panel)), np.round(sup_cells[i_cara],2)))
            gen_tot.append(gen_panel)
    else:
        i_cara = caras[0]
        cell_per_panel = consumo_total/gen_faces[i_cara]
        cell_faces[i_cara] = np.ceil(cell_per_panel)
        sup_cells [i_cara] = cell_faces[i_cara] * area_cell *0.0001 # m^2
        print('{:^16}: {:^15}'.format('Panel', head[i_cara][-7:]))
        print('{:^16} : {:^15}'.format("Generación (kWh)", np.round(consumo_total/3600000,2)))
        print('{:^16} : {:^15}'.format('Cantidad celdas', int(np.ceil(cell_per_panel))))
        print('{:^16} : {:^15} m^2'.format('Superficie total', np.round(sup_cells[i_cara],2)))
        report.write("\n")
        report.write('{:^16}: {:^15}'.format('Panel', head[i_cara][-7:]))
        report.write("\n")
        report.write('{:^16} : {:^15}'.format("Generación (kWh)", np.round(consumo_total/3600000,2)))
        report.write("\n")
        report.write('{:^16} : {:^15}'.format('Cantidad celdas', int(np.ceil(cell_per_panel))))
        report.write("\n")
        report.write('{:^16} : {:^15} m^2'.format('Superficie total', np.round(sup_cells[i_cara],2)))
    return gen_tot, cell_faces, sup_cells

def balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot, report, file, graph = None, case = ""):
    tiempo_cons = to_tiempo(fechas_cons)
    tiempo_gen = to_tiempo(fechas_gen)
    if n>1:
        E_tot_bis = np.interp(tiempo_cons, tiempo_gen, E_tot)
        balance = E_tot_bis- E_cons
        tiempo_bis = tiempo_cons
    else:
        E_tot_bis = []
        E_cons_bis = np.interp(tiempo_gen, tiempo_cons, E_cons)
        E_cons_bis = np.array(E_cons_bis)
        balance = E_tot - E_cons_bis
        tiempo_bis = tiempo_gen
        E_cons = E_cons_bis
    print(f"Balance energético = {np.round(balance[-1]/1000,2)} kJ\t{np.round(balance[-1]/3600000,8)} kWh")
    report.write("\n")
    report.write(f"Balance energético = {np.round(balance[-1]/1000,2)} kJ\t{np.round(balance[-1]/3600000,8)} kWh")
    neg_balance = [i for i,j in zip(balance, tiempo_bis) if i<0]
    tiempo_neg_balance = [j for i,j in zip(balance, tiempo_bis) if i<0]
    batery_cap = np.sum(energy(tiempo_neg_balance, neg_balance))
    print(f"Energía total extra que necesito: {-np.round(batery_cap/1000,2)} kJ\t{-np.round(batery_cap/3600000,8)} kWh (Batería)")
    report.write("\n")
    report.write(f"Energía mínima necesaria en la batería: {-np.round(batery_cap/1000,2)} kJ\t{-np.round(batery_cap/3600000,8)} kWh\n            (Energía acumulada para valores negativos del balance)")
    fechas_bis = to_fechas(tiempo_bis, fechas_gen[0])
    if graph != None:
        plt.figure(graph+ f" {case}", figsize = (15,5))
        plt.plot(fechas_bis, balance, "-o", ms = 2, label = "Balance energético")
        plt.grid()
        plt.title(f"{case} Balance energético")
        plt.xlabel("Tiempo [seg]")
        plt.ylabel("Energía [J]")
        if case != None:
            plt.plot(fechas_bis, E_cons, "-o", ms = 2, label = "Consumo integrado")
        plt.legend()
        plt.savefig(file+f"/{case}/balance.png")
    return tiempo_bis, balance

def balance_case(name, P, file, report):
    print(f"\n----------------------{name}-------(P = {np.round(P,2)})--------------------")
    report.write("\n")
    report.write(f"\n----------------------{name}-------(P = {np.round(P,2)})--------------------")
    try:
        os.mkdir(file+f"/{name}")
    except FileExistsError:
        pass

def add_deg(deg, Imp, Vmp):
    for frac in deg:
        if frac[0].lower() == "i":
            Imp*= float(deg[frac])
        else:
            Vmp*= float(deg[frac])
        Pmp = Imp*Vmp
    return Imp, Vmp, Pmp

def add_degT(Imp, Vmp, T, T_Imp, T_Vmp):
    return (Imp+T_Imp*T) * (Vmp+T_Vmp * T)


def start_report(file, name):
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
    return report

def calc_distribution_cell(input_cells, head, cells, area_cell, sup_disp = 1, posiciones =1):
    if "Max" in input_cells: 
        if input_cells[0] =="S":
            cells = np.floor(cells/area_cell)
        cell_faces = [0] * len(head)
        j = 0   # iterar caras
        i = 0   # celdas ubicadas
        while i<cells:
            try:
                i_cara = posiciones.index(j)
            except IndexError:
                print("--!-- La cantidad de superficie disponible no alcanza para ubicar todas las celdas.")
            cell_max_cara = np.floor(sup_disp[i_cara]/area_cell)
            if cells-i >= cell_max_cara:
                cell_faces [i_cara] = cell_max_cara
            else: cell_faces [i_cara] = cells-i
            i-=cell_faces[i_cara]
            j+=1
    elif input_cells == "NxF":
        cell_faces = cells
        sup_cells = area_cell*cells
    elif input_cells == "SxF":
            sup_cells = cells
            cell_faces = [np.floor(i/area_cell) for i in cells]
    return cell_faces, sup_cells