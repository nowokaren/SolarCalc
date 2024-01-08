# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:38:22 2023

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 1.14
Fecha: 23/3/23
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

def set_params_mec_orbit(mo):
    R_Tierra = 6371 # km
    name =  mo.get("Nombre") # "Prueba"
    Epoch = mo.get("Epoch") # '21 Dec 2023 00:00:00.000'
    SMA =  R_Tierra + mo.getfloat( "Altura") #702 #669 #7035
    ECC = mo.getfloat( "ECC") # 0.001 # 0.00001839 #0.0012
    INC =  mo.getfloat( "INC") #  98.2 #98
    RAAN = mo.getfloat( "RAAN") # 64.3 # 27.5638 #90
    AOP = mo.getfloat( "AOP") #   0 #90
    TA = mo.getfloat("TA") # 0
    step = mo.getfloat("mecanica orbital", "step") # 0.0001     # 0.5    # Days Calcular
    duration_run = mo.get("duration_run") # 0.25  #30      # Days Calcular
    script_name = mo.get("script_name") #"orbit"
    paneles = ast.literal_eval(mo.get("paneles")) # [[-1 , 1, 1], [-0.66795061,  0.74419228,  0.00445376]]
    return name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, step, duration_run, script_name, paneles

def set_params_celda(cell):
    IVcurve_file 	=	cell.get("IVcurve_file")
    area_cell 	=	cell.getfloat("area_cell")	# cm2
    Voc		=	cell.getfloat("Voc")			# V
    Jsc		=	cell.getfloat("Jsc")*area_cell	# A
    Vmp		=	cell.getfloat("Vmp")			# V
    Jmp		=	cell.getfloat("Jmp")*area_cell		# A
    T_Vmp		=	cell.getfloat("Tc_Vmp")		# V/°C
    T_Imp		=	cell.getfloat("Tc_Imp")*area_cell	# A/°C
    return IVcurve_file, area_cell, Voc, Jsc, Vmp, Jmp, T_Vmp, T_Imp

def set_params_condiciones(cond, file, n_cases, len_panels):
    # n: cantidad de paneles agregados manualmente + cantidad de paneles que maximizan cada caso
    calc_cons = cond.getboolean("calc_cons")
    T		=	cond.getfloat("T")          # °C
    porcentaje_celdas =	ast.literal_eval(cond.get("porcentaje_celdas"))   # %
    consumo_file	=	cond.get("consumo_file")
    case		=	cond.get("case")
    caras       = ast.literal_eval(cond.get("caras"))
    distrib_type = cond.get("distrib_type")
    cells = ast.literal_eval(cond.get("cells"))
    sup_disp = ast.literal_eval(cond.get("sup_disp"))
    if calc_cons:
        power_case = "C-"
        if distrib_type == "N-fix":
            cel = "".join(map(str,cells))
            power_case += distrib_type+f"{cel}celdas"
            distribucion = 1
        elif distrib_type == "S-fix":
            cel = "".join(map(str,cells))
            power_case += distrib_type+f"{cel}m2"
            distribucion = 1
        elif "max" in distrib_type:
            sup = "".join(map(str,sup_disp))
            power_case += distrib_type+f"SupDisp{sup}-Cells{cells}"
            distribucion = 1
    else:
        power_case = "G-"
        if distrib_type == "%":
            power_case += "porcentaje"
            distribucion = porcentaje_celdas
        elif "sup_disp" in distrib_type:
            sup = "".join(map(str,sup_disp))
            power_case += distrib_type + sup
            distribucion = sup_disp
            if case == "-WC":
                n_cases-=1
            elif case == "":
                n_cases-=2
            n = n_cases + len_panels
            while len(distribucion)!= n+6:
                distribucion = list(map(float,ast.literal_eval(input(f"La longitud de la variable sup_disp no coincide con la cantidad de planos de la simulación. Definir nuevamente sup_disp con longitud {n+6}:\n"))))
            if "fix" in distrib_type:
                power_case += "caras" + "".join(map(str,caras))
    if "Power_Calc" not in os.listdir(file):
        os.mkdir(file+"/Power_Calc")
        os.mkdir(file+"/Power_Calc/"+power_case)
        report = open(file+"/Power_Calc/"+power_case+"/"+power_case+"-resumen.txt", "w")
    elif power_case not in os.listdir(file+"/Power_Calc"):
        os.mkdir(file+"/Power_Calc/"+power_case)
        report = open(file+"/Power_Calc/"+power_case+"/"+power_case+"-resumen.txt", "w")
    else:
        report = open(file+"/Power_Calc/"+power_case+"/"+power_case+"-resumen.txt", "w")
    report.write(f"---------- Power_Case = {power_case} ----------")
    report.write("\n")
    if calc_cons:
        print("---- Cálculo del consumo posible según la distribución de celdas por caras definida por la usuaria.")
        report.write("---- Cálculo del consumo posible según la distribución de celdas por caras definida por la usuaria.")
        report.write("\n")
        print("---- Las distribución de celdas es:")
        report.write("---- Las distribución de celdas es:")
        report.write("\n")
        if "fix" in distrib_type:
            if "N" in distrib_type:
                print(f"- CANTIDAD DE CELDAS por plano {cells}. ")
                report.write(f"- CANTIDAD DE CELDAS por plano: {cells}. ")
            elif "S" in distrib_type:
                print(f"- SUPERFICIE DE CELDAS por plano {cells}. ")
                report.write(f"- SUPERFICIE DE CELDAS por plano: {cells}. ")
        elif "max" in distrib_type:
            if distrib_type[0]=="S":
                print(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- SUPERFICIE DE CELDAS: {cells} m^2")
                report.write(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- SUPERFICIE DE CELDAS: {cells} m^2")
            elif distrib_type[0]=="N": 
                print(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- CANTIDAD DE CELDAS: {cells} celdas")
                report.write(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- CANTIDAD DE CELDAS: {cells} celdas")
    else:
       print("---- Cálculo del la distribución de celdas por caras según el perfil de consumo definido por la usuaria.")
       report.write("---- Cálculo del la distribución de celdas por caras según el perfil de consumo definido por la usuaria.")
       if "fix" in distrib_type:
            print(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- CARAS: {caras}")
            report.write(f"- SUPERFICIE DISPONIBLE: {sup_disp}\n- CARAS: {caras}")
    
    if abs(sum(porcentaje_celdas)-1)> 1e-4:
        print("--!-- La suma de los porcentajes de las celdas es distinto de 1!")
        report.write("\n")
        report.write("--!-- La suma de los porcentajes de las celdas es distinto de 1!")
    if len(caras) > 1:
        if len(porcentaje_celdas) != len(caras):
            print("--!-- Indicar correctamente qué porcentaje de celdas del total de en cada cara")
            report.write("\n")
            report.write("--!-- Indicar correctamente qué porcentaje de celdas del total de en cada cara")
    return calc_cons, T, porcentaje_celdas, consumo_file, case, caras, cells, distrib_type, sup_disp, power_case, report, distribucion

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

def set_title_nadir(Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV0, BAV1, BAV2, ACT, BCV0, BCV1, BCV2):
    for i  in [BAV0, BAV1, BAV2, BCV0, BCV1, BCV2,]:
       if i<1 and i>0:
           i=str(i)[2:]
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA,2)}-ECC{str(ECC)[2:]}-INC{np.round(INC,1)}-RAAN{np.round(RAAN,2)}-AOP{np.round(AOP,2)}-TA{round(TA,2)}_Nad{ARB[0]}{BAV0}{BAV1}{BAV2}-{ACT[0]}{BCV0}{BCV1}{BCV2}"
    return title

def set_title_tak(Epoch, SMA, ECC, INC, RAAN, AOP, TA, DCM11, DCM12, DCM13, DCM21, DCM22, DCM23, DCM31, DCM32, DCM33, AVX, AVY, AVZ):
    DCM = cp.copy(DCM11)
    for i in [DCM11, DCM12, DCM13, DCM21, DCM22, DCM23, DCM31, DCM32, DCM33]:      
        if i - round(i) != 0.0:
            i = round(i - round(i),2)
        DCM += i
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA,2)}-ECC{str(ECC)[2:]}-INC{np.round(INC,1)}-RAAN{np.round(RAAN,2)}-AOP{np.round(AOP,2)}-TA{round(TA,2)}_3AxKin{DCM}-AngVel{AVX}{AVY}{AVZ}"
    return title


def att_NadirPointing(nadir):
    ARB =  nadir.get("AttitudeReferenceBody") # "Earth"    # AttitudeReferenceBody
    BAV =  ast.literal_eval(nadir.get( "BodyAlignmentVector")) #[1,0,0]     # BodyAlignmentVectorX, BodyAlignmentVectorY, BodyAlignmentVectorZ
    ACT =  nadir.get("AttitudeConstraintType") #"Velocity"  # AttitudeConstraintType
    BCV =  ast.literal_eval(nadir.get("BodyConstraintVector")) #[0,0,1]     # BodyConstraintVectorX, BodyConstraintVectorY, BodyConstraintVectorZ     
    return ARB, BAV, ACT, BCV

def att_ThreeAxisKinematic(tak):
    DCM = np.zeros((3,3))
    DCM1    =   ast.literal_eval(tak.get("DCM1"))/np.linalg.norm(ast.literal_eval(tak.get("DCM1")))         # Matriz de cosenos directores (Ortonormal)
    DCM2    =   ast.literal_eval(tak.get("DCM2"))/np.linalg.norm(ast.literal_eval(tak.get("DCM2"))) 
    DCM3    =   ast.literal_eval(tak.get("DCM3"))/np.linalg.norm(ast.literal_eval(tak.get("DCM3"))) 
    DCM [0,:] = DCM1
    DCM [1,:] = DCM2
    DCM [2,:] = DCM3
    if not np.allclose(np.transpose(DCM), np.linalg.inv(DCM)):
        print("La matriz de cosenos directores indicada no es ortogonal. Ingresar una matriz válida:")
        DCM1    =   np.linalg.norm(ast.literal_eval(input("DCM1 = (formato = [0,0,0])")))         # Matriz de cosenos directores (Ortonormal)
        DCM2    =   np.linalg.norm(ast.literal_eval(input("DCM2 = (formato = [0,0,0])"))) 
        DCM3    =   np.linalg.norm(ast.literal_eval(input("DCM3 = (formato = [0,0,0])")))
        att_ThreeAxisKinematic(tak)
    AVX = tak.getfloat("AngularVelocityX")        # Velocidades angulares (deg/sec)
    AVY = tak.getfloat("AngularVelocityY") 
    AVZ = tak.getfloat("AngularVelocityZ")
    return DCM1, DCM2, DCM3, AVX, AVY, AVZ


def curvaIV(Vx, file_name, graph = False):
    curvaIV = pd.read_excel(file_name, sheet_name='Hoja1' )
    
    I = curvaIV["I"]
    V = curvaIV["V"]
    P = curvaIV["P"]
    Ix = np.interp(Vx, V, I)
    Px = np.interp(Vx, V, P)
    if graph:
        plt.figure(graph)
        plt.plot(V, I, "-o", ms = 2, label = "Corriente")
        plt.plot(V, P, "-o", ms = 2, label = "Potencia")
        plt.legend()
        plt.grid()
        plt.plot(Vx,Ix, "or")
        plt.plot(Vx,Px, "or")
    return Ix, Px


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
    j=0
    for i in f:
        for d in datos:
            if d in i:
                if d in dic and d == "Panel ":
                    j+=1
                    d +=str(j)
                if d == ("Periodo" or "Duración del eclipse más largo") :
                    dic[d] = i.split("\t")[1]
                else:
                    dic[d] = i.split("\t")[-1]
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
        plt.savefig(file+"consumo.png")
    return fechas, duracion, potencia, eclip, consumo

def pot_faces(fechas, cos, Px, head, area_cell, file, Potencia, graph = None, case = ""): 
    P_faces = np.array([Px*cosi for cosi in cos])
    Potencia[f"{case[:3]}-fechas"] = fechas
    Potencia[f"{case[:3]}-tiempo"] = to_tiempo(fechas)
    for i, cara in enumerate(head):
        Potencia[f"{case[:3]}-{cara[3:]}-W"] = P_faces[:,i]
    if graph != None:
        for j, i in enumerate(head):
            plt.figure(graph + f"{case}", figsize = (20,5))
            plt.title(f"{case} Potencia de las caras")
            plt.plot(fechas, P_faces[:,j], "-o", ms = 2, label = i[3:])
            plt.xlabel("Fecha")
            plt.ylabel(f"Potencia [W/celda]  (Area celda = {area_cell} cm2)")
            plt.legend()
            plt.grid()
            plt.savefig(file+f"/{case}/pot_faces.png")
    return P_faces, Potencia

def pot_total(fechas, P_faces,  fechas_cons, P_cons, cell_faces, area_cell, periodo, report, file, Potencia, graph = None, case = ""):
    tiempo = to_tiempo(fechas)
    P_faces = np.array(P_faces)
    P_total = [np.sum(np.array(P_faces[i,:])*cell_faces) for i in range(len(P_faces[:,0]))]
    P_gen_media = np.trapz(P_total, tiempo)/tiempo[-1]
    print (f"Potencia media generada = {np.round(P_gen_media,2)} W/s")
    P_gen_media_orb = np.trapz(P_total, tiempo)/periodo
    print (f"Potencia media orbital generada = {np.round(P_gen_media_orb,2)} W/órbita")
    report.write("\n")
    report.write(f"Potencia media generada = {np.round(P_gen_media,2)} W/s")
    report.write("\n")
    report.write(f"Potencia media orbital generada = {np.round(P_gen_media_orb,2)} W/órbita")
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
    Potencia[f"{case[:3]}-Total-W"] = P_total
    return P_total, Potencia

def energy(fechas, potencia, report = None, graph = None, case = None,):
    '''case: (str) Caso particular que se está calculando. Se imprimirá en la consola.
    Si case distinto de None, se creará un DataFrame y la salida seran 2 variables (una tupla). En ese caso debería quedar: E, Consumo = energia(..)'''
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
        Consumo = pd.DataFrame()
        Consumo[ "Tiempo-s"] = tiempo
        Consumo["Tiempo-fechas"] = tiempo
        Consumo["Potencia-W"] = potencia
        Consumo["Energía-J"] = E
        Consumo["Energía-kWh"] = np.array(E)/3600000
        if type(fechas[0]) == datetime.datetime:
            Consumo["seg"] = tiempo
            Consumo["fechas"] = fechas
        else: 
            Consumo["seg"] = fechas
        return np.array(E), Consumo
    else:
        return np.array(E)
    
def energy_faces(fechas, P_faces, head, area_cell, report, file, Energia, graph = None, case = ""):
    E_faces = []
    E_faces_tot = []
    tiempo = to_tiempo(fechas)
    print(f"Energía total generada en la simulación de {np.round(tiempo[-1],2)} seg: ({case})")
    report.write(f"\nEnergía total por celda generada en la simulación de {np.round(tiempo[-1],2)} seg: ({case})")
    for j, i in enumerate(head):
        E = energy(fechas, P_faces[:,j]) # Jolues
        E_faces.append(np.array(E))
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
    Energia[f"{case[:3]}-fechas"] = fechas
    Energia[f"{case[:3]}-seg"] = tiempo
    for i, cara in enumerate(head):
        Energia[f"{case[:3]}-{cara[3:]}-J"] = E_faces[i]
        Energia[f"{case[:3]}-{cara[3:]}-kWh"] = E_faces[i]/3600000
    return np.array(E_faces), np.array(E_faces_tot), Energia

def energy_total(fechas, E_faces, cell_faces, report, file, Energia, graph = None, case = ""):
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
    Energia[f"{case[:3]}-Total-J"] = E_tot
    Energia[f"{case[:3]}-Total-kWh"] = np.array(E_tot)/3600000
    return np.array(E_tot), Energia

def to_periodic(new, n, fechas0, potencia0, start_date,  tiempo_max, report, graph = None):
    '''tiempo_max: tiempo máximo que se desea que dure la nueva función periodica.'''
    tiempo0 = np.array(to_tiempo(fechas0))
    tiempo_aux = cp.copy(tiempo0)
    t0 = time.process_time()
    n = int(np.ceil(n))
    for i in range(1,n):
        last = tiempo0 + tiempo_aux[-1]
        tiempo_aux = np.concatenate((tiempo_aux, last[1:]))
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

def orden(head, E_faces_tot, report, print_ = True): 
    E_faces_tot_ordenado = sorted(E_faces_tot, reverse=True) 
    posiciones = [] 
    for i in E_faces_tot: 
        posiciones.append(E_faces_tot_ordenado.index(i)) 
    if print_:
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

def distribution_cells(consumo_total, gen_faces, head, area_cell, report, distrib_type, caras = 1, distribucion = 1, case = ""):
    '''Si len(caras) != 1: hay que elegir valor para distribución distinto de 1'''
    gen_tot = []
    if distrib_type == "%":
        for i in caras:
            if abs(gen_faces[i])< 1e-4:
                print(f"--!-- La cara {i} elegida para poner paneles tiene radiación nula.")
                report.write("\n")
                report.write(f"--!-- La cara {i} elegida para poner paneles tiene radiación nula.")
        report.write("\n")
        print('{:^10} {:^5} {:^15} {:^15} {:^15}'.format('Panel', "(%)", "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        report.write("\n")
        report.write('{:^10} {:^5} {:^15} {:^15} {:^15}'.format('Panel', "(%)", "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        cell_faces = np.zeros(len(head))
        sup_cells = np.zeros(len(head))
        for j, i in enumerate(distribucion):
            gen_panel = consumo_total*i
            i_cara = caras[j]
            cell_per_panel = gen_panel/gen_faces[i_cara]
            cell_faces[i_cara] = np.ceil(cell_per_panel)
            sup_cells [i_cara] = cell_faces[i_cara] * area_cell *0.0001 # m^2
            print('{:^10} {:^5} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], i*100, np.round(gen_panel/3600000,2), int(np.ceil(cell_per_panel)), np.round(sup_cells[i_cara],2)))
            report.write("\n")
            report.write('{:^10} {:^5} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], i*100, np.round(gen_panel/3600000,2), int(np.ceil(cell_per_panel)), np.round(sup_cells[i_cara],2)))
            gen_tot.append(gen_panel)
    elif "sup_disp" in distrib_type:
        report.write("\n")
        print('{:^10} {:^15} {:^15} {:^15}'.format('Panel', "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        report.write("\n")
        report.write('{:^10} {:^15} {:^15} {:^15}'.format('Panel', "Generación (kWh)", 'Cantidad celdas', 'Superficie total'))
        try:
            posiciones = orden(head, gen_faces*distribucion, report, False)
        except ValueError:
            print(f"La logitud de sup_disp ({len(distribucion)}) no coincide con la cantidad de paneles en la simulación ({len(gen_faces)}). ")
            raise SystemExit
        cells = sum([np.floor(i) for i in np.array(distribucion)/(area_cell/10000)])
        # cells = np.floor(np.sum(np.array(distribucion))/(area_cell/10000))
        cell_faces = [0] * len(head)
        j = 0   # iterar caras
        i = 0   # celdas ubicadas
        gen = 0
        while i<cells and gen<consumo_total:
            if "max" in distrib_type:
                if j not in posiciones:
                    print("--!-- El espacio disponible en las caras no alcanza.")
                    report.write("--!-- El espacio disponible en las caras no alcanza.")
                    break
                i_cara = posiciones.index(j)
            elif "fix" in distrib_type:
                if j == len(caras):
                    print("--!-- El espacio que hay en las caras seleccionadas no alcanza.")
                    report.write("--!-- El espacio que hay en las caras seleccionadas no alcanza.")
                    break
                i_cara = caras[j]
            cell_max_cara = int(distribucion[i_cara]/(area_cell/10000))
            gen_panel_max = cell_max_cara * gen_faces[i_cara]
            if consumo_total-gen>gen_panel_max:
                cell_faces [i_cara] = cell_max_cara
                gen_panel = gen_panel_max
            else: 
                gen_panel = consumo_total-gen
                cell_faces [i_cara] = np.ceil(gen_panel/gen_faces[i_cara])
            sup_cells = cell_faces[i_cara] * (area_cell/10000)
            i+=cell_faces [i_cara]
            gen += gen_panel
            j+=1
            print('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,2), cell_faces[i_cara], np.round(sup_cells,2)))
            report.write("\n")
            report.write('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,2), cell_faces[i_cara], np.round(sup_cells,2)))
        print("\n")
        report.write("\n")
    return cell_faces, sup_cells

def balance_energy(n, fechas_cons, E_cons, fechas_gen, E_tot, report, file, Energia, graph = None, case = ""):
    tiempo_cons = to_tiempo(fechas_cons)
    tiempo_gen = to_tiempo(fechas_gen)
    # if n<1:
    #     E_tot_bis = np.interp(tiempo_cons, tiempo_gen, E_tot)
    #     balance = E_tot_bis- E_cons
    #     tiempo_bis = tiempo_cons
    # else:
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
    if tiempo_neg_balance != []:
        batery_cap = np.sum(energy(tiempo_neg_balance, neg_balance))
        print(f"Energía mínima necesaria en la batería: {-np.round(batery_cap/1000,2)} kJ = {-np.round(batery_cap/3600000,4)} kWh\n            (Energía acumulada para valores negativos del balance)")
        report.write("\n")
        report.write(f"Energía mínima necesaria en la batería: {-np.round(batery_cap/1000,2)} kJ = {-np.round(batery_cap/3600000,4)} kWh\n            (Energía acumulada para valores negativos del balance)")
    else:
        print("El balance nunca es negativo.")
        report.write("\n")
        report.write("El balance nunca es negativo.")
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
    Energia[f"{case[:3]}-Balance-kJ"] = balance/1000
    Energia[f"{case[:3]}-Balance-kWh"] = balance/3600000
    Energia[f"{case[:3]}-Balance-time-s"] = tiempo_bis
    Energia[f"{case[:3]}-Balance_fecha"] = fechas_bis
    return tiempo_bis, balance, Energia

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

def add_degT(Imp, Vmp, T, T_Imp, T_Vmp, T_ref = 28):
    return (Imp+T_Imp*(T-T_ref)) * (Vmp+T_Vmp * (T-T_ref))


def start_report(file, name, power_case):
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

def calc_distribution_cell(E_gen_tot, distrib_type, head, cells, area_cell, report, sup_disp = 1, posiciones = 1):
    
    print('{:^10} {:^15} {:^15} {:^15}'.format('Panel', "Generación (kWh)", 'N_celdas', 'Superficie (m^2)'))
    report.write("\n")
    report.write('{:^10} {:^15} {:^15} {:^15}'.format('Panel', "Generación (kWh)", 'N_celdas', 'Superficie (m^2)'))
    
    if "max" in distrib_type: 
        if len(sup_disp) != len(head):
            sup_disp = ast.literal_eval(input(f"El vector de superficie disponible por panel tiene una longitud distinta a la de la cantidad de paneles. Ingresar una lista de longitud {len(head)}:  "))
        if distrib_type[0] =="S":
            cells = np.floor(cells/(area_cell/10000))
        cell_faces = [0] * len(head)
        posiciones = orden(head, E_gen_tot*sup_disp, report, False)
        j = 0   # iterar caras
        i = 0   # celdas ubicadas
        while i<cells:
            try:
                i_cara = posiciones.index(j)
            except ValueError:
                print("--!-- La cantidad de superficie disponible no alcanza para ubicar todas las celdas.")
                break
            cell_max_cara = np.floor(sup_disp[i_cara]/(area_cell/10000))
            if cells-i >= cell_max_cara:
                cell_faces [i_cara] = cell_max_cara
            else: 
                cell_faces [i_cara] = cells-i
            sup_cells = cell_faces[i_cara] * (area_cell/10000)
            gen_panel = cell_faces[i_cara] * E_gen_tot[i_cara]
            i+=cell_faces[i_cara]
            j+=1
            print('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,2), cell_faces[i_cara], np.round(sup_cells,2)))
            report.write("\n")
            report.write('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,2), cell_faces[i_cara], np.round(sup_cells,2)))
        print("\n")
        report.write("\n")
    elif "fix" in distrib_type:
        if len(cells) != len(head):
            cells = ast.literal_eval(input(f"cell_faces la dimensión de cell_faces no coincide con la cantidad de paneles con cosenos calculados ({len(head)}). Ingrese una lista de dimensión {len(head)} de la cantidad de celdas por cara:\n"))
        if distrib_type[0] =="S":
            cells = [np.floor(i) for i in np.array(cells)/(area_cell/10000)]
        cell_faces = cells
        sup_cells = area_cell*np.array(cells)/10000
        for i_cara in range(len(cells)):
            sup_cell = sup_cells[i_cara]
            gen_panel = cell_faces[i_cara] * E_gen_tot[i_cara]
            if cells[i_cara] > 0:
                print('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,4), cell_faces[i_cara], np.round(sup_cell,2)))
                report.write("\n")
                report.write('{:^10} {:^15} {:^15} {:^15}'.format(head[i_cara][-7:], np.round(gen_panel/3600000,4), cell_faces[i_cara], np.round(sup_cell,2)))
        print("\n")
        report.write("\n")
    else:
        print("distrib_type INVALIDO")
        distrib_type = input("Ingrese un valor válido para distrib_type:\n")
        cell_faces, sup_cells = calc_distribution_cell(distrib_type, head, cells, area_cell, sup_disp, posiciones)
    return cell_faces, sup_cells