# -*- coding: utf-8 -*-
"""
@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 1.20
Fecha: 20/3/2023
"""

from load_gmat import *
import os
import shutil
import time
import numpy as np
import ast
import datetime
import copy as cp

def set_title_nadir(Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV0, BAV1, BAV2, ACT, BCV0, BCV1, BCV2, model_file):
    for i  in [BAV0, BAV1, BAV2, BCV0, BCV1, BCV2,]:
       if i<1 and i>0:
           i=str(i)[2:]
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA,2)}-ECC{str(ECC)[2:]}-INC{np.round(INC,1)}-RAAN{np.round(RAAN,2)}-AOP{np.round(AOP,2)}-TA{round(TA,2)}_Nad{ARB[0]}{BAV0}{BAV1}{BAV2}-{ACT[0]}{BCV0}{BCV1}{BCV2}"
    return title

def set_title_tak(Epoch, SMA, ECC, INC, RAAN, AOP, TA, DCM11, DCM12, DCM13, DCM21, DCM22, DCM23, DCM31, DCM32, DCM33, AVX, AVY, AVZ, model_file):
    DCM = cp.copy(DCM11)
    for i in [DCM11, DCM12, DCM13, DCM21, DCM22, DCM23, DCM31, DCM32, DCM33]:      
        if i - round(i) != 0.0:
            i = round(i - round(i),2)
        DCM += i
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA,2)}-ECC{str(ECC)[2:]}-INC{np.round(INC,1)}-RAAN{np.round(RAAN,2)}-AOP{np.round(AOP,2)}-TA{round(TA,2)}_3AxKin{DCM}-AngVel{AVX}{AVY}{AVZ}"
    return title

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
    step = mo.getfloat("step") # 0.0001     # 0.5    # Days Calcular
    duration_run = mo.getfloat("duration_run") # 0.25  #30      # Days Calcular
    stepWC = mo.getfloat("stepWC") # 0.0001     # 0.5    # Days Calcular
    duration_runWC = mo.getfloat("duration_runWC") # 0.25  #30      # Days Calcular
    script_name = mo.get("script_name") #"orbit"
    paneles = ast.literal_eval(mo.get("paneles")) # [[-1 , 1, 1], [-0.66795061,  0.74419228,  0.00445376]]
    if abs(mo.getfloat("start_periodWC")) < 0.01:
        start_periodWC = mo.getfloat("start_periodWC")
        end_periodWC = mo.getfloat("end_periodWC")
    else:    
        start_periodWC = datetime.datetime.strptime(mo.get("start_periodWC"), '%d %b %Y %H:%M:%S.%f')
        end_periodWC = datetime.datetime.strptime(mo.get("end_periodWC"), '%d %b %Y %H:%M:%S.%f')
    if abs(mo.getfloat("start_periodWC")) < 0.01:
        startWC = mo.getfloat("startWC")
        endWC = mo.getfloat("endWC")
    else:
        startWC = datetime.datetime.strptime(mo.get("startWC"), '%d %b %Y %H:%M:%S.%f')
        endWC = datetime.datetime.strptime(mo.get("endWC"), '%d %b %Y %H:%M:%S.%f')
    return name, Epoch, SMA, ECC, INC, RAAN, AOP, TA, step, duration_run, stepWC, duration_runWC, script_name, paneles, start_periodWC, end_periodWC, startWC, endWC
   
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

def set_params_op_avanz(section_ini):
    oa = section_ini
    force_overwrite = oa.getboolean("force_overwrite") # True
    corte_subrun = oa.getfloat("corte_subrun") # 1 # días
    scale_actitud = oa.getfloat("scale_actitud") #700
    scale_sol = oa.getfloat("scale_sol") #2000
    scale_panel = oa.getfloat("scale_panel") #2000
    lines = oa.getboolean("lines") # True  #Agregar lineas al gráfico de la simulación larga.
    model_file = "\'"+os.getcwd()+"\\"+ oa.get("model") + "\'"
    return force_overwrite, corte_subrun, scale_actitud, scale_sol, scale_panel, lines, model_file

def set_prop(objeto, script, prop, valor):
    '''Modifica el valor de una propiedad en el script de GMAT según un elemento (satélite, Tierra, Sol, etc.), una propiedad (SMA, Epoch, etc.) y su valor correspondiente.
    - objeto: Objeto (satélite) a modificar.
    - script: Nombre del script de base que se desea utilizar para la simulación.
    - prop: Propiedad del objeto que se desea modificar.
    - valor: Valor que se desea asignar a dicha propiedad.'''
    for j, i in enumerate(script):
        if f"GMAT {objeto}.{prop} =" in i:
            script[j] = f"GMAT {objeto}.{prop} = {valor};"
    return script

def set_duration(script, duration):
    '''Modifica la duración de la simulación en el script de GMAT (en días).
    - script: Nombre del script de base que se desea utilizar para la simulación.
    - duration: Duración de la simulación.'''
    for j, i in enumerate(script):
        if "While Sat.ElapsedDays" in i:
            script[j] = f"While Sat.ElapsedDays < {duration};"
    return script

def set_step(script, step):
    '''Modifica el step de la simulación en el script de GMAT (en días)
    - script: Nombre del script de base que se desea utilizar para la simulación.
    - step: Paso de sampleo o toma de datos. GMAT escribirá los datos obtenidos por la simulación cada un step (en días).'''
    for j, i in enumerate(script):
        if "Propagate DefaultProp(Sat)" in i:
            script[j] = "   Propagate DefaultProp(Sat) {Sat.ElapsedDays = "+str(step)+ "};"
    return script

def set_script(script, file, values, step, duration, ini_path, case = 1, force_overwrite = False):
    '''Modifica la configuración definida para la simulación (parámetros orbitales, de sampleo, etc.) en el script de GMAT.
    - script: Nombre del script de base que se desea utilizar para la simulación.
    - file: Path de la carpeta de la corrida actual. "Runs/{title}" 
    - values: Lista [Epoch, SMA, ECC, INC, RAAN, AOP, TA, Attitude, AttitudeReferenceBody, BodyAlignmentVectorX, BodyAlignmentVectorY, BodyAlignmentVectorZ, AttitudeConstraintType, BodyConstraintVectorX, BodyConstraintVectorY, BodyConstraintVectorZ]
    - step: Paso de sampleo o toma de datos. GMAT escribirá los datos obtenidos por la simulación cada un step (en días).
    - duration: Duración de la simulación.
    - case: Nombre de la simulación a la que refiere (optativo)'''
    print("Configurando el script de GMAT.")
    f = open("Scripts_GMAT/"+script+".script", "r")
    aux = f.read()
    script_aux = aux.split("\n")
    props = ["Epoch", "SMA", "ECC", "INC", "RAAN", "AOP", "TA", "AttitudeReferenceBody", "BodyAlignmentVectorX", "BodyAlignmentVectorY", "BodyAlignmentVectorZ", "AttitudeConstraintType", "BodyConstraintVectorX", "BodyConstraintVectorY", "BodyConstraintVectorZ", "ModelFile"]
    for val, prop in zip(values, props):
        script_aux = set_prop("Sat", script_aux, prop, val)
    name = file.split("/")[-1]
    if case == 1:
        if name not in os.listdir("../Runs"):   # Si no está la carpeta de esta corrida, crearla
            os.mkdir(file)
        else: 
            if force_overwrite == False:
                respuesta = input("---- Ya existe una simulacion con esta orbita y actitud. Esta segura de que quiere hacer esta corrida? Se borraran los resultados anteriores. S: Si, N: No\n")
                if respuesta in ["s", "S"]:
                    shutil.rmtree(file)
                    os.mkdir(file)
                else:
                    raise SystemExit
            else:
                shutil.rmtree(file)
                os.mkdir(file)
    if "angulos" not in os.listdir(file):
        os.mkdir(file+"/angulos")
    if case != 1:
        name = name+"-"+case
        
    script_aux = set_prop("EclipseLocator1", script_aux, "Filename", f"\'Eclipse-{name}.txt\'")
    script_aux = set_prop("ReportFile1", script_aux, "Filename",f"\'{name}.txt\'")
    script_aux = set_step(script_aux, step)   # Days
    script_aux = set_duration(script_aux, duration)  # Days
    f = open("Scripts_GMAT/"+script+".script", "w")
    f.write("\n".join(script_aux))
    f.close()
    if case == 1:
        shutil.copy("Scripts_GMAT/"+script+".script", file + "/"+script+"-"+name.split("_")[0]+".script" )
        name = name.split("_")[0]
        shutil.copy("../"+ini_path[0], file+"/"+name+"-"+ini_path[0])


def print_report(file, name, values, period, num_orbs, axes, duration_run, step, name_plano_max, duracion_eclip_max = 1, case = 1):
    '''Genera un reporte con las propiedades de la simulación.
    - file: Path de la carpeta de la corrida actual. "Runs/{title}"
    - name: Nombre elegido para la corrida actual.
    - values: Lista [Epoch, SMA, ECC, INC, RAAN, AOP, TA, Attitude, AttitudeReferenceBody, BodyAlignmentVectorX, BodyAlignmentVectorY, BodyAlignmentVectorZ, AttitudeConstraintType, BodyConstraintVectorX, BodyConstraintVectorY, BodyConstraintVectorZ]
    - period: periodo de la órbita calculado por GMAT para la primera simulación
    - num_orb: número total de órbitas simuladas.'''
    props = ["Epoch", "SMA", "ECC", "INC", "RAAN", "AOP", "TA", "AttitudeReferenceBody", "BodyAlignmentVectorX", "BodyAlignmentVectorY", "BodyAlignmentVectorZ", "AttitudeConstraintType", "BodyConstraintVectorX", "BodyConstraintVectorY", "BodyConstraintVectorZ"]
    
    # Archivo .txt con los datos de la orbita

    if case == 1:
        report = open(file+"/"+name+"-resumen.txt", "w")
        report.write("Parámetros de la simulación:\n\n")
        for val, prop in zip(values, props):
            report.write(f"- {prop} =\t{val}\n")
        report.write(f"- Duración =\t{duration_run}\tdías\n")
        report.write(f"- Step =\t{step}\tdías\n")
        report.write("\nResultados:\n")
        report.write(f"- Periodo =\t{period}\tsegundos\n")
        report.write(f"- Cantidad de órbitas simuladas =\t{num_orbs}\n")
        report.write("- Paneles solares:")
        if len(axes) != 6:
            for j,i in enumerate(axes[6:]):
                if len(i)>1:
                    report.write(f"\n\tPanel {j+1}:\t{np.round(list(i[1]),4)}")
                else:
                    report.write(f"\n\tPanel {j+1}:\t{np.round(i,4)}")
            if name_plano_max != []:
                report.write(f"- Dirección óptima -  caso simulación 1 o larga =\t{np.round(axes[-1],4)}\t= Panel\t{len(axes[6:])}")
        else:
            report.write(" None")
        
    else:
        report = open(file+"/"+name+"-resumen.txt", "a")
        report.write(f"\n- Duración del eclipse más largo =\t{duracion_eclip_max}\tseg")
        report.write(f"\n\nParámetros de la simulación {case}:\n\n")
        for val, prop in zip(values, props):
            report.write(f"- {prop} =\t{val}\n")
        report.write(f"- Duración =\t{duration_run}\tdías\n")
        report.write(f"- Step =\t{step}\tdías\n")
        report.write(f"\nResultados {case}:\n")
        report.write(f"- Cantidad de órbitas simuladas =\t{num_orbs}\n")
        for i in range(len(name_plano_max)):            
            report.write(f"- Dirección óptima -  caso {name_plano_max[len(name_plano_max)-1-i]} =\t{np.round(axes[len(axes)-1-i],4)}\t= Panel\t{len(axes)-6-i}\n")
    report.close()


def run_gmat(script, file, case = 1):
    '''Ejecuta el script de GMAT y mueve los resultados desde OUTPUT_PATH hasta la carpeta donde se ubican todos los resultados de la simulación.
    - script: Nombre del script de base que se desea utilizar para la simulación.
    - file: Path de la carpeta de la corrida actual. "Runs/{title}'''
    print("Iniciando simulacion de GMAT.")
    gmat.LoadScript("Scripts_GMAT/"+script+".script")
    t0 = time.process_time()
    gmat.RunScript()
    t1 = time.process_time()

    # Muevo 
    txt = [i for i in os.listdir("../Runs") if i[-4:]==".txt"] # Muevo todos los .txt sueltos en Runs, a la carpeta de la corrida
    for i in txt:
        shutil.move("../Runs/"+i, file+"/"+i)
    print(f"Fin de la simulacion de GMAT  (Duracion = {t1-t0} segundos)")
    if case != 1:
        title = file.split("/")[-1]
        os.remove(file+"/Eclipse-"+title+"-"+case+".txt")