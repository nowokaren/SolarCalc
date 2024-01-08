# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from load_gmat import *
import os
import shutil
import time
import numpy as np

def set_title(Epoch, SMA, ECC, INC, RAAN, AOP, TA, ARB, BAV0, BAV1, BAV2, ACT, BCV0, BCV1, BCV2, model_file, Attitude):
    for i  in [BAV0, BAV1, BAV2, BCV0, BCV1, BCV2,]:
       if i<1 and i>0:
           i=str(i)[2:]
    title = f"Ep{Epoch[3:6]}{Epoch[12:14]}-SMA{np.round(SMA)}-ECC{str(ECC)[2:]}-INC{np.round(INC)}-RAAN{np.round(RAAN)}-AOP{np.round(AOP)}-TA{TA}_{Attitude[:3]}{ARB[0]}{BAV0}{BAV1}{BAV2}-{ACT[0]}{BCV0}{BCV1}{BCV2}"
    return title

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

def set_script(script, file, values, step, duration, case = 1, force_overwrite = False):
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
        if name not in os.listdir("Runs"):   # Si no está la carpeta de esta corrida, crearla
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


def print_report(file, name, values, period, num_orbs, paneles, duration_run, step, case = 1):
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
            report.write(f"- {prop} = {val}\n")
        report.write(f"- Duración = {duration_run} días\n")
        report.write(f"- Step = {step} días\n")
        report.write("\nResultados:\n")
        report.write(f"- Periodo = {period} segundos\n")
        report.write(f"- Cantidad de órbitas simuladas = {num_orbs}\n")
        report.write("- Paneles solares:")
        if paneles != []:
            for j,i in enumerate(paneles):
                report.write(f"\n\tPanel {j+1}: {i}")
        else:
            report.write(" None")
    else:
        report = open(file+"/"+name+"-settings.txt", "a")
        # for i in range(1+len(values)+4+len(paneles)+2):
        #     report.readline()
        report.write(f"\n\nParámetros de la simulación {case}:\n\n")
        for val, prop in zip(values, props):
            report.write(f"- {prop} = {val}\n")
        report.write(f"- Duración = {duration_run} días\n")
        report.write(f"- Step = {step} días\n")
        report.write(f"\nResultados {case}:\n")
        report.write(f"- Cantidad de órbitas simuladas = {num_orbs}\n")
        
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
    txt = [i for i in os.listdir("Runs") if i[-4:]==".txt"] # Muevo todos los .txt sueltos en Runs, a la carpeta de la corrida
    for i in txt:
        shutil.move("Runs/"+i, file+"/"+i)
    print(f"Fin de la simulacion de GMAT  (Duracion = {t1-t0} segundos)")
    if case != 1:
        title = file.split("/")[-1]
        os.remove(file+"/Eclipse-"+title+"-"+case+".txt")