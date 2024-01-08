# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 12:29:43 2022

@author: Karen Nowogrodzki
Contact: nowo.karen@gmail.com
Version: 1.20
Fecha: 20/3/2023
"""

'''Para poder graficar los vectores fueron definidos como dos puntos en el espacio, es decir:
vector = [[x1,y1,z1],[x2,y2,z2]]'''
    
import pandas as pd
import math as mt
import numpy as np
import datetime
import matplotlib.pyplot as plt
import copy as cp
import matplotlib.lines as mlines
from scipy.optimize import minimize
import os


def rot3DX(x, a):
    a = a*np.pi/180
    R = np.array([[1,0,0],[0,mt.cos(a), -mt.sin(a)],[0,mt.sin(a), mt.cos(a)]], dtype = object)
    return np.matmul(x,R)

def rot3DY(x, a):
    a = a*np.pi/180
    R = np.array([[mt.cos(a),0, mt.sin(a)],[0,1,0],[-mt.sin(a),0, mt.cos(a)]], dtype = object)
    return np.matmul(x,R)

def rot3DZ(x, a):
    a = a*np.pi/180
    R = np.array([[mt.cos(a), -mt.sin(a),0],[mt.sin(a), mt.cos(a),0],[0,0,1]], dtype = object)
    return np.matmul(x,R)  

def rotEuler(x,w3,w2,w1):
    '''Rota un vector con los ángulos de Euler ingresados con la rotación 321.
    - x: (array 1,3) vector a ser rotado
    - w1, w2, w3: (floats) Ángulos de Euler generados por GMAT que representan la actitud del satélite en cada instante de tiempo.
    - x_rot: (array 1,3) vector rotado'''
    
    w1 *= np.pi/180
    w2 *= np.pi/180
    w3 *= np.pi/180
    c1 = mt.cos(w1)
    c2 = mt.cos(w2)
    c3 = mt.cos(w3)
    s1 = mt.sin(w1)
    s2 = mt.sin(w2)
    s3 = mt.sin(w3)
  
    R321 = np.array([[c2*c3, c2*s3, -s2], [-c1*s3+s1*s2*c3, c1*c3+s1*s2*s3, s1*c2], [s1*s3+c1*s2*c3, -s1*c3+c1*s2*s3, c1*c2]], dtype = object)
    x_rot = np.matmul(x,R321)
    return x_rot 

# Cálculo el ángulo entre dos vectores u y v
def ang_vec(u, v):
    '''Calcula en ángulo y el coseno entre dos vectores.
    - u, v: Vectores cuyo ángulo se quiere calcular.
    - ang: Angulo en grados
    - cos: coseno de dicho ángulo'''
    cos = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle
    ang = np.arccos(np.clip(cos, -1, 1))*180/np.pi
    return ang, cos

# Funtamentals of astrodynamics and applications (David A. Vallado) Sec.5.1.1, Pag.277
def time_params(UTC):
    '''Recibe un string de un Epoch y numeros enteros referentes a la fecha.
    - UTC: "21 Mar 2000 11:59:28.000"
    - yr, mon, day, hs, mins, sec: (ints) Año, Mes, Día, Hora, Minuto, Segundo. '''
    months = {"Jan":1 , "Feb":2 , "Mar":3 , "Apr":4 , "May":5 , "Jun":6 , "Jul":7 , "Aug":8 , "Sep":9 , "Oct":10 , "Nov":11 , "Dec":12}
    day = int(UTC[0:2])
    mon = months[UTC[3:6]]
    yr = int(UTC[7:11])
    hs = int(UTC[12:14])
    mins = int(UTC[15:17])
    sec = int(UTC[18:20])
    if sec==60:
        sec=59
    return yr, mon, day, hs, mins, sec

def sun_vector(yr, mon, day, hs, mins, sec):
    '''Calcula el vector entre la posición del sol y la posición del satélite. Tiene un error de 2.4°"
    - yr, mon, day, hs, mins, sec: (ints) Año, Mes, Día, Hora, Minuto, Segundo.
    - sun_vect: (list 1,3) Vector Tierra → Sol en unidades astronómicas.'''
    JD = 367*yr - np.floor((7*(yr+np.floor((mon+9)/12)))/4) + np.floor(275*mon/9) + day + 1721013.5 + ((((sec/60)+mins)/60)+hs)/24  # Julian Date
    T = (JD -2451545)/(36525)                                                             # Number of Julian Centuries
    lambda_sun = 280.46 + 36000.771 * T                                                   # Mean Longitude of the Sun
    M_sun = (357.5291092 + 35999.05034 * T)  *np.pi/180                                   # Mean Anomaly for the Sun
    lambda_ecl = (lambda_sun + 1.914666471 + 0.019994643 * np.sin(2*M_sun)) * np.pi/180   # Mean Logitude of the ecleptic
    epsilon = (23.439291 -46.836 * T - 0.000183 * T**2  ) * np.pi/180                     # Obliquity of the ecleptic
    r = 1.000140612   - 0.016708617 * np.cos(M_sun) - 0.000139589 * np.cos(2*M_sun)       # Distance of the Earth to the Sun (Unidades Astronomicas)
    #r=1    # Normalizado                                                                              # Vector unitario
    sun_vect = [r * np.cos(lambda_ecl), r * np.sin(lambda_ecl)*np.cos(epsilon), r * np.sin(lambda_ecl)*np.sin(epsilon) ]
    return sun_vect

def data_eclip(file, case = 1, raise_exit = True):
    '''Lee el archivo de eclipses de la simulación deseada y genera listas de tiempos de inicio, finalización y duración de cada eclipse.
    - file: Path de la carpeta de la corrida actual. "Runs/{title}"
    - case: Nombre de la simulación a la que refiere (optativo)
    - start: (list) Lista de datetimes de las fechas de los inicios de cada eclipse.
    - stop: (list) Lista de datetimes de las fechas de los finales de cada eclipse.
    - duration: (list) Lista de la duración de cada eclipse. '''
    name = file.split("/")[-1]  # No es el name de la corrida
    if case != 1:
        name = name +"-"+case
    path = file+"/Eclipse-"+name+".txt"
    start = [] 
    stop = []
    duration = []
    with open(path) as f:
        for j, i in enumerate(f):
            if j>2:
                if i == "\n":
                    break
                start.append(datetime.datetime(*time_params(i[:20]))) 
                stop.append(datetime.datetime(*time_params(i[28:48])))
                duration.append(float(i[58:66])) # seg
    if duration == []:
        print(" ---- El periodo seleccionado no tiene eclipses")
        if raise_exit:
            raise SystemExit
    return start, stop, duration

def data_orbit(file, scale = 1000, case = 1):
    '''Lee el archivo de la órbita de la simulación deseada y entrega listas con la información obtenida por la simulación (tiempos, posiciones, ángulos, parámetros orbitales keplerianos, período, numero de órbita, posición del sol)
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}
    - scale: (int) Escala en la que hay que dibujar los vectores para que se vean.
    - case: (str) Nombre de la simulación a la que refiere (optativo)
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - position: (list) Lista de vectores (listas) que indican la posición del satélite en cada step de la órbita simulada.
    - sun_vect: (list) Lista de vectores (listas) que indican la posición del sol en cada step de la órbita simulada.
    - euler_angs: (list) Lista de vectores (listas) que indican los ángulos de Euler de la actitud del satélite en cada step de la órbita simulada.
    - t_orb_relativo: (pandas - Series) Lista de tiempo de la órbita transcurrido desde el Epoch inicial en segundos.
    - period: (float) Periodo del satélite en la órbita en segundos.
    - num_orb: (array) Lista que indica el número de órbita simulada al que pertenece cada step de la órbita simulada.
    - raan, aop, ta, inc, ecc, sma: (pandas - Series) Listas de cada uno de los 6 parámetros keplerianos orbitales para cada step de la órbita simulada.
    - X_orb, Y_orb, Z_orb: (array) Listas con las componentes x, y, z de cada step de la órbita simulada, respectivamente.'''
    name = file.split("/")[-1]
    if case != 1:
        name = name + "-" + case
    ORBIT = pd.read_csv(file+"/"+name+".txt", sep = "\t")
    
    # Set parametros de la orbita
    t_orb = [datetime.datetime(*time_params(i)) for i in ORBIT["Sat.UTCGregorian"]] # t_orb absoluto
    X_orb = np.array([ORBIT["Sat.EarthMJ2000Eq.X"]], dtype = object)
    Y_orb = np.array([ORBIT["Sat.EarthMJ2000Eq.Y"]], dtype = object)
    Z_orb = np.array([ORBIT["Sat.EarthMJ2000Eq.Z"]], dtype = object)
    X_sun = np.array([ORBIT["Sun.EarthMJ2000Eq.X"]], dtype = object)
    Y_sun = np.array([ORBIT["Sun.EarthMJ2000Eq.Y"]], dtype = object)
    Z_sun = np.array([ORBIT["Sun.EarthMJ2000Eq.Z"]], dtype = object)

    position = [[X_orb[0,i],Y_orb[0,i],Z_orb[0,i]] for i in range(len(X_orb[0]))] # Lista de posiciones del satélite en la órbita
    sun_vect = [[X_sun[0,i],Y_sun[0,i],Z_sun[0,i]] for i in range(len(X_sun[0]))] # Lista de posiciones del sol 
    euler_angs = [[ORBIT["Sat.EulerAngle1"].values[i],ORBIT["Sat.EulerAngle2"].values[i],ORBIT["Sat.EulerAngle3"].values[i]] for i in range(len(ORBIT["Sat.EulerAngle1"].values))] # Actitud del satélite en Angulos de Euler
    t_orb_relativo = ORBIT["Sat.ElapsedSecs"]  # t_orb relativo al inicio en segundos
    period = np.mean(ORBIT["Sat.Earth.OrbitPeriod"].values)
    num_orb = num_period(t_orb_relativo, period)
    raan = ORBIT["Sat.EarthMJ2000Eq.RAAN"]
    aop = ORBIT["Sat.EarthMJ2000Eq.AOP"]
    ta = ORBIT["Sat.Earth.TA"]
    inc = ORBIT["Sat.EarthMJ2000Eq.INC"]
    ecc = ORBIT["Sat.Earth.ECC"]
    sma = ORBIT["Sat.Earth.SMA"]
    return t_orb, position, sun_vect, euler_angs, t_orb_relativo, period, num_orb, raan, aop, ta, inc, ecc, sma, X_orb, Y_orb, Z_orb 

def getpos_eclip(t_orb, position, start, stop, duration, num_orb, file, ax1, axes = False, graph = False, color = "red", case = 1):
    '''Con los datos de data_eclip y data_orbit genera listas de los index y posición XYZ de todos los puntos en los que el satélite se encuentra eclipsado.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - position: (list) Lista de vectores (listas) que indican la posición del satélite en cada step de la órbita simulada.
    - start: (list) Lista de datetimes de las fechas de los inicios de cada eclipse.
    - stop: (list) Lista de datetimes de las fechas de los finales de cada eclipse.
    - duration: (list) Lista de la duración de cada eclipse.
    - num_orb: (array) Lista que indica el número de órbita simulada al que pertenece cada step de la órbita simulada.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.
    - graph: (bool)Indica si se desean resaltar en el gráfico las posiciones de los steps en los que el satélite está eclipsado o no. (True o False respectivamente)
    - color: (str) Indicar color con el que se quiere graficar el eclipse
    - case: (str) Nombre de la simulación a la que refiere (optativo)
    - eclipI: (array) Lista de los index de los steps en los que el satélite se encuentra eclipsado.
    - eclip: (array) Array de arrays que contienen información de cada eclipse → [start (datetime), end (datetime), duration (float), [índices (int)], [X (list), Y (list), Z (list)] (array), núm de eclipse (int)]'''
    eclipI = []
    position = np.array(position)
    X_orb = position[:,0]
    Y_orb = position[:,1]
    Z_orb = position[:,2]
    h = 0
    d = 0
    eclip = []
    title = file.split("/")[-1]
    for i,j,k in zip(start,stop,duration):
        if j-i < datetime.timedelta( minutes = 2 ):
            d += k
            if h==0:
                on = i
            h=1
            continue
        else:
            on = i
            h = 0
            d += k
            eclipIndex = [m for m in range(len(t_orb)) if t_orb[m]>on and t_orb[m]<j]
            if eclipIndex == []:
                continue
            X_eclip = np.array([[X_orb[m] for m in range(len(t_orb)) if t_orb[m]>on and t_orb[m]<j]], dtype = object)
            Y_eclip = np.array([[Y_orb[m] for m in range(len(t_orb)) if t_orb[m]>on and t_orb[m]<j]], dtype = object)
            Z_eclip = np.array([[Z_orb[m] for m in range(len(t_orb)) if t_orb[m]>on and t_orb[m]<j]], dtype = object)
            eclipI += eclipIndex
            eclip.append([on, j, d, eclipIndex, X_eclip, Y_eclip, Z_eclip, num_orb[eclipIndex[0]]])
            d = 0
            if graph == True:
                ax1.plot_wireframe(X_eclip, Y_eclip, Z_eclip, linewidth=4, color = color, label = "Eclipse")
                figname = "/actitud.png"
                if case != 1:
                    figname = f"/actitud-{case}.png"
                plt.title(title)
                plt.savefig(file+figname)
    eclip = np.array(eclip, dtype = object)
    info_leyenda = [("lightblue", "X+"), ("blue", "Y+"), ("darkblue","Z+"),("y", "Sol"), ("yellow", "Vernal Equinox"), ("red", "Eclipse")]
    if len(axes)>5:
        for j, i in enumerate(axes[6:]):
            info_leyenda.append(((0.1+ 5/(j+6), 0.2, 0.5), f"Panel {j+1}"))
    ax1.legend([mlines.Line2D([0], [0], color=clave[0], lw=2) for clave in info_leyenda],[clave[1] for clave in info_leyenda], loc='upper left') 
    return eclipI, eclip

def graph_earth(ax1):
    '''grafica la Tierra en la figura actitud.png'''
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    R_e = 6378
    x = R_e * np.outer(np.cos(u), np.sin(v))
    y = R_e * np.outer(np.sin(u), np.sin(v))
    z = (R_e*0.8) * np.outer(np.ones(np.size(u)), np.cos(v))  # Achatamiento
    ax1.plot_surface(x, y, z, color='b')
    return x, y, z

def sist_ref(ax1 = 1,  graph = False, graph_equinox = False, color = "grey"):
    '''Genera los ejes unitarios del sistema de referencia ECI (y sus versiones negativas) → x = Vernal equinox
    - scale: (int) Escala en la que hay que dibujar los vectores para que se vean.
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.
    - graph: (bool) Indica si se desean graficar los ejes de coordenadas en el gráfico las posiciones de los steps. (True o False respectivamente)
    - graph_equinox: (bool) Indica si se desea graficar la dirección de x (equinox Vernal)
    - color: (str) color con el que se van a graficar los ejes.
    - i_s, j_s, k_s: (array de 2 array de 3 floats c/u) Versores positivos del sistema de referencia ECI.
    - negi_s, negj_s, negk_s: (array de 2 array de 3 floats c/u) Versores negativos del sistema de referencia ECI.'''
    i_s = np.array([[0,0,0],[1,0,0]], dtype = object) 
    j_s = np.array([[0,0,0],[0,1,0]], dtype = object)
    k_s = np.array([[0,0,0],[0,0,1]], dtype = object)
    negi_s = np.array([[0,0,0],[-1,0,0]], dtype = object)
    negj_s = np.array([[0,0,0],[0,-1,0]], dtype = object)
    negk_s = np.array([[0,0,0],[0,0,-1]], dtype = object)
    if graph_equinox == True:
        i_sun = cp.copy(i_s)*1000
        ax1.plot_wireframe(np.array([i_sun[:,0]], dtype = object), np.array([i_sun[:,1]], dtype = object), np.array([i_sun[:,2]], dtype = object), color = "yellow", label = "Vernal Equinox")
    if graph == True:
        for l in [i_s,j_s,k_s]:
            ax1.plot_wireframe(np.array([l[:,0]], dtype = object), np.array([l[:,1]], dtype = object), np.array([l[:,2]], dtype = object), color = color)
    return i_s, j_s, k_s, negi_s, negj_s, negk_s

def get_sunvect(t_orb, ax1, orb = False, scale = False, graph = False, color = "y"):
    '''Calcula y grafica el vector entre la posición del sol y la posición del satélite
    - t_orb: (list) Lista de parámetros temporales de los datetimes de cada step de la órbita simulada.
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.
    - orb: (array 1,3) Vector posición del satélite. (Solo ingresar si graph=True)
    - scale: (int) Escala para graficar los vectores (Solo ingresar si graph=True)
    - color: (str) Color en el que se quiere graficar el vector sol.
    - s =  sun_vect: (list 1,3) Vector Tierra → Sol en unidades astronómicas. '''
    s = sun_vector(*t_orb)
    if graph == True:
        sun = np.add(orb,np.array([[0,0,0],s])*scale)         # Translación a la posición del satélite en el instante o posición [m]
        ax1.plot_wireframe(np.array([sun[:,0]]), np.array([sun[:,1]]), np.array([sun[:,2]]), color = color)    
        sun = np.add(orb,np.array([[0,0,0],s])*scale)         # Translación a la posición del satélite en el instante o posición [m]
        ax1.plot_wireframe(np.array([sun[:,0]]), np.array([sun[:,1]]), np.array([sun[:,2]]), color = color)
    return s

def add_eclip(cos, eclipIndex, eclip_stamp):
    '''Anula todos los cosenos de los puntos de la órbita donde el satélite se encuentra eclipsado y anula cualquier coseno negativo (agrega las sombras).
    - index: (list) Lista de index de cada step de la órbita simulada.
    - cos: (array) Lista de listas de cosenos de cada eje (len=6) para cada step de la órbita simulada.
    - eclipIndex (list) Lista de indices de los steps en los que el satélite esta eclipsado.
    - cos: (array) Idem cos anterior con las sombras y eclipses.'''
    cos = np.array(cos)
    n = len(cos[0])
    for h, i in enumerate(cos):
        if h in eclipIndex:
            cos[h] = np.ones(n)*eclip_stamp
    cos = np.array(cos)                
    return cos

def add_shadow(cos):
    '''Anula todos los cosenos de los puntos de la órbita donde el satélite se encuentra eclipsado y anula cualquier coseno negativo (agrega las sombras).
    - cos: (array) Lista de listas de cosenos de cada eje (len=6) para cada step de la órbita simulada.
    - cos: (array) Idem cos anterior con las sombras.'''
    cos = cos.astype(float)
    cos[np.isnan(cos)]=0
    cos [cos < 0] = 0
    return cos

def graph_cos(t_orb, cos, file, dist_sun, case = 1, line = False, period = [], eclip_dur = []):
    '''Grafica cosenos y ángulos vs tiempo para cada cara en un mismo grafico.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - cos: (array) Lista de listas de cosenos de cada eje (len=6) para cada step de la órbita simulada.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - dist_sun: (list) Lista de distancia al sol del satélite en cada step de la órbita simulada.
    - case: (str) Nombre de la simulación a la que refiere (optativo)
    - line: (bool) Indica si se quieren graficar lineas entre los puntos o no. (True o False respectivamente)'''
    print(f"Creando graficos (case = {case}).")
    cos = cos.astype(float)
    cos[np.isnan(cos)]=0
    fig, ax = plt.subplots(2, figsize=(15, 10))
    ax[0].set_title(file.split("/")[-1])
    #fig.suptitle(file.split("/")[-1])
    label = ["Cara_+X", "Cara_+Y", "Cara_+Z", "Cara_-X", "Cara_-Y", "Cara_-Z"]
    if line == True:
        draw = "-o"
    else: 
        draw = "o"
    for i in range(len(cos[0])):
        if i>5:
            label.append(f"Panel_{i-5}")
        ang = np.arccos(list(cos[:,i]))*180/np.pi
        ax[0].plot(t_orb, cos[:,i], draw, label = label[i],markersize = 4, linewidth = 0.5)
        ax[1].plot(t_orb, ang, draw, label = label[i], markersize = 4,  linewidth = 0.5)
    ax[0].set_ylabel("Cosenos", fontsize=13)
    ax[1].set_ylabel("Angulos", fontsize=13)
    ax[1].set_xlabel("Date", fontsize=13)
    ax[0].legend(loc='upper left')
    ax[1].legend(loc='upper right')
    ax[0].grid()
    ax[1].grid()
    title = "Sol-NormalSat"
    if case != 1:
        title+="-"+case
    if case == "periodoWC":
        ax[1].set_title(f"Periodo = {str(period)} s    Eclipse = {str(eclip_dur)} s")
    plt.savefig(file+f"/{title}.png")
    
def report_cos(file, t_orb, dist_sun, cos, case = 1):
    '''Genera un reporte .txt con los cosenos de cada cara y panel desplegable, la distancia al sol y la fecha de cada punto de la órbita.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - cos: (array) Lista de listas de cosenos de cada eje (len=6) para cada step de la órbita simulada.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - dist_sun: (list) Lista de distancia al sol del satélite en cada step de la órbita simulada.
    - case: (str) Nombre de la simulación a la que refiere (optativo)'''
    name_cos = "cosenos"
    if case != 1:
        name_cos += "-"+case
    report = open(file+f"/angulos/{name_cos}"+".txt", "w") # Archivo .txt con los datos de la orbita
    head = "Fecha\tDistancia"
    label = ["Cara_+X", "Cara_+Y", "Cara_+Z", "Cara_-X", "Cara_-Y", "Cara_-Z"]
    for i in range(len(cos[0])):
        if i>5:
            label.append(f"Panel_{i-5}")
    for i in label:
        head += f"\tcos{i}"
    head+= "\n"
    report.write(head)
    for date, dist, c in zip(t_orb, dist_sun, cos):
        line = f"{date}\t{dist}"
        for i in range(len(label)):
            line+= f"\t{c[i]}"
        line+= "\n"
        report.write(line)
    report.close()
    
def report_ang(file, t_orb, dist_sun, cos, case = 1):
    '''Genera un reporte .txt con los angulos de cada cara y panel desplegable, la distancia al sol y la fecha de cada punto de la órbita.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - cos: (array) Lista de listas de cosenos de cada eje (len=6) para cada step de la órbita simulada.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - dist_sun: (list) Lista de distancia al sol del satélite en cada step de la órbita simulada.
    - case: (str) Nombre de la simulación a la que refiere (optativo)'''
    name_cos = "angulos"
    if case != 1:
        name_cos += "-"+case
    report = open(file+f"/angulos/{name_cos}"+".txt", "w") # Archivo .txt con los datos de la orbita
    head = "Fecha\tDistancia"
    label = ["Cara_+X", "Cara_+Y", "Cara_+Z", "Cara_-X", "Cara_-Y", "Cara_-Z"]
    for i in range(len(cos[0])):
        if i>5:
            label.append(f"Panel_{i-5}")
    for i in label:
        head += f"\tang{i}"
    head+= "\n"
    report.write(head)
    for date, dist, c in zip(t_orb, dist_sun, cos):
        line = f"{date}\t{dist}"
        for i in range(len(label)):
            ang = np.arccos(c[i])*180/np.pi
            line+= f"\t{ang}"
        line+= "\n"
        report.write(line)
    report.close()

def cos_dir(u, v):
    '''Calcula el vector de cosenos directores
    - u:(array 3,3) 3 ejes del sistema de referencia
    - v: (array 1,3) vector a calcular sus cosenos
    - panel_cos: (array 1,3) vector de cosenos directores del vector normal al panel desplegable.'''
    panel_cos = np.zeros(3)
    for i in range(3):
        panel_cos[i] = np.dot(u[i],v)/(np.linalg.norm(u[i])*np.linalg.norm(v))
    return panel_cos

def get_axbody(g, l, euler_angs, orb, scale_actitud = 700, scale_panel=1500, ax1 = 1, graph = False, ax_colors = ["lightblue", "blue", "darkblue"], labels = False):
    '''Genera los vectores rotados con los ángulos de Euler y opcionalmente lo grafica en la posición del satélite correspondiente.
    - g: (int) Indica si se grafica el eje en el sistema body. Si es 0,1 o 2 → Grafica los con los colores de ax_colors. Si es 3, 4 o 5 no los grafica (refiere a las caras negativas del satélite) y si es 6 o más refiere a los paneles desplegables y los grafica.
    - l: (array 2,3) vector unitario a modificar
    - euler_angs: (list 1,3) Vector de ángulos de Euler de la rotación de la actitud.
    - orb: (list 1,3) Vector de posición en la orbita.
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.
    - graph: (bool) Indica si se desean graficar los ejes de coordenadas en el gráfico las posiciones de los steps. (True o False respectivamente)
    - ax_body: (array 2,3) Lista del vector que indica el origen del vector rotado y el de final.
    - ax_body_rot: (lista 1,3) Vector de la actitud del satélite con origen en 0.'''
    ax_body = np.array([l[0] , rotEuler(l[1]-l[0], *euler_angs)], dtype = object)  # Rotación body según actitud en orbit[m]
    ax_body_rot = ax_body[1]
    if graph == True:                                     # Grafico el sistema de coordenadas rotado según la attitude (solo los positivos g<3)
        label = ["X+", "Y+", "Z+"]    
        if g>5:
            ax_body = np.add(ax_body*scale_panel, orb)                        # Translación body según orbita
            ax1.plot_wireframe(np.array([ax_body[:,0]], dtype = object), np.array([ax_body[:,1]], dtype = object), np.array([ax_body[:,2]], dtype = object), color = (0.1+ 5/g, 0.2, 0.5), label = f"Panel {g-5}")
        elif g<3:
            ax_body = np.add(ax_body*scale_actitud, orb)                        # Translación body según orbita
            ax1.plot_wireframe(np.array([ax_body[:,0]], dtype = object), np.array([ax_body[:,1]], dtype = object), np.array([ax_body[:,2]], dtype = object), color = ax_colors[g], label = label[g]) 
    return ax_body, ax_body_rot

def add_panels(axes, paneles):
    '''Agrega a los ejes del sistema de referencia, el vector normal a cada panel desplegable.
    - axes: (list) Lista de vectores (arrays 2,3) con origen y final de cada eje del sistema de referencia.
    - paneles: (list) Lista de vectores (list 1,3) de las normales a los paneles desplegables.
    - axes: Idem anterior pero con los ejes de los paneles agregados.'''
    if paneles != []:
        for panel in paneles:
            axes.append(np.array([[0,0,0],panel], dtype = object))
    return axes

def paneles_desplegables(paneles, bodysun_cos):
    '''Calcula el coseno entre el vector normal del panel y los del sistema body, después suma los ángulos panel-body y body-sun y entrega sus cosenos en panelsun_cos para calcularlos una vez realizada la ejecución del código cosenos.
    - paneles: (list) Lista de vectores (list 1,3) de las normales a los paneles desplegables.
    - bodysun_cos: (list) Lista de los cosenos obtenidos para las 6 caras del satélite.'''
    sist_body = np.eye(3)
    panelsun_cos = []
    for i in paneles:
        panel_cos = cos_dir(sist_body, i) # Cosenos del vector del panel respecto del sistema body[]  
        for k in bodysun_cos:
            panelsun_cos.append(np.dot(panel_cos,k[:3]))    # Cosenos del vector del panel respecto del sol
    return np.array(panelsun_cos)


def coseno_P_sol(P, cosenos):
  cosenos_P_sol = []
  for cara in cosenos:
    cosenos_P_sol.append(np.dot(P, cara))
  return cosenos_P_sol



def cos_sun (euler_angs, position, axes, t_orb, scale_sol, scale_actitud, scale_panel, ax1 = 1, sun_vect = False, calc_sun = False, graph = False):
    '''Calcula el coseno del ángulo entre el vector sol (Tierra-Sol) y la normal de cada cara de un satélite y su distancia al sol en un punto específico de su órbita
    - position: (list) Lista de vectores (listas) que indican la posición del satélite en cada step de la órbita simulada.
    - euler_angs: (list) Lista de vectores (listas) que indican los ángulos de Euler de la actitud del satélite en cada step de la órbita simulada.
    - axes: (list) Lista de vectores (arrays 2,3) con origen y final de cada eje del sistema de referencia.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - t_orb_relativo: (pandas - Series) Lista de tiempo de la órbita transcurrido desde el Epoch inicial en segundos.
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.
    - scale (int): Cuando agrandar los vectores para graficarlos.
    - sun_vect: (array 1,3) Vector posición del sol para cada step de la órbita simulada. Necesario ingresar solo si se indicó calc_sun = False
    - calc_sun: (bool) Indica si se desea calcular la posición del sol con la función get_sunvect o ingresarlo como sun_vect. (True o False respectivamente).
    - graph: (bool) Indica si se desean graficar los vectores de la actitud del satélite.
    - index: (list) Lista de los index de cada step de la órbita simulada.
    - cos: (list) Lista de vectores (array 1,# con #=6+ numero de paneles desplegables) con los cosenos de cada eje de sistema body y paneles desplegables.
    - dist_sun: (list) Lista de distancias entre el satélite y el sol para cada step de la órbita simulada.'''
    e = 0
    cos = []
    index = []
    dist_sun = []
    if len(position) < 15:
        dt = 1
    elif len(position) < 30: 
        dt = 2
    elif len(position) < 100:
        dt = 5
    else: 
        dt = int(len(euler_angs)/50)
    for m, o in enumerate(position):   # Recorro todos los punto de la orbita
        labels = m == dt-1   
        e += 1                                    # Grafico cada dt datos
        if e==dt:                       
            e=0
            graph = True
        else:
            graph = False
        cosm = []
        for g, ax in enumerate(axes):             # Recorro todas las caras del satelite de interés
            l = cp.copy(ax)

        # Ejes de actitud del satélite
            ax_body, aux = get_axbody(g, l, euler_angs[m], position[m], scale_actitud, scale_panel, ax1 = ax1, graph = graph, labels = labels)
        # Genero vector al sol para la posición position[m]
            if calc_sun == True:
                sun, r = get_sunvect(t_orb[m], position[m], scale_sol, graph = graph)
            else:
                sun = np.array(sun_vect[m])-np.array(position[m])
                r = np.linalg.norm(sun)
                sun = sun/r
                if graph == True and g == len(axes)-1:
                    s = np.add(position[m],np.array([[0,0,0],sun], dtype = object)*scale_sol)         # Translación a la posición del satélite en el instante o posición [m]
                    ax1.plot_wireframe(np.array([s[:,0]], dtype = object), np.array([s[:,1]], dtype = object), np.array([s[:,2]], dtype = object), color = "y", label = "Sol")
                # Array de cosenos entre vector sol y normal del satelite para el punto position[m]
                cosm.append(ang_vec(aux, sun)[1])
        dist_sun.append(r)
        cos.append(np.array([*cosm], dtype = object))
        index.append(m)
    info_leyenda = [("lightblue", "X+"), ("blue", "Y+"), ("darkblue","Z+"),("y", "Sol"), ("yellow", "Vernal Equinox")]
    if len(axes)>5:
        for j, i in enumerate(axes[6:]):
            info_leyenda.append(((0.1+ 5/(j+6), 0.2, 0.5), f"Panel {j+1}"))
    ax1.legend([mlines.Line2D([0], [0], color=clave[0], lw=2) for clave in info_leyenda],[clave[1] for clave in info_leyenda], loc='upper left')
    l = 0
    cos = np.array(cos)
    return index, cos, dist_sun


  
def num_period(t_orb_relativo, period):
    '''Entrega una lista que indica a qué numero de periodo pertenece cada punto de la órbita.
    - t_orb_relativo: (pandas - Series) Lista de tiempo de la órbita transcurrido desde el Epoch inicial en segundos.
    - period: (float) Periodo del satélite en la órbita en segundos.
    - num_orb: (array) Lista de números enteros que indican el número de período al que pertenece cada step de la órbita simulada.'''
    f = 1
    num_orb = []
    for j, i in enumerate(t_orb_relativo):
        while i > f*period:
            f+=1
        num_orb.append([j,f])
    num_orb = np.array(num_orb)
    return num_orb 




def cosenos(euler_angs, position, t_orb, X_orb, Y_orb, Z_orb, paneles, file, scale_sol, scale_actitud, scale_panel, calc_sun = False, sun_vect = False, graph = False, case = 1):
    '''Calcula los cosenos entre la cara del satelite y el sol, y dibuja ilustrativamente la orientacion del satelite
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - position: (list) Lista de vectores (listas) que indican la posición del satélite en cada step de la órbita simulada.
    - sun_vect: (list) Lista de vectores (listas) que indican la posición del sol en cada step de la órbita simulada.
    - euler_angs: (list) Lista de vectores (listas) que indican los ángulos de Euler de la actitud del satélite en cada step de la órbita simulada.
    - t_orb_relativo: (pandas - Series) Lista de tiempo de la órbita transcurrido desde el Epoch inicial en segundos.
    - X_orb, Y_orb, Z_orb: (array) Listas con las componentes x, y, z de cada step de la órbita simulada, respectivamente.
    - paneles: (list) Lista de vectores (list 1,3) de las normales a los paneles desplegables.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - sun_vect: (array 1,3) Vector posición del sol para cada step de la órbita simulada. Necesario ingresar solo si se indicó calc_sun = False
    - calc_sun: (bool) Indica si se desea calcular la posición del sol con la función get_sunvect o ingresarlo como sun_vect. (True o False respectivamente).
    - graph: (bool) Indica si se desean graficar los vectores de la actitud del satélite.
    - case: (str) Nombre de la simulación a la que refiere (optativo)
    - index: (list) Lista de los index de cada step de la órbita simulada.
    - cos: (list) Lista de vectores (array 1,# con #=6+ numero de paneles desplegables) con los cosenos de cada eje de sistema body y paneles desplegables.
    - dist_sun: (list) Lista de distancias entre el satélite y el sol para cada step de la órbita simulada.
    - ax1: (matplotlib.axes._subplots.Axes3DSubplot) Plot donde graficar la actitud y la órbita.'''
    fig = plt.figure(figsize = (8,8))
    ax1 = fig.add_subplot(111, projection='3d')
    # Sistema de coordenadas
    axes = list(sist_ref(ax1 = ax1, graph = True, graph_equinox = True))

    # Paneles desplegables (Vectores en el sistema de coordenadas body)
    axes = add_panels(axes, [v/np.linalg.norm(v) for v in paneles])
    if graph == True: 
        ax1.set_xlabel("X (km)")
        ax1.set_ylabel("Y (km)")
        ax1.set_zlabel("Z (km)")
        # Grafico la órbita y Configuro la figura
        if case != 1:
            ax1.plot_wireframe(X_orb, Y_orb, Z_orb, color = "grey")
        ax1.set_xlim(np.min(X_orb)-1000, np.max(X_orb)+1000)
        ax1.set_ylim(np.min(Y_orb)-1000, np.max(Y_orb)+1000)
        ax1.set_zlim(np.min(Z_orb)-1000, np.max(Z_orb)+1000)
    # Calculo los cosenos 
    index, cos, dist_sun = cos_sun (euler_angs, position, axes, t_orb, scale_sol, scale_actitud, scale_panel, ax1 = ax1, calc_sun = False, sun_vect=sun_vect, graph = graph)
    title = file.split("/")[-1]
    if graph == True:
        figname = "/actitud.png"
        if case != 1:
            figname = f"/actitud-{case}.png"
        # ax1.suptitle(title)
        # ax1.text(-1800,0,15500, f"\nPaneles = {paneles}", fontsize = 13)
        ax1.set_title(title)
        # ax1.text(-10,-10,0,  f"{paneles}")
        plt.savefig(file+figname)
    if case == 1:
        return index, cos, dist_sun, ax1, axes
    else:
        return index, cos, dist_sun, ax1

def eclip_duration(start, duration, t_orb,  file):
    '''Grafica la duración de cada eclipse en función de la fecha.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - start: (list) Lista de datetimes de las fechas de los inicios de cada eclipse.
    - duration: (list) Lista de la duración de cada eclipse.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"'''
    plt.figure(figsize=(10,2))
    st = []
    dur = []
    j = 0
    for date in t_orb[1:-1]:
        if date>start[0] and j==0:
            for i in start:
                st.append(i)
                dur.append(duration[start.index(i)])
            j=1
        if date>start[-1] or date<start[0]:
            st.append(date)
            dur.append(0)
    plt.plot(st,dur, "-o", markersize=3)
    plt.tick_params(labelsize=7)
    plt.xlabel("Fecha", fontsize=7)
    plt.ylabel("Duracion (s)", fontsize=7)
    plt.grid()
    plt.title(file.split("/")[-1])
    plt.savefig(file+"/duracion_eclip.png")
    
def select_timedelta(file, start, end, graph = False, case = 1):
    '''Genera un reporte y un gráfico de los cosenos de un período de tiempo seleccionado.
    - start: (datetime) Timestamp de inicio del período de tiempo.
    - end: (datetime) Timestamp de final del período de tiempo.
    - file: (str) Path de la carpeta de la corrida actual. "Runs/{title}"
    - graph: (bool) Indica si se desea graficar.'''
    fecha = []
    cosen = []
    dist = []
    if case != 1:
        case="-"+case
    else:
        case= ""
    with open(file+f"/angulos/cosenos{case}.txt", "r") as cosenos:
        for line in cosenos:
            if "Fecha" in line:
                label = []
                for i in line.split("\t")[2:]:
                    label.append(i[3:])
            elif datetime.datetime.strptime(line.split("\t")[0], "%Y-%m-%d %H:%M:%S")>start and datetime.datetime.strptime(line.split("\t")[0], "%Y-%m-%d %H:%M:%S")<end:
                fecha.append(datetime.datetime.strptime(line.split("\t")[0], "%Y-%m-%d %H:%M:%S"))
                c = line.split("\t")[2:]
                c[-1]= c[-1][:-1]
                for j,i in enumerate(c):
                    c[j]=float(i)
                cosen.append(c)
                dist.append(float(line.split("\t")[1]))
    plt.figure(figsize=(10,5))
    cosen = np.array(cosen)
    if graph == True:
        for j,i in enumerate(label):
            plt.plot(fecha, cosen[:,j], label = i)
        plt.legend()
        plt.xlabel("Fecha")
        plt.ylabel("cosenos Sol-CaraSat")
        plt.savefig(file+"/periodoWC")
    return fecha, cosen, dist

def subrunWC_params(start, stop, duration, startWC, endWC, t_orb, step, duration_run, period):
    '''Con las fechas de inicio, finalización y duración de los eclipses selecciona el eclipse más largo y define los parámetros de la órbita medidos por la simulación anual para el Epoch de la nueva corrida. El Epoch también lo elije tal que incluya el eclipse más largo del año.
    - start: (list) Lista de datetimes de las fechas de los inicios de cada eclipse.
    - stop: (list) Lista de datetimes de las fechas de los finales de cada eclipse.
    - duration: (list) Lista de la duración de cada eclipse.
    - t_orb: (list) Lista de datetimes de cada step de la órbita simulada.
    - period: (float) Periodo del satélite en la órbita en segundos.'''
    if abs(duration_run) < 0.01:
        duration_run =   0.5+(period/(24*60*60))   # Days Calcular
    if abs(startWC)  < 0.01: 
        startWC = start[duration.index(max(duration))]
        endWC = stop[duration.index(max(duration))]
    #ep = startWC.strftime("%d %b %Y %H:%M:%S.000")
    #step =     30 /(24*60*60)        # Dias Calcular
    i = 0
    while t_orb[i]<startWC:
        i+=1
        if i == len(t_orb):
            print("---- El eclipse mas largo es el ultimo del periodo de tiempo seleccionado.")
            break
    if i !=0:
        i-=1
    start_run = t_orb[i]
    Epoch =  "\'"+start_run.strftime("%d %b %Y %H:%M:%S.000")+"\'"
    if i == 0:
        print(" ---- El eclipse mas largo es el primero del período de tiempo seleccionado.")

    return startWC, endWC, start_run, Epoch, i, duration_run, step

def plano_max(file, cos1, case = 1, graph = True):
    '''cos = Matriz cuyas columnas son el coseno de cada cara y paneles desplegables (+X, +Y, +Z, -X, -Y, -Z, P1, ...)'''
    cos = cos1.astype(float)
    cos[cos<0] = 0
    cos[np.isnan(cos)]=0
    
    ro = np.linspace(0,90, 2000)
    fi = np.linspace(0,180, 2000)
    
    cos0 = np.sum(cos[:,0]-cos[:,3])
    cos1 = np.sum(cos[:,1]-cos[:,4])
    cos2 = np.sum(cos[:,2]-cos[:,5])
    
    R, F = np.meshgrid(ro,fi)
    
    Z = np.sin(R*np.pi/180)*np.cos(F*np.pi/180)*cos0 + np.sin(R*np.pi/180)*np.sin(F*np.pi/180)*cos1 + np.cos(R*np.pi/180)*cos2
    
    i , j = np.where(Z == np.max(Z))
    r = R[0,:][j[0]]
    f = F[:,0][i[0]]
    x = np.sin(r*np.pi/180)*np.cos(f*np.pi/180)
    y = np.sin(r*np.pi/180)*np.sin(f*np.pi/180)
    z = np.cos(r*np.pi/180)
    
    Z_max = np.max(Z)/len(cos[:,0])
    
    print(f"---- La normal al plano que maximiza el coseno es [{np.round(x,3)}, {np.round(y,3)}, {np.round(z,3)}].")
    if graph:
        plt.figure(f"Plano Max {case}")
        plt.pcolormesh(R, F, Z)
        plt.colorbar()
        plt.xlabel('θ (theta)')
        plt.ylabel('φ (phi)')
        plt.plot(r,f,"*r")
        plt.title(f'Normal del plano de máximo coseno = [{np.round(x,3)}, {np.round(y,3)}, {np.round(z,3)}] - {case}')
        if "planos_optimizacion" not in os.listdir(file):
            os.mkdir(file+"/planos_optimizacion")
        plt.savefig(file+f"/planos_optimizacion/planos-{case}.png")
    return [x,y,z], Z_max

def optim_iter(cos):
    cos = cos.astype(float)
    cos[np.isnan(cos)] = 0
    cos[cos<0]=0
    cos0 = np.sum(cos[:,0]-cos[:,3])
    cos1 = np.sum(cos[:,1]-cos[:,4])
    cos2 = np.sum(cos[:,2]-cos[:,5])
    def sum_cos(v, cos0 = cos0, cos1 = cos1, cos2 = cos2):
        z = -(v[0]*cos0 + v[1]*cos1 + v[2]*cos2)
        return z
    Vx = minimize(sum_cos, [-1,1,0])
    n_max = Vx.x/np.linalg.norm(Vx.x)
    print(f"---- La normal al plano que maximiza el coseno es [{n_max}]")
    return n_max

def add_plano_max(file, cos, axes, paneles, name_plano_max, graph = True, case = 1):
    n_max, Z_max = plano_max(file, cos, graph = graph, case = case)
    axes.append(n_max)
    paneles.append(n_max)
    index_plano = len(axes)-6
    print(f'Panel {index_plano} = Óptimo para el caso {case}.')
    cos_max = coseno_P_sol(n_max, cos[:,:3])
    cos = np.c_[cos, cos_max]
    name_plano_max.append(case)
    return cos, axes, paneles, n_max, Z_max, index_plano, name_plano_max