# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:26:40 2022

@author: Ignacio Fernandez Mora y Jorge Ruiz

"""

import numpy as np
import matplotlib.pyplot as plt


def AGnormal(cantidad, numeroDeCiudades):
    ''' Esta función crea una población inicial con individuos con un genoma aleatorio
    cantidad = Cantidad de individuos para la población
    numeroDeCiudades = cantidad de ciudades 
    return
    Array con la población generada
    '''
    
    #Se genera un arreglo con los índices de las ciudades, correpondiente a cada una de ellas
    indiceCiudades = np.linspace(0, numeroDeCiudades-1,numeroDeCiudades).astype(int)
    
    poblacionInicial = []
    
    #Se inicia un ciclo para generar un individuo con un genoma aleatorio en cada iteración
    i = 0
    while i < cantidad:
        #Se le añade a la población inicial una permutación aleatoria del arreglo indiceCiudades
        poblacionInicial.append(np.random.permutation(indiceCiudades).tolist())
        i += 1
    return np.array(poblacionInicial)

def Optimizacion(cromosoma, CoordenadasCiudades):
    ''' Calcula el valor de ajuste en cada individuo
    cromosoma = Cromosoma del indivuo
    CoordenadasCiudades =  Coordenas de la ciudades generadas
    return
    Se calcula el inverso de la distancia
    '''
    
    #Variable distancia, donde se guarda el acumulado de las distancias euclidianas entre nodos adyacentes
    distancia = 0
    
    #se comienza el ciclo desde -1 para obtener el último
    for i in range(-1,len(cromosoma)-1):
        #Se calcula la distancia euclidiana y se suma a la variable distancia
        distancia += np.sqrt(np.sum(np.power(CoordenadasCiudades[cromosoma[i]]-CoordenadasCiudades[cromosoma[i+1]],2)))

    return 1/distancia


def Mutacion(cromosoma, probabilidad):
    ''' Función que muta el cromosoma con cierta probabilidad
    cromosoma = Cromosoma del indivuo
    Probabilidad =  Probabilidad en una mutación
    return
    Cromosoma mutado si cumple
    '''
    #Se evalua si el dato aleatorio es menor a la probabilidad y se realiza la mutación 
    if np.random.random() < probabilidad:
        x1 = 0
        x2 = 0
        #En el caso de que los índices sean iguales, se continúa generando índices aleatorios hasta que sean diferentes
        while x1 == x2:
            x1, x2 = np.random.randint(0, len(cromosoma),2)

        #Se intercambian de lugar los elementos en los índices aleatorios obtenidos
        cromosoma[x1], cromosoma[x2] = cromosoma[x2], cromosoma[x1]
    
    return cromosoma


def AGModificada(cantidad, CoordenadasCiudades):
    ''' Función generadora de una población con un individuo ventajoso
    cantidad = cantidad de individuos por población
    CoordenadasCiudades =  Coordenas de la ciudades generadas
    return
    Array con la población generada
    '''
    
    #Se ingresa a una variable la longitud del arrglo de las coordenadas de las ciudades y se inicia la lista para contener la población inicial
    numeroDeCiudades = len(CoordenadasCiudades)
    poblacionInicial = []
    
    #Se inicia el ciclo para crear los individuos de la población
    i = 0
    while i < cantidad:
        
        #Se genera un arreglo con los índices de las ciudades 
        indiceCiudades = np.linspace(0, numeroDeCiudades-1,numeroDeCiudades).astype(int).tolist()
        nodoInicial = indiceCiudades.pop(np.random.randint(0,numeroDeCiudades))
        
        individuoActual = []
        individuoActual.append(nodoInicial)
        
        #Se define la variable nodoPasado para guardar el valor del último índice ingresado para calcular en la siguiente iteración su ciudad más cercano
        nodoPasado = nodoInicial
        while len(indiceCiudades) > 0:
            #Se ordenan los indices dependiendo de la distancia entre las ciudades en orden ascendente
            indiceCiudades.sort(key=lambda x: np.sqrt(np.sum(np.power(CoordenadasCiudades[x]-CoordenadasCiudades[nodoPasado],2))), reverse=False)
            nodoPasado = indiceCiudades.pop(0)
            individuoActual.append(nodoPasado)
        
        #Si la población no es la primera, se realizan entre 3 y 10 mutaciones (aleatoriamente) en el individuo
        if i != 0:
            for j in range(0,np.random.randint(3,11)):
                #Se realiza el mismo procedimiento que en mutación
                x1 = 0
                x2 = 0
                while x1 == x2:
                    x1, x2 = np.random.randint(0, len(individuoActual),2)

                individuoActual[x1], individuoActual[x2] = individuoActual[x2], individuoActual[x1]
        
        poblacionInicial.append(individuoActual)
        i += 1
    return np.array(poblacionInicial)

#Se procede con el inicio del programa

def calculoCordenadas(numeroCiudad):
    ''' Se hace el calculo de las coordenadas para las ciudades'''
    cordenadaX = 0.1 * ((9+13 * numeroCiudad**2) % 200 )
    cordenadaY = 0.1 * ((7 + 1327 * numeroCiudad) % 200)
    resultado = [cordenadaX, cordenadaY]
    return resultado

def cordenadasGlobales(numeroCiudades):
    ''' Crea un arreglo con las coordenadas de las ciudades'''
    conjuntoCordenadas = []
    for ciudad in range(numeroCiudades):
        conjuntoCordenadas += [calculoCordenadas(ciudad)]
    return conjuntoCordenadas

CoordenadasCiudades = np.array(cordenadasGlobales(100))


#Se definien los parámetros para la simulación
Poblacion = 100
probabilidadMutacion = 0.2
Iteraciones = 300

#Se crea una lista para guardar los mejores individuos 
# En resultados se guardan los mejores individuos y los promedio de la función de ajuste de cada iteración
mejorCamino = []
promedioAjuste = []


#Se define la población inicial
#poblacion = AGnormal(Poblacion,len(CoordenadasCiudades)) #Para la simulación normal
poblacion = AGModificada(Poblacion,CoordenadasCiudades) #Para la simulación con un individuo con ventaja


n = Iteraciones
while n >= 0:
    #Se define la variable para acumular la suma de los valores de ajuste para calcular su promedio
    caminos = poblacion.tolist()
    promedioGeneracional = 0
    #Se itera sobre la población para calcular el valor de ajuste de cada individuo y calcular el promedio generacional
    for i in range(0,len(caminos)):
        #Se guarda en la variable camino el camino en el índice actual y su valor de ajuste y se suma el valor de ajuste al acumulador
        caminos[i] = (caminos[i],Optimizacion(caminos[i],CoordenadasCiudades))
        promedioGeneracional += caminos[i][1]
    
    #Se divide la suma total de valores de ajuste entre el total de individuos para calcular el promedio generacional
    promedioGeneracional /= len(caminos)
    
    #Se ordenan los individuos por su función de ajuste    
    caminos.sort(key=lambda x: x[1], reverse=True)

    #Se guarda el valor de ajuste del mejor individuo y se guarda el promedio generacional en la lista de resultados
    mejorValorAjuste = caminos[0][1]
    promedioAjuste.append([promedioGeneracional,mejorValorAjuste])
    
    #Si no ha habido un mejor camino, solamente se añade el mejor camino actual a la lista
    if mejorCamino == []:
        mejorCamino.append((caminos[0][0], mejorValorAjuste))
    elif mejorValorAjuste > mejorCamino[len(mejorCamino)-1][1]:
        mejorCamino.append((caminos[0][0], mejorValorAjuste))
    
    #Se aplican las mutaciones a cada individuo de la población, dependiendo de la probabilidad
    for i in range(0,len(poblacion)):
        poblacion[i] = Mutacion(poblacion[i],probabilidadMutacion)
    
    n -= 1

resultados = np.array(promedioAjuste)

#Se mapean los índices en el cromosoma de los mejores caminos a las respectivas coordenadas de cada ciudad
coordenadasMejorCamino = []
for camino in mejorCamino:
    caminoAlter = camino[0].copy()
    caminoAlter.append(caminoAlter[0])
    coordenadasMejorCamino.append(CoordenadasCiudades[caminoAlter])
#Se muestra al usuario la distancia recorrida
print('Distancia recorrida fue de {0:.2f}'.format(1/mejorCamino[len(mejorCamino)-1][1]))


#Se grafican los resultados
global ax
fig, ax = plt.subplots()
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Resultados de la Simulación')

#Se muestra como primero de los mejores caminos encontrados (el peor de los mejores) como gráfico inicial
ax.plot(coordenadasMejorCamino[0][:,0], coordenadasMejorCamino[0][:,1])
ax.scatter(CoordenadasCiudades[:,0],CoordenadasCiudades[:,1],s=50,color='red')
ax.scatter(coordenadasMejorCamino[0][0,0],coordenadasMejorCamino[0][0,1],s=50,color='green')



#Gráficas de los valores de ajuste promedio y mejor por generación
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Generación')
ax1.set_ylabel('Valor promedio de la función de optimización')
ax1.plot(np.linspace(0, Iteraciones,Iteraciones+1).astype(int),resultados[:,0])
ax1.set_title('Resultados de la Simulación por Generación')

fig2, ax2 = plt.subplots()
ax2.set_xlabel('Generación')
ax2.set_ylabel('Mejor valor de la función de optimización')
ax2.plot(np.linspace(0, Iteraciones,Iteraciones+1).astype(int),resultados[:,1])
ax2.set_title('Resultados del Mejor Individuo por Generación')

#Se muestran las gráficas
plt.show()
