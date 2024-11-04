from bacteria import bacteria
from chemiotaxis import chemiotaxis

import numpy

poblacion = []
path = "C:\secuenciasBFOA\multifasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1                                              #numero de gaps a insertar 
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)                           #mejor bacteria   
tempBacteria = bacteria(path)                        #bacteria temporal para validaciones
original = bacteria(path)                           #bacteria original sin gaps
globalNFE = 0                                       #numero de evaluaciones de la funcion objetivo

dAttr= 0.1 #0.1
wAttr= 0.2 #0.2
hRep=dAttr
wRep= 10    #10


def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return
      

for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))


for iteracion in range(iteraciones):
    # Ajuste dinámico de parámetros de atracción y repulsión
    factor = iteracion / iteraciones  # Factor de progresión (entre 0 y 1)
    
    dAttr_iter = dAttr * (1 - factor)    # Disminuir atracción con iteraciones
    wAttr_iter = wAttr * (1 - factor)
    hRep_iter = hRep * (1 + factor)      # Aumentar repulsión con iteraciones
    wRep_iter = wRep * (1 + factor)

    # Para cada bacteria en la población, realizar el tumbo-nado y evaluación
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        bacteria.autoEvalua()

    # Aplicar quimiotaxis con los parámetros ajustados
    chemio.doChemioTaxis(poblacion, dAttr_iter, wAttr_iter, hRep_iter, wRep_iter)

    # Actualizar el conteo global de evaluaciones de la función objetivo
    globalNFE += chemio.parcialNFE
    best = max(poblacion, key=lambda x: x.fitness)
    
    # Si se encuentra una mejor bacteria, clonarla
    if (veryBest == None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)

    print("Iteración: ", iteracion, "Interacción: ", veryBest.interaction, "Fitness: ", veryBest.fitness, " NFE:", globalNFE)
    
    # Ajustar el número de bacterias eliminadas dinámicamente
    eliminar_porcentaje = 0.5 - (0.3 * factor)  # Disminuir el porcentaje de eliminación con el tiempo
    chemio.eliminarClonar(path, poblacion)
    
    # Insertar bacterias aleatorias para mejorar la diversidad
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
    
    print("Población: ", len(poblacion))


veryBest.showGenome()
validaSecuencias(path, veryBest)