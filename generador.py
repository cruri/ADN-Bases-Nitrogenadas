"""
generador.py
============
Módulo para la generación sintética de secuencias de ADN
utilizando modelos estocásticos (Cadenas de Markov).
"""

import random
import config
from logger_setup import obtener_logger

logger = obtener_logger(__name__)

# =============================================================================
# ENTRENAMIENTO DEL MODELO
# =============================================================================

def construir_matriz_transicion(conteo_kmers, verbose=True):
    """
    Construye la matriz de transición de Markov a partir de los conteos de k-mers.
    
    Parameters
    ----------
    conteo_kmers : dict
        Diccionario con frecuencias absolutas de k-mers (ej. {'AA': 1849, 'AC': 1641})
    verbose : bool
        Si True, imprime un resumen de la matriz.
        
    Returns
    -------
    dict
        Matriz de probabilidades (ej. {'A': {'A': 0.3, 'C': 0.2, ...}})
    """
    matriz_frecuencias = {}
    
    # 1. Agrupar las frecuencias absolutas
    for kmer, frecuencia in conteo_kmers.items():
        # kmer[:-1] toma todo excepto la última letra (Estado actual)
        # kmer[-1] toma solo la última letra (Siguiente estado)
        estado_actual = kmer[:-1]
        siguiente_estado = kmer[-1]
        
        if estado_actual not in matriz_frecuencias:
            matriz_frecuencias[estado_actual] = {}
            
        matriz_frecuencias[estado_actual][siguiente_estado] = frecuencia

    # 2. Convertir frecuencias a probabilidades condicionadas (de 0 a 1)
    matriz_probabilidades = {}
    
    for estado, transiciones in matriz_frecuencias.items():
        total_transiciones = sum(transiciones.values())
        matriz_probabilidades[estado] = {}
        
        for siguiente, freq in transiciones.items():
            prob = freq / total_transiciones
            matriz_probabilidades[estado][siguiente] = prob
            
    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print("MATRIZ DE TRANSICIÓN DE MARKOV CONSTRUIDA")
        print("="*60)
        for estado in sorted(matriz_probabilidades.keys()):
            # Formatear la salida para que sea fácil de leer
            probs_str = ", ".join([f"{k}: {v:.3f}" for k, v in matriz_probabilidades[estado].items()])
            print(f"Estado '{estado}' -> [{probs_str}]")
        print("="*60 + "\n")
        
    logger.info("Matriz de transición construida exitosamente.")
    return matriz_probabilidades


# =============================================================================
# GENERACIÓN SINTÉTICA
# =============================================================================

def generar_cadena_markov(matriz, longitud_deseada, semilla=None):
    """
    Genera una secuencia sintética "caminando" por la matriz de transición.
    
    Parameters
    ----------
    matriz : dict
        La matriz de probabilidades generada previamente.
    longitud_deseada : int
        Número total de bases que tendrá la secuencia artificial.
    semilla : str, optional
        Secuencia inicial para forzar el arranque. Si es None, arranca al azar.
        
    Returns
    -------
    list
        Lista con las bases generadas artificialmente.
    """
    cadena_sintetica = []
    estados_posibles = list(matriz.keys())
    
    # 1. Elegir el punto de partida
    if semilla and semilla in estados_posibles:
        estado_actual = semilla
    else:
        # Arranque aleatorio
        estado_actual = random.choice(estados_posibles)
        
    cadena_sintetica.append(estado_actual)
    longitud_estado = len(estado_actual)
    
    # 2. Generar el resto de la cadena
    while len("".join(cadena_sintetica)) < longitud_deseada:
        
        if estado_actual in matriz:
            # Extraer las opciones y las probabilidades (el "dado cargado")
            opciones = list(matriz[estado_actual].keys())
            probabilidades = list(matriz[estado_actual].values())
            
            # Lanzar el dado usando los pesos estadísticos
            siguiente_estado = random.choices(opciones, weights=probabilidades, k=1)[0]
            cadena_sintetica.append(siguiente_estado)
            
            # Actualizar el estado para el siguiente paso (memoria deslizante)
            estado_actual = ("".join(cadena_sintetica))[-longitud_estado:]
            
        else:
            # Mecanismo de seguridad: Si llegamos a un estado sin salida
            logger.warning(f"Callejón sin salida en el estado '{estado_actual}'. Reiniciando aleatoriamente.")
            estado_actual = random.choice(estados_posibles)

    # 3. Limpiar y devolver la lista del tamaño exacto
    cadena_final = list("".join(cadena_sintetica)[:longitud_deseada])
    logger.info(f"Cadena sintética generada. Longitud final: {len(cadena_final)} bases.")
    
    return cadena_final

def guardar_secuencia_sintetica(cadena, nombre_archivo="secuencia_sintetica.fasta"):
    """
    Guarda la secuencia generada en un archivo físico, formateada en bloques.
    
    Parameters
    ----------
    cadena : list o str
        La secuencia sintética generada.
    nombre_archivo : str
        Nombre del archivo de salida.
    """
    # Si viene como lista, la unimos en un solo string
    if isinstance(cadena, list):
        cadena = "".join(cadena)
        
    ruta_salida = config.OUTPUT_DIR / nombre_archivo
    
    try:
        with open(ruta_salida, 'w', encoding='utf-8') as file:
            # Escribir el encabezado típico de un archivo FASTA
            file.write(">Secuencia_Sintetica_Markov_Orden1\n")
            
            # Escribir la secuencia en bloques de 80 caracteres por línea (estándar bioinformático)
            bloque = 80
            for i in range(0, len(cadena), bloque):
                file.write(cadena[i : i+bloque] + "\n")
                
        logger.info(f"Secuencia sintética guardada exitosamente en: {ruta_salida}")
        return True
    
    except IOError as e:
        logger.error(f"Error al guardar la secuencia sintética: {e}")
        return False