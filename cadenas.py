"""
cadenas.py
==========
Módulo para procesamiento de secuencias de ADN/ARN.

Funcionalidades:
- Lectura de archivos de texto
- Limpieza de secuencias
- Generación de caminatas aleatorias
- Regresión lineal
- Rotación de datos
"""

import string
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

import config
from logger_setup import obtener_logger

logger = obtener_logger(__name__)


# =============================================================================
# LECTURA Y VALIDACIÓN DE ARCHIVOS
# =============================================================================

def leer_archivo(ruta_archivo):
    """
    Lee un archivo de texto de forma segura.
    
    Parameters
    ----------
    ruta_archivo : str o Path
        Ruta al archivo a leer
    
    Returns
    -------
    str
        Contenido del archivo
    
    Raises
    ------
    FileNotFoundError
        Si el archivo no existe
    IOError
        Si hay error al leer el archivo
    """
    try:
        logger.info(f"Leyendo archivo: {ruta_archivo}")
        with open(ruta_archivo, 'r', encoding='utf-8') as file:
            contenido = file.read()
        logger.info(f"Archivo leído exitosamente. Tamaño: {len(contenido)} caracteres")
        return contenido
    
    except FileNotFoundError:
        logger.error(f"Archivo no encontrado: {ruta_archivo}")
        raise FileNotFoundError(f"No se encontró el archivo: {ruta_archivo}")
    
    except IOError as e:
        logger.error(f"Error al leer archivo {ruta_archivo}: {e}")
        raise IOError(f"No se pudo leer el archivo: {ruta_archivo}")


def validar_extension(ruta_archivo):
    """
    Valida que la extensión del archivo sea soportada.
    
    Parameters
    ----------
    ruta_archivo : str o Path
        Ruta al archivo
    
    Returns
    -------
    bool
        True si la extensión es válida
    
    Raises
    ------
    ValueError
        Si la extensión no es soportada
    """
    from pathlib import Path
    extension = Path(ruta_archivo).suffix.lower() # Se extrae la extensión del archivo
    
    if extension not in config.EXTENSIONES_SOPORTADAS:
        raise ValueError(
            f"Extensión no soportada: {extension}. "
            f"Use: {config.EXTENSIONES_SOPORTADAS}"
        )
    return True


# =============================================================================
# LIMPIEZA DE SECUENCIAS
# =============================================================================

def limpieza(cadena_sucia):
    """
    Limpia una cadena removiendo caracteres no deseados.
    
    Mantiene solo las bases nitrogenadas válidas (A, C, G, T).
    
    Parameters
    ----------
    cadena_sucia : str
        Cadena con posibles caracteres inválidos
    
    Returns
    -------
    list
        Lista con solo las bases válidas
    
    Raises
    ------
    ValueError
        Si la cadena limpia está vacía
    """
    # Convertir a mayúsculas para normalizar
    cadena_sucia = cadena_sucia.upper()
    
    # Crear tabla de traducción para eliminar caracteres no deseados
    tabla_traduccion = str.maketrans("", "", config.CARACTERES_A_ELIMINAR)
    
    # Aplicar limpieza
    datos_limpios = cadena_sucia.translate(tabla_traduccion)
    
    # Convertir a lista de caracteres
    cadena_limpia = list(datos_limpios)
    
    if not cadena_limpia:
        logger.error("La cadena limpia está vacía después de la limpieza")
        raise ValueError("La cadena limpia está vacía. Verifica el archivo de entrada.")
    
    # Validar que solo contenga bases válidas
    bases_invalidas = set(cadena_limpia) - config.BASES_VALIDAS
    if bases_invalidas:
        logger.warning(f"Se encontraron bases inválidas: {bases_invalidas}")
        logger.warning("Se eliminarán las bases inválidas")
        cadena_limpia = [base for base in cadena_limpia if base in config.BASES_VALIDAS]
    
    if config.VERBOSE:
        print("\n" + "="*60)
        print("CADENA LIMPIA")
        print("="*60)
        print(f"Longitud: {len(cadena_limpia)} bases")
        print(f"Composición: A={cadena_limpia.count('A')}, C={cadena_limpia.count('C')}, "
              f"G={cadena_limpia.count('G')}, T={cadena_limpia.count('T')}")
        print(f"Primeros 100 caracteres: {''.join(cadena_limpia[:100])}")
        print("="*60 + "\n")
    
    logger.info(f"Cadena limpia generada: {len(cadena_limpia)} bases válidas")
    return cadena_limpia


# =============================================================================
# CAMINATA ALEATORIA
# =============================================================================

def caminata_aleatoria(cadena_limpia, verbose=True):
    """
    Genera una caminata aleatoria basada en las bases nitrogenadas.
    
    Regla:
    - A o T: -1 (baja)
    - C o G: +1 (sube)
    
    Parameters
    ----------
    cadena_limpia : list
        Lista con las bases válidas
    verbose : bool
        Si True, muestra información detallada
    
    Returns
    -------
    np.ndarray
        Array con la posición en cada paso
    """
    if not cadena_limpia:
        raise ValueError("La cadena no puede estar vacía")
    
    caminata = [0]  # Posición inicial
    posicion = 0
    
    for base in cadena_limpia:
        if base in ['A', 'T']:
            posicion -= 1
        elif base in ['C', 'G']:
            posicion += 1
        else:
            logger.warning(f"Base desconocida encontrada: {base}")
            continue
        
        caminata.append(posicion)
    
    caminata_array = np.array(caminata)
    
    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print("CAMINATA ALEATORIA")
        print("="*60)
        print(f"Longitud: {len(caminata_array)} pasos")
        print(f"Posición inicial: {caminata_array[0]}")
        print(f"Posición final: {caminata_array[-1]}")
        print(f"Mínimo: {caminata_array.min()}")
        print(f"Máximo: {caminata_array.max()}")
        print(f"Media: {caminata_array.mean():.4f}")
        print(f"Desviación estándar: {caminata_array.std():.4f}")
        print("="*60 + "\n")
    
    logger.info(f"Caminata aleatoria generada: {len(caminata_array)} pasos")
    return caminata_array


# =============================================================================
# REGRESIÓN LINEAL
# =============================================================================

def regresion_lineal(x, y, verbose=True):
    """
    Realiza una regresión lineal sobre la caminata aleatoria.
    
    Parameters
    ----------
    x : np.ndarray
        Índices de los pasos (1D o 2D)
    y : np.ndarray
        Posiciones de la caminata
    verbose : bool
        Si True, muestra estadísticas
    
    Returns
    -------
    tuple
        (modelo_entrenado, pendiente)
    
    Raises
    ------
    ValueError
        Si los datos son inválidos
    """
    # Validar datos
    if x.size == 0 or y.size == 0:
        raise ValueError("Los datos no pueden estar vacíos")
    
    if len(x.flatten()) != len(y):
        raise ValueError(f"Dimensiones incompatibles: x={len(x.flatten())}, y={len(y)}")
    
    x_flat = x.flatten()
    
    # Análisis de correlación
    corr_pearson, p_value = pearsonr(x_flat, y)
    
    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        x, y,
        train_size=config.TRAIN_SIZE,
        random_state=config.RANDOM_STATE,
        shuffle=config.SHUFFLE
    )
    
    # Crear y entrenar modelo
    modelo = LinearRegression(fit_intercept=False)
    modelo.fit(X_train, y_train)
    
    # Calcular métricas
    r2_score = modelo.score(x, y)
    pendiente = modelo.coef_[0]
    
    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print("ANÁLISIS DE REGRESIÓN LINEAL")
        print("="*60)
        print(f"Coeficiente de correlación Pearson: {corr_pearson:.6f}")
        print(f"P-value: {p_value:.6e}")
        print(f"Pendiente (coeficiente): {pendiente:.6f}")
        print(f"Coeficiente de determinación R²: {r2_score:.6f}")
        print(f"Tamaño de entrenamiento: {len(X_train)}")
        print(f"Tamaño de test: {len(X_test)}")
        print("="*60 + "\n")
    
    logger.info(f"Regresión lineal completada. Pendiente: {pendiente:.6f}, R²: {r2_score:.6f}")
    
    return modelo, pendiente


# =============================================================================
# ROTACIÓN DE DATOS
# =============================================================================

def rotacion(x, y, theta, usar_radianes=True, verbose=True):
    """
    Aplica una rotación en 2D y analiza los datos rotados.
    
    Parameters
    ----------
    x : np.ndarray
        Índices (1D o 2D)
    y : np.ndarray
        Posiciones de caminata
    theta : float
        Ángulo de rotación (en radianes o grados)
    usar_radianes : bool
        Si True, theta está en radianes; si False, en grados
    verbose : bool
        Si True, muestra información detallada
    
    Returns
    -------
    dict
        Diccionario con datos rotados y modelo
    """
    if not usar_radianes:
        theta = np.radians(theta)
    
    x_flat = x.flatten()
    
    # Matriz de rotación
    matriz_rotacion = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])
    
    # Combinar datos
    datos_originales = np.vstack([x_flat, y])
    
    # Aplicar rotación
    datos_rotados = matriz_rotacion @ datos_originales
    
    x_rotado = datos_rotados[0, :]
    y_rotado = datos_rotados[1, :]
    
    # Regresión en datos rotados
    x_rotado_2d = x_rotado.reshape(-1, 1)
    modelo_rotado = LinearRegression(fit_intercept=False)
    modelo_rotado.fit(x_rotado_2d, y_rotado)
    
    pendiente_rotada = modelo_rotado.coef_[0]
    r2_rotado = modelo_rotado.score(x_rotado_2d, y_rotado)
    
    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print("ANÁLISIS DE DATOS ROTADOS")
        print("="*60)
        print(f"Ángulo de rotación: {theta:.4f} radianes ({np.degrees(theta):.2f}°)")
        print(f"Pendiente rotada: {pendiente_rotada:.6f}")
        print(f"R² rotado: {r2_rotado:.6f}")
        print("="*60 + "\n")
    
    logger.info(f"Rotación completada con ángulo {np.degrees(theta):.2f}°")
    
    return {
        'x_rotado': x_rotado,
        'y_rotado': y_rotado,
        'modelo': modelo_rotado,
        'pendiente': pendiente_rotada,
        'r2': r2_rotado,
        'theta': theta
    }


# =============================================================================
# VISUALIZACIÓN
# =============================================================================

def graficar_caminata(x, y, y_pred, pendiente, nombre_archivo=None):
    """
    Grafica la caminata aleatoria con su regresión lineal.
    
    Parameters
    ----------
    x : np.ndarray
        Índices
    y : np.ndarray
        Posiciones
    y_pred : np.ndarray
        Predicciones del modelo
    pendiente : float
        Valor de la pendiente
    nombre_archivo : str, optional
        Si se proporciona, guarda la gráfica
    """
    plt.figure(figsize=config.FIGSIZE_CAMINATA, dpi=config.DPI)
    plt.plot(x, y, label='Caminata Aleatoria', linewidth=1.5, alpha=0.8)
    plt.plot(x, y_pred, label=f'Regresión Lineal (m={pendiente:.6f})', 
             linestyle='--', linewidth=2, color='red')
    plt.xlabel('Número de Elemento en la Cadena', fontsize=12)
    plt.ylabel('Posición (+1/-1)', fontsize=12)
    plt.title('Caminata Aleatoria de Secuencia de ADN', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if nombre_archivo and config.GUARDAR_GRAFICAS:
        ruta_salida = config.DIRECTORIO_GRAFICAS / nombre_archivo
        plt.savefig(ruta_salida, dpi=config.DPI, bbox_inches='tight')
        logger.info(f"Gráfica guardada: {ruta_salida}")
    
    plt.show()


def graficar_rotacion(x_rotado, y_rotado, y_pred_rotado, theta, nombre_archivo=None):
    """
    Grafica los datos rotados con la regresión lineal.
    
    Parameters
    ----------
    x_rotado : np.ndarray
        Datos X rotados
    y_rotado : np.ndarray
        Datos Y rotados
    y_pred_rotado : np.ndarray
        Predicciones rotadas
    theta : float
        Ángulo de rotación (en radianes)
    nombre_archivo : str, optional
        Si se proporciona, guarda la gráfica
    """
    plt.figure(figsize=config.FIGSIZE_ROTACION, dpi=config.DPI)
    plt.scatter(x_rotado, y_rotado, label='Datos Rotados', alpha=0.6, s=30)
    plt.plot(x_rotado, y_pred_rotado, label='Regresión Lineal', 
             linestyle='-', linewidth=2, color='red')
    plt.xlabel('X Rotado', fontsize=12)
    plt.ylabel('Y Rotado', fontsize=12)
    plt.title(f'Caminata Aleatoria Rotada ({np.degrees(theta):.2f}°)', 
              fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if nombre_archivo and config.GUARDAR_GRAFICAS:
        ruta_salida = config.DIRECTORIO_GRAFICAS / nombre_archivo
        plt.savefig(ruta_salida, dpi=config.DPI, bbox_inches='tight')
        logger.info(f"Gráfica guardada: {ruta_salida}")
    
    plt.show()

# =============================================================================
# ÁNALISIS K-MERS
# =============================================================================

def analizar_kmers(cadena_limpia, k=2, verbose=True):
    """
    Calcula la frecuencia y probabilidad de k-mers (subsecuencias) en la cadena.
    
    Parameters
    ----------
    cadena_limpia : list
        Lista con las bases válidas (A, C, G, T).
    k : int
        Tamaño del k-mer (2 para dímeros, 3 para codones, etc.).
    verbose : bool
        Si True, imprime las estadísticas en consola.
        
    Returns
    -------
    tuple
        (conteo_kmers, probabilidades)
    """
    if not cadena_limpia or len(cadena_limpia) < k:
        raise ValueError(f"La cadena es demasiado corta para analizar k-mers de tamaño {k}")
    
    # Diccionario
    conteo_kmers = {}
    
    N = len(cadena_limpia)
    total_ventanas = N - k + 1

    for i in range(total_ventanas):
        ventana = cadena_limpia[i:i+k]
        # Convertimos a cadena (str)
        kmer = "".join(ventana)
        conteo_kmers[kmer] = conteo_kmers.get(kmer,0) + 1

    probabilidades = {}
    for kmer, conteo in conteo_kmers.items():
        probabilidades[kmer] = conteo / total_ventanas # Frecuencia Relativa

    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print(f"ANÁLISIS DE K-MERS (k={k})")
        print("="*60)
        print(f"Total de subsecuencias encontradas: {total_ventanas}")
        print(f"K-mers únicos identificados: {len(conteo_kmers)}")
        print("\nTop 5 k-mers más frecuentes: ")

        kmers_ordenados = sorted(probabilidades.items(), key=lambda item: item[1], reverse=True)

        for kmer, prob in kmers_ordenados[:5]:
            conteo_real = conteo_kmers[kmer]
            print(f"- {kmer}: {conteo_real} apariciones (Prob: {prob:.4f})")
        print("\n" + "="*60)

    logger.info(f"Análisis de {k}-mers completado. {len(conteo_kmers)} combinaciones únicas encontradas.")

    return conteo_kmers, probabilidades


# =============================================================================
# ÁNALISIS FOURIER
# =============================================================================

def analizar_fourier(y_rotado, verbose=True):
    """
    Calcula la Transformada Rápida de Fourier (FFT) y el Espectro de Potencia
    de la caminata aleatoria rotada.
    
    Parameters
    ----------
    y_rotado : np.ndarray
        Arreglo numérico con las posiciones de la caminata sin tendencia.
    verbose : bool
        Si True, muestra la frecuencia dominante en consola.
        
    Returns
    -------
    tuple
        (frecuencias_filtradas, espectro_potencia_filtrado)
    """
    N = len(y_rotado)
    if N == 0:
        raise ValueError("El arreglo de datos rotados no puede estar vacío.")
    
    fft_valores = np.fft.fft(y_rotado)
    espectro_potencia = np.abs(fft_valores) ** 2

    frecuencias = np.fft.fftfreq(N, d=1.0)
    
    mitad = N // 2
    frecuencias_filtradas = frecuencias[:mitad]
    espectro_potencia_filtrado = espectro_potencia[:mitad]

    if verbose and config.VERBOSE:
        # Excluimos el índice 0 (frecuencia cero / componente DC) para buscar el pico real
        idx_pico = np.argmax(espectro_potencia_filtrado[1:]) + 1
        frec_dominante = frecuencias_filtradas[idx_pico]
        potencia_max = espectro_potencia_filtrado[idx_pico]
        
        # Calcular el periodo (Periodo = 1 / Frecuencia)
        periodo = 1 / frec_dominante if frec_dominante != 0 else float('inf')
        
        print("\n" + "="*60)
        print("ANÁLISIS DE FRECUENCIAS (FOURIER)")
        print("="*60)
        print(f"Frecuencia dominante: {frec_dominante:.4f}")
        print(f"Periodo equivalente: {periodo:.2f} bases")
        print(f"Potencia del pico: {potencia_max:.2f}")
        print("="*60 + "\n")
        
    logger.info(f"Análisis de Fourier completado. Mitad de espectro: {len(frecuencias_filtradas)} puntos.")
    
    return frecuencias_filtradas, espectro_potencia_filtrado


def graficar_espectro(frecuencias, espectro, nombre_archivo=None):
    """
    Genera el gráfico del Espectro de Potencia (Periodograma) y lo guarda en disco
    sin mostrarlo en pantalla.
    
    Parameters
    ----------
    frecuencias : np.ndarray
        Arreglo con las frecuencias calculadas.
    espectro : np.ndarray
        Arreglo con la potencia correspondiente a cada frecuencia.
    nombre_archivo : str, optional
        Nombre del archivo con el que se guardará la gráfica.
    """
    # Usamos la configuración de tamaño de la rotación o una personalizada
    plt.figure(figsize=config.FIGSIZE_ROTACION, dpi=config.DPI)
    
    # Graficamos la señal en el dominio de la frecuencia
    plt.plot(frecuencias, espectro, color='darkviolet', linewidth=1.5, label='Espectro de Potencia')
    
    # Configuraciones estéticas de la gráfica
    plt.xlabel('Frecuencia ($f$)', fontsize=12)
    plt.ylabel('Potencia ($|X(f)|^2$)', fontsize=12)
    plt.title('Espectro de Potencia de la Secuencia de ADN (FFT)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    # Guardar automáticamente en la ruta preconfigurada
    if nombre_archivo and config.GUARDAR_GRAFICAS:
        ruta_salida = config.DIRECTORIO_GRAFICAS / nombre_archivo
        plt.savefig(ruta_salida, dpi=config.DPI, bbox_inches='tight')
        logger.info(f"Gráfica de Fourier guardada en: {ruta_salida}")
    
    # Cerrar la figura explícitamente para no consumir memoria ni mostrar ventana
    plt.close()


def calcular_hurst_dfa(caminata, min_ventana=10, max_ventanas_divisor=4, num_puntos=20, verbose=True):
    """
    Calcula el Exponente de Hurst usando Análisis de Fluctuaciones sin Tendencia (DFA).
    
    Parameters
    ----------
    caminata : np.ndarray
        Arreglo numérico con la caminata aleatoria.
    min_ventana : int
        Tamaño mínimo de la ventana (caja).
    max_ventanas_divisor : int
        Fracción del tamaño total para la ventana máxima (ej. N/4).
    num_puntos : int
        Cantidad de tamaños de ventana a evaluar.
    verbose : bool
        Si True, imprime el resultado y la interpretación en consola.
        
    Returns
    -------
    tuple
        (Hurst, log_n, log_F, modelo_hurst)
    """
    N = len(caminata)
    if N < min_ventana * 2:
        raise ValueError("La caminata es demasiado corta para el análisis DFA.")

    # 1. Generar tamaños de ventana espaciados logarítmicamente
    max_ventana = N // max_ventanas_divisor
    # Usamos np.unique para evitar tamaños duplicados al redondear a enteros
    ventanas = np.unique(
        np.logspace(np.log10(min_ventana), np.log10(max_ventana), num_puntos).astype(int)
    )
    
    fluctuaciones = []
    
    # 2. Iterar sobre cada tamaño de caja (n)
    for n in ventanas:
        num_ventanas = N // n
        
        # Recortar la caminata para que sea un múltiplo exacto de 'n'
        caminata_recortada = caminata[:num_ventanas * n]
        
        # Dividir eficientemente en fragmentos usando reshape
        fragmentos = caminata_recortada.reshape(num_ventanas, n)
        
        x_local = np.arange(n)
        suma_varianzas = 0
        
        # 3. Eliminar tendencia local en cada fragmento
        for fragmento in fragmentos:
            # np.polyfit ajusta un polinomio (grado 1 = línea recta)
            coeficientes = np.polyfit(x_local, fragmento, 1)
            tendencia = np.polyval(coeficientes, x_local)
            
            # Restar la tendencia y calcular la varianza de los residuos
            residuos = fragmento - tendencia
            suma_varianzas += np.var(residuos)
            
        # 4. Calcular la fluctuación cuadrática media (F_n)
        promedio = suma_varianzas / num_ventanas
        fluctuaciones.append(np.sqrt(promedio))
        
    # 5. Transformación log-log para la ley de potencia
    log_n = np.log10(ventanas)
    log_F = np.log10(fluctuaciones)
    
    # 6. Regresión lineal global para obtener la pendiente (Exponente de Hurst)
    log_n_2d = log_n.reshape(-1, 1)
    modelo_hurst = LinearRegression(fit_intercept=True)
    modelo_hurst.fit(log_n_2d, log_F)
    
    H = modelo_hurst.coef_[0]
    r2 = modelo_hurst.score(log_n_2d, log_F)
    
    # 7. Reporte e Interpretación
    if verbose and config.VERBOSE:
        print("\n" + "="*60)
        print("ANÁLISIS DE EXPONENTE DE HURST (MÉTODO DFA)")
        print("="*60)
        print(f"Exponente de Hurst (H): {H:.4f}")
        print(f"Coeficiente R² de la ley de potencia: {r2:.4f}")
        print("-" * 60)
        
        # Interpretación física
        if H < 0.45:
            estado = "Antipersistente (Reversión constante a la media)"
        elif 0.45 <= H <= 0.55:
            estado = "Caminata Aleatoria Pura (Movimiento Browniano, sin memoria)"
        else:
            estado = "Persistente (Correlaciones positivas a largo plazo, fractal)"
            
        print(f"Comportamiento detectado: {estado}")
        print("="*60 + "\n")
        
    logger.info(f"Análisis DFA completado. Exponente de Hurst: {H:.4f}")
    
    return H, log_n, log_F, modelo_hurst

def graficar_hurst(log_n, log_F, modelo_hurst, H, nombre_archivo=None):
    """
    Genera el gráfico log-log de la fluctuación vs tamaño de ventana 
    para ilustrar el Exponente de Hurst.
    """
    plt.figure(figsize=config.FIGSIZE_ROTACION, dpi=config.DPI)
    
    # Puntos reales medidos por el algoritmo
    plt.scatter(log_n, log_F, color='teal', label='Fluctuaciones DFA ($F(n)$)', alpha=0.7, s=40)
    
    # Línea de regresión (La pendiente es H)
    y_pred = modelo_hurst.predict(log_n.reshape(-1, 1))
    plt.plot(log_n, y_pred, color='red', linestyle='--', linewidth=2.5, 
             label=f'Ajuste de Ley de Potencia (Pendiente $H = {H:.4f}$)')
    
    # Estética
    plt.xlabel('$\log_{10}(n)$ (Tamaño de ventana)', fontsize=12)
    plt.ylabel('$\log_{10}(F(n))$ (Fluctuación)', fontsize=12)
    plt.title('Análisis DFA: Correlaciones a Largo Alcance en ADN', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    if nombre_archivo and config.GUARDAR_GRAFICAS:
        ruta_salida = config.DIRECTORIO_GRAFICAS / nombre_archivo
        plt.savefig(ruta_salida, dpi=config.DPI, bbox_inches='tight')
        logger.info(f"Gráfica de Hurst guardada en: {ruta_salida}")
    
    plt.close()