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
