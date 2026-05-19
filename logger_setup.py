"""
logger_setup.py
===============
Configuración centralizada de logging para el proyecto.
"""

import logging
import sys
from config import NIVEL_LOG, FORMATO_LOG, ARCHIVO_LOG, LOGS_DIR

def obtener_logger(nombre_modulo):
    """
    Obtiene un logger configurado para un módulo específico.
    
    Parameters
    ----------
    nombre_modulo : str
        Nombre del módulo (típicamente __name__)
    
    Returns
    -------
    logging.Logger
        Logger configurado listo para usar
    """
    logger = logging.getLogger(nombre_modulo)
    
    # No agregar handlers si ya los tiene (evitar duplicados)
    if logger.handlers:
        return logger
    
    # Establecer nivel
    logger.setLevel(getattr(logging, NIVEL_LOG))
    
    # Crear formateador
    formatter = logging.Formatter(FORMATO_LOG)
    
    # Handler para consola
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, NIVEL_LOG))
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Handler para archivo
    try:
        LOGS_DIR.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(ARCHIVO_LOG)
        file_handler.setLevel(logging.DEBUG)  # Archivo siempre en DEBUG
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except Exception as e:
        logger.warning(f"No se pudo crear handler de archivo: {e}")
    
    return logger
