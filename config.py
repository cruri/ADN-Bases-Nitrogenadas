"""
config.py
=========
Configuración centralizada del proyecto de secuencias de ADN.
Permite cambiar parámetros sin modificar el código principal.
"""

import os
from pathlib import Path

# ============================================================================
# RUTAS DEL PROYECTO
# ============================================================================
PROJECT_ROOT = Path(__file__).parent.absolute()
DATA_DIR = PROJECT_ROOT / "data"
INPUT_DIR = DATA_DIR / "input"
OUTPUT_DIR = DATA_DIR / "output"
LOGS_DIR = PROJECT_ROOT / "logs"

# Crear directorios si no existen
for directory in [INPUT_DIR, OUTPUT_DIR, LOGS_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# ============================================================================
# CONFIGURACIÓN DE LECTURA DE ARCHIVOS
# ============================================================================
ARCHIVOS_ENTRADA = {
    "homo_sapiens": INPUT_DIR / "Homo_sapiens.txt",
    "secuencia": INPUT_DIR / "sequence.fasta",
    "citocromo_c": INPUT_DIR / "citocromo_c_Homo_sapiens.fasta"
    # Agregar más referencias aquí según sea necesario
}

# Extensiones de archivo soportadas
EXTENSIONES_SOPORTADAS = {".txt", ".fasta", ".fa"}

# ============================================================================
# CONFIGURACIÓN DE LIMPIEZA
# ============================================================================
# Bases nitrogenadas válidas
BASES_VALIDAS = {'A', 'C', 'G', 'T'}

# Caracteres a eliminar durante la limpieza
CARACTERES_A_ELIMINAR = r"!@#$%^&*()-_=+[]{};:',<.>?/\|`~\n\r\t "

# ============================================================================
# CONFIGURACIÓN DE REGRESIÓN LINEAL
# ============================================================================
TRAIN_SIZE = 0.8  # 80% para entrenamiento, 20% para test
RANDOM_STATE = 1234  # Estado fijo para reproducibilidad
SHUFFLE = True  # Revolver datos durante el split

# ============================================================================
# CONFIGURACIÓN DE ROTACIÓN
# ============================================================================
USAR_RADIANES = True  # Si False, usa grados directamente

# ============================================================================
# CONFIGURACIÓN DE VISUALIZACIÓN
# ============================================================================
FIGSIZE_CAMINATA = (12, 7)
FIGSIZE_ROTACION = (10, 8)
DPI = 100
GUARDAR_GRAFICAS = True
DIRECTORIO_GRAFICAS = OUTPUT_DIR / "graficas"

# Crear directorio de gráficas
DIRECTORIO_GRAFICAS.mkdir(parents=True, exist_ok=True)

# ============================================================================
# CONFIGURACIÓN DE LOGGING
# ============================================================================
NIVEL_LOG = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL
FORMATO_LOG = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
ARCHIVO_LOG = LOGS_DIR / "cadenas.log"

# ============================================================================
# CONFIGURACIÓN DE VERBOSIDAD
# ============================================================================
VERBOSE = True  # Imprimir información detallada en consola
MOSTRAR_ESTADISTICAS = True  # Mostrar estadísticas después de procesar

# ============================================================================
# VALIDACIÓN DE CONFIGURACIÓN
# ============================================================================
def validar_configuracion():
    """Valida que la configuración sea válida."""
    assert 0 < TRAIN_SIZE < 1, "TRAIN_SIZE debe estar entre 0 y 1"
    assert len(BASES_VALIDAS) == 4, "Debe haber exactamente 4 bases válidas"
    return True

if __name__ == "__main__":
    print("✓ Configuración cargada correctamente")
    print(f"  Proyecto root: {PROJECT_ROOT}")
    print(f"  Input dir: {INPUT_DIR}")
    print(f"  Output dir: {OUTPUT_DIR}")
