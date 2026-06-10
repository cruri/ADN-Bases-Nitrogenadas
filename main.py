"""
main.py
=======
Script principal para el análisis de secuencias de ADN.

Flujo:
1. Leer archivo de entrada
2. Limpiar cadena
3. Generar caminata aleatoria
4. Realizar regresión lineal
5. Aplicar rotación
6. Generar gráficas
"""

import sys
import argparse
import numpy as np
from pathlib import Path

import config
from logger_setup import obtener_logger
from cadenas import (
    leer_archivo,
    validar_extension,
    limpieza,
    analizar_kmers, 
    caminata_aleatoria,
    regresion_lineal,
    rotacion,
    graficar_caminata,
    graficar_rotacion
)

logger = obtener_logger(__name__)


def main(ruta_entrada, verbose=True):
    """
    Función principal que orquesta todo el pipeline.
    
    Parameters
    ----------
    ruta_entrada : str o Path
        Ruta al archivo de entrada
    verbose : bool
        Si True, muestra información detallada
    
    Raises
    ------
    Exception
        Captura y registra cualquier error durante el procesamiento
    """
    try:
        logger.info("="*70)
        logger.info("INICIANDO ANÁLISIS DE SECUENCIA DE ADN")
        logger.info("="*70)
        
        # =====================================================================
        # PASO 1: VALIDACIÓN Y LECTURA
        # =====================================================================
        logger.info("\n[PASO 1/5] Validando y leyendo archivo...")
        
        ruta_entrada = Path(ruta_entrada)
        validar_extension(ruta_entrada) # True -> Todo Ok | False -> Mala extension
        contenido = leer_archivo(ruta_entrada)
        
        # =====================================================================
        # PASO 2: LIMPIEZA
        # =====================================================================
        logger.info("\n[PASO 2/5] Limpiando cadena...")
        
        cadena_limpia = limpieza(contenido)
        
        # =====================================================================
        # ANÁLISIS DE K-MERS
        # =====================================================================
        logger.info("\n[EXTRA] Analizando k-mers (dimeros)...")
        conteno_kmers, probabilidades_kmers = analizar_kmers(cadena_limpia, k=2, verbose=verbose)
        # =====================================================================
        # PASO 3: CAMINATA ALEATORIA
        # =====================================================================
        logger.info("\n[PASO 3/5] Generando caminata aleatoria...")
        
        caminata = caminata_aleatoria(cadena_limpia, verbose=True)
        
        # Preparar datos para regresión
        x = np.arange(len(caminata)).reshape(-1, 1)
        
        # =====================================================================
        # PASO 4: REGRESIÓN LINEAL
        # =====================================================================
        logger.info("\n[PASO 4/5] Realizando regresión lineal...")
        
        modelo, pendiente = regresion_lineal(x, caminata, verbose=True)
        y_pred = modelo.predict(x)
        
        # Graficar caminata
        if config.GUARDAR_GRAFICAS:
            nombre_grafica = f"caminata_aleatoria_{Path(ruta_entrada).stem}.png"
            graficar_caminata(x, caminata, y_pred, pendiente, nombre_grafica)
        else:
            graficar_caminata(x, caminata, y_pred, pendiente)
        
        # =====================================================================
        # PASO 5: ROTACIÓN
        # =====================================================================
        logger.info("\n[PASO 5/5] Aplicando rotación...")
        
        # Usar el valor absoluto de la pendiente como ángulo
        angulo_rotacion = abs(pendiente)
        resultado_rotacion = rotacion(
            x, caminata, angulo_rotacion,
            usar_radianes=config.USAR_RADIANES,
            verbose=True
        )
        
        # Graficar datos rotados
        y_pred_rotado = resultado_rotacion['modelo'].predict(
            resultado_rotacion['x_rotado'].reshape(-1, 1)
        )
        
        if config.GUARDAR_GRAFICAS:
            nombre_grafica = f"rotacion_{Path(ruta_entrada).stem}.png"
            graficar_rotacion(
                resultado_rotacion['x_rotado'],
                resultado_rotacion['y_rotado'],
                y_pred_rotado,
                resultado_rotacion['theta'],
                nombre_grafica
            )
        else:
            graficar_rotacion(
                resultado_rotacion['x_rotado'],
                resultado_rotacion['y_rotado'],
                y_pred_rotado,
                resultado_rotacion['theta']
            )
        
        # =====================================================================
        # RESUMEN FINAL
        # =====================================================================
        logger.info("\n" + "="*70)
        logger.info("ANÁLISIS COMPLETADO EXITOSAMENTE")
        logger.info("="*70)
        
        if verbose:
            print("\n" + "="*60)
            print("RESUMEN DEL ANÁLISIS")
            print("="*60)
            print(f"Archivo procesado: {ruta_entrada.name}")
            print(f"Longitud de cadena: {len(cadena_limpia)} bases")
            print(f"Pasos en caminata: {len(caminata)}")
            print(f"\nRegresión Lineal:")
            print(f"  - Pendiente: {pendiente:.6f}")
            print(f"  - R²: {modelo.score(x, caminata):.6f}")
            print(f"\nRotación (ángulo: {np.degrees(resultado_rotacion['theta']):.2f}°):")
            print(f"  - Pendiente rotada: {resultado_rotacion['pendiente']:.6f}")
            print(f"  - R² rotado: {resultado_rotacion['r2']:.6f}")
            print("="*60 + "\n")
        
        logger.info(f"Los logs se encuentran en: {config.ARCHIVO_LOG}")
        logger.info(f"Las gráficas se guardaron en: {config.DIRECTORIO_GRAFICAS}")
        
        return {
            'exito': True,
            'cadena_limpia': cadena_limpia,
            'conteo kmers' : conteno_kmers,
            'probabilidaes kmers' : probabilidades_kmers,
            'caminata': caminata,
            'modelo': modelo,
            'pendiente': pendiente,
            'resultado_rotacion': resultado_rotacion
        }
    
    except Exception as e:
        logger.error(f"Error durante el procesamiento: {str(e)}", exc_info=True)
        print(f"\n❌ Error: {str(e)}", file=sys.stderr)
        return {'exito': False, 'error': str(e)}


def crear_argumentos():
    """Crea y retorna el parser de argumentos."""
    parser = argparse.ArgumentParser(
        description='Análisis de secuencias de ADN con caminata aleatoria',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python main.py -i data/input/Homo_sapiens.txt
  python main.py -i data/input/prueba.txt -v
  python main.py -i data/input/secuencia.txt --no-graficas
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Ruta al archivo de entrada (txt o fasta)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=True,
        help='Modo verbose (por defecto activado)'
    )
    
    parser.add_argument(
        '--no-graficas',
        action='store_true',
        help='No guardar gráficas en disco'
    )
    
    return parser


if __name__ == "__main__":
    # Validar configuración
    try:
        config.validar_configuracion()
    except AssertionError as e:
        logger.error(f"Error en configuración: {e}")
        sys.exit(1)
    
    # Parsear argumentos
    parser = crear_argumentos()
    args = parser.parse_args()
    
    # Ajustar configuración según argumentos
    if args.no_graficas:
        config.GUARDAR_GRAFICAS = False

    if args.input in config.ARCHIVOS_ENTRADA:
        ruta_a_analizar = config.ARCHIVOS_ENTRADA[args.input]
    else:
        ruta_a_analizar = args.input
    
    # Ejecutar análisis
    resultado = main(ruta_a_analizar, verbose=args.verbose)
    
    # Salir con código apropiado
    sys.exit(0 if resultado['exito'] else 1)
