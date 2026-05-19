# DNA Sequences Analysis - Generador de Secuencias Sintéticas

## Descripción del Proyecto

Este proyecto implementa un análisis estadístico avanzado de secuencias de ADN con el objetivo de generar secuencias sintéticas que sean estadísticamente indistinguibles de las secuencias naturales.

### Objetivo Principal
Replicar secuencias de bases nitrogenadas (A, C, G, T) de forma que pasen pruebas estadísticas rigurosas, siendo imposible distinguir entre una secuencia sintética y una natural basándose en características estadísticas.

## Características Principales

✅ **Lectura robusta de archivos** - Manejo automático de errores  
✅ **Limpieza de datos** - Eliminación de caracteres no deseados  
✅ **Caminata aleatoria** - Análisis del patrón de comportamiento  
✅ **Regresión lineal** - Determinación de tendencias  
✅ **Rotación de datos** - Análisis de distribución  
✅ **Logging completo** - Seguimiento detallado de ejecución  
✅ **Configuración centralizada** - Fácil personalización  
✅ **Gráficas automáticas** - Visualización de resultados  

## Estructura del Proyecto

```
dna_sequences_project/
├── main.py                  # Script principal de ejecución
├── config.py               # Configuración centralizada
├── cadenas.py              # Módulo con funciones de análisis
├── logger_setup.py         # Configuración de logging
├── requirements.txt        # Dependencias de Python
├── .gitignore             # Archivos a ignorar en Git
├── README.md              # Este archivo
│
├── data/
│   ├── input/             # Archivos de entrada (.txt, .fasta)
│   └── output/            # Resultados y gráficas
│
├── logs/
│   └── cadenas.log        # Archivo de registro detallado
│
└── tests/                 # (Próximo paso) Pruebas unitarias
    └── test_cadenas.py
```

## Instalación

### Requisitos Previos
- Python 3.8 o superior
- pip (gestor de paquetes de Python)
- Git

### Pasos de Instalación

1. **Clonar el repositorio**
```bash
git clone https://github.com/tu-usuario/dna-sequences-analysis.git
cd dna_sequences_project
```

2. **Crear entorno virtual (recomendado)**
```bash
python -m venv venv
source venv/bin/activate  # En Linux/macOS
# o
venv\Scripts\activate     # En Windows
```

3. **Instalar dependencias**
```bash
pip install -r requirements.txt
```

## Uso

### Uso Básico

```bash
python main.py -i data/input/secuencia.txt
```

### Opciones Disponibles

```bash
# Ayuda
python main.py --help

# Con modo verbose (detallado)
python main.py -i data/input/Homo_sapiens.txt -v

# Sin guardar gráficas
python main.py -i data/input/secuencia.txt --no-graficas
```

### Estructura de Archivo de Entrada

El archivo debe ser un archivo de texto (.txt) que contenga bases nitrogenadas. Ejemplo:

```
ATGCATGCATGC...
```

El programa automáticamente:
- Elimina espacios en blanco
- Ignora caracteres especiales
- Convierte a mayúsculas
- Valida que solo contenga A, C, G, T

## Pipeline de Análisis

### 1. **Lectura del Archivo**
- Valida la extensión
- Lee el contenido
- Maneja errores de forma segura

### 2. **Limpieza de Cadena**
- Elimina caracteres no deseados
- Normaliza a mayúsculas
- Valida bases nitrogenadas
- Genera estadísticas de composición

### 3. **Caminata Aleatoria**
Genera una secuencia de pasos basada en la regla:
- **A o T**: -1 (paso hacia abajo)
- **C o G**: +1 (paso hacia arriba)

**Objetivo**: Modelar el comportamiento estadístico de la secuencia

### 4. **Regresión Lineal**
- Calcula correlación de Pearson
- Ajusta un modelo lineal
- Calcula R² (coeficiente de determinación)
- Determina la pendiente (tendencia)

### 5. **Rotación de Datos**
- Aplica una rotación 2D con el ángulo |pendiente|
- Analiza la distribución después de la rotación
- Proporciona insights sobre la orientación de la tendencia

## Configuración

Edita `config.py` para personalizar:

```python
# Tamaño de entrenamiento para regresión
TRAIN_SIZE = 0.8

# Estado aleatorio para reproducibilidad
RANDOM_STATE = 1234

# Guardar gráficas automáticamente
GUARDAR_GRAFICAS = True

# Modo verbose
VERBOSE = True
```

## Resultados

Después de ejecutar el análisis, encontrarás:

- **Logs detallados**: `logs/cadenas.log`
- **Gráficas**: `data/output/graficas/`
  - `caminata_aleatoria_[nombre].png`
  - `rotacion_[nombre].png`

## Ejemplos de Salida

```
======================================================================
ANÁLISIS COMPLETADO EXITOSAMENTE
======================================================================

RESUMEN DEL ANÁLISIS
============================================================
Archivo procesado: Homo_sapiens.txt
Longitud de cadena: 318 bases

Regresión Lineal:
  - Pendiente: 0.123456
  - R²: 0.456789

Rotación (ángulo: 7.07°):
  - Pendiente rotada: -0.987654
  - R² rotado: 0.654321
============================================================
```

## Próximas Fases

- [ ] Soporte para archivos FASTA
- [ ] Generación de secuencias sintéticas
- [ ] Pruebas estadísticas avanzadas
- [ ] Interfaz gráfica (GUI)
- [ ] Análisis paralelo para archivos grandes
- [ ] Base de datos de resultados

## Contribución

Para contribuir al proyecto:

1. Fork el repositorio
2. Crea una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. Commit tus cambios (`git commit -m 'Add some AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

## Autores

- Fajardo Baltazar Gustavo Alexis - Servicio Social / Proyecto Académico

## Licencia


## Contacto

- Email: ga.fajardobaltazar@ugto.mx
- GitHub: [@cruri](https://github.com/cruri)

## Agradecimientos

- Dr.  Jóse de Jesús Bernal Alvarado

---

**Última actualización**: 2026
**Versión**: 1.0.0
