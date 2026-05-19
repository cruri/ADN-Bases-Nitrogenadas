# 🚀 Inicio Rápido

## 1️⃣ Configuración Inicial (Primera vez)

```bash
# Clonar el repositorio
git clone https://github.com/tu-usuario/dna-sequences-analysis.git
cd dna_sequences_project

# Crear entorno virtual
python3 -m venv venv
source venv/bin/activate

# Instalar dependencias
pip install -r requirements.txt
```

## 2️⃣ Preparar Datos

Coloca tu archivo de texto en:
```
data/input/mi_secuencia.txt
```

El archivo debe contener bases nitrogenadas (A, C, G, T).

## 3️⃣ Ejecutar Análisis

```bash
# Análisis básico
python main.py -i data/input/mi_secuencia.txt

# Con información detallada
python main.py -i data/input/mi_secuencia.txt -v

# Sin guardar gráficas
python main.py -i data/input/mi_secuencia.txt --no-graficas
```

## 4️⃣ Consultar Resultados

- **Logs**: `logs/cadenas.log`
- **Gráficas**: `data/output/graficas/`

## 📊 Ejemplo Completo

```bash
# 1. Activar entorno
source venv/bin/activate

# 2. Ejecutar
python main.py -i data/input/Homo_sapiens.txt

# 3. Ver resultados en consola y gráficas generadas
```

## 🆘 Solución de Problemas

### "No se encontró el módulo numpy"
```bash
pip install -r requirements.txt
```

### "Archivo no encontrado"
Verifica que el archivo existe en `data/input/` con el nombre exacto

### "Error de permisos en Linux"
```bash
chmod +x main.py
```

## 📝 Ejemplo de Archivo de Entrada

Crea `data/input/prueba.txt`:
```
ATGCATGCATGCGGGCCCAAATTT
ATGCATGCATGCGGGCCCAAATTT
```

Luego ejecuta:
```bash
python main.py -i data/input/prueba.txt
```

¡Listo! 🎉
