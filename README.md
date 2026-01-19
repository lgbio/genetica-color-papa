# Pasos y scripts del trabajo con AGROSAVIA de genética del color en papa.

## LOG
jan/18: r0.35 : Done paper revision about outliers.

Jun/22: r0.34 : Improved Heritability and Genomic Selection Full. Plots Spanish.

May/29: r0.33 : Working on heritabiliy
May/08: r0.32 : Cloning to KLinux for Paper revision

Mar/28: r0.31 : Filtrado dos traits en Plots de variación de marcadores. Cambiado background a blanco
Mar/06: r0.3  : Creados scrips selección genómica variando número de marcadores
Ene/01: r0.2  : Corregido Heredabilidad. MultiGWAS con 3 Tools (No Shesis).
Oct/26: r0.1  : Version Inicial de los pasos, revisados hasta el 11

## Indicaciones
- Están ya revisados los pasos hasta el 11 (Los que tienen "x" todavia no)

- Los scripts de cada proceso toman las entradas del directorio "inputs" y dejan los resultados en el directorio "outputs"

- Todos los scripts se ejecutan con solo llamarlos y el mismo toma las entradas y guarda las salidas.

- Los procesos largos como el 08 de MultiGWAS por ahora no se ejecuta directamente, solo está los resultados y los archivos de configuración para lanzar MultiGWAS.

- Traté de colocar la estructura de "inputs" "outputs" para que no les toque meterse al código, que puede ser complicado entenderlo. Si necesitan algo de uno de estos scripts me dicen.

- Finalmente, el ambiente de desarrollo es en linux (ubuntu 22), R versión 4.4, y para reproducirlo he creado un ambiente de Conda (miniconda) que  instala el R y las librerías necesarias con la que están trabajando actualmente los scripts. Para esto necesitan el archivo "renv40.yml" que es está en el siguiente repo de github con todos los procesos y el cual deben descargar:

1. Descargar desde github el repositorio del proyecto (desde bash):
git clone https://github.com/lgarreta/genetica-color-procesos.git

2. Instalar Conda
    2.1. Descargar el instalador (desde bash):
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    2.2. Instalar Miniconda o Anaconda (desde bash):
bash Miniconda3-latest-Linux-x86_64.sh

    2.3. Verificar la instalación (desde bash): 
conda --version

3. Instalar r-curl via conda (desde bash):
conda install -n renv40 -c conda-forge r-curl

4. Crear el entorno a partir del archivo .yml (desde bash)
conda env create --name renv40 --file renv40.yml

5. Activar el nuevo entorno
conda activate renv40

6. Ubicarse en algún directorio de los diferentes procesos (desde bash):
    cd 05-Crear-RHSTable-Correlaciones-entre-fenotipos

    6.1. Ejecutar el Script (Creación tablas para correlaciones):
    R CMD 01-create-tables-correlations-TwoCodesColorTraits.R

    6.2. Ejecutar el Script (Creación tablas para correlaciones):
    R CMD 02-create-plots-correlations-TwoCodesColorTraits.R

Y así para los otros Scripts






