La distribución de carpetas tiene los siguientes requerimientos (una vez instalados GMAT y Python 3.7):

1) Dentro de load_gmat.py, la variable "StartUp" debe coincidir con la ubicación (path) de api_startup_file.txt.
2) Dentro de api_startup_file.txt, la variable "OUTPUT_PATH" debe coincidir con la ubicación de la carpeta "Runs", donde se archivarán las carpetas con los resultados de cada órbita simulada.
3) Dentro de api_startup_file.txt, la variable "ROOT_PATH" debe coincidir con la ubicación de la carpeta "GMAT" obtenida en la descompresión del archivo descargado desde la página de GMAT.