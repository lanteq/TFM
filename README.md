# TFM 

### Función "StratifiedCV_AnalysisFS"

Se trata de una función de aprendizaje máquina que integra cuatro métodos de selección de biomarcadores para datos de expresión génica. Ha sido implementada para el Trabajo de Fin de Máster, para el Máster Universitario de Ciencia de Datos e Ingeniería de Computadores por la Universidad de Granada. Se da permiso para cualquier tipo de proyecto en el que necesite usarse parcial o completamente el código, agradeciendo el reconocimiento de su autoría

### Guía y observaciones

Una vez que se han preparado los datos (hemos obtenido los valores de expresión de gen/matriz de expresión, se ha realizado un tratamiento del efecto batch,..) usaremos la función que hemos denominado `StratifiedCV_AnalysisFS` en la que se aplicará una validación cruzada estratificada con `L` folds con el objetivo de realizar:

+ una extracción de genes basándonos en el LFC (magnitud de la diferencia de expresión de cada gen en la muestra tumoral respecto a la sana).
+ una selección de características con el método elegido - RF, mRMR, IG/Semántica (autores Qi & Tang), RL/Semántica (nueva propuesta).
+ Opción de balancear las clases con dos técnicas diferentes (downsampling y upsampling) pre y post selección de características.
+ Normalización min-max.
+ Optimización de parámetros y predicción con el clasificador elegido.

Se unen las predicciones de cada fold y se realiza una evaluación según F1-score.

##### Parámetros de entrada
+ `X`: Una matriz de expresión o data.frame que contiene los genes en las columnas y las muestras en las filas.
+ `Y`: Un vector o factor que contiene las etiquetas para cada una de las muestras en el objeto de datos.
+ `ngenes`: número de genes a seleccionar (huella génica).
+ `fs`: método de selección de genes. Los valores deben ser uno de los cuatro siguientes "rf", "mrmr", "semantic", "semanticRL".
+ `sub_fs`: Su valor por defecto es "backward", es un parámetro que solo tiene sentido si se selecciona el método "semanticRL" ya que el primer bloque ofrece dos variantes (selección hacía atrás y hacia delante).
+ `clasif`: algoritmo de clasificación. Posibles valores: "knn" o "svm".
+ `L`: número de folds (validación cruzada estratificada).
+ `lfc` y `cov` parámetros para el filtrado de genes
+ `balanced`: lógica. Si `TRUE` se realiza balanceamiento. Por defecto `FALSE`.
+ `after`: lógica. Si `TRUE` se realiza el balanceamiento posterior a la selección. Por defecto `FALSE`.
+ `down`: lógica. Si `TRUE` se realiza balanceamiento downsampling. Por defecto `FALSE`, balanceamiento upsampling.

