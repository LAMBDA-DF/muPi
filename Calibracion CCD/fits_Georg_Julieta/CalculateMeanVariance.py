#%%
#30/3/2024: este codigo no me funciona porque no tengo suficiente memoria para correrlo
from astropy.io import fits
import numpy as np
import os
from astropy.nddata import block_reduce
# Directorio que contiene los archivos FITS
#directorio = r'C:\Users\georg\OneDrive\Escritorio\Labo 6\codigos labo 6\im_sintetizada'
directorio = r'F:\Imagenes_con_camara'
#directorio = r'F:\filtros_camara'
#directorio = r'D:\imagenes de javier\fits'


# Contador para llevar la cuenta del número total de imágenes
total_imagenes = 0

datos_imagenes = []

i=0
# Itera sobre todos los archivos FITS en el directorio``
for archivo in os.listdir(directorio):
    i+=1
    print(i)
    if archivo.endswith('.fits'):
        ruta_archivo = os.path.join(directorio, archivo)
        with fits.open(ruta_archivo) as hdulist:
            datos_imagen = hdulist[0].data.astype(np.float64)-256 # hay que restar 256 si no son de las matrices reducidas
            #datos_imagen=datos_imagen[:1500, :2400]
            #datos_imagen = block_reduce(datos_imagen, (8, 8),np.mean)
            datos_imagenes.append(datos_imagen)
            total_imagenes += 1
            if total_imagenes == 1:
                promedio_datos = datos_imagen.copy()  # Inicializar con los datos de la primera imagen
            else:
                promedio_datos += (datos_imagen - promedio_datos) / total_imagenes

# Calcular la varianza del conjunto de datos de todas las imágenes
#varianza=[]
#for i in datos_imagen:
 #   x=(i-promedio_datos)**2
 #   varianza.append(x)
carga=datos_imagen.flatten()
##print(varianza, "lista varianza")
#varianza_datos=np.sum(varianza)/(total_imagenes-1)

varianza_datos = np.sum((datos_imagen - promedio_datos) ** 2 for datos_imagen in datos_imagenes) / (total_imagenes-1)

# Convertir los datos de la varianza a uint16
#varianza_final = varianza_acumulada.astype(np.uint16)

# Guardar la imagen FITS con los datos de la varianza
#nombre_archivo_varianza = 'varianza_imagen.fits'
#fits.writeto(nombre_archivo_varianza, varianza_final, overwrite=True)

#print(f"Se ha creado la imagen de varianza '{nombre_archivo_varianza}' con los datos de varianza del conjunto de {total_imagenes} imágenes.")

print(promedio_datos.shape, "dimension del promedio")
print(varianza_datos.shape, "dimension de la varianza")
print(promedio_datos, "matriz promedio")
print(varianza_datos, "matriz varianza")
hay_negativos = np.any(promedio_datos < 0)

# Imprimir el resultado
if hay_negativos:
    print("Sí, hay valores negativos en promedio_datos.")
else:
    print("No, no hay valores negativos en promedio_datos.")
#-------------------------------
# Guardar la imagen FITS con los datos de la varianza
#nombre_archivo_varianza = 'varianza_imagen.fits'
#fits.writeto(nombre_archivo_varianza, varianza_acumulada, overwrite=True)

# Guardar la imagen FITS con los datos del promedio
#nombre_archivo_promedio = 'promedio_imagen.fits'
#fits.writeto(nombre_archivo_promedio, promedio_datos, overwrite=True)

#print(f"Se han creado las imágenes de promedio '{nombre_archivo_promedio}' y varianza '{nombre_archivo_varianza}'.")

import matplotlib.pyplot as plt

# Leer los datos del archivo FITS de promedio y varianza
#promedio_datos = fits.getdata('promedio_imagen.fits')
#varianza_datos = fits.getdata('varianza_imagen.fits')

# Convertir los datos de la imagen de 2D a 1D para facilitar el trazado


promedio_datos = block_reduce(promedio_datos, (8, 8), np.mean)
varianza_datos = block_reduce(varianza_datos, (8, 8), np.mean)
promedio_flat = promedio_datos.flatten()
varianza_flat = varianza_datos.flatten()
ganancia=varianza_flat/promedio_flat
index =np.argwhere(ganancia<1.2)
ganancia_reducida=ganancia[index]
#------------
'''''
umbral = 1000
elementos_mayor_umbral = varianza_flat[varianza_flat > umbral]
print(elementos_mayor_umbral)
print(len(elementos_mayor_umbral))
print(varianza_flat)
print(len(varianza_flat))
'''''

from scipy.optimize import curve_fit

def modelo_lineal(x, g, b):
    return g * x + b
#------

parametros, covarianza = curve_fit(modelo_lineal, promedio_flat, varianza_flat)
# Obtén los parámetros óptimos
a_optimo = parametros[0]
b_optimo = parametros[1]
'''''
hist, bins = np.histogram(ganancia_reducida, bins=500)

# Encontrar el centro de cada bin
bin_centers = (bins[1:] + bins[:-1]) / 2

# Definir una función gaussiana
def gaussiana(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / stddev) ** 2)

# Ajustar la gaussiana al histograma
parametros_optimos_gauss, covarianza_gauss = curve_fit(gaussiana, bin_centers, hist, p0=[1, np.mean(ganancia_reducida), np.std(ganancia_reducida)])

# Imprimir los parámetros óptimos
print("Parámetros óptimos de la gaussiana:")
print("Amplitud:", parametros_optimos_gauss[0])
print("Media:", parametros_optimos_gauss[1])
print("Desviación estándar:", parametros_optimos_gauss[2])
print("Matriz de covarianza:", covarianza_gauss)
'''''
#----




# Imprime el valor del parámetro 'a' óptimo
print("Ganancia por curve_fit:", a_optimo)

# Calcula la varianza ajustada utilizando el modelo lineal y los parámetros óptimos
varianza_ajustada = modelo_lineal(promedio_flat, a_optimo, b_optimo)

print("la ganancia es", np.var(varianza_flat,ddof=1)/np.mean(promedio_flat))
print("matriz de covarianza",covarianza)


#----------------------------------
plt.figure(0)
#plt.figure(figsize=(8, 6))
plt.scatter(promedio_flat, ganancia, s=1, alpha=0.5)
#plt.ylim(0, 1000)
plt.xlabel('Esperanza')
plt.ylabel('Ganancia')
plt.grid(True)

# Graficar la varianza vs el promedio
plt.figure(1)
#plt.figure(figsize=(8, 6))
plt.scatter(promedio_flat, varianza_flat, s=1, alpha=0.5, label='Datos')
plt.plot(promedio_flat, varianza_ajustada, color='red', label='Ajuste lineal')
#plt.scatter(promedio_flat, ganancia, s=1, alpha=0.5)
#plt.ylim(0, 500)
plt.xlabel('Esperanza [ADU]',fontsize=14)
plt.ylabel('Varianza [ADU]$^2$',fontsize=14)
#plt.title('Varianza vs Esperanza (azul)',fontsize=16  )
plt.grid(True)
plt.legend(fontsize=14)


# Crear el mapa de colores

plt.figure(2)
plt.hist2d(promedio_flat, varianza_flat, bins=100, cmap='inferno')
#plt.ylim(0, 1000)
plt.xlabel('Esperanza')
plt.ylabel('Varianza')


# Agregar barra de colores
plt.colorbar(label='Intensidad')

#%%

# Graficar el histograma de la ganancia
plt.figure(3)

plt.hist(ganancia_reducida, bins=500, label='Datos', alpha=0.5, color='b')

# Graficar el ajuste gaussiano
#plt.plot(bin_centers, gaussiana(bin_centers, *parametros_optimos_gauss), color='r', label='Ajuste gaussiano')

# Configurar etiquetas y título
plt.xlabel('Ganancia [ADU/e$^-$]',fontsize=14)
plt.ylabel('Frecuencia',fontsize=14)
#plt.title('Histograma de ganancia',fontsize=16)
#plt.ylim(0,1.2e5)
plt.grid(True)
#plt.legend()
#plt.hist(ganancia,bins=300)

plt.xlim(0,1)

#plt.legend()

plt.figure(4)
plt.hist(carga,bins=50)
plt.title('Histograma de carga (m11)')
plt.ylabel('Frecuencia')
plt.xlabel('Carga')
#plt.xlim(0,800)
plt.grid(True)

#plt.figure(5)
#plt.title("Residuos")
#plt.plot(promedio_flat,residuos)
plt.show()
''''
esperanza_deseada = 500
margen = 1  # Puedes ajustar el tamaño del margen según tus necesidades

# Definir el rango de esperanza
esperanza_minima = esperanza_deseada - margen
esperanza_maxima = esperanza_deseada + margen

# Lista para almacenar las varianzas dentro del rango de esperanza
varianzas_esperanza_deseada = []

# Iterar sobre los pares de esperanza y varianza
for esperanza, varianza in zip(promedio_flat, varianza_flat):
    # Verificar si la esperanza está dentro del rango deseado
    if esperanza_minima <= esperanza <= esperanza_maxima:
        # Agregar la varianza a la lista
        varianzas_esperanza_deseada.append(varianza)

# Imprimir la lista de varianzas dentro del rango de esperanza deseado
print("Varianzas dentro del rango de esperanza", esperanza_minima, "-", esperanza_maxima, ":", np.var(varianzas_esperanza_deseada,ddof=1))
'''''

ganancia_promedio=np.mean(ganancia)
error_ganancias=np.var(ganancia,ddof=1)
error_ganancia_promedio=error_ganancias/(np.sqrt(len(promedio_datos-1)))
print("Esta es la ganancia promedio= ",ganancia_promedio," y el error es",error_ganancia_promedio)
ganancia_e=promedio_flat/varianza_flat
'''''
with open('Ganancia_e_exposicion_1s_numero_de_imagenes_60_componente_M22.txt', 'w') as file:
    # Escribe la cabecera
    file.write('Ganancia\n')
    
    # Escribe los datos
    for v in zip(ganancia_e):
        # f-string para formatear la cadena con los valores de v y p
        file.write(f'{v}\n')

residuos=varianza_ajustada-varianza_flat
plt.figure(8)
plt.scatter(promedio_flat,residuos, s=1, alpha=0.5)
plt.show()
'''''
