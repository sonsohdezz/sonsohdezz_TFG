# SCRIPT PARA CREAR ARCHIVO DE ENTRADA DE PHYLOSEQ A PARTIR DE DATOS EGGNOG


# __________________________________________________________________________________________________________________________________________


## MÓDULOS ##


# Importamos los módulos que se van a usar.
import sys
import time
import os
import csv
import re
import pandas as pd
import collections
import getopt
import shutil



## AYUDA Y COMPROBACIÓN DE ARGUMENTOS ##


# Función de ayuda.
def help():
	print("\nSINOPSIS\n")
	print("\tpython script.py -i <inputdir> -o <outputdir>\n")
	print("\nDESCRIPCIÓN\n")
	print("\t Este programa sirve para generar un archivo de entrada a partir de las taxonomías para diferentes muestras genómicas obtenidas con eggNOG-mapper.")
	print("\t El archivo resultado consiste en una tabla los OTUs encontrados, la taxonomía de cada uno, y el conteo de estos para cada muestra.")
	print("\t Las opciones de comando necesarios se describen a continuación:\n")
	print("\t [inputdir] \t Path a la carpeta que contiene las taxonomías obtenidas para las distintas muestras.")
	print("\t [outputdir] \t Path a la carpeta donde se quieran almacenar los archivos resultado.\n")

	time.sleep(2)
	sys.exit()


# Función de comprobación de argumentos.
def main(argv):

	inputdir = ''
	outputdir = ''

	try:
		opts, args = getopt.getopt(argv, "hi:o:")
	except getopt.GetoptError:
		help()
		sys.exit()

	# Verificamos que se han introducido todos los argumentos necesarios.
	if len(opts) < 1 :
		help()
		sys.exit()

	# Para cada parámetro, comprobamos el argumento introducido y almacenamos la información si es correcta.
	for opt, arg in opts:

		if opt in ("-h"):
			help()
			sys.exit()

		elif opt in ("-i"):
			inputdir = arg
		
		elif opt in ("-o"):
			outputdir = arg

		else:
			help()
			sys.exit()

	return inputdir, outputdir


# __________________________________________________________________________________________________________________________________________


## EXTRACCIÓN DE DATOS DE LOS ARCHIVOS DE TAXONOMÍA DE LAS MUESTRAS ##



# Función de generación de una lista con los path de los archivos contenidos en el directorio dado.
def files(inputdir):

	taxonomic_files = os.listdir(inputdir) 				# Recorremos el directorio.
	files_paths = []									# Creamos una lista para almacenar los paths de los archivos.	

	# Comprobamos la naturaleza de los archivos y añadimos el path de cada uno a la lista creada justo antes.
	for file in taxonomic_files:
		file_path = os.path.join(inputdir, file)
		if os.path.isfile(file_path) and file.endswith('contigs.csv'):
			files_paths.append(file_path)

	return files_paths


# Función de generación de un diccionario con la columna 'contig_ID' y la columna 'taxonomy' a partir de los archivos input.
def OTUs(files_paths):

	# Creamos un diccionario vacío, donde se almacenarán los distintos OTUs que se vayan encontrando.
	OTUs_dic = {}

	# Creamos una lista, donde se almacenarán los nombres de los distintos OTUs.
	OTUs_list = []

	n = 0 		# Inicializamos una variable para ir nombrando a los OTUs.

	# Recorremos la lista de paths de los archivos contenidos en el directorio dado.
	for input_file in files_paths:

		# Leemos los datos del archivo de taxonomía como un dataframe.
		data = pd.read_csv(input_file, delimiter=",")

		# Recorremos los elementos del dataframe.
		for i in range(len(data)):
			
			s = (data.iloc[i]['TAXONOMY'])           	# Almacenamos los datos de la columna de taxonomías en una variable 'taxas'.
			taxa = s

			if taxa in OTUs_dic.values():
				continue
			else:
				n += 1
				OTUs_list.append("OTU" + str(n))			
				OTUs_dic[("OTU" + str(n))] = taxa
				
	return OTUs_dic, OTUs_list


# Función de generación de una lista con los conteos de cada OTU para cada muestra.
def count(files_paths, OTUs_dic):

	# Creamos un diccionario vacío donde irán los conteos de cada OTU para las distintas muestras, por orden.
	count_dic = {}				

	# Recorremos la lista de paths de los archivos contenidos en el directorio dado.
	for input_file in files_paths:

		file_name = input_file.split("/")[-1].split("_lineage_contigs.csv")[0]		# Sacamos el nombre del archivo.
		
		# Leemos los datos del archivo de taxonomía como un dataframe.
		data = pd.read_csv(input_file, delimiter=",")
		
		# Inicializamos una lista donde irán los conteos para cada OTU.
		count_list = []				

		# Recorremos el archivo de la muestra y contamos el número de veces que aparece cada OTU.
		for OTU in OTUs_dic.values():		

			# Realizamos un try except ya que algunos OTUs no están en todas las muestras, y se debe almacenar un conteo de 0.
			try:
				count = data['TAXONOMY'].value_counts()[OTU]
				count_list.append(count)
			except KeyError:
				count_list.append(0)

		# En el diccionario asignamos la lista de conteos al nombre del archivo.
		count_dic[file_name] = count_list
				
	return count_dic



## EJECUCIÓN Y CREACIÓN DEL DATAFRAME GENÉRICO DE OTUs ##


# Ejecutamos las distintas funciones.
inputdir, outputdir = main(sys.argv[1:])
files_paths = files(inputdir)
OTUs_dic, OTUs_list = OTUs(files_paths)
count_dic = count(files_paths, OTUs_dic)

# Generamos un DATAFRAME de OTUs a partir del diccionario con las taxonomias únicas para todas las muestras.
df = pd.DataFrame.from_dict([OTUs_dic])
dfT = df.T									# Invertimos columnas y filas para obtener los contigs en las filas.
dfT.index.name = ('OTU')					# Le damos nombre a la primera columna.

# Añadimos los diccionarios con los conteos en las muestras como columnas.
for name, counts in count_dic.items():
	dfT[name] = counts

# Hacemos algunos cambios en el formato a nivel de nombre.
dfT = dfT.rename(columns = {0:'TAXONOMY'})		
dfT.index = dfT.index.str.replace('OTU', '', regex=True)
dfT['TAXONOMY'] = dfT['TAXONOMY'].str.replace(r'\[|\]|\(|\)|', '', regex=True).replace(',',';', regex=True)


# Guardamos el dataframe generado en un archivo .csv
dfT.to_csv('phyloseq_dataframe.tsv')

# Definimos las rutas para mover el archivo multifasta creado al directorio destino.
src_folder = os.getcwd()
dst_folder = outputdir

loc_1 = src_folder + "/" + 'phyloseq_dataframe.tsv'
loc_2 = dst_folder + "/" + 'phyloseq_dataframe.tsv'

shutil.move(loc_1, loc_2)

print("\n(1/4) Archivo genérico para phyloseq creado y guardado como 'phyloseq_dataframe.tsv'.")



# _____________________________________________________________________________________________________________________________


## EXTRACCIÓN DE DATOS PARA CREACIÓN DE 3 TABLAS DE INPUT ##


# Función que genera una lista con las muestras.
def samples(files_paths):

	samples_dic = {}
	n = 0

	for input_file in files_paths:
		
		file_name = input_file.split("/")[-1].split("_lineage_contigs.csv")[0]		# Sacamos el nombre del archivo.
		n += 1
		samples_dic[file_name] = 'S' + str(n)

	return samples_dic


# Función que genera un diccionario con los OTUs y su taxonomía en diccionario.
def taxa_dic(OTUs_dic, OTUs_list):

	# Creamos un diccionario de referencia con la que comparar los rangos que extraigamos de dic_ranks.
	ref_dic = {'SK':'Domain', 'P':'Phylum', 'C':'Class', 'O':'Order', 'F': 'Family', 'G':'Genus', 'SG':'Species'}	

	# Creamos un diccionario nuevo.
	taxonomy_dic = {}
		
	for OTU, taxonomy in OTUs_dic.items():		# Recorremos los elementos del diccionario de OTUs (OTU y taxonomía en forma de lista).

		# Creamos un diccionario nuevo.
		taxa_dic = {}

		# Recorremos la taxonomía (en formato de string) y la convertimos en una lista.
		taxonomy = taxonomy.split('(')[1].split(')')[0].replace("'", "")
		taxonomy = list(taxonomy.split(", "))

		for element in taxonomy:				# Recorreamos cada uno de los elementos de la lista taxonómica.

			element = element.replace(",","") 
			rank = element.split('_')[0]		# Definimos cuál es el rango. 
			taxa = element.split('_')[1]		# Definimos cuál es el nombre.

			# Detectamos el rango en nuestro diccionario ref_dic.
			if rank in ref_dic:
				new_rank = ref_dic[rank]		# Reemplazamos el rango de iniciales por un rango descriptivo.
				taxa_dic[new_rank] = taxa 		# Asignamos el rango a su nombre en el taxa como un diccionario.

		taxonomy_dic[OTU] = taxa_dic

	return taxonomy_dic



## EJECUCIÓN Y CREACIÓN DE LOS 3 DATAFRAMES INPUT ##


# PRIMER DATAFRAME #

samples_dic = samples(files_paths)

# Generamos un DATAFRAME de muestras a partir del diccionario creado en la función anterior y lo guardamos en un archivo .csv
df_samples = pd.DataFrame.from_dict([samples_dic])
df_samples = df_samples.T
df_samples = df_samples.rename(columns = {0:'SAMPLES'})

# Guardamos el archivo resultado en el directorio destino.
# df_samples.to_csv('phyloseq_samples.csv')
df_samples.to_csv('phyloseq_samples.tsv', sep="\t")

# Definimos las rutas para mover el archivo multifasta creado al directorio destino.
src_folder = os.getcwd()
dst_folder = outputdir

loc_1 = src_folder + "/" + 'phyloseq_samples.tsv'
loc_2 = dst_folder + "/" + 'phyloseq_samples.tsv'

shutil.move(loc_1, loc_2)

print("(2/4) Archivo de muestras creado y guardado en el directorio destino como 'phyloseq_samples.tsv'.")


# SEGUNDO DATAFRAME #

# Generamos un DATAFRAME de OTUs a partir del dataframe genérico y lo guardamos en un archivo .csv
df_OTUs = pd.DataFrame.from_dict([OTUs_dic])
df_OTUs = df_OTUs.T								

for name, counts in count_dic.items():
	df_OTUs[name] = counts

df_OTUs = df_OTUs.drop(df_OTUs.columns[[0]], axis=1)

# Guardamos el archivo resultado en el directorio destino.
# df_OTUs.to_csv('phyloseq_OTUs.csv')
df_OTUs.to_csv('phyloseq_OTUs.tsv', sep="\t")

# Definimos las rutas para mover el archivo multifasta creado al directorio destino.
src_folder = os.getcwd()
dst_folder = outputdir

loc_1 = src_folder + "/" + 'phyloseq_OTUs.tsv'
loc_2 = dst_folder + "/" + 'phyloseq_OTUs.tsv'

shutil.move(loc_1, loc_2)

print("(3/4) Archivo de OTUs creado y guardado en el directorio destino como 'phyloseq_OTUs.tsv'.")


# TERCER DATAFRAME #

taxonomy_dic = taxa_dic(OTUs_dic, OTUs_list)

# Generamos un DATAFRAME con las taxonomías acordes a cada OTU a partir del diccionario de OTUs y lo guardamos en un archivo .csv
df_taxa = pd.DataFrame.from_dict(taxonomy_dic)
df_taxa = df_taxa.T
df_taxa = df_taxa.where(pd.notnull(df_taxa), "Unknown")

# Guardamos el archivo resultado en el directorio destino.
# df_taxa.to_csv('phyloseq_taxas.csv')
df_taxa.to_csv('phyloseq_taxas.tsv', sep="\t")

# Definimos las rutas para mover el archivo multifasta creado al directorio destino.
src_folder = os.getcwd()
dst_folder = outputdir

loc_1 = src_folder + "/" + 'phyloseq_taxas.tsv'
loc_2 = dst_folder + "/" + 'phyloseq_taxas.tsv'

shutil.move(loc_1, loc_2)

print("(4/4) Archivo de taxonomías creado y guardado en el directorio destino como 'phyloseq_taxas.tsv'.\n")





	









