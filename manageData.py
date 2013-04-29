"""
Yet another attempt to quickly and consistently handel the data.
20121001 RAT
""" 
import numpy as np

def load_SerisMatrix(dataPath):
	"""Simple method to get sample ids, feature ids and 
	the data matrix out of a simple series matrix that 
	is used in the Geo database.
	"""
	fin = open(dataPath)
	dataread = False
	for line in fin:
		line = line.rstrip('\n').lstrip()
		tmp = line.split('\t')
		if tmp[0] == '!Sample_data_row_count':
			n = int(tmp[1].rstrip('"').lstrip('"'))
		if line=='!series_matrix_table_end':break
		if dataread==True:
			info = line.rstrip('\n').lstrip().split('\t')
			featureID.append(info[0].lstrip('"').rstrip('"'))
			# ok we could have missing values 
			tmp = np.array(info[1:],dtype='|S15')
			tmp[tmp==''] = 'nan'
			data[i,:] = np.array(tmp,dtype=float) 
			i = i+1

		if line=='!series_matrix_table_begin':
			dataread = True
			line = fin.next()
			sampleID = line.rstrip().lstrip().split('\t')[1:]
			m = len(sampleID)
			data = np.zeros((n,m))
			featureID = []
			i = 0

	# remove the annoying quotes from sample ids
	for i in range(len(sampleID)):
		sampleID[i] = sampleID[i].lstrip('"').rstrip('"')

	return data, np.array(featureID,dtype=str), np.array(sampleID,dtype=str)

def convert_serisMat2SRF(dataPath):
	"""Convert a seris matrix to standard rat format
	with a labled data matrix of samples (rows) vs features	(cols)
	and files for the sampleData and sampleMeta.  
	NOTE: the only sample info will be the ids, other
	info must be added in a proper format.
	NOTE: Features may not be in same order as
	the platform feature information, most of the 
	code in this section does a check for this;
	however, we do not correct the order here; 
	furthermore, one must take care not to write 
	custom code that does not check the order.
	"""
	
	# load the data
	data,featureID,sampleID = manageData.load_SerisMatrix(dataPath)
	# get the dir 
	tmp = dataPath[-1::-1]
	m = tmp.find('/')
	outDir = tmp[m:][-1:0:-1]
	fout = open(outDir+'/data.dat','w')
	fout.write('nan\t')
	# we get a matrix of features vs samples, we want samples vs features
	data = data.T
	featureID.tofile(fout,sep='\t')
	fout.write('\n')
	n = len(sampleID)
	for i in range(n):
		fout.write(sampleID[i]+'\t')
		data[i].tofile(fout,sep='\t')
		fout.write('\n')

	fout.close()
	fout = open(outDir+'/sampleData.dat','w')
	fout.write('sample_ID\n')
	n = len(sampleID)
	for i in range(n):
		fout.write(sampleID[i]+'\n')

	fout.close()
	fout = open(outDir+'/sampleMeta.dat','w')
	fout.write('sample_ID\tstr\tGEO Identification tag')
	fout.close()
			
		

def load_platformData(platform,platformPath='',orderCheck=True):
	"""Get the feature data from a platform into a dictionary"""
	defaultPlatformPath ='/Users/rtasseff/DevMod/DB/platforms'
	
	if platformPath == '':
		platformPath = defaultPlatformPath

	featureData,featureMeta = loadDataMetaFiles(platformPath+'/'+platform+'/featureData.dat',platformPath+'/'+platform+'/featureMeta.dat',orderCheck,True)
	return featureData,featureMeta

def load_exp(dataDir,platformPath='',orderCheck=True):
	"""Basic load method for typical platform data, like the affy gene chips.
	Assumes the information is in a rat-standard format, we have
	a data, meta, sampleData, sampleMeta and featureData file avalible.
	Currently we assume featureData has a header and no featureMeta is avalible.
	We also assume standard names and that all data (except for feature data) 
	is in a single folder (dataDir).
	If order check is true we look for feature_ID and sample_ID to compare to the
	columns and rows of the data matrix.  Ensures order is right.
	We return
	meta dictonary a genral discription of the experiment, some keys are standard but not all
	data matrix of the primary data
	featureData dictonary for info on the features (typically probe meta data)
	sampleData dictonary for info on the samples (time, class labels, ect..)
	"""
	defaultPlatformPath ='/Users/rtasseff/DevMod/DB/platforms'
	
	if platformPath == '':
		platformPath = defaultPlatformPath

	# get the meta data (it at least has the platform info)
	meta =  loadMetaFile(dataDir+'/meta.dat')

	# load up the feature data from the correct platform
	platform = meta['platform']

	
	featureData = loadDataMetaFiles(platformPath+'/'+platform+'/featureData.dat',platformPath+'/'+platform+'/featureMeta.dat',orderCheck)
	
	sampleData = loadDataMetaFiles(dataDir+'/sampleData.dat',dataDir+'/sampleMeta.dat',orderCheck)

	data = np.loadtxt(dataDir+'/data.dat',dtype=str,delimiter='\t')
	
	if orderCheck:
		# check the sample order
		tmp = data[1:,0]
		try:
			_doOrderCheck(tmp,sampleData['sample_ID'])
		except ValueError:
			raise ValueError('Order of the samples in data and sampleData do not match.')
		tmp = data[0,1:]
		try:
			_doOrderCheck(tmp,featureData['feature_ID'])
		except ValueError:
			raise ValueError('Order of the features in data and featureData do not match.')
	
	# made it here, good to go
	data = np.array(data[1:,1:],dtype=float)

	return meta, data, sampleData, featureData


	
	
def _doOrderCheck(list1,list2,nMax = 10000):
	# list1 is current list2 is master
	n = min(len(list2),len(list1))
	n = min(n,nMax)
	for i in range(n):
		if list1[i]!=list2[i]:
			print i
			print list1[i]
			print list2[i]
			raise ValueError('List do not match!')
			


def loadDataMetaFiles(dataPath,metaPath,orderCheck=True,getMeta=False):
	"""Loads up the data meta combo files 
	and creates a proper dictionary for the data.
	Used for sampleData and featureData.
	"""
	# lets get the sample data now...
	sampleMeta = np.loadtxt(metaPath,dtype=str,delimiter='\t')
	sampleDataMatrix = np.loadtxt(dataPath,dtype=str,delimiter='\t')
	
		

	# the sampleMeta is not like meta, the entries are not genric 
	# we have specific entries for each sample
	# a label, a data type and a discription, here we need the first two
	sampleData = {}
	if getMeta: metaDic = {}
	if sampleMeta.ndim==2:
		n = len(sampleMeta)
		for i in range(n):
			# run a check for consistance
			if orderCheck and sampleMeta[i,0]!=sampleDataMatrix[0,i]:
					raise ValueError('inconsistency in meta/data files for col '+str(i)+' in '+dataPath)
			if sampleMeta[i,1]=='float':
				sampleData[sampleMeta[i,0]] = np.array(sampleDataMatrix[1:,i],dtype=float)
			elif sampleMeta[i,1]=='int':
				sampleData[sampleMeta[i,0]] = np.array(sampleDataMatrix[1:,i],dtype=int)
			elif sampleMeta[i,1]=='str':
				sampleData[sampleMeta[i,0]] = sampleDataMatrix[1:,i]
			else:
				raise ValueError('column '+str(i)+' form data is of unknown datatype: '+sampleMeta[i,1]+', in '+dataPath)
			if getMeta:
				metaDic[sampleMeta[i,0]] = sampleMeta[i,2]
	else:
		if orderCheck and sampleMeta[0]!=sampleDataMatrix[0]:
					raise ValueError('inconsistency in meta/data files for col '+str(i)+' in '+dataPath)
		if sampleMeta[1]=='float':
			sampleData[sampleMeta[0]] = np.array(sampleDataMatrix[1:],dtype=float)
		elif sampleMeta[1]=='int':
			sampleData[sampleMeta[0]] = np.array(sampleDataMatrix[1:],dtype=int)
		elif sampleMeta[1]=='str':
			sampleData[sampleMeta[0]] = sampleDataMatrix[1:]
		else:
			raise ValueError('column '+str(i)+' form data is of unknown datatype: '+sampleMeta[1]+', in '+dataPath)
		if getMeta:
			metaDic[sampleMeta[0]] = sampleMeta[2]
		

	if getMeta:
		return sampleData, metaDic
	else:
		return sampleData
	



def loadMetaFile(fileName,metaType='single'):
	""" Make a dictionary out of a rat-standard 
	meta file.  Give me the path and I give you the 
	dic.
	""" 
	meta = {}
	if metaType=='single':
		fin = open(fileName)
	
		for line in fin:
			line = line.lstrip().rstrip()
			if line.startswith('>>'):
				data = line.split('\t')
				meta[data[1]] = data[2]
	elif metaType=='list':
		fin = np.loadtxt(fileName,dtype=str,delimiter='\t')
		n = len(fin)
		for i in range(n):
			meta[fin[i,0]]=fin[i,1:]
	return(meta)
