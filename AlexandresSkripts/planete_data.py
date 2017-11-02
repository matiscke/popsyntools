import os
import numpy
import tables

def nameToNumber(name):
	try:
		seppos = name.find( "SIM" )
		if seppos >= 0:
			name = name[ seppos + 3 : ]
		return name.lstrip( "0" )
	except:
		return name

def areSame( *names ):
	curNumber = None
	for name in names:
		number = nameToNumber( name )
		if curNumber is None:
			curNumber = number
		elif curNumber != number:
			return False
	return True

def listFromFile(path):
	entries = []
	for line in open( path, "r" ):
		entries.append( line.strip() )
	return NameList( entries )

class NameList(object):
	def __init__(self, names):
		self.raw = names
		self.names = None
		self.numbers = None
		self.iters = None
		self.asNumbers = False
		self.cur = None

	def _buildNames(self):
		cdNSigni = 0
		simNSigni = 0

		cdNumbers = []
		simNumbers = []

		for entry in self.raw:
			if entry.startswith( "CD_" ):
				seppos = entry.find( "_SIM" )
				cdNumber = entry[ 3 : seppos ]
				simNumber = entry[ seppos + 4 : ]
			elif entry.startswith( "SIM" ):
				cdNumber = None
				simNumber = entry[ 3 : ]
			else:
				cdNumber = None
				simNumber = entry

			cdNumbers.append( cdNumber )
			simNumbers.append( simNumber )

			nSigni = 0 if cdNumber is None else len( cdNumber ) - cdNumber.rfind( "0" ) - 1
			if cdNSigni < nSigni:
				cdNSigni = nSigni
			nSigni = len( simNumber ) - simNumber.rfind( "0" ) - 1
			if simNSigni < nSigni:
				simNSigni = nSigni

		self.names = []
		for i in range( len( self ) ):
			padSimNumber = simNumbers[ i ]
			if len( simNumbers[ i ] ) < simNSigni:
				for i in range( simNSigni - len( simNumbers[ i ] ) ):
					padSimNumber = "0" + padSimNumber

			if cdNSigni == 0:
				self.names.append( "SIM{}".format( padSimNumber[ -simNSigni : ] ) )
			else:
				self.names.append( "CD_{}_SIM{}".format( cdNumbers[ i ][ -cdNSigni : ], padSimNumber[ -simNSigni : ] ) )

	def _buildNumbers(self):
		self.numbers = [ nameToNumber( raw ) for raw in self.raw ]

	def _buildIters(self):
		if self.asNumbers:
			if self.numbers is None:
				self._buildNumbers()
			self.iters = self.numbers
		else:
			if self.names is None:
				self._buildNames()
			self.iters = self.names

	def __len__(self):
		return len( self.raw )

	def __iter__(self):
		self.cur = -1
		self.max = len( self )
		return self

	def __next__(self):
		self.cur += 1
		if self.cur >= self.max:
			raise StopIteration

		return self[ self.cur ]

	def __getitem__(self, key):
		if self.iters is None:
			self._buildIters()

		return self.iters[ key ]

	def next(self):
		if self.cur is None:
			iter(self)

		return self.__next__()

	def __contains__(self, item):
		if item in self.raw:
			return True

		number = nameToNumber( item )
		if self.numbers is None:
			self._buildNumbers()

		return number in self.numbers

def getMin( entries ):
	ret = entries[0]
	for j in range( 1, len( entries ) ):
		if ret is None or (not entries[j] is None and entries[j] < ret):
			ret = entries[j]
	return ret

def getCommonList( *lists ):
	for lst in lists:
		lst.asNumbers = True

	listslen = len( lists )
	entries = [ lst.next() for lst in lists ]
	uniq = []

	while True:
		if areSame( *entries ):
			cur = entries[0]
			uniq.append( cur )
		else:
			cur = getMin( entries )

		for j in range( listslen ):
			if entries[j] == cur:
				try:
					entries[j] = lists[j].next()
				except StopIteration:
					return NameList( uniq )

def getMergedList( *lists ):
	for lst in lists:
		lst.asNumbers = True

	listslen = len( lists )
	entries = [ lst.next() for lst in lists ]
	merged = []

	while True:
		cur = getMin( entries )

		if cur is None:
			return NameList( merged )

		merged.append( cur )

		for j in range( listslen ):
			if entries[j] == cur:
				try:
					entries[j] = lists[j].next()
				except StopIteration:
					entries[j] = None

class FilterFinished(object):
	def getName(self):
		return "finished"

	def __call__(self, system):
		if not system.line is None and len( system.line ) > 200:
			# The fortran driver/manager system that writes to simulation_list.dat
			return "TER" in system.line

		try:
			# The python driver that send planete's output to the "sortie" file
			sortie = open( system.path + "/sortie", "r" )
			gotLine = False
			for line in sortie:
				if "DISK DISAPPEARED" in line or "TERMINE" in line:
					gotLine = True
					break
			sortie.close()
			return gotLine
		except FileNotFoundError:
			pass

		try:
			os.stat( system.path + "/evolution_001.outputdat" )
			return True
		except FileNotFoundError:
			pass

		try:
			os.stat( system.path + "/tenddisk" )
			return True
		except FileNotFoundError:
			return False

class FilterMinTime(object):
	def __init__(self, minTime = 1.e4):
		self.minTime = minTime

	def getName(self):
		return "mintime"

	def __call__(self, system):
		oldLast = system.onlyLast( True )

		try:
			lastTime = system[ 1 ][ -1, 1 ]
		except:
			lastTime = 0.

		system.onlyLast( oldLast )

		return lastTime > self.minTime

class FilterMassJump(object):
	def __init__(self, factor = 1.5):
		self.factor = factor

	def getName(self):
		return "massjump"

	def __call__(self, system):
		oldLast = system.onlyLast( False )

		for planete in system:
			mass = planete[ :, 4 ]
			for i in range( len( mass ) - 1 ):
				if mass[ i + 1 ] / mass[ i ] > self.factor and mass[ i ] > 1.:
					system.onlyLast( oldLast )
					return False

		system.onlyLast( oldLast )
		return True

filters = [ FilterFinished(), FilterMinTime(), FilterMassJump() ]

class PlaneteSystem(object):
	'''
	Base class for planetary system (i.e. one planete simulation)
	'''

	def __init__(self, pop, name):
		self.pop = pop
		self.name = name
		self.line = None
		self.filters = {}
		self.last = False

	def __len__(self):
		raise NotImplementedError

	def __getitem__(self, key):
		raise NotImplementedError

	def __iter__(self):
		self.cur = 0
		self.max = len( self )
		return self

	def __next__(self):
		self.cur += 1
		if self.cur > self.max:
			raise StopIteration

		return self[ self.cur ]

	def next(self):
		return self.__next__()

	def onlyLast(self, v):
		old = self.last
		self.last = v
		return old

	def setFilter(self, name, val):
		self.filters[ name ] = val

	def getFilter(self, name):
		if not name in self.filters:
			if self.pop is None:
				return None
			for filt in self.pop.filters:
				if filt.getName() == name:
					self.filters[ name ] = filt( self )

		return self.filters[ name ]

	def listOtherData(self):
		raise NotImplementedError

	def getOtherData(self, name):
		raise NotImplementedError

	def listEMPSData(self):
		raise NotImplementedError

	def getEMPSData(self, name):
		raise NotImplementedError

	def getDiskLifetime(self):
		raise NotImplementedError

	def getText(self, name):
		raise NotImplementedError

class PlaneteSystemFiles(PlaneteSystem):
	'''
	Planetary system whose data comes directly from planete output.
	'''

	def __init__(self, pop, name, path, line):
		PlaneteSystem.__init__(self, pop, name)
		self.path = path
		self.line = line
		self.nplanets = -1

	def __len__(self):
		if self.nplanets >= 0:
			return self.nplanets

		self.nplanets = 0
		self.data = []

		for n in range( 1000 ):
			try:
				os.stat( "{}/tracks_{:03d}.outputdat".format( self.path, n + 1 ) )
				self.nplanets += 1
				self.data.append( None )
			except FileNotFoundError:
				break

		return self.nplanets

	def __getitem__( self, key ):
		nmax = len( self )
		if key < 1 or key > nmax:
			raise Exception( "Invalid planet number requested for system {}: {} but there are only {}".format( self.name, key, nmax ) )

		if not self.data[ key - 1 ] is None:
			return self.data[ key - 1 ]
		elif self.last:
			return numpy.loadtxt( "{}/resultat_{:03d}.outputdat".format( self.path, key ), ndmin = 2 )
		else:
			#self.data[ key - 1 ] = numpy.loadtxt( "{}/tracks_{:03d}.outputdat".format( self.path, key ), ndmin = 2 )
			self.data[ key - 1 ] = numpy.genfromtxt( "{}/tracks_{:03d}.outputdat".format( self.path, key ), invalid_raise = False )
			return self.data[ key - 1 ]

	def listOtherData(self):
		ret = []

		try:
			for entry in sorted( os.listdir( self.path ) ):
				if entry.endswith( ".outputdat" ) and not entry.startswith( "evolution_" ) and not entry.startswith( "resultat_" ) and not entry.startswith( "tracks_" ):
					ret.append( entry[ : -10 ] )

		except:
			pass

		return ret

	def getOtherData(self, name):
		filepath = "{}/{}.outputdat".format( self.path, name )
		res = os.stat( filepath )
		if res.st_size == 0:
			return None
		else:
			try:
				return numpy.loadtxt( filepath, ndmin = 2 )
			except:
				return None

	def listEMPSData(self):
		ret = []

		try:
			for entry in sorted( os.listdir( self.path + "/emps_output" ) ):
				if entry.endswith( ".out" ):
					ret.append( entry[ : -4 ] )
		except:
			pass

		return ret

	def getEMPSData(self, name):
		filepath = "{}/emps_output/{}.out".format( self.path, name )
		try:
			os.stat( filepath )
			return numpy.loadtxt( filepath, delimiter = ",", ndmin = 2 )
		except FileNotFoundError:
			return None

	def getDiskLifetime(self):
		for line in open( self.path + "/sortie", "r" ):
			if line.startswith( " setting disklifetime : " ):
				return float( line.split()[3] )
		return None

	def getText(self, name):
		if "/" in name:
			raise Exception( "Malformed file name: {}".format( name ) )

		try:
			return open( self.path + "/" + name, "r" ).read()
		except FileNotFoundError:
			return None

	def convertToTable( self, handle, group ):
		handle.set_node_attr( group, "n", len( self ) )

		if not self.line is None:
			handle.set_node_attr( group, "line", self.line )

		disklifetime = self.getDiskLifetime()
		if not disklifetime is None:
			handle.set_node_attr( group, "disklifetime", disklifetime )

		for filtName in self.filters:
			handle.set_node_attr( group, "filter_{}".format( filtName ), self.filters[ filtName ] )

		for i in range( len( self ) ):
			handle.create_array( group, "planet_{:03d}".format( i + 1 ), self[ i + 1 ] )

		for other in self.listOtherData():
			otherData = self.getOtherData( other )
			if otherData is not None and len( otherData ):
				handle.create_array( group, "other_{}".format( other ), otherData )

		for emps in self.listEMPSData():
			empsData = self.getEMPSData( emps )
			if empsData is not None and len( empsData ):
				handle.create_array( group, "emps_{}".format( emps ), empsData )

		for text in [ "sortie", "erreur" ]:
			textCont = self.getText( text )
			if not textCont is None:
				handle.set_node_attr( group, "text_{}".format( text ), textCont )


class PlaneteSystemTable(PlaneteSystem):
	'''
	Planetary system whose data comes from aggregated HDF5 file.
	'''

	def __init__(self, pop, key):
		PlaneteSystem.__init__(self, pop, key)

		try:
			self.group = self.pop.handle.root._f_get_child( self.name )
			self.accessor = self.group._f_get_child
		except AttributeError:
			self.group = self.pop.handle.root._f_getChild( self.name )
			self.accessor = self.group._f_getChild

		self.nplanets = -1
		self.other_data = []
		self.emps_data = []
		self.text = {}
		self.disklifetime = None

		for attrName in self.group._v_attrs._f_list():
			if attrName == "line":
				self.line = self.group._v_attrs[ attrName ]
			elif attrName == "n":
				self.nplanets = self.group._v_attrs[ attrName ]
			elif attrName == "disklifetime":
				self.disklifetime = self.group._v_attrs[ attrName ]
			elif attrName.startswith( "filter_" ):
				self.setFilter( attrName[ 7 : ], self.group._v_attrs[ attrName ] )
			elif attrName.startswith( "text_" ):
				self.text[ attrName[ 5 : ] ] = self.group._v_attrs[ attrName ]

		try:
			nodes_itr = self.group._f_iter_nodes()
		except AttributeError:
			nodes_itr = self.group._f_iterNodes()

		for child in nodes_itr:
			childName = child._v_name
			if childName.startswith( "other_" ):
				self.other_data.append( childName[ 6 : ] )
			if childName.startswith( "emps_" ):
				self.emps_data.append( childName[ 5 : ] )

	def __len__(self):
		if self.nplanets >= 0:
			return self.nplanets

		self.nplanets = 0

		for n in range( 1000 ):
			try:
				self.accessor( "{:03d}".format( n + 1 ) )
				self.nplanets += 1
			except:
				break

		return self.nplanets

	def __getitem__( self, key ):
		nmax = len( self )
		if key < 1 or key > nmax:
			raise Exception( "Invalid planet number requested for system {}: {} but there are only {}".format( self.name, key, nmax ) )

		try:
			array = self.accessor( "planet_{:03d}".format( key ) )
		except:
			array = self.accessor( "{:03d}".format( key ) )

		if self.last:
			return array[ slice( array.nrows - 1, array.nrows, None ), ... ]
		else:
			return array[ ..., ... ]

	def listOtherData(self):
		return self.other_data

	def getOtherData(self, name):
		try:
			return self.accessor( "other_{}".format( name ) )[ ..., ... ]
		except tables.exceptions.NoSuchNodeError:
			return None

	def listEMPSData(self):
		return self.emps_data

	def getEMPSData(self, name):
		try:
			return self.accessor( "emps_{}".format( name ) )[ ..., ... ]
		except tables.exceptions.NoSuchNodeError:
			return None

	def getText(self, name):
		if name in self.text:
			return self.text[ name ]
		else:
			return None

	def getDiskLifetime(self):
		return self.disklifetime

class PlaneteSystems(object):
	'''
	Base class for collections of planetary system, such as single simulations, populations, etc.
	'''
	def __init__(self):
		self.reqFilters = True
		self.filters = filters
		self.iterfilters = []

	def __iter__(self):
		self.cur = 0
		self.max = len( self )
		return self

	def __next__(self):
		while True:
			self.cur += 1
			if self.cur > self.max:
				raise StopIteration

			ret = self[ self.cur - 1 ]
			good = True
			for name in self.iterfilters:
				if not ret.getFilter( name ):
					good = False
					break
			if good:
				return ret

	def next(self):
		return self.__next__()

	def setFilters(self, filters):
		self.filters = filters

	def setIterFilters(self, filters):
		self.iterfilters = filters

	def getNameList(self):
		raise NotImplementedError

class PlaneteSystemsSingle( PlaneteSystems ):
	'''
	List of independent systems.
	'''
	def __init__( self, systems ):
		PlaneteSystems.__init__(self)
		self.systems = systems

	def __len__(self):
		return len( self.systems )

	def __getitem__( self, key ):
		return PlaneteSystemFiles( None, key, self.systems[ key ], None )


class PlaneteSystemsPopulation( PlaneteSystems ):
	'''
	Base class for population.
	'''

	def __init__(self, base, select):
		PlaneteSystems.__init__(self)
		self.base = base
		self.setSelect(select)
		self.systems = None

	def setSelect(self, select):
		self.done = False

		if select is None:
			self.select = None
		elif type(select) is NameList:
			self.select = select
		elif len( select ) > 0:
			self.select = NameList( select )
		else:
			self.select = None

	def _build(self):
		raise NotImplementedError

	def __len__(self):
		if not self.done:
			self.done = True
			self._build()

		return len( self.systems )

	def __getitem__( self, key ):
		if not self.done:
			self.done = True
			self._build()

		return self.systems[ key ]

class PlaneteSystemsPopulationFiles( PlaneteSystemsPopulation ):
	def __init__(self, base, select):
		PlaneteSystemsPopulation.__init__(self, base, select)

	def _build(self):
		self.systems = []

		names, entries, lines = self.getNameList()

		for i in range( len( names ) ):
			entry = entries[ i ]
			name = names[ i ]
			line = lines[ i ]
			if self.select is None or name in self.select:
				system = PlaneteSystemFiles( self, name, self.base + "/" + entry, line )
				self.systems.append( system )

	def getNameList(self):
		try:
			handle = open( self.base + "/simulation_list.dat" )
		except FileNotFoundError:
			handle = None

		entries = []
		lines = []

		if handle is None:
			for entry in sorted( os.listdir( self.base ) ):
				if entry.startswith( "CD_" ):
					entries.append( entry )
					lines.append( None )
		else:
			for line in handle:
				entry =  "{}_{}".format( line[ 0 : 17 ], line[ 119 : 136 ] )
				try:
					os.stat( self.base + "/" + entry )
					entries.append( entry )
					lines.append( line.strip( "\n" ) )
				except FileNotFoundError:
					pass

		names = NameList( entries )

		return names, entries, lines

class PlaneteSystemsTable( PlaneteSystemsPopulation ):
	def __init__(self, base, select):
		PlaneteSystemsPopulation.__init__(self, base, select)
		self.handle = None

	def __del__(self):
		if not self.systems is None:
			del self.systems
		if not self.handle is None:
			self.handle.close()

	def _build(self):
		self.systems = []

		names, entries, lines = self.getNameList()
		for name in names:
			if self.select is None or name in self.select:
				self.systems.append( name )

	def getNameList(self):
		entries = []
		lines = []

		# Really, thank you pytables for breaking backward compatibility
		if self.handle is None:
			try:
				pytables_open = tables.open_file
			except AttributeError:
				pytables_open = tables.openFile

			self.handle = pytables_open( self.base, "r" )

		try:
			group_itr = self.handle.walk_groups()
		except AttributeError:
			group_itr = self.handle.walkGroups()

		for group in group_itr:
			groupName = group._v_pathname
			if groupName.count( "/" ) == 1 and len( groupName ) > 1:
				entries.append( groupName[ 1 : ] )
				lines.append( None )

		names = NameList( entries )

		return names, entries, lines

	def __getitem__( self, key ):
		if not self.done:
			self.done = True
			self._build()

		return PlaneteSystemTable( self, self.systems[ key ] )

class PlaneteSystemsMix(object):
	'''
	Base class for collection of collections of planetary system
	'''
	def __init__(self, colls):
		self.colls = colls

	def __iter__(self):
		self.idx = 0
		self.cur = self.colls[ self.idx ]
		iter( self.cur )
		self.max = len( self.colls )
		return self

	def __next__(self):
		try:
			return self.cur.next()
		except StopIteration:
			while True:
				self.idx += 1
				if self.idx >= self.max:
					raise StopIteration
				self.cur = self.colls[ self.idx ]
				iter( self.cur )
				try:
					return self.cur.next()
				except StopIteration:
					pass

	def next(self):
		return self.__next__()

class MultiSystems(object):
	'''
	Comibination of various populations.
	'''
	def __init__(self, systems):
		self.systems = [ getSystems( system ) for system in systems ]
		self.choice = 'all'
		self.filters = []
		self.select = None
		self.init = None
		self.entries = None
		self.nmin = None

	def getCommonList(self):
		allNames = []
		for system in self.systems:
			names, entries, lines = system.getNameList()
			allNames.append( names )

		return getCommonList( *allNames )

	def setCommon(self):
		self.nmin = len( self.systems )
		self.select = self.getCommonList()

		for system in self.systems:
			system.setSelect( self.select )

	def getMergedList(self):
		allNames = []
		for system in self.systems:
			names, entries, lines = system.getNameList()
			allNames.append( names )

		return getMergedList( *allNames )

	def setMerged(self):
		self.nmin = 1
		self.select = self.getMergedList()

		for system in self.systems:
			system.setSelect( None )

	def which(self, choice):
		self.choice = choice

	def setIterFilters(self, filters):
		self.filters = filters

	def __len__(self):
		return len( self.systems )

	def __iter__(self):
		if self.choice == 'all' or self.choice == 'any':
			if self.choice == 'all':
				self.setCommon()
			else:
				self.setMerged()

			for system in self.systems:
				system.setIterFilters( [] )
				iter( system )
			self.init = True
			return self
		else:
			system = self.systems[self.choice]
			system.setIterFilters( self.filters )
			return iter( system )

	def __next__(self):
		if self.init:
			self.init = False
			self.entries = []
			for system in self.systems:
				self.entries.append( system.next() )

		while True:
			cur = getMin( [ None if system is None else system.name for system in self.entries ] )
			if cur is None:
				raise StopIteration

			ret = []
			for i in range( len( self.systems ) ):
				if not self.entries[i] is None and self.entries[i].name == cur:
					ret.append( self.entries[i] )
					try:
						self.entries[i] = self.systems[i].next()
					except StopIteration:
						self.entries[i] = None
				else:
					ret.append( None )

			n = 0
			for i in range( len( ret ) ):
				if not ret[i] is None:
					good = True
					for name in self.filters:
						if not ret[i].getFilter( name ):
							good = False
					if good:
						n += 1
					else:
						ret[i] = None

			if n >= self.nmin:
				return ret

	def next(self):
		return self.__next__()

def convertToTable( systems, filename ):
	# Really, thank you pytables for breaking backward compatibility
	try:
		pytables_open = tables.open_file
	except AttributeError:
		pytables_open = tables.openFile

	handle = pytables_open( filename, "w" )
	root = handle.root

	for system in systems:
		for filt in filters:
			system.getFilter( filt.getName() )

		group = handle.create_group( root, system.name )
		system.convertToTable( handle, group )

	handle.close()

def getSystems( desc, select = None ):
	items = desc.split( "," )
	first = True
	allSingle = True
	paths = []
	pops = []
	for item in items:
		subitems = item.split( ":" )
		itempath = subitems[ 0 ]
		if os.path.isfile( itempath ):
			thisType = "table"
		elif os.path.isdir( itempath ):
			if os.path.isfile( "{}/tracks_001.outputdat".format( itempath ) ):
				thisType = "single"
			else:
				thisType = "files"
		else:
			raise Exception( "Broken path {}".format( itempath ) )

		if thisType == "table":
			allSingle = False
			sel = select if len( subitems ) == 1 else subitems[ 1 : ]
			pops.append( PlaneteSystemsTable( itempath, sel ) )
		elif thisType == "files":
			allSingle = False
			sel = select if len( subitems ) == 1 else subitems[ 1 : ]
			pops.append( PlaneteSystemsPopulationFiles( itempath, sel ) )
		else:
			paths.append( itempath )
			pops.append( PlaneteSystemsSingle( [ itempath ] ) )

	if len( pops ) == 1:
		return pops[ 0 ]
	elif allSingle:
		return PlaneteSystemsSingle( paths )
	else:
		return PlaneteSystemsMix( pops )

