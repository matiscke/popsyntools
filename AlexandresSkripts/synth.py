interactive = False

import planete_data
import sys
import math
import argparse
import numpy
import matplotlib
if not interactive:
	matplotlib.use( "pdf" )
import matplotlib.pyplot

parser = argparse.ArgumentParser()

parser.add_argument( "systems", nargs = "+", help = "planete systems to process" )
parser.add_argument( "-p", "--panels", help = "Panels to show" )
parser.add_argument( "-c", "--c", help = "Color scale" )
parser.add_argument( "-x", "--x", help = "Y axis" )
parser.add_argument( "-y", "--y", help = "X axis" )
parser.add_argument( "-s", "--select-file", help = "file containing list of systems to plot" )
parser.add_argument( "-o", "--output-file", default = "test.pdf", help = "output file" )
parser.add_argument( "-t", "--tracks", action = "store_true", help = "output file" )
parser.add_argument( "-u", "--publish", action = "store_true", help = "output file" )

args = parser.parse_args()

old_cmap = True
putColorBar = True
defaultColor = "coref" # "coref", "coll-n", "coll-m", "envf"

def makeTexStr( str ):
	if publish:
		return "$\\mathrm{" + str.replace( " ", "\ " ) + "}$"
	else:
		return "$\\mathdefault{" + str.replace( " ", "\ " ) + "}$"

def getConf( item, what, ax = None ):
	if ax.endswith( "-pref" ):
		key = ax[ :-5 ]
		axwhat = key + what
	elif ax.endswith( "-suff" ):
		key = ax[ :-5 ]
		axwhat = what + key
	else:
		key = what
		axwhat = what

	if axwhat in item:
		return item[ axwhat ]

	if type( item[ key ] ) == type( "" ):
		if what in quantities[ item[ key ] ]:
			return quantities[ item[ key ] ][ what ]
	else:
		for iteritem in item[ key ]:
			if what in quantities[ iteritem ]:
				return quantities[ iteritem ][ what ]

	return None

def setScale( item, axis, which ):
	scale = getConf( item, "scale", which + "-pref" )
	if scale == None:
		return

	kwargs = {}

	if scale == "symlog":
		linthresh = getConf( item, "linthresh", which + "-suff" )
		if linthresh != None:
			kwargs[ "linthresh" + which[0] ] = linthresh

	if which[0] == "x":
		axis.set_xscale( scale, **kwargs )
		if publish:
			if item[ "x" ] == "a":
				axis.axvline( 20., lw = 0.1, linestyle = "--", color = "k" )
	else:
		axis.set_yscale( scale, **kwargs )

	if scale == "linear" or scale == None:
		scilimits = getConf( item, "scilimits", which + "-suff" )
		if scilimits != None:
			axis.ticklabel_format( scilimits = scilimits, axis = which[0] )

def checkInteresting( xdata, ydata, status ):
	interesting = False

	#if status != 0 and xdata > 10. and xdata < 20. and ydata < 0.1:
	#	print( status )
	#	interesting = True

	if status != 0:
		return interesting

	#if xdata > 20. and ydata > 1.e2:
	#	interesting = True
	#if ydata > 1.e3:
	#	interesting = False

	#if ydata > 190.:
	#	print( xdata, ydata )
	#	interesting = True

	return interesting

def prepareScatter( plot, nplan, data, colls ):
	xname = plot[ "x" ]
	yname = plot[ "y" ]
	cname = plot[ "c" ] if "c" in plot else defaultColor

	if type( yname ) == list:
		yname = yname[0]

	xquant = quantities[ xname ]
	yquant = quantities[ yname ]
	cquant = quantities[ cname ] if cname in quantities else None

	xindex = xquant[ "column" ]
	yindex = yquant[ "column" ]
	cindex = cquant[ "column" ] if not cquant is None else None

	single = len( args.systems ) == 1
	interesting = False

	if not "xdata" in plot:
		plot[ "xdata" ] = []
		plot[ "ydata" ] = []
		plot[ "cdata" ] = []
		plot[ "edata" ] = []
		plot[ "sdata" ] = []

	if comp:
		common = None
		same = True
		sameFate = True
		for item in data:
			if item is None:
				same = False
				break
			elif common == None:
				common = item[ -1, 67 ]
			elif common != item[ -1, 67 ]:
				sameFate = False
				same = False
				break
		exists = common == 0
		if same:
			common = None
			for item in data:
				if common == None:
					common = item[ -1, xindex ]
				elif ( 2. * abs( common - item[ -1, xindex ] ) / abs( common + item[ -1, xindex ] ) ) > 0.1:
					same = False
					break
		if same:
			common = None
			for item in data:
				if common == None:
					common = item[ -1, yindex ]
				elif ( 2. * abs( common - item[ -1, yindex ] ) / abs( common + item[ -1, yindex ] ) ) > 0.1:
					same = False
					break
		if same:
			if exists:
				plot[ "xdata" ].append( data[ 0 ][ -1, xindex ] )
				plot[ "ydata" ].append( data[ 0 ][ -1, yindex ] )
				plot[ "cdata" ].append( data[ 0 ][ -1, 2 ] / data[ 0 ][ -1, 4 ] if old_cmap else data[ 0 ][ -1, 3 ] / data[ 0 ][ -1, 4 ] )
				plot[ "edata" ].append( -1 )
				plot[ "sdata" ].append( 5. )
				interesting = checkInteresting( data[ 0 ][ -1, xindex ], data[ 0 ][ -1, yindex ], 0 )
		else:
			c = 0
			xl = []
			yl = []
			for item in data:
				if not item is None and item[ -1, 67 ] == 0:
					plot[ "xdata" ].append( item[ -1, xindex ] )
					plot[ "ydata" ].append( item[ -1, yindex ] )
					plot[ "cdata" ].append( item[ -1, 2 ] / item[ -1, 4 ] if old_cmap else item[ -1, 3 ] / item[ -1, 4 ] )
					plot[ "edata" ].append( c )
					plot[ "sdata" ].append( 5. ) #if sameFate else 40. )
					xl.append( item[ -1, xindex ] )
					yl.append( item[ -1, yindex ] )
					interesting = interesting or checkInteresting( item[ -1, xindex ], item[ -1, yindex ], item[ -1, 67 ] )
				c += 1
			#if sameFate:
			#	plot[ "a" ].plot( xl, yl, "g", lw = 0.3 )
	else:
		nsys = 0
		for item in data:
			if single or item[ -1, 67 ] == 0:
				xdata = item[ -1, xindex ]
				ydata = item[ -1, yindex ]

				interesting = interesting or checkInteresting( xdata, ydata, item[ -1, 67 ] )
				skip = False

				if item[ -1, 67 ] == 0:
					if cname == "coll-n":
						if len( colls[ nsys ] ) == 0:
							cdata = -1
						else:
							cdata = colls[ nsys ][ nplan ][ "n" ] if nplan in colls[ nsys ] else 0
					elif cname == "coll-m":
						cdata = colls[ nsys ][ nplan ][ "m" ] if nplan in colls[ nsys ] else 0.
					elif cname == "coref":
						cdata = item[ -1, 2 ] / item[ -1, 4 ]
					elif cname == "envf":
						cdata = item[ -1, 3 ] / item[ -1, 4 ]
					else:
						cdata = item[ -1, cindex ]
						#if item[ -1, cindex ] <= 0.:
						#	cdata = -10.
						#else:
						#	cdata = math.log( item[ -1, cindex ] )
				else:
					if item[ -1, 67 ] == 2:
						xdata = 1000.
					cdata = -1 #float( "nan" )
					if nplan in colls[ nsys ]:
						#cdata = colls[ nsys ][ nplan ][ "regime" ] / 8.
						if item[ -1, 67 ] < 2:
							xdata = colls[ nsys ][ nplan ][ "a" ]
					skip = ( cname == "coll-n" or cname == "coll-m" )

				if not skip:
					plot[ "xdata" ].append( xdata )
					plot[ "ydata" ].append( ydata )
					plot[ "cdata" ].append( cdata )
					plot[ "edata" ].append( None )
					plot[ "sdata" ].append( 5. )
			nsys += 1

	return interesting

def plotTrack( plot, data, suffix, narg, nsys, npla, naxi, multi ):
	diffs = [ "marker", "linestyle", "color" ]

	xname = plot[ "x" ]
	ynames = plot[ "y" + suffix ]
	if type( ynames ) == type( "" ):
		ynames = [ ynames ]

	if len( ynames ) > 1:
		multi.append( "column" )

	axis = plot[ "a" + suffix ]

	xquant = quantities[ xname ]
	xindex = xquant[ "column" ]

	kwargs = {
		"color": "black",
		"linestyle": "-",
		"marker": None,
	}

	if "planet" in multi:
		planDiff = diffs.pop()
	if "arg" in multi:
		arguDiff = diffs.pop()
	if "system" in multi:
		systDiff = diffs.pop()
	if "axis" in multi:
		axisDiff = diffs.pop()
	if "column" in multi:
		thisDiff = diffs.pop()

	ncol = 0
	for yname in ynames:
		yquant = quantities[ yname ]
		yindex = yquant[ "column" ]
		xdata = data[ :, xindex ]
		ydata = data[ :, yindex ]

		if track:
			if yquant[ "all" ]:
				if npla == 0:
					# Select everything
					select = numpy.ones_like( data[ :, 67 ], dtype = bool )
				else:
					continue
			else:
				select = data[ :, 67 ] == 0 # emps_status
				#select = numpy.logical_and( data[ :, 0 ] >= 480, data[ :, 0 ] <= 550 )

			xdata = xdata.compress( select )
			ydata = ydata.compress( select )

			# Differentiation
			if "planet" in multi:
				kwargs[ planDiff ] = diffVals[ planDiff ][ npla % len( diffVals[ planDiff ] ) ]
			if "arg" in multi:
				kwargs[ arguDiff ] = diffVals[ arguDiff ][ narg % len( diffVals[ arguDiff ] ) ]
			if "system" in multi:
				kwargs[ systDiff ] = diffVals[ systDiff ][ nsys % len( diffVals[ systDiff ] ) ]
			if "axis" in multi:
				kwargs[ axisDiff ] = diffVals[ axisDiff ][ naxi % len( diffVals[ axisDiff ] ) ]
			if "column" in multi:
				kwargs[ thisDiff ] = diffVals[ thisDiff ][ ncol % len( diffVals[ thisDiff ] ) ]

		ncol += 1

		if "compute" in yquant:
			if yquant[ "compute" ] == "rcoll":
				ydata = ( 3. / 4. / math.pi * data[ :, 4 ].compress( select ) * 5.97424e+27 / ydata ) ** ( 1. / 3. )
			elif yquant[ "compute" ] == "rratio":
				ydata = data[ :, 14 ].compress( select ) / data[ :, 9 ].compress( select )
			elif yquant[ "compute" ] == "rcapcore":
				ydata = data[ :, 27 ].compress( select ) / data[ :, 9 ].compress( select )
			elif yquant[ "compute" ] == "rcaptot":
				ydata = data[ :, 27 ].compress( select ) / data[ :, 14 ].compress( select )
			elif yquant[ "compute" ] == "diff":
				old = ydata.copy()
				old[ 1 : ] = old[ : -1 ]
				old[ 0 ] = 0.
				ydata = ydata - old
			elif yquant[ "compute" ] == "+dtall":
				if npla == 0:
					time = data[ :, 1 ]
					old = time.copy()
					old[ 1 : ] = old[ : -1 ]
					old[ 0 ] = 0.
					dt = time - old
					axis.plot( data[ :, xindex ][ 1 : ], dt[ 1 : ], color = "lime", linestyle = kwargs[ "linestyle" ], lw = 0.3 )
			elif yquant[ "compute" ] == "ediff":
				ydata = (data[ :, 52 ] - data[ :, 53 ]).compress( select )
			elif yquant[ "compute" ] == "mdotgasratio":
				ydata = (data[ :, 8 ] / data[ :, 104 ]).compress( select )
			elif yquant[ "compute" ] == "lcont":
				if data[ 600, 108 ] == 0.:
					ydata = (data[ :, 107 ] + data[ :, 108 ]).compress( select )
				else:
					ydata = (data[ :, 107 ] + data[ :, 108 ] + data[ :, 55 ]).compress( select )
			elif yquant[ "compute" ] == "lcomp":
				ydata = (data[ :, 5 ] + data[ :, 73 ] - data[ :, 72 ] - data[ :, 79 ]).compress( select )
			else:
				raise Exception( "OOps" )

		if "test" in yquant:
			ye = data[ :, 62 ].compress( select )
			bottom = ydata * ( 1. - ye )
			top = ydata * ( 1. + ye )
			del kwargs[ "marker" ]
			axis.fill_between( xdata, bottom, top, alpha = 0.7, **kwargs )
		else:
			axis.plot( xdata, ydata, lw = 0.3, **kwargs )

def finaliseScatter( fig, plot, i ):
	if comp:
		edge = []
		line = []
		for v in plot[ "edata" ]:
			if v == -1:
				edge.append( "blue" )
				line.append( 0.6 )
			elif v == 0:
				edge.append( "black" )
				line.append( 0.6 )
			elif v == 1:
				edge.append( "red" )
				line.append( 0.6 )
			else:
				edge.append( v )
				line.append( 0. )
	else:
		edge = None
		line = None

	colorType = defaultColor if not "c" in plot or plot[ "c" ] is None else plot[ "c" ]

	if colorType == "coll-n":
		cmap = matplotlib.colors.ListedColormap( [ "b", "g", "y", "r" ] )
		cmap.set_under( "black" )
		norm = matplotlib.colors.BoundaryNorm( [ -0.5, 0.5, 1.5, 2.5, 3.5 ], ncolors = cmap.N )
		colors = plot[ "cdata" ]
		vmin = None
		vmax = None
	else:
		norm = matplotlib.colors.LogNorm()
		if colorType == "coll-m":
			vmin = 0.01
			vmax = 100.
		elif colorType == "coref" or colorType == "envf":
			vmin = 0.
			vmax = 1.
			norm = None
		elif colorType == "a":
			vmin = 0.01
			vmax = 100.
		elif colorType == "mplan":
			vmin = 0.01
			vmax = 1000.
		else:
			vmin = None
			vmax = None

		if popColor or comp:
			cmap = None
			if comp:
				colors = edge
			else:
				colors = "black" if i == 0 else "red"
			edge = None
			line = None
		else:
			cmap = matplotlib.cm.get_cmap( "jet" if old_cmap else None )
			cmap.set_under( "black" )
			cmap.set_bad( "black" )
			colors = plot[ "cdata" ]

	map = plot[ "a" ].scatter( plot[ "xdata" ], plot[ "ydata" ], s = plot[ "sdata" ], c = colors, edgecolors = edge, cmap = cmap, norm = norm, vmin = vmin, vmax = vmax, linewidths = line )
	if putColorBar and not popColor and not comp:
		cbar = fig.colorbar( map )
		if colorType == "coll-n":
			cbar.set_ticks( [ 0., 1., 2., 3. ] )
			cbar.set_label( makeTexStr( "Number of collisions" ) )
		elif colorType == "coll-m":
			cbar.set_label( makeTexStr( "Mass of the largest impactor [M_\\oplus]" ) )
		elif colorType == "coref":
			cbar.set_label( makeTexStr( "Core mass fraction" ) )
		elif colorType == "envf":
			cbar.set_label( makeTexStr( "Envelope mass fraction" ) )

	plot[ "xdata" ] = []
	plot[ "ydata" ] = []
	plot[ "cdata" ] = []
	plot[ "edata" ] = []
	plot[ "sdata" ] = []

def getData( systems, key, onlyLast ):
	ret = []

	for system in systems:
		if onlyLast:
			if system is None:
				missing = True
				data = None
			else:
				system.onlyLast( True )
				try:
					data = system[ key ]
					missing = False
				except Exception:
					data = None
					missing = True
				system.onlyLast( False )
				if not missing and len( data[ -1, : ] ) > phaseCol and not data[ -1, phaseCol ] in phaseType:
					data = None # Reload from full track
		else:
			missing = system is None
			data = None

		if data is None and not missing:
			data = system[ key ]
			if len( data[ -1, : ] ) > phaseCol:
				listPhaseType = [ t for t in phaseType ] # There is certainly a better way to copy.
				firstType = listPhaseType.pop()
				mask = data[ :, phaseCol ] == firstType
				for nextType in listPhaseType:
					mask = numpy.logical_or( mask, data[ :, phaseCol ] == nextType )
				data = data.compress( mask, axis = 0 )

		ret.append( data )

	return ret

def getCollList( system ):
	ret = {}
	colldata = system.getEMPSData( "collhist" )

	if colldata is None: #or len( colldata[ 0, : ] ) < 34:
		return ret

	for k in range( len( colldata[ :, 0 ] ) ):
		regime = int( colldata[ k, 32 ] ) if len( colldata[ k, ... ] ) > 32 else -1
		i1 = int( colldata[ k, 3 ] )
		i2 = int( colldata[ k, 4 ] )

		if i1 in ret:
			n1 = ret[ i1 ][ "n" ] + 1
			m1 = ret[ i1 ][ "m" ]
		else:
			n1 = 1
			m1 = 0.

		if i2 in ret:
			n2 = ret[ i2 ][ "n" ] + 1
			m2 = ret[ i2 ][ "m" ]
		else:
			n2 = 1
			m2 = 0.

		d1 = math.sqrt( colldata[ k, 5 ] ** 2 + colldata[ k, 6 ] ** 2 + colldata[ k, 7 ] ** 2 )
		d2 = math.sqrt( colldata[ k, 8 ] ** 2 + colldata[ k, 9 ] ** 2 + colldata[ k, 10 ] ** 2 )

		m1 = max( m1, colldata[ k, 24 ] * 1.9884138333164502e+33 / 5.97424e+27 )

		if regime > 2:
			import collscaling

			conf = collscaling.Conf()
			collscaling.conf_unit_merc( conf )
			collscaling.conf_model( conf, 2 )

			big = collscaling.Body( mass = colldata[ k, 23 ], radius = colldata[ k, 17 ], pos_x = colldata[ k, 5 ], pos_y = colldata[ k, 6 ], pos_z = colldata[ k, 7 ], vel_x = colldata[ k, 11 ], vel_y = colldata[ k, 12 ], vel_z = colldata[ k, 13 ] )
			small = collscaling.Body( mass = colldata[ k, 24 ], radius = colldata[ k, 18 ], pos_x = colldata[ k, 8 ], pos_y = colldata[ k, 9 ], pos_z = colldata[ k, 10 ], vel_x = colldata[ k, 14 ], vel_y = colldata[ k, 15 ], vel_z = colldata[ k, 16 ] )
			res, reg = collscaling.scale( conf, big, small, 2, 1 )

			rmap = { 1: 0, 2: 4, 3: 2, 4: 2 }
			r1 = rmap[ reg ]
			r2 = rmap[ reg ]
		else:
			#if not i1 in ret:
			r1 = 6
			#if not i2 in ret:
			r2 = 8

		ret[ i1 ] = { "regime": r1, "a": d1, "m": m1, "n": n1 }
		ret[ i2 ] = { "regime": r2, "a": d2, "m": m2, "n": n2 }

	return ret

track = args.tracks
publish = not track or args.publish

symEarth = "{\oplus}" if publish else "{Earth}"
symJup = "{\\jupiter}" if publish else "{Jupiter}"
symSun = "{\\astrosun}" if publish else "{Sun}"
symYear = "{yr}"

quantities = {
	"n": { "column": 0, "scale": "linear", "label": "n", "scilimits": ( 0, 3 ), "all": True },
	"t": { "column": 1, "scale": "log", "min": 1.e3, "max": 1.e7, "pmin": 1.e3, "pmax": 1.e7, "label": "Time [" + symYear + "]", "all": True },
	"dtall": { "column": 1, "scale": "log", "label": "timestep [" + symYear + "]", "compute": "diff", "all": True },
	"mcore": { "column": 2, "scale": "log", "label": "Core mass [M_" + symEarth + "]", "pmin": 1.e-2, "pmax": 5.e2, "all": False },
	"menv": { "column": 3, "scale": "log", "label": "Envelope mass [M_" + symEarth + "]", "pmin": 1.e-2, "pmax": 1.e4, "all": False },
	"mplan": { "column": 4, "scale": "log", "label": "Planet mass [M_" + symEarth + "]", "pmin": 1.e-2, "pmax": 1.e4, "all": False },
	"ltot": { "column": 5, "scale": "log", "label": "L_{tot} [L_" + symJup + "]", "all": False },
	"lcore": { "column": 6, "scale": "log", "label": "L_{core} [L_" + symJup + "]", "all": False },
	"mdotcore": { "column": 7, "scale": "log", "min": 1.e-8, "max": 1.e-2, "label": "Mdot core [M_" + symEarth + "/" + symYear + "]", "all": False },
	"mdotgasl": { "column": 8, "scale": "log", "min": 1.e-8, "max": 1.e-2, "label": "Mdot gas [M_" + symEarth + "/" + symYear + "]", "all": False },
	"rcore": { "column": 9, "scale": "log", "label": "Core radius [R_" + symJup + "]", "all": False },
	"pcore": { "column": 11, "scale": "log", "label": "Central pressure [bar]", "all": False },
	"tcore": { "column": 12, "scale": "log", "label": "Central temperature [K]", "all": False },
	"rhocen": { "column": 13, "scale": "linear", "label": "Central density [g/cm^3]", "all": False },
	"rtot": { "column": 14, "scale": "log", "label": "Outer radius [R_" + symJup + "]", "all": False },
	"rratio": { "column": 14, "scale": "log", "label": "Rtot/Rcore", "compute": "rratio", "all": False },
	"a": { "column": 18, "scale": "log", "label": "Semi-major axis [AU]", "pmin": 1.e-1 if track else 3.e-3, "pmax": 1.e3 if track else 1.5e3, "all": False },
	"ae": { "column": 18, "scale": "log", "label": "Semi-major axis [AU]", "test": True, "all": False },
	"mdiskg": { "column": 20, "scale": "linear", "label": "Gas disc mass [M_" + symSun + "]", "all": True },
	"mdiskp": { "column": 21, "scale": "linear", "label": "Planetesimal disc mass [M_" + symEarth + "]", "all": True },
	"rroche": { "column": 22, "scale": "log", "label": "Hills radius [R_" + symJup + "]", "all": True },
	"racc": { "column": 23, "scale": "log", "label": "Bondi radius [R_" + symJup + "]", "all": True },
	"dt": { "column": 24, "scale": "log", "min": 1.e-1, "max": 1.e5, "label": "Timestep [" + symYear + "]", "compute": "+dtall", "all": False },
	"sigmamean": { "column": 25, "scale": "log", "label": "Mean solid surface density [g/cm^2]", "all": False },
	"rfeed": { "column": 26, "scale": "log", "label": "Feeding zone width [AU]", "all": False },
	"rcap": { "column": 27, "scale": "log", "label": "Capture radius [R_" + symJup + "]", "all": False },
	"rcapcore": { "column": 27, "scale": "log", "label": "Rcap/Rcore", "compute": "rcapcore", "all": False },
	"rcaptot": { "column": 27, "scale": "log", "label": "Rcap/Rtot", "compute": "rcaptot", "all": False },
	"tneb": { "column": 29, "scale": "log", "label": "Tneb [K]", "all": False },
	"pneb": { "column": 30, "scale": "log", "label": "Pneb [bar]", "all": False },
	"type_mig": { "column": 31, "scale": "linear", "label": "Migration type", "all": False },
	"mejetot": { "column": 32, "scale": "linear", "label": "Ejected solid mass [M_" + symEarth + "]", "all": True },
	"macctot": { "column": 33, "scale": "linear", "label": "Accreted solid mass [M_" + symEarth + "]", "all": True },
	"miso": { "column": 34, "scale": "linear", "label": "Isolation mass [M_" + symEarth + "]", "all": False },
	"sigmagas": { "column": 35, "scale": "log", "label": "Gas surface density at planet location [g/cm^2]", "all": False },
	"dtpl": { "column": 39, "scale": "log", "label": "Time step lim [" + symYear + "]", "all": False },
	"rhocore": { "column": 41, "scale": "linear", "label": "Core density [g/cm^3]", "all": False },
	"mgazacc": { "column": 42, "scale": "linear", "label": "Accreted gas mass [M_" + symSun + "]", "all": False },
	"mgazevap": { "column": 43, "scale": "linear", "label": "Evaporated gas mass [M_" + symSun + "]", "all": False },
	"pout": { "column": 44, "scale": "log", "label": "Pout [bar]", "all": False },
	"tout": { "column": 45, "scale": "log", "label": "Tout [K]", "all": False },
	"lcont": { "column": 51, "scale": "log", "label": "L_{env,cont}", "all": False },
	"enew": { "column": 52, "scale": "log", "label": "E_{new}", "all": False },
	"eold": { "column": 53, "scale": "log", "label": "E_{old}", "all": False },
	"ediff": { "column": 53, "scale": "log", "label": "E_{diff}", "compute": "ediff", "all": False },
	"kenerg": { "column": 54, "scale": "log", "label": "kenerg", "all": False },
	"kenergdiff": { "column": 54, "scale": "linear", "min": -0.1, "max": 0.2, "compute": "diff", "label": "kenerg", "all": False },
	"lacccore": { "column": 55, "scale": "log", "label": "L_{core,acc}", "all": False },
	"ep": { "column": 56, "scale": "log", "label": "e_p", "all": False },
	"ip": { "column": 57, "scale": "log", "label": "i_p", "all": False },
	"e": { "column": 62, "scale": "log", "label": "Eccentricity", "pmin": 1.e-8, "pmax": 1., "all": False },
	"i": { "column": 63, "scale": "log", "label": "Inclination", "pmin": 1.e-15, "pmax": 10., "all": False },
	"typemig": { "column": 66, "scale": "linear", "min": 0., "max": 20., "label": "Type of migration", "all": False },
	"tmig": { "column": 69, "scale": "symlog", "linthresh": 1.e4, "label": "Migration timescale [" + symYear + "]", "all": False },
	"tmige": { "column": 70, "scale": "log", "label": "Ecc. damp. timescale [" + symYear + "]", "pmin": 1.e-1, "pmax": 1.e10, "all": False },
	"mdotgas": { "column": 73, "scale": "log", "label": "Gas acc rate [M_" + symEarth + "/" + symYear + "]", "all": False },
	"lactual": { "column": 100, "scale": "log", "label": "L_{actual}", "all": False },
	"corrlesti": { "column": 101, "scale": "log", "label": "corrLesti", "all": False },
	"correesti": { "column": 102, "scale": "log", "label": "corrEesti", "all": False },
	"mdotgasmax": { "column": 104, "scale": "log", "label": "Max gas acc rate [M_" + symEarth + "/" + symYear + "]", "all": False },
	"mdotgasratio": { "column": 8, "scale": "log", "min": 1.e-4, "max": 1.1, "label": "Mdotgas/Mdotgasmax", "compute": "mdotgasratio", "all": False },
	"dtmode": { "column": 78, "scale": "linear", "label": "Timestep limitation factor", "all": False },
	"lcomp": { "column": 6, "scale": "log", "label": "L_{sum}", "compute": "lcomp", "all": False },
	"lcontcore": { "column": 108, "scale": "log", "label": "L_{core,cont}", "all": False },
	"lcontenv": { "column": 107, "scale": "log", "label": "L_{env,cont}", "all": False },
	"lcontsum": { "column": 108, "scale": "log", "label": "L_{cont}", "compute": "lcont", "all": False },
}

panels = args.panels
if panels is None:
	panels = "usual" if track else "ma"

if panels == "debug":
	plots = [
		{ "x": "t", "y": "menv" },
		{ "x": "t", "y": "lcore" },
		#{ "x": "t", "y": "a" },
		#{ "x": "t", "y": "rtot", "y2": "rcoll", "ymin": 1.e8, "ymax": 1.e12, "y2min": 1.e8, "y2max": 1.e12 },
		#{ "x": "t", "y": "rhocen" },
		{ "x": "t", "y": "mdotgasl" },
		{ "x": "t", "y": "mdotcore" },
		#{ "x": "n", "y": "t", "yscale": "linear", "ymin": None, "scilimitsy": ( 0, 0 ) },
		#{ "x": "t", "y": "tmig" },
	]
elif panels == "usual":
	plots = [
		{ "x": "a", "y": "mplan" },
		{ "x": "t", "y": "mplan" },
		{ "x": "t", "y": "a" },
		{ "x": "n", "y": "t", "yscale": "linear", "ymin": None, "scilimitsy": ( 0, 0 ) },
		#{ "x": "t", "y": "menv" },
	]
elif panels == "dtmode":
	plots = [
		{ "x": "t", "y": "mplan" },
		{ "x": "n", "y": "dt" },
		{ "x": "n", "y": "dtmode" },
		{ "x": "n", "y": "t", "yscale": "linear", "ymin": None, "scilimitsy": ( 0, 0 ) },
	]
elif panels == "mass":
	plots = [
		{ "x": "t", "y": "mcore" },
		{ "x": "t", "y": "menv" },
		{ "x": "t", "y": "mplan" },
		{ "x": "n", "y": "t", "yscale": "linear", "ymin": None, "scilimitsy": ( 0, 0 ) },
	]
elif panels == "lumi":
	plots = [
		#{ "x": "t", "y": [ "lcomp", "lactual" ] },
		#{ "x": "t", "y": "corrlesti" },
		{ "x": "t", "y": [ "ltot", "lcore" ] },
		{ "x": "t", "y": "lacccore" },
		#{ "x": "t", "y": "lcontcore" },
		#{ "x": "t", "y": "lcontenv" },
	]
elif panels == "rad":
	plots = [
		{ "x": "t", "y": "rcore" },
		{ "x": "t", "y": "rtot" },
		{ "x": "t", "y": "rcapcore" },
		{ "x": "t", "y": "rcaptot" },
	]
elif panels == "rratio":
	plots = [
		{ "x": "t", "y": "rratio" },
	]
elif panels == "haha":
	plots = [
		{ "x": "n", "y": [ "menv", "lcore" ], "xscale": "linear", "xmin": 500, "xmax": 1250 },
	]
elif panels == "damp":
	plots = [
		{ "x": "a", "y": "mplan" },
		{ "x": "t", "y": "mplan" },
		{ "x": "t", "y": "a" },
		{ "x": "t", "y": "tmige", "ymin": 1.e-1, "ymax": 1.e10 },
	]
elif panels == "mig":
	plots = [
		{ "x": "a", "y": "mplan" },
		{ "x": "t", "y": "mplan" },
		{ "x": "t", "y": "a" },
		{ "x": "t", "y": "tmig" },
	]
elif panels == "pube":
	plots = [
		{ "x": "a", "y": "mplan" },
		{ "x": "t", "y": "mplan" },
		{ "x": "t", "y": "a" },
		{ "x": "t", "y": "e" },
	]
elif panels == "da":
	plots = [
		{ "x": "t", "y": "a" },
		{ "x": "t", "y": "tmige", "ymin": 1.e-1, "ymax": 1.e10 },
	]
elif panels == "at":
	plots = [
		{ "x": "t", "y": "a", },
	]
elif panels == "aet":
	plots = [
		{ "x": "t", "y": "ae", }, #"ymin": 1e-1, "ymax": 2e1 },
	]
elif panels == "orb":
	plots = [
		{ "x": "t", "y": "a", },
		{ "x": "t", "y": "e", },
		{ "x": "t", "y": "i", },
	]
elif panels == "mt":
	plots = [
		{ "x": "t", "y": "mplan", },
	]
elif panels == "mct":
	plots = [
		{ "x": "t", "y": "mcore", },
	]
elif panels == "met":
	plots = [
		{ "x": "t", "y": "menv", },
	]
elif panels == "men":
	plots = [
		{ "x": "n", "y": "menv", },
	]
elif panels == "ma":
	plots = [
		{ "x": "a", "y": "mplan" },
	]
elif panels == "mca":
	plots = [
		{ "x": "a", "y": "mcore" },
	]
elif panels == "ea":
	plots = [
		{ "x": "a", "y": "e" },
	]
elif panels == "rcap":
	plots = [
		{ "x": "t", "y": "rcap" },
	]
elif panels == "bound":
	plots = [
		{ "x": "t", "y": "pneb" },
		{ "x": "t", "y": "tneb" },
		{ "x": "t", "y": "pout" },
		{ "x": "t", "y": "tout" },
	]
elif panels == "sigma":
	plots = [
		{ "x": "t", "y": "sigmagas" },
		{ "x": "t", "y": "sigmamean" },
	]
elif panels == "typemig":
	plots = [
		{ "x": "t", "y": "typemig", "xscale": "linear", "xmin": 0., "xmax": 1.e6 },
	]
elif panels == "mdotcore":
	plots = [
		{ "x": "t", "y": "mdotcore" },
	]
elif panels == "mdotgas":
	plots = [
		{ "x": "t", "y": "mdotgasl", "y2": "mdotgasmax", "yscale": "log", "y2scale": "log", "ymin": 1.e-8, "ymax": 1., "y2min": 1.e-8, "y2max": 1. },
	]
elif panels == "mdisk":
	plots = [
		{ "x": "t", "y": "mdiskg", "xscale": None },
		{ "x": "t", "y": "mdiskp", "xscale": None },
	]
elif panels == "mdiskp2":
	plots = [
		{ "x": "t", "y": "mejetot" },
		{ "x": "t", "y": "macctot" },
	]
elif panels == "p2":
	plots = [
		{ "x": "t", "y": "ep" },
		{ "x": "t", "y": "ip" },
		{ "x": "t", "y": "mplan" },
		{ "x": "n", "y": "t" },
	]
elif panels == "rcoll":
	plots = [
		{ "x": "t", "y": [ "rcoll5", "rcoll4", "rcoll2", "rcap", "rcoll3", "rcoll1" ], "ymin": 2.e-2, "ymax": 2.e2, "y2min": 1.e-2, "y2max": 1.e4, "y2": [ "mplan", "menv", "mcore" ],
   "xi": "t", "yi": [ "rcoll5", "rcoll4", "rcoll2", "rcap", "rcoll3", "rcoll1" ], "xiscale": "linear", "ximin": 570.e3, "ximax": 580.e4, "yimin": 1., "yimax": 2.e2, "xipmin": 572.e3, "xipmax": 582.e3, "yipmin": 1., "yipmax": 2.e2, "xilabel": "", "yilabel": "" }
	]
else:
	plots = [
		{ "x": args.x, "y": args.y.split( "," ), "c": args.c },
	]

# Differentiation
diffVals = {
	"linestyle": [ "-", (0., (20., 15.)), (0., (20., 15., 8., 15.)), (0., (8., 6., 2., 6.)), ":", (0., (5., 3.)) ],
	"color": [ "k", "r", "g", "b", "m", "y", "c", "0.7", "0.5", "0.3" ],
	#           1    2    3    4    5    6    7    8      9      10
	#"color": [ "none", "none", "none", "none", "none", "none", "none", "none", "none", "k" ]
	"marker": [ ".", None ],
}

if publish:
	matplotlib.rc( "text", usetex = True )
	matplotlib.rc( "text.latex", preamble = "\\usepackage{wasysym}" )

two = len( plots ) == 2

if two:
	figsize = ( 6, 3 )
elif publish:
	figsize = ( 6 if putColorBar else 5, 4 )
else:
	figsize = None

fig = matplotlib.pyplot.figure( figsize = figsize )

hasTwin = False
nplots =  len( plots )
nx = int( math.sqrt( nplots ) )
ny = int( math.ceil( nplots / nx ) )
for c in range( nplots ):
	plots[ c ][ "a" ] = fig.add_subplot( nx, ny, c + 1 )
	if "y2" in plots[ c ]:
		plots[ c ][ "a2" ] = plots[ c ][ "a" ].twinx()
		hasTwin = True
	if "yi" in plots[ c ]:
		plots[ c ][ "ai" ] = fig.add_axes( [ 0.18, 0.65, 0.3, 0.3 ] )

if two:
	fig.subplots_adjust( left = 0.12, right = 0.85 if hasTwin else 0.97, bottom = 0.20, top = 0.93, hspace = 0.25, wspace = 0.27 )
elif publish:
	fig.subplots_adjust( left = 0.12, right = 0.85 if hasTwin else 0.98, bottom = 0.11, top = 0.97, hspace = 0.25, wspace = 0.27 )
else:
	fig.subplots_adjust( left = 0.10, right = 0.85 if hasTwin else 0.97, bottom = 0.10, top = 0.95, hspace = 0.25, wspace = 0.27 )

select = None
if args.select_file:
	select = []
	for line in open( args.select_file ):
		select.append( line.strip() )

phaseCol = 77
phaseType = [1, 2]

overallMulti = []
if len( args.systems ) > 1:
	overallMulti.append( "arg" )

planList = []

this = planete_data.MultiSystems( args.systems )

comp = not track and len( args.systems ) == 2
popColor = not track and len( args.systems ) > 1
any = True

if comp:
	if any:
		iterList = [ "any" ]
	else:
		iterList = [ "all" ]
else:
	iterList = range( len( args.systems ) )

#for systems in args.systems:
#	this = planete_data.getSystems( systems, select )
#	if len( this ) == 0:
#		print( "Note: {}, no corresponding system found.".format( systems ) )

if not track:
	this.setIterFilters( [ "finished", "mintime", "massjump" ] )

narg = 0

for item in iterList:
	this.which( item )

	argMulti = overallMulti
	if len( this ) > 1:
		argMulti.append( "arg" )

	nsys = 0
	for system in this:
		if type( system ) is list:
			sysl = system
		else:
			sysl = [ system ]

		if not track:
			for sys in sysl:
				if not sys is None:
					sys.onlyLast( True )

		systMulti = argMulti

		if len( planList ) > 0:
			if len( planList ) > 1:
				systMulti.append( "planet" )
		else:
			for sys in sysl:
				if not sys is None and len( sys ) > 1:
					systMulti.append( "planet" )
					break

		nplans = 0
		colls = []
		for sys in sysl:
			if not sys is None:
				nplans = max( nplans, len( sys ) )
				if not track:
					colls.append( getCollList( sys ) )

		for nplan in range( nplans ):
			if len( planList ) > 0 and not nplan + 1 in planList:
				continue

			for c in range( nplots ):
				interesting = False

				if not track:
					interesting = prepareScatter( plots[ c ], nplan + 1, getData( sysl, nplan + 1, True ), colls )

				if track or interesting:
					if interesting:
						print( sysl[0].name, nplan + 1 )

					plotMulti = systMulti
					if "y2" in plots[ c ]:
						plotMulti.append( "axis" )

					nitem = 0
					items = getData( sysl, nplan + 1, False )

					for data in items:
						plotTrack( plots[ c ], data, "", nitem if comp else narg, nsys, nplan, 0, plotMulti )
						if "y2" in plots[ c ]:
							plotTrack( plots[ c ], data, "2", nitem if comp else narg, nsys, nplan, 1, plotMulti )
						if "yi" in plots[ c ]:
							plotTrack( plots[ c ], data, "i", nitem if comp else narg, nsys, nplan, 0, plotMulti )
						nitem += 1

		nsys += 1

	if not track:
		for c in range( nplots ):
			finaliseScatter( fig, plots[ c ], narg )

	narg += 1

for c in range( nplots ):
	setScale( plots[ c ], plots[ c ][ "a" ], "x" )
	setScale( plots[ c ], plots[ c ][ "a" ], "y" )
	if "y2" in plots[ c ]:
		setScale( plots[ c ], plots[ c ][ "a2" ], "x" )
		setScale( plots[ c ], plots[ c ][ "a2" ], "y2" )
	if "yi" in plots[ c ]:
		setScale( plots[ c ], plots[ c ][ "ai" ], "xi" )
		setScale( plots[ c ], plots[ c ][ "ai" ], "yi" )

	xmin = getConf( plots[ c ], "pmin" if publish else "min", "x-pref" )
	xmax = getConf( plots[ c ], "pmax" if publish else "max", "x-pref" )
	if xmin != None and xmax != None:
		plots[ c ][ "a" ].set_xlim( xmin, xmax )
		if "y2" in plots[ c ]:
			plots[ c ][ "a2" ].set_xlim( xmin, xmax )

	ymin = getConf( plots[ c ], "pmin" if publish else "min", "y-pref" )
	ymax = getConf( plots[ c ], "pmax" if publish else "max", "y-pref" )
	if ymin != None and ymax != None:
		plots[ c ][ "a" ].set_ylim( ymin, ymax )

	if "y2" in plots[ c ]:
		y2min = getConf( plots[ c ], "pmin" if publish else "min", "y2-pref" )
		y2max = getConf( plots[ c ], "pmax" if publish else "max", "y2-pref" )
		if y2min != None and y2max != None:
			plots[ c ][ "a2" ].set_ylim( y2min, y2max )

	if "yi" in plots[ c ]:
		ximin = getConf( plots[ c ], "pmin" if publish else "min", "xi-pref" )
		ximax = getConf( plots[ c ], "pmax" if publish else "max", "xi-pref" )
		if ximin != None and ximax != None:
			plots[ c ][ "ai" ].set_xlim( ximin, ximax )

		yimin = getConf( plots[ c ], "pmin" if publish else "min", "yi-pref" )
		yimax = getConf( plots[ c ], "pmax" if publish else "max", "yi-pref" )
		if yimin != None and yimax != None:
			plots[ c ][ "ai" ].set_ylim( yimin, yimax )

	xlabel = getConf( plots[ c ], "label", "x-pref" )
	if xlabel != None:
		plots[ c ][ "a" ].set_xlabel( makeTexStr( xlabel ) )

	ylabel = getConf( plots[ c ], "label", "y-pref" )
	if ylabel != None:
		plots[ c ][ "a" ].set_ylabel( makeTexStr( ylabel ) )

	if "y2" in plots[ c ]:
		ylabel2 = getConf( plots[ c ], "label", "y2-pref" )
		if ylabel2 != None:
			plots[ c ][ "a2" ].set_ylabel( makeTexStr( ylabel2 ) )

	if "yi" in plots[ c ]:
		xilabel = getConf( plots[ c ], "label", "xi-pref" )
		if xilabel != None:
			plots[ c ][ "ai" ].set_xlabel( makeTexStr( xilabel ) )

		yilabel = getConf( plots[ c ], "label", "yi-pref" )
		if yilabel != None:
			plots[ c ][ "ai" ].set_ylabel( makeTexStr( yilabel ) )

	#plots[ c ][ "a" ].set_xlim( 1.e4, 8.e6 )
	#plots[ c ][ "a2" ].set_xlim( 1.e4, 8.e6 )

if interactive:
	matplotlib.pyplot.ion()
	fig.show()
	matplotlib.pyplot.pause(-1)
else:
	fig.savefig( args.output_file, orientation = "landscape" )
