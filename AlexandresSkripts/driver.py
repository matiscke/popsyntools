#!/usr/bin/env python

# Uncomplicated planete driver
# The number of the run is read from the "SGE_TASK_ID" environment variable
# This number is then matched with the corresonding line of simulation list file.

# The evironment variable is set by the Sun Grid Engine when starting array jobs,
# using qsub [options] -t start:end <program>
# This allows a one to one map with task in the Sun Grid Engine and simulation.
# A maximum number of concurrently running simulations can be set by adding "-tc <num>"

from __future__ import print_function

import os
import os.path
import sys
import shutil
import socket
import signal
import tempfile
import subprocess

# ----------
# PARAMETERS
# ----------
# Debugging of the driver itself
debug = False
# Run planete with a single call; needs a code which has r1785 included
oneshot = True
# Use a local temporary directory when running planete
tmp = True
# Synchronise the content of the local temporary directory with main storage after the dump file has been written
# Only if tmp is True
sync = True
# Write planete's output to "sortie" and "erreur" files.
# Only if tmp is False; always the case otherwise.
outputToFile = True
# Whether planete's dump file works correctly; needs a code which has r1828 and r1947 included
canRestart = True
# Path to the list of simulations
listPath = "simulation_list.dat"
# Path to the planete executable
codePath = "planete"
# Do not change unless you have a good reason
codeName = os.path.basename( codePath )
# Directory containing the template layout of simulations
templatePath = "template"
# File to update with parameters from simulation_list.dat
# Note: the driver DOES NOT set the initial value of the "lecture" parameter, so it MUST already have the correct value
# ("non" for oneshot = False, "all" for oneshot = True).
fileUpdate = [ "controles", "donnees", "donnees_MP" ]
# Additional files to remove when a simulation finished
fileRemove = [ "fort.10", "fort.11", "fort.12", "fort.13", "fort.14", "fort.15", "dumpfile.dat", "dumpfile_old.dat", codeName ]
# Remove files 
removeTemplateItems = True
# --------------
# END PARAMETERS
# --------------

# Signal handling
class ChildFinished( BaseException ):
	pass

def child_handler( signum, frame ):
	raise ChildFinished

# First, get our number
if os.getenv( "SGE_TASK_ID" ) and os.getenv( "SGE_TASK_ID" ) != "undefined":
	number = int( os.getenv( "SGE_TASK_ID" ) )
else:
	raise Exception( "Unable to determine driver number." )

# Now get the corresponding line from the list
listFile = open( listPath, "r" )
listID = 0
for listItem in listFile:
	listID += 1
	if listID == number:
		line = listItem[:-1]
		break
else:
	raise Exception( "Could not get line number {} from {}".format( number, listPath ) )

listFile.close()

if debug:
	print( "Got line from list: {}".format( line ) )

# Parse data
data = {
	"CDname": line[ 0 : 17 ],
	"CDnumber": line[ 3 : 17 ],
	"fgp": line[ 20 : 34 ],
	"diskM": line[ 37 : 51 ],
	"a_in": line[ 54 : 68 ],
	"a_out": line[ 71 : 85 ],
	"expo": line[ 88 : 102 ],
	"windM": line[ 105 : 119 ],
	"simName": line[ 119 : 136 ],
	"simNumber": line[ 122 : 136 ],
	"a_start": line[ 139 : 153 ],
	"t_start": line[ 156 : 170 ],
}

# And we need to invert fgp to fpg
# fgp: gas to dust ratio (fraction gaz poussiere in French)
# fpg: dust to gas ratio (fraction poussiere gaz in French)
data[ "fpg" ] = "{:22.16E}".format( 1. / float( data[ "fgp" ] ) )

if debug:
	print( "Got data from line:", data )

# Good, we can now create the directory
dirName = data[ "CDname" ] + "_" + data[ "simName" ]
baseDir = os.getcwd()

if debug:
	print( "Simulation directory:", dirName )

# The final directory where results should be
simDir = baseDir + "/" + dirName

# Do not restart a simulation that is already finished
try:
	if tmp and not sync:
		os.stat( simDir )
	else:
		os.stat( simDir + "/end" )
	finished = True
except:
	finished = False

if finished:
	print( "Simulation already finished" )
	sys.exit( 0 )

# Determine the working directory
if tmp:
	runDir = tempfile.gettempdir() + "/" + dirName
else:
	runDir = simDir

# Create the directory if it doesn't exist yet
try:
	os.stat( runDir )
	simDirExist = True
except:
	# If we are to run in a temporary directory and the final one already exists,
	# then copy it instead of the template (in the case of restart).
	if tmp:
		try:
			os.stat( simDir )
			simDirExist = True
		except:
			simDirExist = False
	else:
		simDirExist = False
	if simDirExist:
		shutil.copytree( simDir, runDir )
	else:
		shutil.copytree( templatePath, runDir )

try:
	os.stat( runDir + "/" + codeName )
except:
	shutil.copy( codePath, runDir + "/" + codeName )

# Check for dump file
try:
	os.stat( runDir + "/dumpfile.dat" )
	hasRestart = True
except:
	hasRestart = False

doRestart = canRestart and hasRestart

if debug and doRestart:
	print( "Restarting" )

# And now we can go inside the directory
os.chdir( runDir )

# Update contents of files
for name in fileUpdate:
	contents = open( name, "r" ).read()
	contents = contents.format( **data )
	open( name, "w" ).write( contents )

# Flush output so that it appears at the correct location
if tmp:
	print( "Switching to temporary directory" )
	print( "Host:",  socket.gethostname() )
	print( "Directory:", runDir )
	# Close the standard output and error and send them to files in the temporary directory.
	# In the case of the cluster, this allows to be completely independent of /home0.
	sys.stdout.close()
	sys.stderr.close()
	os.close( 1 )
	os.close( 2 )
	sys.stdout = open( "sortie", "a" )
	sys.stderr = open( "erreur", "a" )
	sub_out = sys.stdout
	sub_err = sys.stderr
elif outputToFile:
	sub_out = open( "sortie", "a" )
	sub_err = open( "erreur", "a" )
else:
	# Flush the output so that is gets in the correct order
	sys.stderr.flush()
	sys.stdout.flush()
	sub_out = None
	sub_err = None

if oneshot or doRestart:
	res = -1 # No initialisation step
else:
	res = subprocess.call( [ "./" + codeName ], stdout = sub_out, stderr = sub_err )

if res == 0 or doRestart:
	update = open( "controles", "r+" )
	update.seek( 0 )
	update.write( "oui" )
	update.close()

if res == 0 or res == -1:
	# A big mess just to be able to sync directory contents while waiting for a subprocess to complete
	watch = False
	if tmp and sync:
		try:
			import pyinotify
			class ProcessDumpEnd( pyinotify.ProcessEvent ):
				def process_IN_MOVED_TO( self, event ):
					if event.name == "dumpfile.dat":
						try:
							# shutil.copytree() does not work if the target directory already exists
							os.stat( simDir )
							shutil.rmtree( simDir )
						except:
							pass
						try:
							# Do not break on error
							shutil.copytree( runDir, simDir )
						except:
							pass

				def process_default( self, event ):
					pass
			watch = True
		except:
			pass

	if watch:
		# Watch for changes to the dump file
		wm = pyinotify.WatchManager()
		notifier = pyinotify.Notifier( wm, default_proc_fun = ProcessDumpEnd() )
		wm.add_watch( runDir, pyinotify.IN_MOVED_TO )

		# Add signal handler so that an exception is raised when a child terminates
		# This needs to be done *after* the pyinotify stuff, as it spawns processes.
		signal.signal( signal.SIGCHLD, child_handler )

		# Prepare
		proc = None
		res = None

		while res is None:
			try:
				if proc is None:
					# Now that we are catching the signal, start the process
					proc = subprocess.Popen( [ "./" + codeName ], stdout = sub_out, stderr = sub_err )

				notifier.loop()
			except ChildFinished:
				# A child process finished; check whether it is our one
				res = proc.poll()
	else:
		res = subprocess.call( [ "./" + codeName ], stdout = sub_out, stderr = sub_err )

if debug:
	print( "Program exit code:", res )

successful = res == 0

# Create an "end" file; this will serve in case the driver is called again to not restart the simulation.
if sub_out is None:
	# Status is unknown
	open( "end", "w" ).write( "UNKNOWN STATUS" )

	if debug:
		print( "Unknown status" )
else:
	# Copy the last line of "sortie" to "end"
	for out_line in open( "sortie" ):
		# This is just to get the last line...
		pass

	open( "end", "w" ).write( out_line )

	if debug:
		print( "Last line:", out_line )

	if successful:
		# Well, you would expect that planete exits with a non-zero status in case of error.
		# This is mostly not the case, so you need to check that the output state that the code finished successfully.
		successful = "TERMINE" in out_line
		if not successful:
			ret = 1

if debug:
	print( "Successfull:", successful )

if successful:
	# Simulation finished without error
	# Remove temporary files
	for name in fileRemove:
		try:
			os.remove( name )
		except:
			pass

	# Clean up from the common template stuff, unless new files were created inside 
	# (the case of emps_output for instance).
	if removeTemplateItems:
		baseTemplate = templatePath
		if baseTemplate[ 0 ] != "/":
			baseTemplate = baseDir + "/" + baseTemplate
		try:
			for itemName in os.listdir( baseTemplate ):
				if itemName in fileUpdate:
					continue
				itemPath = baseTemplate + "/" + itemName
				try:
					if os.path.isdir( itemPath ):
						for subitemName in os.listdir( itemPath ):
							try:
								os.remove( itemName + "/" + subitemName )
							except:
								pass
						os.rmdir( itemName )
					else:
						os.remove( itemName )
				except:
					pass
		except:
			pass

if tmp:
	# Move directory away from temporary storage.
	# In the case the original directory is unavailable (storage crash, NFS error),
	# it is moved directly to /tmp so that the results are not lost once the driver exists
	# (the Sun Grid Engines deleted the temporary directory at job's end).
	back = False
	try:
		os.chdir( baseDir )
		back = True
	except:
		if debug:
			print( "chdir() to", baseDir, "failed. Using temporary directory." )

		baseDir = "/tmp"
		simDir = baseDir + "/" + dirName
		os.chdir( baseDir )

	if back:
		sys.stdout.close()
		sys.stderr.close()
		sys.stdout = open( os.getenv( "SGE_STDOUT_PATH" ), "a" )
		sys.stderr = open( os.getenv( "SGE_STDERR_PATH" ), "a" )
		print( "Finished; back to directory" )
		try:
			os.stat( simDir )
			shutil.rmtree( simDir )
		except:
			pass

	shutil.copytree( runDir, simDir )
	shutil.rmtree( runDir )
else:
	# Close these files in case they are open.
	if sub_out:
		sub_out.close()
	if sub_err:
		sub_err.close()

print( "Exit code:", res )
