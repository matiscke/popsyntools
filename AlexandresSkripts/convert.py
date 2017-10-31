import sys
import planete_data

systems = planete_data.getSystems( sys.argv[ 1 ] )
planete_data.convertToTable( systems, sys.argv[ 2 ] )
