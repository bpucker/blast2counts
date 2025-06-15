### Boas Pucker ###
### pucker@uni-bonn.de ###
### v0.11 ###


__usage__ = """
					python3 analyze_susp_fastqs.py
					--in <FASTQ.GZ_FILE>
					--out <FASTA_FILE>
					"""

__version__ = "v0.1"

import gzip, os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	oputput_file = arguments[ arguments.index('--out')+1 ]

	only_N_counter = 0
	other_counter = 0

	with open( oputput_file, "w" ) as out:
		with gzip.open( input_file, 'rt', encoding='utf-8' ) as f:
			line = f.readline()
			while line:
				seq = f.readline().strip()
				if len( seq ) == seq.count("N"):
					only_N_counter += 1
				else:
					other_counter += 1
					out.write( '>' + line + seq + "\n" )
				x = f.readline()
				qual = f.readline()
				line =f.readline()

	sys.stdout.write( "only N:" + str( only_N_counter ) + "\n" )
	sys.stdout.write( "total reads: " + str( only_N_counter+other_counter ) + "\n" )
	sys.stdout.write( str( only_N_counter / ( only_N_counter+other_counter ) ) + "\n" )
	sys.stdout.flush()

if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
