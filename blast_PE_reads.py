### Boas Pucker ###
### pucker@uni-bonn.de ###

__version__ = "v0.11"

__usage__ = """
					python3 blast_PE_reads.py (""" + __version__ + """)
					--ref <REFERENCE_GENOME_SEQ_FILE>
					--gff <ANNOTATION_FILE>
					--reads1 <R1_READS_IN_FASTA_FILE>
					--reads2 <R2_READS_IN_FASTA_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--cpus <NUM_THREADS_FOR_BLAST>[10]
					--wordsize <BLASTN_WORD_SIZE>[10]
					--minsim <MIN_BLAST_HIT_SIMILARITY_%>[80]
					--minlen <MIN_BLAST_HIT_LENGTH>[25]
					--maxeval <MAX_EVALUE_CUTOFF>[0.001]
					--minscore <MIN_BLAST_SCORE>[30]
					"""

import os, sys, subprocess



def load_BLAST_results( blast1_result_file, blast2_result_file, min_sim, min_len, max_evalue, min_score ):
	"""! @brief load BLAST results passing filters """
	data = {}
	best_hits1 = {}
	failed_reads = {}
	# --- load BLAST hits from first file --- #
	with open( blast1_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				best_hits1[ parts[0] ]	#check if there is already a best hit for this read (ignore line)
			except KeyError:
				try:
					failed_reads[ parts[0] ]	#check if the best hit was already too bad (ignore line and following hits)
				except KeyError:
					good_hit_status = False
					if float( parts[2] ) > min_sim:
						if float( parts[3] ) > min_len:
							if float( parts[-2] ) < max_evalue:
								if float( parts[-1] ) > min_score:
									good_hit_status = True
					if good_hit_status:
						data.update( { parts[0]: { 	'chr1': parts[1], 'sstart1': parts[6], 'send1': parts[7], 'evalue1': parts[-2], 'score1': parts[-1],
																	'chr2': None, 'sstart2': None, 'send2': None, 'evalue2': None, 'score2': None } } )
						best_hits1.update( { parts[0]: None } )
					else:
						failed_reads.update( { parts[0]: None } )
			line = f.readline()
	best_hits1 = {}	#overwrite dictionary
	
	# --- load BLAST hits from second file --- #
	best_hits2 = {}
	with open( blast2_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				best_hits2[ parts[0] ]	#check if there is already a best hit for this read (ignore line)
			except KeyError:
				try:
					failed_reads[ parts[0] ]	#check if the best hit was already too bad (ignore line and following hits)
				except KeyError:
					good_hit_status = False
					if float( parts[2] ) > min_sim:
						if float( parts[3] ) > min_len:
							if float( parts[-2] ) < max_evalue:
								if float( parts[-1] ) > min_score:
									good_hit_status = True
					if good_hit_status:	#only work with good hits, because bad hits have been counted based on R1
						try:
							data[ parts[0] ]['chr2'] = parts[1]
							data[ parts[0] ]['sstart2'] = parts[6]
							data[ parts[0] ]['send2'] = parts[7]
							data[ parts[0] ]['evalue2'] = parts[-2]
							data[ parts[0] ]['score2'] = parts[-1]
						except KeyError:
							data.update( { parts[0]: { 	'chr2': parts[1], 'sstart2': parts[6], 'send2': parts[7], 'evalue2': parts[-2], 'score2': parts[-1],
																		'chr1': None, 'sstart1': None, 'send1': None, 'evalue1': None, 'score1': None } } )
							del failed_reads[ parts[0] ]	#remove read ID from failed reads
						best_hits2.update( { parts[0]: None } )
			line = f.readline()
	return data, failed_reads


def analyze_distribution_across_chromosomes( data ):
	"""! @brief anayses the distribution of BLAST hits across the different chromosomes """
	
	hits_per_chromosome = {}
	for key in sorted( list( data.keys() ) ):
		try:
			hits_per_chromosome[ data[key]['chr1'] ] += 1
		except KeyError:
			hits_per_chromosome.update( { data[key]['chr1']: 1 } )
			
	sys.stdout.write( str( hits_per_chromosome ) + "\n" )
	sys.stdout.flush()


def  load_gene_positions_from_GFF( gff_file ):
	"""! @brief load gene positions from GFF """
	
	gene_positions_per_chr = {}
	genes = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					genes.append( parts[-1] )
					try:
						gene_positions_per_chr[ parts[0] ].append( { 'id': parts[-1], 'start': int( parts[3] ), 'end': int( parts[4] ) } )
					except KeyError:
						gene_positions_per_chr.update( { parts[0]: [ { 'id': parts[-1], 'start': int( parts[3] ), 'end': int( parts[4] ) } ] } )
			line = f.readline()
	return gene_positions_per_chr, genes


def feature_counting( data, gene_positions, genes ):
	"""! @brief count mapped reads per gene """
	
	counts_per_gene = {}
	for gene in genes:
		counts_per_gene.update( { gene: 0 } )
	
	div_PE_mapping = 0	#number of reads not mapping to same chromosome/contig
	for key in list( data.keys() ):	#run this analysis for each BLAST hit
		if data[ key ]['chr1'] == data[ key ]['chr2']:	#this also excludes read pairs if only one mate is mapped
			try:
				put_genes = gene_positions[ data[ key ]['chr1'] ]	#get positions of all genes on the relevant chromosome
			except KeyError:
				sys.stdout.write( "WARNING: mapping to contig without annotated genes: " + str( data[ key ] ) + "\n" )
				sys.stdout.flush()
				put_genes = []
			qstart1, qend1, qstart2, qend2 = int( data[ key ]['sstart1'] ), int( data[ key ]['send1'] ), int( data[ key ]['sstart2'] ), int( data[ key ]['send2'] )
			for put in put_genes:
				if put['start'] < qend1:
					if put['end'] > qstart1:
						if put['start'] < qend2:
							if put['end'] > qstart2:
								counts_per_gene[ put['id'] ] += 1
		else:
			div_PE_mapping += 1
	return counts_per_gene


def main( arguments ):
	"""! @brief run everything """
	
	reference_file = arguments[ arguments.index('--ref') + 1 ]
	gff_file = arguments[ arguments.index('--gff') + 1 ]
	read1_fasta_input_file = arguments[ arguments.index('--reads1') + 1 ]
	read2_fasta_input_file = arguments[ arguments.index('--reads2') + 1 ]
	output_folder = arguments[ arguments.index('--out') + 1 ]
	
	if output_folder[-1] != '/':
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--cpus' in arguments:
		try:
			cpus = int( arguments[ arguments.index('--cpus') + 1 ] )
		except:
			cpus = 10
	else:
		cpus = 10
	
	if '--minsim' in arguments:
		try:
			min_sim = int( arguments[ arguments.index('--minsim') + 1 ] )
		except:
			min_sim = 80
	else:
		min_sim = 80
	
	if '--minlen' in arguments:
		try:
			min_len = int( arguments[ arguments.index('--minlen') + 1 ] )
		except:
			min_len = 25
	else:
		min_len = 25
	
	if '--maxeval' in arguments:
		try:
			max_evalue = float( arguments[ arguments.index('--maxeval') + 1 ] )
		except:
			max_evalue = 0.001
	else:
		max_evalue = 0.001
	
	if '--minscore' in arguments:
		try:
			min_score = int( arguments[ arguments.index('--minscore') + 1 ] )
		except:
			min_score = 30
	else:
		min_score = 30
	
	
	if '--wordsize' in arguments:
		try:
			wordsize = int( arguments[ arguments.index('--wordsize') + 1 ] )
		except:
			wordsize = 10
	else:
		wordsize = 10

	# --- run BLAST --- #
	blastdb = output_folder + "blastdb"
	cmd = "makeblastdb -in " + reference_file + " -out " + blastdb + " -dbtype nucl"
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()
	
	blast1_result_file = output_folder + "blast1_results.txt"
	if not os.path.isfile( blast1_result_file ):
		cmd = "blastn -query " + read1_fasta_input_file + " -db " + blastdb + " -out " + blast1_result_file + " -outfmt 6 -evalue 0.01 -num_threads " + str( cpus ) + " -word_size " + str( wordsize )
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	blast2_result_file = output_folder + "blast2_results.txt"
	if not os.path.isfile( blast2_result_file ):
		cmd = "blastn -query " + read2_fasta_input_file + " -db " + blastdb + " -out " + blast2_result_file + " -outfmt 6 -evalue 0.01 -num_threads " + str( cpus ) + " -word_size " + str( wordsize )
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	
	# --- load all data passing cutoff --- #
	data, failed_reads = load_BLAST_results( blast1_result_file, blast2_result_file, min_sim, min_len, max_evalue, min_score )
	sys.stdout.write( "number of mapped reads (at least one): "+ str( len( data.keys() ) ) + "\n" )	#at least one read resulted in hit
	sys.stdout.write( "number of unspecifically mapped or unmapped reads: " + str( len( failed_reads.keys() ) ) + "\n" )
	sys.stdout.flush()

	# --- generate output file with clean hits --- #
	clean_blast_hit_file = output_folder + "clean_blast_hits.txt"
	with open( clean_blast_hit_file, "w" ) as out:
		for key in sorted( list( data.keys() ) ):
			info = data[ key ]
			out.write( "\t".join( [ key, info['chr1'], info['sstart1'], info['send1'], info['evalue1'], info['score1'], info['sstart2'], info['send2'], info['evalue2'], info['score2'] ] ) + "\n" )
			
	# --- analyze clean hits over chromosomes --- #
	analyze_distribution_across_chromosomes( data )
	
	# -- load gene positions --- #
	gene_positions, genes = load_gene_positions_from_GFF( gff_file )
	sys.stdout.write( "number of identified genes: " + str( len( genes ) ) + "\n" )
	sys.stdout.flush()
	
	# --- count hits per gene --- #
	counts_per_gene_file = output_folder + "counts_per_gene.txt"
	reads_per_gene = feature_counting( data, gene_positions, genes )
	with open( counts_per_gene_file, "w" ) as out:
		for key in list( reads_per_gene.keys() ):
			if reads_per_gene[ key ] > 0:
				sys.stdout.write( key + "\t" + str( reads_per_gene[ key ] ) + "\n" )
				sys.stdout.flush()
			out.write( key + "\t" + str( reads_per_gene[ key ] ) + "\n" )


if '--ref' in sys.argv and '--gff' in sys.argv and '--reads1' in sys.argv and '--reads2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
