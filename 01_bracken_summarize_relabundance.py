#!/usr/bin/env python

'''Summarize taxonomic relative abundance from kraken2/bracken 
for plots & analysis.

This script takes the following input:
- directory containing bracken tsv files

Taxa with relative abundance values are returned with following output:
- specified output file & path containting a single tsv file with 
taxa as rows, samples as columns, and relative abundance as cell values.

'''

import argparse, csv, glob, os
from collections import defaultdict
import pandas as pd


def parse_bracken_tsvs(input_directory, minimum_abundance, verbose):
	# generate dictionary of unique taxa & samples
	# initialize dictionary & list objects
	os.chdir(input_directory)				# change to input directory
	file_list = glob.glob("*.G.bracken")		# generate list of genus *.bracken files
	#file_list = glob.glob("*.S.bracken")		# generate list of species *.bracken files
	unique_taxa_list = []					# initialize taxa list
	MG_IDs = []
	id_mgid_taxon_reads_relab = defaultdict(list)
	mgid_taxa_list = defaultdict(list)


	print('Parsing bracken files...')

	for file in file_list:

		# get list of mg_ids
		# choose filename delimeter:
		#mg_id = file.rstrip().split('_')[0]
		mg_id = file.rstrip().split('.')[0]

		if mg_id not in MG_IDs:
			MG_IDs.append(mg_id)

		# parse file
		with open(file, 'r') as F:
			next(F)
			for line in F:
				X = line.rstrip().split('\t')
				taxon = X[0]	# species (or other tax level) name
				reads = X[5] 	# bracken revised assigned reads
				rel_ab = X[6]	# fraction of total reads assigned (relative abundance)

				# add taxa to unique_taxa_list if not listed
				if taxon not in unique_taxa_list:
					unique_taxa_list.append(taxon)

				# test if rel_ab value is greater than minimum
				if float(rel_ab) > float(minimum_abundance):

					# write mg_id, taxa, reads, relab to dictionary
					id_mgid_taxon_reads_relab[f"{mg_id}_{taxon}"] = {'mg_id' : mg_id,
						'taxon' : taxon, 'reads' : reads, 'rel_ab' : rel_ab}

					# add taxa to list of taxa for an MG
					mgid_taxa_list[mg_id].append(taxon)

	return MG_IDs, unique_taxa_list, id_mgid_taxon_reads_relab, mgid_taxa_list


def generate_relab_matrix(unique_taxa_list, MG_IDs, 
		id_mgid_taxon_reads_relab, mgid_taxa_list, reads_please, verbose):
	print('Generating matrices...')
	
	# initialize objects
	matrix = defaultdict(list)
	r_matrix = defaultdict(list)

	# output counts if verbose
	if verbose:
		print(f"{len(MG_IDs)} MG IDs | {len(unique_taxa_list)} unique taxa parsed")

	# create relative abundance matrix
	for MG in MG_IDs:

		#add MG taxa rel_ab values to matrix
		if verbose:
			print(f"   Updating matrix with {MG} relative abundance values...")
		for taxon in unique_taxa_list:
			if taxon in mgid_taxa_list[MG]:
				matrix[MG].append(id_mgid_taxon_reads_relab[f"{MG}_{taxon}"]['rel_ab'])
				# if verbose:
				# 	unique_key=f"{MG}_{taxon}"
				# 	print(f"	 {taxon} (relative abundance {id_mgid_taxon_reads_relab[unique_key]['rel_ab']}) added to matrix")
			else:
				matrix[MG].append(0)
				# if verbose:
				# 	print(f"	 {taxon} (relative abundance 0) added to matrix")		

		# add MG taxa read count values to r_matrix
		#	THIS OPTION NOT YET WORKING
		if reads_please:
			if verbose:
				print(f"   Updating matrix with {MG} read count values...")
			for taxon in unique_taxa_list:
				if taxon in mgid_taxa_list[MG]:
					r_matrix[MG].append(id_mgid_taxon_reads_relab[f"{MG}_{taxon}"]['reads'])
					if verbose:
						unique_key=f"{MG}_{taxon}"
						print("unique key", unique_key)
						print(f"	 {taxon} (read count {id_mgid_taxon_reads_relab[unique_key]['reads']}) added to matrix")
				else:
					r_matrix[MG].append(0)
					if verbose:
						print(f"	 {taxon} (read count 0) added to matrix")		

	print('')
	print('Converting completed matrices to dataframe...')

	relab_matrix = pd.DataFrame(matrix, index = unique_taxa_list)
	relab_matrix = relab_matrix.loc[(relab_matrix!=0).any(1)]
	relab_matrix.sort_index(inplace=True)
	relab_matrix.sort_index(axis=1, inplace=True)

	if reads_please:
		reads_matrix = pd.DataFrame(r_matrix, index = unique_taxa_list)
		reads_matrix = reads_matrix.loc[(reads_matrix!=0).any(1)]
		reads_matrix.sort_index(inplace=True)
		reads_matrix.sort_index(axis=1, inplace=True)	

	return relab_matrix, reads_matrix

def main():
	# configure argparse arguments & pass to functions.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-b', '--bracken_dir',
		help = 'Please specify directory path containing bracken output files.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-s', '--suffix',
		help = 'Please specify bracken file suffix (e.g. *.bracken).',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify output file path (& optional prefix).',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-m', '--minimum_abundance',
		help = 'Specify minimum abundance to filter taxa to include. Default 0.001',
		metavar = '',
		default=0.001
		)
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	parser.add_argument(
		'-r', '--reads_please',
		help = 'Toggle to produce reads_matrix from bracken output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	if args['verbose']:
		print('')
		print(f"Bracken directory: {args['bracken_dir']} | suffix: {args['suffix']}")
		print('')

	MG_IDs, unique_taxa_list, id_mgid_taxon_reads_relab, mgid_taxa_list= parse_bracken_tsvs(args['bracken_dir'], 
		args['minimum_abundance'], args['verbose'])
	relab_matrix, reads_matrix= generate_relab_matrix(unique_taxa_list, MG_IDs, id_mgid_taxon_reads_relab,
		mgid_taxa_list, args['bracken_dir'], args['verbose'])

#write output tsv files
	relab_matrix.to_csv(f"{args['output']}/relab_matrix.tsv", sep='\t')

	if args['reads_please']:
		reads_matrix.to_csv(f"{args['output']}/reads_matrix.tsv", sep='\t')

	# print('MG_IDs:')
	# print(MG_IDs)
	# print('')
	# print('unique_taxa_list:')
	# print(unique_taxa_list)
	# print('')
	# print('id_mgid_taxon_reads_relab:')
	# print(id_mgid_taxon_reads_relab)
	# print('')
	# print('mgid_taxa_list:')
	# print(mgid_taxa_list)
	# print('')
	# print('relab_matrix:')
	# print(relab_matrix)
	# print('')
	# print('reads_matrix:')
	# print(reads_matrix)
	# print('')

	if args['verbose']:
		#print('')
		#print('bracken dictionary:')
		#print(bracken_dict)
		# print('')
		# print('relab dictionary:')
		# print(relab_matrix)
		# print('')
		# print('Unique taxa list:')
		# print(unique_taxa_list)
		# print(len(unique_taxa_list), 'unique taxa detected.')
		print('')
		print(len(unique_taxa_list),'Taxa with bracken relative abundance values were identified.')

	print('Looks like everything completed!')
	print(f"Output files written to: {args['output']}")

if __name__ == "__main__":
	main()