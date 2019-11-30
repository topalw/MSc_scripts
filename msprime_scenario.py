#!/usr/bin/env python3
import msprime
import math

def simple_sim():
	names = []
	for i in range(1,11):
		names.append('main' + str(i))
	for i in range(1,11):
		names.append('CR' + str(i))
	for i in range(1,11):
		names.append('CY' + str(i))
	Nmain = 10000
	NCR = 4000
	NCY = 5000
	gen = 3		
	t1 = 5500	#2 x 5500 = 11000 yago
	t2 = 6000	#= 12000 years ago
	t3 = 19500	# 39000 years ago
	t4 = 20000	# 40000 
	growth_CY = 0.004605 	#growth reflects only starting size - final size in exp growth
	growth_CR = 0.004158
	m = 1e-4	# to force mainland to island
	r = 1e-7	# 1 per gen per chr
	length = 1e8	#chr 
	mu = 5e-8
	population_configurations = [
		msprime.PopulationConfiguration( sample_size = 20, initial_size = Nmain ),
		msprime.PopulationConfiguration( sample_size = 20, initial_size = NCR ),
		msprime.PopulationConfiguration( sample_size = 20, initial_size = NCY )
		]
	migration_matrix = [
	[ 0, 0, 0],
	[ 0, 0, 0],
	[ 0, 0, 0],
	]
	demographic_events = [
	#CY is growing 
	msprime.PopulationParametersChange(
		time = t1, initial_size = NCY, growth_rate = growth_CY, population_id = 2),
	msprime.MigrationRateChange(
		time = t1, rate = m, matrix_index = (0, 2)),
	msprime.MigrationRateChange(
		time = t1, rate = m, matrix_index = (2, 0)),
	#CY merges with mainland
	msprime.MassMigration(
		time = t2, source = 2, destination = 0, proportion = 1.0),
	msprime.MigrationRateChange(
		time = t2, rate = 0),
	#CR is growing
	msprime.PopulationParametersChange(
		time = t3, initial_size = NCR, growth_rate = growth_CR, population_id = 1),
	msprime.MigrationRateChange(
		time = t3, rate = m, matrix_index = (0, 1)),
	msprime.MigrationRateChange(
		time = t3, rate = m, matrix_index = (1, 0)),
	#CR merges with mainland
	msprime.MassMigration(
		time = t4, source = 1, destination = 0, proportion = 1.0),
	msprime.MigrationRateChange(
		time = t4, rate = 0)
	]
#	dd = msprime.DemographyDebugger(
#		population_configurations = population_configurations,
#		migration_matrix = migration_matrix,
#		demographic_events = demographic_events
#	)
#	dd.print_history()
	tree_sequence = msprime.simulate(
		population_configurations = population_configurations,
		migration_matrix = migration_matrix,
		demographic_events = demographic_events,
		length = length,
		mutation_rate = mu,
		recombination_rate = r
		)
	with open('out.vcf', 'w') as vcf_file:
		tree_sequence.write_vcf(vcf_file, ploidy = 2)
simple_sim()
