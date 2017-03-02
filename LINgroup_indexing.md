# The method of LINgroup indexing
## Introduction
LINgroup indexing is one of the methods adopted to reduce the number of ANI calculations for the discovery of the subject genome in the database with the highest ANI with the newly submitted genome. It takes the advantage of LINs already assigned to the genomes in the database, and LINgroup showing the clustering of genomically similar genomes at different thresholds corresponding to the different positions of LIN.## Pseudo code
```Python
# Variables used
previous_route = ""
N = Number of LIN positions
n = The current position
LIN_dictionary = {LIN:Genome_ID}
similarity = {Genome_ID:ANI}

# FUNCTIONS
def determine_current_step(previous_route, n, similarity):
	ANI_storage = {}
	subset_genomes = A subset of genomes fetched from database that all have the previous route
	unique_numbers = unique numbers assigned to the subset of genomes at position n
	for each_unique_number in unique_numbers:
		potential_route = previous_route + each_unique_number
		next_unique_numbers = unique numbers assigned at position n+1 for the genomes having potential_route
		for each_next_unique_numbers in next_unique_numbers:
			potential_subject_LIN = potential_route + each_next_unique_number + "0"*(N-n-1)
			potential_subject_Genome_ID = LIN_dictionary[potential_subject_LIN]
			if potential_subject_Genome_ID not in similarity:
				similarity[potential_subject_Genome_ID] = Calculated ANI
			if potential_route not in ANI_storage:
				ANI_storage[potential_route] = similarity[potential_subject_Genome_ID]
			else:
				ANI_storage[potential_route].append(similarity[potential_subject_Genome_ID])
		LIN_ANI_max_storage[potential_route] = max(ANI_storage[potential_route])
	update_route = max(LIN_ANI_max_storage)
	return update_route, n+1
	
def LINgroup_indexing():
	previous_route = ""
	N = Number of LIN positions
	n = The current position, initiates with 0
	LIN_dictionary = {LIN:Genome_ID}
	similarity = {} # {Genome_ID:ANI}
	while n < N:
		previous_route, n = determine_current_step(previous_route, n, similarity)
	# Now n = N, and the previous_route is N-1 in length
	ANI_storage = {}
	subset_genomes = A subset of genomes fetched from database that all have the previous route
	unique_numbers = unique numbers assigned to the subset of genomes at position n, which is the last position at this point
	for each_unique_numbers in unique_numbers:
		potential_subject_LIN = previous_route + each_unique_numbers
		potential_subject_Genome_ID = LIN_dictionary[potential_subject_LIN]
		ANI_storage[potential_subject_Genome_ID] = similarity[potential_subject_Genome_ID] # They are already calculated at position n-1
	best_Genome_ID = max(ANI_storage)
	return best_Genome_ID, LIN, ANI
		
```