import numpy as np
import util
import graph_tool as gt
from tqdm import tqdm


class evolution():
	def __init__(self, pop, f_landscape):
		'''
		pop: a population object, including the info. of all genetic variants and num. of individuals
		f_landscape: a fitness landscape/graph object, including the info. of fitness and neighborhood structure
		'''
		self.population = pop
		self.fitness_landscape = f_landscape
		#self.generation = 0

		pop_f = [ self.fitness_landscape.graph.vertex_properties["fitness"][i] for i in self.population.pop ]
		avg = np.mean(pop_f)
		var = np.var(pop_f)
		self.population.avg_f = avg  #determines the average population fitness
		self.population.var_f = var


	def __str__(self):
		return "\n> The population has been evolved for {} generations. \n>> The average fitness = {} \n>> The starting genotype is {}".format(self.population.generation, self.population.avg_f, np.unique(self.population.pop))

	__repr__ = __str__


	def update_status(self):
		#avg, var = util.weighted_variance(self.f_landscape.fitness, self.pop.individual)
		pop_f = [ self.fitness_landscape.graph.vertex_properties["fitness"][i] for i in self.population.pop ]
		avg = np.mean(pop_f)
		var = np.var(pop_f)
		self.population.avg_f = avg  #determines the average population fitness
		self.population.var_f = var
		self.population.generation += 1


	def plot_allel_freq(self):
		pass

	def evolve(self, mutation_rate=0.01, num_gen=1000, replicates=3, saving_dir=""):
		print('\n> Loading settings for evolution simulation using Fisher-Wright Model')
		print('\n> This simulatio will run for {} generations with mutation rate = {} and repeat for {} times'.format(num_gen, mutation_rate, replicates))

		# make alias to the PropertyMap for 'variant_index to fitness' mapping
		# Note that PropertyMap is under graph-tool library, which is similar to python dictionary but not the same  
		f_dic = self.fitness_landscape.graph.vertex_properties["fitness"]

		# use set to store unque elements, plz don't use list... 
		diversity_explored = set()
		# What is this for??
		substitutions = {0:np.unique(self.population.pop)}

		max_fit = max(f_dic.get_array())
		num_offspring = len(self.population.pop)
		pop_size = num_offspring # just a alias
		# Allocate the momory to store all the population trace, which is represented by a matrix of dim = (pop_size, num_gen, replicates)
		pop_record = np.zeros((num_offspring, num_gen, replicates))
		# Copy from https://gitlab.com/devinbendixsen/innovation_ribozyme_fls/blob/master/Scripts/RiboEvolve.py
		DATA_sims = {'final_fitness':[],'final_diversity':[],'unique_sequences_explored':[],'mut_ben_total':[],'mut_del_total':[],'mut_nut_total':[],'sub_ben_total':[],'sub_del_total':[],'sub_nut_total':[]}
		gen_fit_plot = np.zeros((num_gen, replicates))



		print('\n> Start evolving on the landscape')
		# for loop for replicates, this loop could be parallelized
		for i in range(replicates):
			mut_effect_record = np.zeros(6) # mut_ben, mut_del, mut_nut, sub_ben, sub_del, sub_nut
			# for loop for serial evoluation, this loop cannot be parallelized 
			for j in tqdm(range(num_gen)):
				# while loop for individual-based simulation
				offspring_count = 0
				while offspring_count < num_offspring:
					# Step1: randomly choose a member in the parent population
					parent = np.random.choice(self.population.pop)

					if np.random.rand() <= (f_dic[parent]/max_fit): # determines if the parent is 'fit' enough to reproduce
						offspring_count += 1
						if np.random.rand() <= mutation_rate: # determines if the offspring will have a mutation
							# To do: make this choosing process a biased one...
							offspring_genotype = int(np.random.choice([k for k in self.fitness_landscape.graph.vertex(parent).all_neighbors()])) #randomly chooses one of the mut_options
							pop_record[offspring_count-1,j,i] = offspring_genotype
							if (f_dic[offspring_genotype] - f_dic[parent]) >= (1/num_gen): #determines if the mutation is beneficial
								mut_effect_record[0] += 1
							elif (f_dic[parent] - f_dic[offspring_genotype]) >= (1/num_gen): #determines if the mutation is deleterious
								mut_effect_record[1] += 1
							else: #concludes that the mutation is neutral
								mut_effect_record[2] += 1
						else:
							pop_record[offspring_count-1,j,i] = parent #if offspring is not mutated, copys the parent genotype into offspring

				self.population.pop = pop_record[:,j,i]  #makes the offspring population the parent population for the next generation
				self.update_status()
				x = np.unique(self.population.pop) # determines the number of unique genotypes in the new population

				for seq in x:
					if seq not in diversity_explored: # adds each unique sequence to the diversity explored for this replicate
						diversity_explored.add(seq)
					if seq not in substitutions.values() and (np.count_nonzero(self.population.pop == seq) >= (0.9*num_offspring)): #determines if there is a new substitution mutation that is 90% of the population
						substitutions[len(substitutions)] = seq

				gen_fit_plot[j,i] = self.population.avg_f #determines the average population fitness

			DATA_sims['final_fitness'].append(self.population.avg_f)
			DATA_sims['final_diversity'].append(len(np.unique(self.population.pop)))
			DATA_sims['unique_sequences_explored'].append(len(diversity_explored))



			# save to files
			#np.save(self.fitness_landscape.name[:4] + "_gen_fit", gen_fit_plot)

			# reset the population
			self.population.pop = self.population.ori_pop

		return gen_fit_plot



			# Problematic, not editted yet
		"""
			DATA_sims['mut_total'].append(mut_effect_record)
			for seq in substitutions:
				if seq != 0:
					if f_dic[substitutions[seq-1]] - f_dic[substitutions[seq]] <= (1/gen): #determines if the substitution is beneficial
						mut_effect_record[3] += 1
					elif f_dic[substitutions[seq-1]] - f_dic[substitutions[seq]] >= (1/gen): #determines if the substitution is deleterious
						mut_effect_record[4] += 1
					else:#concludes that the substituion is neutral
						mut_effect_record[5] +=1 
			
			DATA_sims['sub_ben_total'].append(sub_ben)
			DATA_sims['sub_del_total'].append(sub_del)
			DATA_sims['sub_nut_total'].append(sub_nut)
		"""
			#if ((i)%10==0):
			#	print('Replicates : ',i)
			

			#return diversity_explored, DATA_sims, gen_fit_plot, pop_record

			