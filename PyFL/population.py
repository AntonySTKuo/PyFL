import numpy as np


class population():
	def __init__(self, g_index_list, num_individuals, mode = 'homogeneous', *args, **kwargs):
		self.pop = None
		self.ori_pop = None
		self.generation = 0
		self.avg_f = None
		self.var_f = None


		# random 1000 / designated starting genotypes
		if mode == 'homogeneous':
			starting_genotype = kwargs.get('starting_genotype_index', None)
			self.pop = [starting_genotype for i in range(num_individuals)]
			self.ori_pop = [starting_genotype for i in range(num_individuals)]
		elif mode == 'rand_homogeneous':
			start_genotype = np.random.choice(g_index_list)
			self.pop = [ start_genotype for i in range(num_individuals)]
			self.ori_pop = [ start_genotype for i in range(num_individuals)]			
		elif mode == 'rand_heterogeneous':
			num_starting_genotypes = kwargs.get('diversity', None)
			proportion_starting_genotypes = kwargs.get('proportion', None)
			start_genotypes = np.random.choice(g_index_list, num_starting_genotypes)
			self.pop = [ g for idx, g in enumerate(start_genotypes) for j in range(proportion_starting_genotypes[idx]) ]
			self.ori_pop = [ g for idx, g in enumerate(start_genotypes) for j in range(proportion_starting_genotypes[idx]) ]
		else:
			raise Exception('mode can only be homogeneous or heterogeneous!!')

	def __str__(self):
		return "It's a population with {} genotypes".format(np.unique(self.pop).size)

	__repr__ = __str__

