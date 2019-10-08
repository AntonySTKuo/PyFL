import os, sys
import numpy as np
import pandas as pd
import graph_tool.all as gt
import matplotlib,pylab as plt
import concurrent.futures
from time import time
from tqdm import tqdm
from math import sqrt
dirname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dirname)
import util

#import statsmodels.api as sm

#gt.openmp_enabled()
#gt.openmp_get_num_threads()
#gt.openmp_set_num_threads(6)

##################################################
def initializeFL(name):
    exists = os.path.isfile(dirname+'/../data/'+name+'.xml.gz')
    if exists:
        print('\n> Graph Loading...')
        t1 = time()
        fl_graph = gt.load_graph(dirname+'/../data/'+name+'.xml.gz')
        t2 = time()
        print('Time usage:', round(t2-t1, 2), 'seconds')

    else:
        ### input
        print('\n> Fitness Graph Constructing...')
        df = pd.read_csv(dirname+'/../data/'+name+'.csv')
        len1 = len(df)
        df = df.dropna(axis=0, subset=['LogMean'])
        df.reset_index(drop=True, inplace=True)
        len2 = len(df)
        genotypes = df['SeqID']
        phenotypes = df['LogMean']
        f_dic = dict(zip(genotypes,phenotypes))
        index_array = list(range(len2))
        index_dic = dict(zip(genotypes,index_array))
        fl_graph = gt.Graph()

        ### create property maps
        eprop_fDiff = fl_graph.new_edge_property('double')
        for i in tqdm(range(len2)):
            mut_neigh = util.neighbors(genotypes[i])
            f1 = phenotypes[i]
            for seq in mut_neigh:
                f2 = f_dic.get(seq, np.nan)
                if (f2 > f1) and (not np.isnan(f2)):
                    e = fl_graph.add_edge(i, index_dic[seq])
                    eprop_fDiff[e] = f2 - f1
            v = fl_graph.vertex(i)
        ### turn property maps into 'internal' properties
        fl_graph.vertex_properties['fitness'] = fl_graph.new_vertex_property('double', vals=phenotypes)
        fl_graph.vertex_properties['genotype'] = fl_graph.new_vertex_property('string', vals=genotypes)
        fl_graph.edge_properties['fDiff'] = eprop_fDiff

        ### store the graph 
        print('\n> Graph Saving...')
        t1 = time()
        fl_graph.save(dirname+'/../data/'+name+'.xml.gz')
        t2 = time()
        print('Time usage:', round(t2-t1, 2), 'seconds')

    return fl_graph

##################################################
def SingleEvolve(graph, population, mutation_rate, num_gen):
    f_dic = graph.vertex_properties["fitness"]
    max_fit = max(f_dic.get_array())
    num_offspring = len(population.pop)
    
    ## Use set to store unique elements, plz don't use list... 
    diversity_explored = set()
    ## What is this for??
    substitutions = {0:np.unique(population.pop)}
    ## Store all the population trace, which is represented by a matrix of dim = (pop_size, num_gen, replicates)
    pop_record = np.zeros((num_offspring, num_gen))
    ## Mean fitness in each generation
    gen_fit_plot = np.zeros((num_gen))
    ## mut_ben, mut_del, mut_nut, sub_ben, sub_del, sub_nut
    mut_effect_record = np.zeros(6) 

    ## Copy from https://gitlab.com/devinbendixsen/innovation_ribozyme_fls/blob/master/Scripts/RiboEvolve.py
    # DATA_sims = {'final_fitness':[],'final_diversity':[],'unique_sequences_explored':[],'mut_ben_total':[],'mut_del_total':[],'mut_nut_total':[],'sub_ben_total':[],'sub_del_total':[],'sub_nut_total':[]}

    for j in tqdm(range(num_gen)):
        ## while loop for individual-based simulation
        offspring_count = 0
        while offspring_count < num_offspring:
            ## Step1: randomly choose a member in the parent population
            parent = np.random.choice(population.pop)
            ## determines if the parent is 'fit' enough to reproduce
            if np.random.rand() <= (f_dic[parent]/max_fit):
                offspring_count += 1
                ## determines if the offspring will have a mutation
                if np.random.rand() <= mutation_rate:
                    ## To do: make this choosing process a biased one...
                    offspring_genotype = int(np.random.choice([k for k in graph.vertex(parent).all_neighbors()])) ## randomly chooses one of the mut_options
                    pop_record[offspring_count-1,j] = offspring_genotype
                    ## determines if the mutation is beneficial
                    if (f_dic[offspring_genotype] - f_dic[parent]) >= (1/num_gen):
                        mut_effect_record[0] += 1
                    ## determines if the mutation is deleterious
                    elif (f_dic[parent] - f_dic[offspring_genotype]) >= (1/num_gen):
                        mut_effect_record[1] += 1
                    ## concludes that the mutation is neutral
                    else:
                        mut_effect_record[2] += 1
                ## if offspring is not mutated, copys the parent genotype into offspring
                else:
                    pop_record[offspring_count-1,j] = parent
        ## makes the offspring population the parent population for the next generation
        population.pop = pop_record[:,j]  
        
        ### update_status
        pop_f = [graph.vertex_properties["fitness"][i] for i in population.pop ]
        avg = np.mean(pop_f)
        var = np.var(pop_f)
        population.avg_f = avg  #determines the average population fitness
        population.var_f = var
        population.generation += 1

        x = np.unique(population.pop) # determines the number of unique genotypes in the new population
        for seq in x:
            if seq not in diversity_explored: # adds each unique sequence to the diversity explored for this replicate
                diversity_explored.add(seq)
            if seq not in substitutions.values() and (np.count_nonzero(population.pop == seq) >= (0.9*num_offspring)): #determines if there is a new substitution mutation that is 90% of the population
                substitutions[len(substitutions)] = seq

        gen_fit_plot[j] = population.avg_f #determines the average population fitness

    # DATA_sims['final_fitness'].append(population.avg_f)
    # DATA_sims['final_diversity'].append(len(np.unique(population.pop)))
    # DATA_sims['unique_sequences_explored'].append(len(diversity_explored))
    
    return gen_fit_plot


def evolve(graph, population, mutation_rate=0.01, num_gen=1000, replicates=3, saving_dir=""):
    print('\n> Loading settings for evolution simulation using Fisher-Wright Model')
    print('\n> This simulatio will run for {} generations with mutation rate = {} and repeat for {} times'.format(num_gen, mutation_rate, replicates))
    print('\n> Start evolving on the landscape')
    # Using 'concurrent.futures' to parallelize the different replicates
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(SingleEvolve, graph, population, mutation_rate, num_gen) for _ in range(replicates)]

        for f in concurrent.futures.as_completed(results):
            print('Done', len(f.result()))

##################################################
class FitnessLandscape():
    def __init__(self, name, graph=0, directed=True):
        self.name = name

        if graph == 0: self.graph = initializeFL(name)
        else: self.graph = graph

        if not directed: self.graph.set_directed(False)

        pks = []
        vls = []
        for v in self.graph.vertices():
            if v.in_degree() == 27:
                pks.append(v)
            if v.out_degree() == 27:
                vls.append(v)
        self.peaks = pks
        self.valleys = vls

        print('\n> Construction of Fitness Graph Completed !!')
        print('# of nodes:', self.graph.num_vertices())
        print('# of edges:', self.graph.num_edges())
        print('# of peaks:', len(self.peaks))
        print('# of valleys:', len(self.valleys))
    '''
    def __str__(self):
        return 'The name of fitness landscape is '+self.name+'.'

    __repr__ = __str__
    '''
    
    def subgraph(self, f_min=2.8, f_max=3.5, draw=False):
        func = lambda v: (self.graph.vp.fitness[v] > f_min) and (self.graph.vp.fitness[v] < f_max)
        subgraph = gt.GraphView(self.graph, vfilt=func, skip_properties=False)
        subgraph.purge_vertices()

        return subgraph
    
    def GraphDraw(self):
        F = self.graph.vp.fitness
        sqrt2 = np.vectorize(sqrt)
        F.a = sqrt2(10**((F.a-min(F.a))/(max(F.a)-min(F.a))+1))*2

        PG = gt.pagerank(self.graph)
        PG.a = PG.a*100

        gt.graph_draw(self.graph, vertex_size=F, vertex_fill_color=PG, output=dirname+'/../result/'+self.name+'.png')
        # output_size = (size, size)
        # vertex_size = v_size
        # edge_pen_width = E_PWDITH
    
    def FitnessDistribute(self, save=False):
        od_list = np.zeros(28)
        for v in self.graph.vertices():
            od_list[v.out_degree()] += 1
        
        plt.bar(range(len(od_list)), od_list)
        if save:
            plt.savefig(dirname+'/../result/'+'FitnessDistribute.png')
        else: plt.show()
        
        return od_list

    '''
    def PageRank():
        gt.pagerank(self.graph)
        return 0
    '''
    