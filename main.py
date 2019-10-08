import numpy as np
import matplotlib.pylab as plt
import PyFL.PyFL as PyFL
import PyFL.evolution as ev
import PyFL.population as pop

import graph_tool.all as gt

print("\n###########################################################")    
print("1. Setting up Fitness landscape")
fl = PyFL.FitnessLandscape('arti_SDR_union25')
#fl.graph.list_properties()

#fl_sub = PyFL.FitnessLandscape(fl.name+'_subgraph', graph=fl.subgraph(f_min=2.7))
#fl_sub.GraphDraw()
#gt.GraphWidget(fl_sub.graph, fl_sub.graph.vp.fitness)

#test = np.load('arti_gen_fit.npy')

print("\n###########################################################")
print("2. Setting up Population")
# SD = 'UAAGGAGGU'
# anti_SD = 'AUUCCUCCA'
genotypes = [x for x in fl.graph.vp.genotype]
starting_genotype = np.random.choice(len(genotypes))
print('Initial Genotype is', genotypes[starting_genotype])
#anti_SD_index = genotypes.index(anti_SD)
arti_pop = pop.population(len(genotypes), 1000, 'homogeneous', starting_genotype_index=starting_genotype)

print("\n###########################################################")
print("3. Start Evolution Simulation")
n = 10
PyFL.evolve(fl.graph, arti_pop, replicates=n)




# for _ in range(5):
# 	print("\n###########################################################")
# 	print("2. Setting up Population")
# 	SD = 'UAAGGAGGU'
# 	anti_SD = 'AUUCCUCCA'
# 	genotypes = [x for x in fl.graph.vp.genotype]
# 	starting_genotype = np.random.choice(len(genotypes))
# 	print('Initial Genotype is', genotypes[starting_genotype])
# 	#anti_SD_index = genotypes.index(anti_SD)
# 	arti_pop = pop.population(len(genotypes), 1000, 'homogeneous', starting_genotype_index=starting_genotype)

# 	print("\n###########################################################")
# 	print("3. Initialize a Evolution Simulation Instance")
# 	arti_evo = ev.evolution(arti_pop, fl)

# 	print("\n###########################################################")
# 	print("4. Start Evolution Simulation")
# 	n = 10
# 	gen_fit_plot = arti_evo.evolve(mutation_rate=0.01, num_gen=1000, replicates=n)
# 	np.save('arti_gen_fit_'+genotypes[starting_genotype]+'x10.npy', gen_fit_plot)



# print("\n###########################################################")
# print("5. Plot results")
# gen_fit_plot = np.load(fl.name[:4]+"_gen_fit.npy")
# for i in range(n):
#     plt.plot(gen_fit_plot[:,i])
# plt.show()




'''
from math import sqrt

g = gt.price_network(1500)
deg = g.degree_property_map("in")
sqrt2 = np.vectorize(sqrt)
deg.a = 4 * (sqrt2(deg.a) * 0.5 + 0.4)
ebet = gt.betweenness(g)[1]
ebet.a /= ebet.a.max() / 10.
eorder = ebet.copy()
eorder.a *= -1
pos = gt.sfdp_layout(g)
control = g.new_edge_property("vector<double>")
for e in g.edges():
    d = sqrt(sum((pos[e.source()].a - pos[e.target()].a) ** 2)) / 5
    control[e] = [0.3, d, 0.7, d]

gt.graph_draw(g, pos=pos, vertex_size=deg, vertex_fill_color=deg, vorder=deg, edge_color=ebet, eorder=eorder, edge_pen_width=ebet, edge_control_points=control, output="graph-draw.pdf")
'''

'''
g = gt.collection.data["netscience"]
g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
g.purge_vertices()
state = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)
t = gt.get_hierarchy_tree(state)[0]
tpos = pos = gt.radial_tree_layout(t, t.vertex(t.num_vertices() - 1), weighted=True)
cts = gt.get_hierarchy_control_points(g, t, tpos)
pos = g.own_property(tpos)
b = state.levels[0].b
shape = b.copy()
shape.a %= 14
gt.graph_draw(g, pos=pos, vertex_fill_color=b, vertex_shape=shape, edge_control_points=cts, edge_color=[0, 0, 0, 0.3], vertex_anchor=0, output="netscience_nested_mdl.pdf")
'''











print()


