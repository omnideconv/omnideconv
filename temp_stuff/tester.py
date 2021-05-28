import anndata
import numpy as np
import pandas as pd
import autogenes as ag
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Interface for autogenes')
parser.add_argument('-sc',help='The path to the single cell .h5ad file')
parser.add_argument('-out',help='The path to the output pickle file',default="autogenes.pickle")
parser.add_argument('-ngen',type=int,help='Number of generations. The higher, the longer it takes', default=2)
parser.add_argument('-mode',help='In standard mode, the number of genes of a selection is allowed to vary arbitrarily. In fixed mode, the number of selected genes is fixed (using nfeatures)', default="standard")
parser.add_argument('-nfeatures',type=int,help='Number of genes to be selected in fixed mode', default=None)
parser.add_argument('-weights',type=int, nargs='+',help='Weights applied to the objectives. For the optimization, only the sign is relevant: 1 means to maximize the respective objective, -1 to minimize it and 0 means to ignore it. The weight supplied here will be the default weight for selection. There must be as many weights as there are objectives', default=[-1,1])
parser.add_argument('-objectives', nargs='+',help='The objectives to maximize or minimize. Must have the same length as weights. The default objectives (correlation, distance) can be referred to using strings. For custom objectives, a function has to be passed. For further details, refer to the respective tutorial.', default=["correlation","distance"])
parser.add_argument('-seed',type=int,help='Seed for random number generators', default=0)
parser.add_argument('-verbose',type=bool,help='If True, output a progress summary of the optimization (the current generation, size of the pareto front, min and max values of all objectives)', default=False)
parser.add_argument('-population_size',type=int,help='Size of every generation (mu parameter)', default=100)
parser.add_argument('-offspring_size',type=int,help='Number of individuals created in every generation (lambda parameter)', default=50)
parser.add_argument('-crossover_pb',type=float,help='Crossover probability', default=0.7)
parser.add_argument('-mutation_pb',type=float,help='Mutation probability', default=0.3)
parser.add_argument('-mutate_flip_pb',type=float,help='Mutation flipping probability (fixed mode)', default=1E-3)
parser.add_argument('-crossover_thres',type=int,help='Crossover threshold (standard mode)', default=1000)
parser.add_argument('-ind_standard_pb',type=float,help='Probability used to generate initial population in standard mode', default=0.1)
parser.add_argument('-plot_weights',help='the path to the single cell .h5ad file', default=0)
parser.add_argument('-plot_objectives',help='the path to the single cell .h5ad file', default=0)
parser.add_argument('-index',help='the path to the single cell .h5ad file', default=0)
parser.add_argument('-close_to',help='the path to the single cell .h5ad file', default=0)

#choices=[0, 1, 2],

args = parser.parse_args()
weights = tuple(args.weights)
objectives = tuple(args.objectives)

adata = sc.read_h5ad(args.sc)
ag.init(adata,celltype_key='label')
ag.optimize(ngen=args.ngen,weights = weights,objectives = objectives,verbose = args.verbose,nfeatures=args.nfeatures,seed=args.seed,mode=args.mode,
            population_size=args.population_size,offspring_size=args.offspring_size,crossover_pb=args.crossover_pb,mutation_pb=args.mutation_pb,
            mutate_flip_pb=args.mutate_flip_pb,crossover_thres=args.crossover_thres,ind_standard_pb=args.ind_standard_pb)
ag.plot(weights=(-1,0))
ag.save(args.out)



bulk_data = pd.read_csv('bulk_small.csv',index_col=0).transpose()
ag.load("savestate.pickle")
coef = ag.deconvolve(bulk_data, model='linear')
print(coef)
res = pd.DataFrame(coef,columns=ag.adata().obs_names,index=bulk_data.index)
print(res)
res.to_csv("out.csv",sep="\t")
