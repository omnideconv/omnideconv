import anndata
import numpy as np
import pandas as pd
import autogenes as ag
import scanpy as sc
import argparse


parser = argparse.ArgumentParser(description='Interface for autogenes - Deconvolution')
parser.add_argument('-bulk',help='The path to the bulk data csv')
parser.add_argument('-input',help='The path to the input pickle file',default="autogenes.pickle")
parser.add_argument('-out',help='The path to the output deconvolution table',default="autogenes_out.tsv")

parser.add_argument('-model',help='Choose a regression model. Available options: NuSVR, non-negative least squares and linear model', default='nusvr')

parser.add_argument('-nu',type=float,help='Parameter for the nusvr deconvolution', default=0.5)
parser.add_argument('-C',type=float,help='Parameter for the nusvr deconvolution', default=0.5)
parser.add_argument('-kernel',help='Parameter for the nusvr deconvolution', default='linear')
parser.add_argument('-degree',type=int,help='Parameter for the nusvr deconvolution', default=3)
parser.add_argument('-gamma',help='Parameter for the nusvr deconvolution', default='scale')
parser.add_argument('-coef0',type=float,help='Parameter for the nusvr deconvolution', default=0.0)
parser.add_argument('-shrinking',type=bool,help='Parameter for the nusvr deconvolution', default=True)
parser.add_argument('-tol',type=float,help='Parameter for the nusvr deconvolution', default=1E-3)
parser.add_argument('-cache_size',type=int,help='Parameter for the nusvr deconvolution', default=200)
parser.add_argument('-verbose',type=bool,help='Parameter for the nusvr deconvolution', default=False)
parser.add_argument('-max_iter',type=int,help='Parameter for the nusvr deconvolution', default=-1)

parser.add_argument('-weights',nargs='+',help='Weights with which to weight the objective values. For example, (-1,2) will minimize the first objective and maximize the the second (with higher weight).', default=None)
parser.add_argument('-index',nargs='+',help='If one int is passed, return pareto[index] If two ints are passed, the first is an objective (0 for the first). The second is the nth element if the solutions have been sorted by the objective in ascending order. For example, (0,1) will return the solution that has the second-lowest value in the first objective. (1,-1) will return the solution with the highest value in the second objective', default=None)
parser.add_argument('-close_to',nargs='+',help='Select the solution whose objective value is closest to a certain value. Assumes (objective,value). For example, (0,100) will select the solution whose value for the first objective is closest to 100.', default=None)
parser.add_argument('-select',type=bool,help='Whether to select a specific solution', default=False)

args = parser.parse_args()

if args.weights is None:
    weights = None
else:
    weights = tuple(args.weights)
if len(args.index) > 1:
    args.index = tuple(args.index)
args.close_to = tuple(args.close_to)


bulk_data = pd.read_csv(args.bulk,index_col=0).transpose()
ag.load(args.input)
if args.select:
    #Copy and key_added parameters are not needed, thus not supported
    ag.select(weights=weights,index=args.index,close_to=args.close_to)
coef = None
if args.model == "nusvr":
    coef = ag.deconvolve(bulk_data, model=args.model, nu=args.nu, C=C, kernel=args.kernel, degree=args.degree,
                         gamma=args.gamma, coef0=args.coef0, shrinking=args.shrinking, tol=args.tol,
                         cache_size=args.cache_size, verbose=args.verbose, max_iter=args.max_iter)
elif args.model == "nnls" or args.model == "linear":
    coef = ag.deconvolve(bulk_data, model=args.model)

if coef == None:
    print("Error, model "+args.model+" not found")
    exit(1)

res = pd.DataFrame(coef,columns=ag.adata().obs_names,index=bulk_data.index)
res.to_csv(args.out,sep="\t")