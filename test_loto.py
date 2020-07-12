from pyloto import pyloto
from pyloto import REC
from pyloto import RGD

import networkx as nx
from pandas import read_csv

df = read_csv("reference.tsv", header=None, sep="\t") #loading file
G = nx.from_pandas_edgelist(df, 0, 1, create_using=nx.DiGraph()) #converting file list into a digraph

#creating a pyloto object
myLoTo =  pyloto(G, nproc = 4)


#######################################################################################################
##
## uncomment this section for a toy use example
##
#######################################################################################################
##basic statics
#print("#TFs","n","#P","#N")
#print(myLoTo.number_TFs, myLoTo.number_genes, myLoTo.true_edges, myLoTo.false_edges)
#
##getting number of graphlets per type
#print("getting number of graphlets per type")
#print(myLoTo.count_graphlets())
#
##looking for how many times a node is forming graphlets
#for node in G.nodes:
#	print(node, myLoTo.count(node))
#
##looking for nodes not forming graphlets
#for node in G.nodes:
#	if(sum(myLoTo.count(node).values()) == 0):
#		print (node,"is not forming graphlets")
#
#print("printing the first ten type 1 graphlets")
#for i in range(10):
#	print (myLoTo.graphlets[1][i])
######################################################################################################	
	
#creating a segond loto object with a shuffled network to use REC and RGD
df = read_csv("inferred.tsv", header=None, sep="\t") #loading file
H = nx.from_pandas_edgelist(df, 0, 1, create_using=nx.DiGraph()) #converting file list into a digraph

#it is not neccessary to create a pyloro object to get the REC
#getting REC
print("computing REC")
global myREC
myREC = REC(myLoTo, H, nproc = 4)
#uncomment the following two lines to print rec for each graphlet
#for rec in myREC:
#	print(rec, myREC[rec])
	
#computing RGD
print("computing RGD for gadE")
myRGD = RGD(myREC,"gadE",nproc = 4)
print(myRGD)

print("computing RGD for multiples nodes")
myRGD = RGD(myREC,list(G.nodes),nproc = 4)#list to parse a list of string instead of nodeviews
for node in myRGD:
	print("RGD for node: ", node,"is",myRGD[node]) 

