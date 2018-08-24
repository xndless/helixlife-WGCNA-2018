#!/usr/bin/python

from sys import argv

script,infile,nodes,outfile = argv

inf = open(infile,"r")
nodes = open(nodes,"r")
outf = open(outfile,"w")

array = []
for line in nodes:
	line = line.strip("\n")
	array.append(line)

for line in inf:
	line = line.strip("\n")
	data = line.split("\t")
	if data[0] in array and data[1] in array:
		out = line+"\n"
		outf.write(out)
	else:
		pass

inf.close()
nodes.close()
outf.close()
