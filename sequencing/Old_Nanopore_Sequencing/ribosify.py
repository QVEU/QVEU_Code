#ribosify.py : to make one RNA.
# Read the file and get the DNA string
import sys
infile=sys.argv[1]
outfile="RNA_"+infile
with open(infile,'r') as IF:
  lines=[line for line in IF]
with open(outfile,'w') as OF:
  for line in lines:
    if line.startswith(">"):
      OF.write(line)
    else:
      line=[char if char!="T" else "U" for char in line]
      OF.write("".join(line))
