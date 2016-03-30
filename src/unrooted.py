#!/usr/bin/env python

from Bio import Nexus
from Bio import Phylo
import sys, getopt
from radialtree import RadialTree

"Just returns the string content of the file"
def read_tree(infile):
    f = open(infile,"r")
    treestr = f.read()
    f.close()
    return(treestr)

def usage():
    return sys.argv[0]+' -i <inputfile> -o <outputfile> -t <pdf|png> [-W <width: default 800> -H <height: default 800>]'

def main(argv):
    goodtypes= ('pdf','png')
    inputfile = None
    outputfile = None
    width = 800
    height = 800
    otype= None
    try:
        opts, args = getopt.getopt(argv,"hi:o:t:W:H:",["ifile=","ofile=","type=","width=","height="])
    except getopt.GetoptError:
        print "Option error"
        print usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-t", "--type"):
            if(arg in goodtypes):
                otype = arg
            else:
                print "Bad type of output (-t)"
                print usage()
                sys.exit(2)
        elif opt in ("-W", "--width"):
            try: 
                width = int(arg)
            except ValueError:
                print "Bad width definition (-W)"
                print usage()
                sys.exit(2)
        elif opt in ("-H", "--height"):
            try: 
                height = int(arg)
            except ValueError:
                print "Bad height definition (-W)"
                print usage()
                sys.exit(2)

    if(inputfile is None or 
       outputfile is None or
       otype is None):
        print "Some arguments are missing"
        print usage()
        sys.exit(2)

    treestr = read_tree(inputfile)
    treeIO = Nexus.Nexus.Nexus("#NEXUS\nBegin trees;\ntree 1 = "+treestr+";\nEnd;")
    radialTree = RadialTree()
    if(otype == "pdf"):
        radialTree.render_pdf(treeIO.trees[0],width,height,outputfile)
    else:
        radialTree.render_png(treeIO.trees[0],width,height,outputfile)
    
if __name__ == "__main__":
   main(sys.argv[1:])
