from Bio import Nexus
from Bio import Phylo
from radialtree import RadialTree

treestr1="(Taxon5,(Taxon3,(Taxon1,Taxon4)),Taxon2);"
treestr2="((2:1,3:1):1,(4:1,5:1):1,(0:1,1:1):1);"
#nexusIO = Nexus.Nexus.Nexus("#NEXUS\nBegin trees;\ntree 1 = (Taxon5:0.5,((Taxon3:0.5,(Taxon1:0.5,Taxon4:0.5):0.5):0.5,Taxon2:0.5):0.5);\nEnd;")
#nexusIO = Nexus.Nexus.Nexus("#NEXUS\nBegin trees;\ntree 1 = "+treestr1+"\nEnd;")
nexusIO = Nexus.Nexus.Nexus("#NEXUS\nBegin trees;\ntree 1 = "+treestr2+"\nEnd;")
radialTree = RadialTree()
#radialTree.render_png(nexusIO.trees[0],800,800, "c.png")
RadialTree.render_png(radialTree,nexusIO.trees[0],800,800,"c.png")
