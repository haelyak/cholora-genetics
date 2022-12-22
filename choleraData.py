import array

def loadMatrix(filename):
    """Load tab separated matrix. Store each row as an array of
    unsigned ints. Return list of arrays."""
    f = open(filename,"r")
    dataL=[]
    while True:
        s = f.readline()
        if s=="":
            break
        dataL.append(array.array('H',map(int,s.split())))
    return dataL


vcN16961geneInfoL = [x.rstrip() for x in open("vcN16961geneInfo.txt","r")]
vcPS15geneInfoL = [x.rstrip() for x in open("vcPS15geneInfo.txt","r")]
vc2740_80geneInfoL = [x.rstrip() for x in open("vc2740_80geneInfo.txt","r")]

vcN16961geneCoordL = [tuple(x.rstrip().split()) for x in open("vcN16961coords.tsv","r")]

vcN16961_vs_vcPS15 = loadMatrix("vcN16961_vs_vcPS15.tsv")
vcN16961_vs_vc2740_80 = loadMatrix("vcN16961_vs_vc2740_80.tsv")

