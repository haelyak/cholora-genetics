from choleraData import *

def printIslands(coordsL,geneInfoL,geneCoordL):
    """For each potential island in coordsL print location and gene
    info."""
    for coords in coordsL:
        print("** Island")
        print("  chrom",geneCoordL[coords[0]][1],)
        print(geneCoordL[coords[0]][2]+"-"+geneCoordL[coords[1]-1][3])
        # print genes from coords[0] to coords[1], not including coords[1]
        for i in range(coords[0],coords[1]):
            print("  "+geneInfoL[i])
        print()

def cholera():
    hasHomologL = hasHomolog(vcN16961_vs_vcPS15,vcN16961_vs_vc2740_80,700)
    coordsL = islands(hasHomologL,12)
    printIslands(coordsL,vcN16961geneInfoL,vcN16961geneCoordL)


def hasHomolog(mat1,mat2,threshold):
    '''Inputs: two matrices of pairwise alignment scores, mat1 and mat2
    and a threshold number
    Outputs: a list of 1's and 0's, 
    If the maximum value in the corresponding row of both mat1 and mat2 is less than threshold, then it returns a 0 in that position, otherwise a 1'''
    homoList = []
    for x in range(len(mat1)):
        if max(mat1[x]) >= threshold or max(mat2[x]) >= threshold:
           homoList += [1]
        else:
            homoList += [0]
    return homoList

def islands(hasHomologL,minSize):
    '''finds islands of genes that are not homologs, and
    returns the indices of those genes'''
    counter = 0
    unsortedIslands = []
    for i in range(len(hasHomologL)):
        if hasHomologL[i] == 0:
            counter += 1
        elif counter >= minSize:
            unsortedIslands += [(counter, i - counter, i)]
            counter = 0
        else:
            counter = 0
    if counter >= minSize:
            unsortedIslands += [(counter, len(hasHomologL) - counter, len(hasHomologL))]
    unsortedIslands.sort()
    unsortedIslands.reverse()
    sortedIslands = []
    for i in range(len(unsortedIslands)):
        sortedIslands += [(unsortedIslands[i][1:])]

    return sortedIslands



