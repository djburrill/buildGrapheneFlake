"""
buildGrapheneFlake.py

Build graphene flake with given radius.
"""

# Imports
import itertools as iter
import numpy as np
import argparse as arg

# Functions
def writeXYZ(fileName,atomCoords):
    '''
    Write XYZ for graphene file.
    '''

    # Variables
    numAtoms = len(atomCoords)

    # Write file
    with open(fileName,'w') as outFile:
        # Write header
        outFile.write(str(numAtoms) + '\n')
        outFile.write('Graphene flake \n')

        # Write atoms
        for atom in atomCoords:
            outStr = '{:<3s} {:>8.2f} {:>8.2f} {:>8.2f}'.format(atom[0],float(atom[1]),float(atom[2]),float(atom[3]))
            outFile.write(outStr+'\n')
            # outFile.write(' '.join(atom) + '\n')

def removeHydrogen(hAtomList,cAtomList):
    '''
    Remove extra hydrogens.

    NOTES
        - Extra hydrogens are defined as not having a carbon attached to them (no carbon with ~1.42 Ang).
    '''

    # Variables
    a = 1.42
    trustRadius = a*1.25
    saveList = []
    finalList = []

    # Iterate through hydrogen atoms
    for index,hAtom in enumerate(hAtomList):
        # Iterate through carbon atoms
        for cAtom in cAtomList:
            # Calculate Cartestian difference
            diffX = np.abs(hAtom[0]-cAtom[0])
            diffY = np.abs(hAtom[1]-cAtom[1])

            # If carbon attached, save
            if (diffX < trustRadius and diffY < trustRadius):
                saveList.append(index)

    # Generate finalList with saved hydrogen atoms
    for index,hAtom in enumerate(hAtomList):
        if (index in saveList):
            finalList.append(hAtom)

    return finalList

def removeDuplicates(atomList):
    '''
    Remove atoms at the same position.

    INPUT
        atomList: List of atoms in XYZ format.
            Ex. [['C',0.0,0.0,0.0],['H',1.0,0.0,0.0]]
    '''

    # Variables
    threshold = 0.5
    finalAtomList = []
    dupIndexList = []

    # Check for duplicates
    for index1 in range(len(atomList)):
        for index2 in range(index1,len(atomList)):
            # Do not count same atoms
            if (index2 == index1):
                continue

            # Calculate Cartestian difference
            diffX = np.abs(float(atomList[index1][1])-float(atomList[index2][1]))
            diffY = np.abs(float(atomList[index1][2])-float(atomList[index2][2]))

            # Check for similarity
            if (diffX <= threshold and diffY <= threshold):
                dupIndexList.append(index2)

    # Store all non-duplicates
    for index in range(len(atomList)):
        if (index not in dupIndexList):
            finalAtomList.append(atomList[index])

    return finalAtomList

def decodeCombination(comboTuple):
    '''
    Decode combination tuple into a translation vector.
    '''

    # Variables
    transLat = [0.0,0.0]
    a = 1.42
    coeff1 = np.sqrt(3)/2.0

    if (len(comboTuple) == 0):
        return transLat

    # Lattice vectors
    lat1 = [a*coeff1,1.5*a]
    lat2 = [2*coeff1*a,0]
    lat3 = [-a*coeff1,1.5*a]
    lat4 = [-2*coeff1*a,0]
    lat5 = [-a*coeff1,-1.5*a]
    lat6 = [a*coeff1,-1.5*a]
    latticeList = [lat1,lat2,lat3,lat4,lat5,lat6]

    # Apply lattice vectors
    for vec in comboTuple:
        transLat[0] += latticeList[int(vec)-1][0]
        transLat[1] += latticeList[int(vec)-1][1]

    return transLat

def buildCyclic(translationVec):
    '''
    Build carbon positions given a translation vector for the center.

    INPUT
        translationVec: 2d translation vector. [x,y]

    OUTPUT
        carbonList: List of atomic positions (2d) of carbon atoms.
    '''

    # Variables
    a = 1.42
    coeff1 = np.sqrt(3)/2.0

    # Carbon atom vectors
    carbon1 = [0.0+translationVec[0],a+translationVec[1]]
    carbon2 = [coeff1*a+translationVec[0],0.5*a+translationVec[1]]
    carbon3 = [coeff1*a+translationVec[0],-0.5*a+translationVec[1]]
    carbon4 = [translationVec[0],-a+translationVec[1]]
    carbon5 = [-coeff1*a+translationVec[0],-0.5*a+translationVec[1]]
    carbon6 = [-coeff1*a+translationVec[0],0.5*a+translationVec[1]]
    carbonList = [carbon1,carbon2,carbon3,carbon4,carbon5,carbon6]

    return carbonList

def buildGraphene(radius):
    '''
    Build graphene sheet.

    NOTES
        - A single hexane has radius 1.
    '''

    # Variables
    coordList = []
    carbonList = []
    hydrogenList = []
    hCoordList = []
    finalCoods = []

    print('Building...')

    # Build up graphene flake
    for rad in range(radius):
        print('\tRadius: '+str(rad+1))

        # Determine combinations of lattice vectors
        combos = iter.combinations_with_replacement('123456',rad)

        for element in combos:
            # Get translation vector
            transVec = decodeCombination(element)

            # Determine carbon coordinates
            carbonCoords = buildCyclic(transVec)
            coordList.extend(carbonCoords)

    # Add Hydrogen atoms
    rad = radius

    # Determine combinations of lattice vectors
    combos = iter.combinations_with_replacement('123456',rad)

    for element in combos:
        # Get translation vector
        transVec = decodeCombination(element)

        # Determine carbon coordinates
        hydrogenCoords = buildCyclic(transVec)
        hCoordList.extend(hydrogenCoords)

    # Remove duplicate carbon atoms
    #coordList = removeDuplicates(coordList)

    # Remove duplicate hydrogen atoms
    #hCoordList = removeDuplicates(hCoordList)

    # Remove extra hydrogens
    hCoordList = removeHydrogen(hCoordList,coordList)

    # Combine to XYZ format
    cFinalCoords = [['C',str(atom[0]),str(atom[1]),'0.0'] for atom in coordList]
    hFinalCoords = [['H',str(atom[0]),str(atom[1]),'0.0'] for atom in hCoordList]

    finalCoords = cFinalCoords + hFinalCoords

    # Remove duplicate carbon atoms
    finalCoords = removeDuplicates(finalCoords)

    return finalCoords

# Main
if (__name__ == '__main__'):
    # Read command line arguments
    parser = arg.ArgumentParser(description='Build finite-sized graphene flakes.')
    parser.add_argument('-r',
                        '--radius',
                        type=int,
                        default=3,
                        help='Radius of graphene sheet.')

    parser.add_argument('-f',
                        '--fileName',
                        type=str,
                        default='graphene_flake.xyz',
                        help='Name of output XYZ file.')

    args = parser.parse_args()

    # Build graphene sheet
    coords = buildGraphene(args.radius)

    # Save coordinates in XYZ
    writeXYZ(args.fileName,coords)
