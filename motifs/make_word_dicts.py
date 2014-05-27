import sys
import os
import openbabel
import glob

OUTPUT_FILE_SUFFIX = '''def get_key(obj, atom, word):
    key = (obj, atom['model'], atom['segi'], atom['chain'], atom['resn'], atom['resi'], word)
    
    return key
    
def get_selection_statement(atom, key, atomNames):
    selectionStatement = None
    if(atom['segi'].strip() == ''):
        selectionStatement = '(object %s and model %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:2]) + list(key[3:6]) + [atomNames])
    else:
        selectionStatement = '(object %s and model %s and segi %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:6]) + [atomNames])
        
    return selectionStatement

def find_aa_word(obj, bondAtoms):
    (key, selectionStatement) = (None, None)
    
    # The bond key is of the form (resn, atom0, atom1) where the atom
    # names are in alphabetical order.
    atomNames = [bondAtoms[0]['name'], bondAtoms[1]['name']]
    atomNames.sort()
    bond_key = (bondAtoms[0]['resn'], atomNames[0], atomNames[1])
    
    if(bond_key in BOND_WORD_DICT):
        word = BOND_WORD_DICT[bond_key]
        # Now form a new key that has the current residue and the word.
        key = get_key(obj, bondAtoms[0], word)
        selectionStatement = get_selection_statement(bondAtoms[0], key, WORDS_DICT[word][key[4]])
            
    return (key, selectionStatement)

def find_suns_atom(obj, selectedAtom):
    (key, selectionStatement) = (None, None)

    atom_key = (selectedAtom['resn'], selectedAtom['name'])
    if(atom_key in ATOMS_DICT):
        key = get_key(obj, selectedAtom, ATOMS_DICT[atom_key])
        selectionStatement = get_selection_statement(selectedAtom, key, atom_key[1])
    
    return (key, selectionStatement)

def __init_plugin__(self):
    pass
'''

################################################################################
def getFileList(dirOfInterestOrFile, query='*', recursive=True):
    '''
    This function will get all files matching the query in the given directory.
    '''
    listOfFiles = []
    # Recursively get all files in this directory and
    # sub-directory.
    if(os.path.isdir(dirOfInterestOrFile)):
        if(recursive):
            for path, dirs, files in os.walk(dirOfInterestOrFile):
                listOfFiles += glob.glob( os.path.join(path, query) )
        else:
            listOfFiles += glob.glob( os.path.join(dirOfInterestOrFile, query) )
    elif(os.path.isfile(dirOfInterestOrFile)):
        for line in open(dirOfInterestOrFile):
            listOfFiles += [line.strip()]

    return listOfFiles

################################################################################
def get_bonds(obConversion, f):
    bonds = []
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, f)
    
    for i in range(mol.NumBonds()):
        bond = mol.GetBond(i)
        res = bond.GetBeginAtom().GetResidue()
        res2 = bond.GetEndAtom().GetResidue()
        bondAtoms = [res.GetAtomID(bond.GetBeginAtom()).strip(), res2.GetAtomID(bond.GetEndAtom()).strip()]
        bondAtoms.sort()
        bonds += [bondAtoms]
    
    return bonds

################################################################################
def write_output(filename, wordsDict, bondWordDict, atomsDict):
    of = open(filename, 'w')
    of.write('WORDS_DICT = ' + str(wordsDict) + '\n\n')
    of.write('BOND_WORD_DICT = ' + str(bondWordDict) + '\n\n')
    of.write('ATOMS_DICT = ' + str(atomsDict) + '\n\n')
    of.write(OUTPUT_FILE_SUFFIX)
    of.close()

def get_atom_key(f):
    key = None
    for l in open(f):
        if(l[0:6] == 'HETATM'):
            key = (l[17:20].strip(), l[12:16].strip())

    return key

################################################################################
if(len(sys.argv) < 3):
    print 'Usage: python make_word_dicts.py <directory> <output file>'
    sys.exit(1)

fl = getFileList(sys.argv[1], '*.pdb')

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('pdb', 'pdb')

wordsDict = {}
bondWordDict = {}
atomsDict = {}

for i, f in enumerate(fl):
    print 'Working on ' + f + ' ' + str(i+1) + ' of ' + str(len(fl)) + '.'
    
    (dir, fn) = os.path.split(f)
    word = os.path.split(dir)[-1]
    resn = os.path.splitext(fn)[0]
    bonds = get_bonds(obConversion, f)
    if(len(bonds) > 0):
        if(word not in wordsDict):
            wordsDict[word] = {}
        for bond in bonds:
            key = (resn, bond[0], bond[1])
            if(resn not in wordsDict[word]):
                wordsDict[word][resn] = '+'.join([bond[0], bond[1]])
            else:
                wordsDict[word][resn] += '+' + '+'.join([bond[0], bond[1]])
            bondWordDict[key] = word
    else:
        key = get_atom_key(f)
        atomsDict[key] = word

write_output(sys.argv[2], wordsDict, bondWordDict, atomsDict)
#open('wordsDict.py', 'w').write('WORDS_DICT = ' + str(wordsDict) + '\n\nBOND_WORD_DICT = ' + str(bondWordDict) + '\n\n')
