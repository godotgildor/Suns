import sys
import os
import openbabel
import glob

OUTPUT_FILE_SUFFIX = '''def find_aa_word(obj, bondAtoms):
    key = None
    selectionStatement = None
    
    # The bond key is of the form (resn, atom0, atom1) where the atom
    # names are in alphabetical order.
    atomNames = [bondAtoms[0]['name'], bondAtoms[1]['name']]
    atomNames.sort()
    bond_key = (bondAtoms[0]['resn'], atomNames[0], atomNames[1])
    
    if(bond_key in BOND_WORD_DICT):
        word = BOND_WORD_DICT[bond_key]
        # Now form a new key that has the current residue and the word.
        key = (obj, bondAtoms[0]['model'], bondAtoms[0]['segi'], bondAtoms[0]['chain'], bondAtoms[0]['resn'], bondAtoms[0]['resi'], word)
        if(bondAtoms[0]['segi'].strip() == ''):
            selectionStatement = '(object %s and model %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:2]) + list(key[3:6]) + [WORDS_DICT[word][key[4]]])
        else:
            selectionStatement = '(object %s and model %s and segi %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:6]) + [WORDS_DICT[word][key[4]]])
            
    return (key, selectionStatement)
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
def write_output(filename, wordsDict, bondWordDict):
    of = open(filename, 'w')
    of.write('WORDS_DICT = ' + str(wordsDict) + '\n\n')
    of.write('BOND_WORD_DICT = ' + str(bondWordDict) + '\n\n')
    of.write(OUTPUT_FILE_SUFFIX)
    of.close()

################################################################################
if(len(sys.argv) < 3):
    print 'Usage: python make_word_dicts.py <directory> <output file>'
    sys.exit(1)

fl = getFileList(sys.argv[1], '*.pdb')

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('pdb', 'pdb')

wordsDict = {}
bondWordDict = {}

for i, f in enumerate(fl):
    print 'Working on ' + f + ' ' + str(i+1) + ' of ' + str(len(fl)) + '.'
    
    (dir, fn) = os.path.split(f)
    word = os.path.split(dir)[-1]
    resn = os.path.splitext(fn)[0]
    bonds = get_bonds(obConversion, f)
    if(word not in wordsDict):
        wordsDict[word] = {}
    for bond in bonds:
        key = (resn, bond[0], bond[1])
        if(resn not in wordsDict[word]):
            wordsDict[word][resn] = '+'.join([bond[0], bond[1]])
        else:
            wordsDict[word][resn] += '+' + '+'.join([bond[0], bond[1]])
        bondWordDict[key] = word

write_output(sys.argv[2], wordsDict, bondWordDict)
open('wordsDict.py', 'w').write('WORDS_DICT = ' + str(wordsDict) + '\n\nBOND_WORD_DICT = ' + str(bondWordDict) + '\n\n')
