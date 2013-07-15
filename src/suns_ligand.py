import pybel

LIGAND_DEFINITIONS = {'phosphate': 'OP(=O)(O)O',
                      'sulfate': 'S(=O)(=O)(O)O',
                      'nitrate': 'N(=O)(=O)O',
                      'selmet': 'C[Se]C',
                      'acetate': 'C(=O)(O)C',
                      'adamantane': 'N[C@@]12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3',
                     }
                     
LIGAND_SELECTION = '_SUNS_LIGAND'

def find_atom_index(frp, ap):
    '''This function will loop through the atom lines in frp
       looking for a line which has the substring given in ap.
       Substring because the getPdbString function of pymol renumbers
       the atom numbers starting at 1, so two selections may have the
       same atoms with different atom numbers.'''
    for i, atom in enumerate(frp):
        if(atom.find(ap) >= 0):
            return i
    return -1
        

def find_atom_indices(fullResiduePdbstr, atomPdbstr):
    '''This function will take a pdb string, and look 
       for which lines each of the given atom strings 
       appear on.'''
    frp = fullResiduePdbstr.split('\n')
    indices = []
    for ap in atomPdbstr:
        index = find_atom_index(frp, ap)
        if(index >= 0):
            indices += [index]
    return indices

def get_pdb_string(fullResiduePdbstr, match):
    '''This function will simply extract the proper
       lines from a PDB string for the atoms that
       matched a motif.'''
    frp = fullResiduePdbstr.split('\n')
    subStr = [frp[i-1] for i in match]
    return '\n'.join(subStr)

def get_atom_names_statement(pdbstr):
    '''This function will generate a list of atom names
       from the given pdb string.'''
    atomNames = ''   
    for l in pdbstr.split('\n'):
        if( (l[0:6] == 'HETATM') or (l[0:4] == 'ATOM') ):
            atomNames += l[12:16].strip() + '+'

    # Remove the final # sign
    return atomNames[0:-1]
    
def look_for_ligand(fullResiduePdbstr, atomPdbstr):
    '''This function will look for the defined ligand
       motifs in this ligand residue.'''
    atomIndices = find_atom_indices(fullResiduePdbstr, atomPdbstr)
    if(len(atomIndices) != 2):
        return
    mol = pybel.readstring('pdb', fullResiduePdbstr)
    for ligand in LIGAND_DEFINITIONS:
        smarts = pybel.Smarts(LIGAND_DEFINITIONS[ligand])
        matches = smarts.findall(mol)
        for match in matches:
            # The Smarts search appears to use 1 based indexing.
            if(((atomIndices[0]+1) in match) and ((atomIndices[1]+1) in match)):
                retVal = {}
                retVal['motif'] = ligand
                retVal['pdbstr'] = get_pdb_string(fullResiduePdbstr, match)
                retVal['atomNames'] = get_atom_names_statement(retVal['pdbstr'])
                return retVal
    
    return {}

def find_ligand_word(cmd, obj, bondAtoms):
    '''This function will take the pymol cmd object, the name
       of the object selected, and the two bonded atoms selected,
       will get the pdb text of the full residue, look for any
       defined ligand motifs and return the first motif found that
       includes the two atoms in question.'''
    # Get the pdb string of the full residue.
    selectStatement = obj + ' and model ' + bondAtoms[0]['model'] + ' and chain ' + bondAtoms[0]['chain'] + ' and resn ' + bondAtoms[0]['resn'] + ' and resi ' + bondAtoms[0]['resi']
    if(bondAtoms[0]['segi'].strip() != ''):
        selectStatement += ' and segi ' + bondAtoms[0][1]
    cmd.select(LIGAND_SELECTION, selectStatement)
    fullResiduePdbstr = cmd.get_pdbstr(LIGAND_SELECTION)
    
    # Get the pdb string of the first atom
    selectStatement2 = selectStatement + ' and name ' + bondAtoms[0]['name']
    cmd.select(LIGAND_SELECTION, selectStatement2)
    atomPdbstr = [cmd.get_pdbstr(LIGAND_SELECTION).split('\n')[0][13:]]
    
    # Get the pdb string of the second atom
    selectStatement2 = selectStatement + ' and name ' + bondAtoms[1]['name']
    cmd.select(LIGAND_SELECTION, selectStatement2)
    atomPdbstr += [cmd.get_pdbstr(LIGAND_SELECTION).split('\n')[0][13:]]
    
    ligandMatch = look_for_ligand(fullResiduePdbstr, atomPdbstr)
    key = None
    if(ligandMatch != {}):
        selectStatement = '(' + selectStatement + ' and name ' + ligandMatch['atomNames'] + ')'
        key = (obj, bondAtoms[0]['model'], bondAtoms[0]['segi'], bondAtoms[0]['chain'], bondAtoms[0]['resn'], bondAtoms[0]['resi'], ligandMatch['motif'])
    else:
        selectStatement = None
                
    return (selectStatement, key)
    