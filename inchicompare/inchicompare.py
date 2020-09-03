import re
import rdkit
import itertools as it
from collections import defaultdict

def split_inchi(inchi):
    try:
        mol = rdkit.Chem.MolFromInchi(inchi)
        r = [rdkit.Chem.MolToInchi(m) for m in rdkit.Chem.rdmolops.GetMolFrags(mol,asMols=True)]
    except:
        r = [inchi]
    return [i for i in r if not i is None]

def split(inchi):
    layers = inchi.split('/') #split into layers with the slash delimiter
    version = layers[0] #save the version header for later

    #deterine how many molcules by splitting the chemical formula
    split_inchis = []
    for i, data in enumerate(layers[1].split(".")):
        numerator = re.findall('^[0-9]+', data)
        if len(numerator) > 0:
            num = int(numerator[0])
            data = data[1:]
            for j in range(num):
                split_inchis.append([data])
        else:
            split_inchis.append([data])

    #add the rest of the layers
    for layer in layers[2:]:
        prefix = layer[0] #save the prefix as it isn't repeated for each molecule
        sep = '.' if prefix in ['m'] else ';'
        data_list = layer[1:].split(sep)

        #check for multipliers
        data_list_expanded = []
        
        for i in range(len(data_list)):
            m = re.search('^([0-9]+)\*(.*)', data_list[i])
            n = int(m.group(1)) if m else 1 
            data = m.group(2) if m else data_list[i]
            for j in range(n):
                data_list_expanded.append(data)
        
        data_list = data_list_expanded
        
        #add the data to the appropriate molecules
        for i, data in enumerate(data_list):
            if len(data) > 0: split_inchis[i].append(prefix+data)

    #join the data into a string
    split_inchi_strs = []
    for split_inchi in split_inchis:
        split_inchi_str = '/'.join([version]+split_inchi)
        split_inchi_strs.append(str(split_inchi_str))
    
    return split_inchi_strs

#I learnt Java first you know
def get_prefixes():
    prefix_dict = {'c': 'connectivity', 
                   'h': 'hydrogen', 
                   'q': 'charge', 
                   'p': 'protonation', 
                   'f': 'fixed H', 
                   'd': 'double bond stereoisomerism', 
                   't': 'tetrahedral stereoisomerism',
                   'i': 'isotopism',
                   'r': 'reconnection'
                  }
    #/m and /s are used as flags specifying absolute and relative stereochemistry and requesting the inverse arrangement of a double bond stereoisomer.
    #I've chosen not to handle these, but that might be a good idea in the future

    return prefix_dict

def parse(inchi):
    inchi = inchi.split('/')
    
    version = inchi[0]
    cf = inchi[1]
    layers = {}
    for layer in inchi[2:]:
        prefix = layer[0]
        data = layer[1:]
        layers[prefix] = data
    
    return version, cf, layers

def filter_layers(inchi, filter):
    version, cf, layers = parse(inchi)

    inchi_layer_order = ['c','h','q','f','d','t','i','r']

    filtered_layers = {}
    for k,v in layers.items():
        if k in filter:
            filtered_layers[k] = v

    ordered_filtered_layers = [version, cf]
    for layer_prefix in inchi_layer_order:
        if layer_prefix in filtered_layers.keys():
            ordered_filtered_layers.append(layer_prefix + filtered_layers[layer_prefix])

    return '/'.join(ordered_filtered_layers)

def compare(inchi1, inchi2):
    #parse inchis
    v1, cf1, inchi1 = parse(inchi1)
    v2, cf2, inchi2 = parse(inchi2)
    
    inchi1['cf'] = cf1
    inchi2['cf'] = cf2
    
    inchi1_keys = set(inchi1.keys())
    inchi2_keys = set(inchi2.keys())
    
    differences = {}
    
    #handle layers common in both inchis
    for key in inchi1_keys & inchi2_keys:
        if inchi1[key] != inchi2[key]:
            differences[key] = [inchi1[key], inchi2[key]]
    
    #handles layers present in one but not the other
    for key in inchi1_keys - inchi2_keys:
        differences[key] = [inchi1[key], None]
    
    for key in inchi2_keys - inchi1_keys:
        differences[key] = [None, inchi2[key]]
    
    return differences

def compare_split(inchi1, inchi2, 
                  filter_compare_types={'...matches...', '...matches a component of...', '...has a component which matches...', '...has a component which matches a component of...'},  # '...has a component which matches a component of...'
                  filter_layers = {'h','q','p','f','i'}):  
    
    results = defaultdict(lambda :defaultdict(dict))
    split_inchi1 = list(split_inchi(inchi1))
    split_inchi2 = list(split_inchi(inchi2))
    
    results['...matches...'][(inchi1,inchi2)] = compare(inchi1,inchi2)
    for i2 in split_inchi2:
        results['...matches a component of...'][(inchi1,i2)] = compare(inchi1,i2)
    for i1 in split_inchi1:
        results['...has a component which matches...'][(i1,inchi2)] = compare(i1,inchi2)
    for i1,i2 in it.product(split_inchi1,split_inchi2):
        results['...has a component which matches a component of...'][(i1,i2)] = compare(i1,i2)
    
    filtered_results = {k:v for k,v in results.items() if k in filter_compare_types}  # remove unauthorised compare_types
    
    # remove compares with unauthorised differences
    for compare_type, compare_data in filtered_results.items():
        for (i1,i2), differences in compare_data.items():
            filtered_results[compare_type] = {(i1,i2):differences for (i1,i2),differences in compare_data.items() if not (set(differences.keys()) - filter_layers)}
    
    filtered_results = {k:v for k,v in filtered_results.items() if v}  # remove empty compare_types
    
    return filtered_results

def join_inchis(inchi_list):
    if not inchi_list:
        return None
    
    r = rdkit.Chem.MolFromInchi(inchi_list[0])
    for i in inchi_list:
        m = rdkit.Chem.MolFromInchi(i)
        if m:
            r = rdkit.Chem.CombineMols(r,m)
    
    try:
        return rdkit.Chem.MolToInchi(r)
    except: return None

def compare_subset(inchi1, inchi2, 
                  filter_layers = {'h','q','p','f','i'}):
    
    split_inchi1 = list(split_inchi(inchi1))
    split_inchi2 = list(split_inchi(inchi2))
    
    matches = []
    for i,((i1,s1),(i2,s2)) in enumerate(it.product(enumerate(split_inchi1), enumerate(split_inchi2))):
        differences = compare(s1,s2)
        if not (set(differences.keys()) - filter_layers):
            matches.append((i,i1,i2,differences))
        
    matches = sorted(matches, key=lambda x:(len(x[3]),x[0]))
    joined_inchi1 = []
    joined_inchi2 = []
    for i,i1,i2,d in matches:
        if not (i1 in joined_inchi1 or i2 in joined_inchi2):
            joined_inchi1.append(i1)
            joined_inchi2.append(i2)
    
    p1 = (len(joined_inchi1),len(split_inchi1))
    p2 = (len(joined_inchi2),len(split_inchi2))
    
    joined_inchi1 = join_inchis([split_inchi1[i] for i in joined_inchi1])
    joined_inchi2 = join_inchis([split_inchi2[i] for i in joined_inchi2])
    
    if (joined_inchi1 is None) or (joined_inchi2 is None):
        return None
    
    differences = compare(joined_inchi1,joined_inchi2)
    
    return joined_inchi1, joined_inchi2, differences, p1, p2

def mol_consistent(mol1,mol2):
    # compare chiral centers
    chiral_compare = {}
    for atom in mol1.GetAtoms():
        k = atom.GetIdx()
        chiral_compare[k] = [atom.GetChiralTag(), None]
    for atom in mol2.GetAtoms():
        k = atom.GetIdx()
        if k in chiral_compare:
            chiral_compare[k][1] = atom.GetChiralTag()
        else:
            chiral_compare[k] = [None, atom.GetChiralTag()]
    
    for k,(a1,a2) in chiral_compare.items():
        if a1==a2:
            continue
        if {a1, a2} & {rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED}:
            continue
        return False, chiral_compare
    
    #compare double bonds
    doublebond_compare = {}
    for bond in mol1.GetBonds():
        k = tuple(sorted([bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]))
        doublebond_compare[k] = [bond.GetStereo(), None]
    for bond in mol2.GetBonds():
        k = tuple(sorted([bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]))
        if k in doublebond_compare:
            doublebond_compare[k][1] = bond.GetStereo()
        else:
            doublebond_compare[k] = [None, bond.GetStereo()]
    
    for (begin_atom, end_atom),(b1,b2) in doublebond_compare.items():
        if b1==b2:
            continue
        if {b1, b2} & {rdkit.Chem.rdchem.BondStereo.STEREONONE, rdkit.Chem.rdchem.BondStereo.STEREOANY}:
            continue
        return False, doublebond_compare
    
    return True, None

def compare_consistent(inchi1, inchi2, filter_layers={'h','f','p','q','i','t','b','m','s'}):
    r = compare_subset(inchi1, inchi2, filter_layers=filter_layers)
    if r:
        joined_inchi1, joined_inchi2, differences, p1, p2 = r
        consistent, evidence = mol_consistent(rdkit.Chem.MolFromInchi(joined_inchi1), rdkit.Chem.MolFromInchi(joined_inchi2))
        if consistent:
            return True, (differences, p1, p2)
        else:
            return consistent, evidence
    else:
        return False, None
    
def strip_inchi(inchi, exclude_inchis):
    def inchi_conn_layer(inchi):
        version, cf, layers = parse(inchi)
        if 'c' in layers:
            return f"InChI=1S/{cf}/c{layers['c']}"
        else:
            return f"InChI=1S/{cf}"
    
    filtered_inchis = []
    for i in split_inchi(inchi):
        c = inchi_conn_layer(i)
        if not c in exclude_inchis:
            filtered_inchis.append(i)
    return join_inchis(filtered_inchis)
