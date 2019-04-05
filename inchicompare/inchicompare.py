#import openbabel as ob
import pybel
import re

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
    inchi1 = parse_inchi(inchi1)
    inchi2 = parse_inchi(inchi2)
    
    inchi1_keys = set(inchi1.keys())
    inchi2_keys = set(inchi2.keys())
    
    print(inchi1_keys)
    print(inchi2_keys)
    print(inchi1_keys.intersection(inchi2_keys))
    differences = {}
    
    #handle layers common in both inchis
    for key in inchi1_keys.intersection(inchi2_keys):
        if inchi1[key] != inchi1[key]:
            differences[key] = [inchi1[key], inchi2[key]]
    
    #handles layers present in one but not the other
    for key in inchi1_keys - inchi2_keys:
        differences[key] = [inchi1[key], None]
                            
    for key in inchi2_keys - inchi1_keys:
        differences[key] = [None, inchi2[key]]
            
    return differences

# These functions require OpenBabel and pybel

#Compare InChIs with OpenBabel molecular fingerprints
def fp_compare(inchi1, inchi2):
    mol1_fp = pybel.readstring('inchi',inchi1).calcfp('fp2')
    mol2_fp = pybel.readstring('inchi',inchi2).calcfp('fp2')
    
    return mol1_fp | mol2_fp

#Convert ot canonical SMILES and back again
def normalise(inchi):
    mol = pybel.readstring('inchi',inchi).write('can')
    return pybel.readstring('can', mol).write('inchi')[:-1]
