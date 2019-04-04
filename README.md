# inchicompare
A few functions for splitting, normalising and comparing InChi strings. Detects exact matches by InChi layer which is useful for detecting ions and stereoisomers. Can split a mixture InChi into single molecule InChi keys.

Functions:
* __split_inchi(inchi)__ where inchi is a str. Returns list of str.

* __get_prefix_dict()__ returns a dict that gives the long title for each of the prefixes (sorry, not very pythonic)

* __parse_inchi(inchi)__ where inchi is a str. Returns dict of all the layers and sublayers. The key of each dict item is the prefix.

* __inchi_str_compare(inchi1, inchi2)__ where inchi1 and inchi2 are str. Returns a dict of all the layers and sublayers that are different beteen the two InChis including when a layer is absent.

The InChi format is explained in [Heller *et al.*, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4486400/)
