# inchicompare
A few functions for splitting, normalising and comparing InChi strings. Detects exact matches by InChi layer which is useful for detecting ions and stereoisomers. Can split a mixture InChi into single molecule InChi keys.

Functions:
* __split(inchi)__ where inchi is a str. Returns list of str.

* __get_prefixes()__ returns a dict that gives the long title for each of the prefixes (sorry, not very pythonic)

* __parse(inchi)__ where inchi is a str. Returns dict of all the layers and sublayers. The key of each dict item is the prefix.

* __compare(inchi1, inchi2)__ where inchi1 and inchi2 are str. Returns a dict of all the layers and sublayers that are different beteen the two InChis including when a layer is absent.

* __normalise(inchi)__ where inchi is as str. returns an InChI string that has been parsed and reformed by OpenBabel.

* __fp_compare(inchi1, inchi2)__ where inchi1 and inchi2 are str. Returns a number between 0 and 1 that represents the similarity between the two molecules (1.0 is identical). Works by using OpenBabel to generate fingerprints for the molecules and caluculate the Tanimoto coefficient between the two ([OpenBabel documentation](https://openbabel.org/wiki/Tutorial:Fingerprints)).

The InChi format is explained in [Heller *et al.*, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4486400/)
