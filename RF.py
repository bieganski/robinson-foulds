#!/usr/bin/env python
# coding: utf-8

# In[197]:


# homeo2.fa

from Bio.Align.Applications import ClustalwCommandline

clust = ClustalwCommandline("clustalw", infile="homeo2.fa")
a = clust()


# In[24]:


# utworzone zostaly pliki .dnd i .aln
from glob import glob
dnd = glob("*dnd")[0]
print(dnd)


# In[22]:


from Bio import Phylo

t1 = Phylo.read('homeo2.dnd', "newick")
print(type(t1))


# In[30]:


from Bio.Align.Applications import MuscleCommandline
muscle = MuscleCommandline(input="homeo2.fa")
out, err = muscle() # out w formacie fasta
print(glob("*"))


# In[49]:


from Bio import AlignIO
from io import StringIO

align_t2 = AlignIO.read(StringIO(out), "fasta")

from Bio.Phylo.TreeConstruction import DistanceCalculator

calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(align_t2)

# dm to nowoobliczona macierz podobieństwa


# In[55]:


# obliczenie dm trwa długo, zapisuję ją zatem na dysk

import pickle
pickle.dump(dm, open("dm.pickle", "wb"))


# In[60]:


# tworzymy drugie drzewo filogenetyczne algorytmem NJ z macierzą dm

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# dm = pickle.load( open( "dm.pickle", "rb" ) )

constructor = DistanceTreeConstructor()
t2 = constructor.nj(dm)

print(type(t2))


# In[61]:


pickle.dump(t1, open("t1.pickle", "wb"))
pickle.dump(t2, open("t2.pickle", "wb"))
# t1 = pickle.load( open( "t1.pickle", "rb" ) )
# t2 = pickle.load( open( "t2.pickle", "rb" ) )


# In[212]:


def napraw_nazwy(name):
    return name.replace("(", "_").replace(")", "_")

terms = t2.get_terminals()
for el in terms:
    el.name = napraw_nazwy(el.name)

s1 = frozenset([cl.name for cl in t1.get_terminals()])
s2 = frozenset([cl.name for cl in t2.get_terminals()])
assert(len(s1.intersection(s2)) == 839)


# In[213]:


m2 = dict() # mapa frozenset(A, B) -> clade z T2, 
           # gdzie type(A) == type(B) == frozenset, gdzie A, B tworzą podział liści drzewa T2

for cl in t2.get_nonterminals():
    termsA = frozenset([str(el) for el in cl.get_terminals()])
    termsB = s2 - termsA
    assert(len(termsA) + len(termsB) == 839)
    m2[frozenset([termsA, termsB])] = cl


# In[214]:


m1 = dict() # mapa frozenset(A, B) -> clade z T1,
           # gdzie type(A) == type(B) == frozenset, gdzie A, B tworzą podział liści drzewa T1

for cl in t1.get_nonterminals():
    termsA = frozenset([str(el) for el in cl.get_terminals()])
    termsB = s1 - termsA
    assert(len(termsA) + len(termsB) == 839)
    m1[frozenset([termsA, termsB])] = cl


# In[216]:


res = 0

def odNonuj(cos):
    return 0 if cos is None else cos

for podzial in m1:
    cl1 = m1[podzial] # istnieje
    cl2 = m2.get(podzial) # moze byc none
    if cl2 is None:
        res += (odNonuj(cl1.branch_length) ** 2)
    else:
        res += ((odNonuj(cl1.branch_length) - odNonuj(cl2.branch_length)) ** 2)

for podzial in m2:
    cl2 = m2[podzial] # istnieje
    cl1 = m1.get(podzial) # moze byc none
    if cl1 is None:
        res += (odNonuj(cl2.branch_length) ** 2)
    else:
        pass # to już było policzone

print("RF(T1, T2) = {}".format(res))

