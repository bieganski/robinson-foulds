from Bio import Phylo
from Bio import AlignIO
from io import StringIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

#tworzę drzewo filogenetyczne za pomocą algorytmu ClustalW
cline = ClustalwCommandline("clustalw", infile="homeo2.fa")
cline()

#tworzę drzewo filogenetyczne za pomocą algorytmu NJ na macierzy podobieństwa obliczonej 
#przy pomocy macierzy BLOSUM62 na podstawie multiuliniowienia otrzymanego przy pomocy algorytmu muscle
clustal_tree = Phylo.read("homeo2.dnd", "newick")

muscle_cline = MuscleCommandline(input="homeo2.fa")
stdout, stderr = muscle_cline()
print(stdout)
aln = AlignIO.read(StringIO(stdout), "fasta")

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
Phylo.draw_ascii(tree)
print(tree.count_terminals())

for clade in tree.get_terminals():
    print(clade.name)
 
# Tworzę listę par: dla każdego poddrzewa mam informację jakie liście znajdują się w danym poddrzewie i jaka 
# jest długość krawędzi łączącej ich z resztą
listMuscle = []
    
for clade in tree.find_clades(): 
    print( clade.name)
    a = clade.get_terminals()
    y = tree.get_path(clade)
    print(clade.branch_length)
    if len(y) > 0:
        #t = tree.distance(a, y[len(y) - 1])
        listMuscle.append((a, clade.branch_length))
    
for i in range(len(listMuscle)):
    print("%%%%%%%%%%%%%%%%%%%%%%")
    print(listMuscle[i][1])
    for clade in listMuscle[i][0]:
        print(clade.name)
    
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# Tworzę taką samą listę par dla drugiego drzewa
listClustal = []
    
for clade in tree.find_clades(): 
    print( clade.name)
    a = clade.get_terminals()
    y = tree.get_path(clade)
    print(clade.branch_length)
    if len(y) > 0:
        #t = tree.distance(a, y[len(y) - 1])
        listClustal.append((a, clade.branch_length))
    
for i in range(len(listClustal)):
    print("%%%%%%%%%%%%%%%%%%%%%%")
    print(listClustal[i][1])
    for clade in listClustal[i][0]:
        print(clade.name)
        
        
wynik = 0
# Przechodzę przez listy i porównuję każdy element z pierwszej listy z każdym z drugiej. 
# Jeśli przecięcie liści jest pełne lub zerowe to znaczy, że natrafiłam na takie same podzbiory i 
# biorę różnicę zapisanych dla nich odległości od reszty drzewa. Jeśli dla jakiegoś elementu z 
# listy nie znalazłam odpowiednika to dodaję do wyniku jego odległość do kwadratu
for i in range(len(listMuscle)):
    found = False
    for j in range(len(listClustal)):
        if found == False:
            interSection = set(listMuscle[i][0]).intersection(listClustal[j][0])
            if len(interSection) == 0 or ((len(interSection) == len(listMuscle[i][0]) and (len(interSection) == len(listlistClustal[j][0])))):
                wynik += (listClustal[j][1] - listMuscle[i][1]) * (listClustal[j][1] - listMuscle[i][1])
                found = True
    if found == False:
        wynik += listMuscle[i][1] * listMuscle[i][1]
            
for i in range(len(listClustal)):
    found = False
    for j in range(len(listMuscle)):
        if found == False:
            interSection = set(listMuscle[j][0]).intersection(listClustal[i][0])
            if len(interSection) == 0 or (len(interSection) == len(listMuscle[j][0]) and len(interSection) == len(listClustal[i][0])):
                found = True
    if found == False:
        wynik += listClustal[i][1] * listClustal[i][1]
        
print(wynik)
