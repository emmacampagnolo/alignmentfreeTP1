
# Alignment free - TP 1

L'objectif du TP est de comparer 5 especes de bactéries entre elles.
Vous trouverez les données en suivant [ce lien](https://we.tl/t-WeWvheBBGX)

## Composer le TP

Vous devez forker ce projet puis compléter ses fonctions.
Le rendu sera le dépot git dans lequel vous aurrez forké.

Le but est d'obtenir toutes les distances paire à paire des différentes bactéries.
Vous pouvez modifier l'affichage final pour obtenir une matrice d'adjacence si vous les souhaitez.

En observant les distances obtenues, que pouvez-vous dire des espèces présentes dans cet échantillon ?


### Remarques/Interprétation des résultats

# OBJECTIF DU TP:
L'objectif de ce TP était de comparer deux génomes rapidement.
On a donc pas voule faire un alignement puisque les alignments 
ont une compléxité très forte. L'idée était de regarder plutôt
la composition en kmers de deux génomes et de la comparer. Ainsi,
il faut déterminer l'ensemble des kmers du génome a et ceux du 
génome b et puis de déterminer l'intersection de ces deux ensembles. 
Un autre moyen utilisé pour résuire la compléxité et ainsi le 
temps d'éxécution de l'algorithe était d'encoder nos kmers en 
entiers pour manipuler des entiers plutôt que des chaînes de 
caractères. 


# REMARQUES:
On a réussi à implémenter la version optimisée de l'algorithme vu 
en cours. 

Nous avons utilisé deux méthodes pour comparer les génomes:

- METHOD 1: 
On construit la liste A des kmers du génome a, on construit ensuite 
le dictionnaire de comptage dico_A (avec la fonction Counter de python). 
Ensuite, on construit progressivement la liste B des kmers du génome
b, mais simultanément, pour chaque kmer lu, on véfie s'il est dans le 
dico_A (le dictionnaire de comptage des kmers du génome a). Si le kmer 
lu dans le génome b est une clé de dico_A, on sait que ce kmer est dans 
l'intersection des deux ensembles. On met à jour dico_A en enlevant le 
kmer lu et on continue de parcourir le génome b. Ainsi, on a construit
la liste inter contenant tous les kmers appartenant à la fois à A et à
B.

- METHOD 2:
On contruit les deux listes A et B avec les kmers des génomes a et b 
(respectiement). Ensuite on trie les deux listes de kmers. Ensuite, 
on peut comparer les éléments des listes triées 1 à 1 afin de déterminer
quels kmers sont identiques entre les deux listes. 

Les deux méthodes renvoient les mêmes résultats. La méthode 1 est la 
plus rapide des deux. En effet, cel n'est pas suprenant. Les deux 
méthodes requièrent la construction des listes des kmers. Par contre,
dans la méthode 1 il ne faut que chercher dans un dictionnaire (opération
en temps constant) et dans la méthode 2 il faut trier des listes qui 
est une étape qui est une opération plus couteuse en temps (de l'ordre
de nlogn). La méthode 1 parmi ces deux méthodes est donc la meilleure.




# RESULTATS: 
On peut voir les résultats dans le fichier image results.png. On a 
la matrice d'ajacence pour les 5 espèces étdiées. On observe que 
certaines espèces sont très simlaires (ex: GCF_000022165.1_ASM2216v1 et
 GCF_008244785.1_ASM824478v1) et d'autres sont très différentes (GCF_000005845.2_ASM584v2
et GCF_000006945.2_ASM694v2). On peut utiliser cette matrice d'adjacence
afin de cosntruire un abre. L'arbre obtenu par l'alogrithme NJ est
aussi inlus (species_tree_NJ.jpg) dans le fichier alinementfreeTP1. 
Dans cet arbre:
A= GCF_000005845.2_ASM584v2
B= GCF_000006945.2_ASM694v2
C= GCF_000008865.2_ASM886v2
D= GCF_000022165.1_ASM2216v1
E= GCF_008244785.1_ASM824478v1



# Alignment free - TP 2


# OBJECTIF DU TP:
L'objectif de ce deuxième TP était d'encore optimiser les fonctions du tp1. 
On voulait, plutôt que devoir garder en mémoire tous les kmers des deux génomes, 
d«en garder que des échantillons que l'on pourrait comparer qui seraient représen-
tatifs des deux populations. 
Pour cela, nous avons codé trois fonctions:
- make_sketch_1: ici, dans cette version, l'idée est de faire l'encodage que nous avons utilisé 
dans le tp1 mais ensuite de ne garder qu'une partie des kmers de chaque génome. Afin que ces échantillons 
puissent âtre comparés d'un génome à un autre, il faut que les kmars gardés en mémoire soient relativement
similaire. L'idée alors est de garder qu'un certain nombre (ici, 1000) de kmers ayant la valeur la plus faible pour chaque génome. 
Par contre, si on fait cela, on risque d'avoir des grands biais sur les kmers conservés puisque les kmers ayant
les valeurs les plus faibles seront des kmers avec beaucoup de nucléotides "A" au début. Il faut donc trouver un 
moyen de "mélanger" les encodages des kmers. On a fdonc utilisé une fonction xorshift qui permet de faire cela.
Concrètement, on initialise une liste vide L et on parcourt les kmers du génome. Tant la liste est de taille inférieure
à la taille de notre sketch (ici, 1000), on ajoute les kmers à notre liste. Ensuite, on regarde si le prochain kmer du génome est 
inférieur au max de la liste. Si oui, on enlève le max et on ajoute le kmer à la liste, sinon on passe au prochain kmer du génome. 
- make_sketch_2: Dans cette deuxième  fonction, on veut encore optimiser la fonction make_sketch_1. Plutôt que garder les 
kmers en mémoire sout forme de liste, on va utilise la structure de "heap", cela permet de plus rapidement obtenir le maximum et de réordonner
le tas. 
- make_sketch_3: Dans cette denrière fonction, l'idée c'est d'utiliser des "casiers" où l'on peut stocker nos kmers. On aura un certain nombre N de casiers
et pour chaque kmer, on le stockera dans le casier kmer%N. Si un autre kmer est déjà présent dans ce casier, on ne garde que celui avec le valeur plus faible.


# RESULTATS: 
Avec les deux premières méthodes, on obtient les mêmes résultats (ce qui est attendu). On a des résultats un peu différents avec la méthode 3.