import random
import math


def alphabet(name):                          #Quelle est la bonne séquence IUPAC pour AA ?
    if (name == "nucleic"):
        return "acgt"
    elif (name == "protein"):
        return "ACDEFGHIKLMNPQRSTVWY"
    elif (name == "iupac_nucleic"):
        return "ACGTURYSWKMBDHVN-"
    elif (name == "iupac_protein"):
        return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"      

    
def randseq(num, alpha):                          #Les methodes built-in sont-elles valides aussi ?
    seqtmp = []
    listalpha = list(alpha)
    for i in range(0, num):
        randnumber = random.randint(0,len(alpha)-1)
        seqtmp.append(listalpha[randnumber])
    seq = "".join(seqtmp)
    return seq


def hamming(seq1, seq2):        #Que fait-on des cas où les deux seq sont de tailles différentes ?
    hammingcount = 0
    listseq1 = list(seq1)
    listseq2 = list(seq2)
    for i,j in zip(listseq1, listseq2):
        if i != j:
            hammingcount+=1
    return hammingcount


def mutate(seq, num_subs):
    listseq = list(seq)
    for i in range(num_subs):
        randnumber = random.randint(0,len(seq)-1)    
        tmpnuc = randseq(1, seq)
        while tmpnuc == listseq[randnumber]:              #On ne veut pas remplacer par le même nuc
             tmpnuc = randseq(1, seq)
        listseq[randnumber] = tmpnuc
    return hamming(seq, listseq)


def experiments(le, su, nb):
    v = []
    for i in range(nb):
        seq = randseq(le,alphabet("nucleic"))
        v.append(mutate(seq,su))
    return v


def mean(data):
    somme = 0.0
    for i in data:
        somme += i
    moy = somme / len(data)
    return moy


def variance(data):
    moy = mean(data)
    sommesq = 0.0
    for i in data:
        sommesq += i**2
    var = (sommesq - len(data)*moy**2)/(len(data)-1)
    return var


def std(data):
    var = variance(data)
    standev = math.sqrt(var)
    return standev


def distanceJC69(means, L):       #Dans la formule, p=moyenne des differents / longueur sequence
    p = means / L                 #C'est la frequence d'apparition d'une subsitution
    d = -3/4*math.log(1.0-p*4/3)
    return d


#Dans generate le=seqlen  nb=nbofexp  xx=nbofsub et dans experiments le=seqlen  su=nbofsub  nb=nbofexp
#On définit le comme la longueur de chaque sequence à utiliser, xx est la liste contenant le nombre
#de substitutions à tester. Pour nb, on fixe un nombre de tentatives dans generate qui correspond au
#nombre de fois que l'on va tester pour chaque valeur de i dans experiments, où i correspond à chaque
#valeur de la liste xx. On obient donc, pour des sequences de longueur 1000, 10 valeurs de Hamming
#distance pour chaque nombre de substitution testé. Ensuite on moyenne ces Hamming distance pour
#pouvoir calculer les JC distances associées à chaque valeur de subsitution.


def generate(le, nb, xx):
    data = []
    hamdist = []
    jcdist = []
    for i in xx:
        v = experiments(le, i, nb)
        m = mean(v)
        hamdist.append(m)
        d = distanceJC69(m, le)
        jcdist.append(d)
    data.append(hamdist)
    data.append(jcdist)
    return data
    

    

#Main

L = 1000
R = 10
x = list(range(100,2200,200))
data = generate(L,R,x)
print(data[0])
print(data[1])
