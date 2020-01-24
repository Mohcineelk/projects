import sqlite3
import numpy as np
from math import *

#Fonction d'Hachage

int1=np.vectorize(int)
def b2h(S):
    n=len(S)
    T=S+'0'*(-n%8)
    m=ceil(n/8)
    B=np.zeros((m,8))
    H=np.zeros(2*m)
    for i in range(m):
        for j in range(8):
            B[i][j]=T[8*i+j]
    for i in range(m):
        h=0
        for j in range(8):
            h=h+B[i][j]*2**j
        H[2*i]=h//16
        H[2*i+1]=h%16
    bi=''
    for i in range(2*m):
        bi=bi+hex(int(H[i]))[2:]
    return bi
def hx(s):
    if '0'<=s<='9': return int(s)
    elif s=='a' or s=='A': return 10
    elif s=='b' or s=='B': return 11
    elif s=='c' or s=='C': return 12
    elif s=='d' or s=='D': return 13
    elif s=='e' or s=='E': return 14
    elif s=='f' or s=='F': return 15
def h2b(H):
    H=H+len(H)%2*'0'
    m=len(H)//2
    B=np.zeros((m,8))
    T=np.zeros(8*m)
    for i in range(m):
        h=16*hx(H[2*i])+hx(H[2*i+1])
        for j in range(8):
            B[i][j]=h%2
            h=h//2
    for i in range(m):
        for j in range(8):
            T[8*i+j]=B[i][j]
    ch=''
    for i in range(len(T)):
        ch=ch+str(int(T[i]))
    return ch
def xor(a,b):
    return int(((not a) and b) or ((not b) and a))
def xors(s,t):
    m=max(len(s),len(t))
    s='0'*(m-len(s))+s
    t='0'*(m-len(t))+t
    X=''
    for i in range(m):
        X=X+str(xor(int(s[i]),int(t[i])))
    return X
def converts_a(s):
    w=len(s)//25
    A=np.zeros((5,5,w))
    for x in range(5):
        for y in range(5):
            for z in range(w):
                A[x][y][z]=s[w*(5*y+x)+z]
    return int1(A)
def converta_s(a):
    w=len(a[0][0])
    s=''
    for y in range(5):
        p=''
        for x in range(5):
            l=''
            for z in range(w):
                l=l+str(a[x][y][z])
            p=p+l
        s=s+p
    return s
def theta(a):
    w=len(a[0][0])
    c=np.zeros((5,w))
    d=np.zeros((5,w))
    A=np.zeros((5,5,w))
    for x in range(5):
        for z in range(w):
            c[x][z]=xor(xor(xor(a[x][0][z],a[x][1][z]),xor(a[x][2][z],a[x][3][z])),a[x][4][z])
    for x in range(5):
        for z in range(w):
            d[x][z]=xor(c[(x-1)%5][z],c[(x+1)%5][(z-1)%w])
    for x in range(5):
        for z in range(w):
            for y in range(5):
                A[x][y][z]=xor(a[x][y][z],d[x][z])
    return int1(A)
def rho(a):
    w=len(a[0][0])
    A=np.zeros((5,5,w))
    for z in range(w):
        A[0][0][z]=a[0][0][z]
    (x,y)=(1,0)
    for t in range(24):
        for z in range(w):    
            A[x][y][z]=a[x][y][(z-((t+1)*(t+2))//2)%w]
        (x,y)=(y,(2*x+3*y)%5)
    return int1(A)
def pi(a):
    w=len(a[0][0])
    A=np.zeros((5,5,w))
    for x in range(5):
        for y in range(5):
            for z in range(w):
                A[y][(2*x+3*y)%5][z]=a[x][y][z]
    return int1(A)
def chsy(a):
    w=len(a[0][0])
    A=np.zeros((5,5,w))
    for x in range(5):
        for y in range(5):
            for z in range(w):
                A[x][y][z]=xor(a[x][y][z],xor(a[(x+1)%5][y][z],1)*a[(x+2)%5][y][z])
    return int1(A)
def rc(t):
    n=t%255
    if n==0: return 1
    R=[1,0,0,0,0,0,0,0]
    for i in range(1,n+1):
        R=[0]+R
        R[0]=xor(R[0],R[8])
        R[4]=xor(R[4],R[8])
        R[5]=xor(R[5],R[8])
        R[6]=xor(R[6],R[8])
        R=R[:8]
    return int(R[0])
def iota(a,i):
    w=len(a[0][0])
    l=int(log2(w))
    A=np.zeros((5,5,w))
    for x in range(5):
        for y in range(5):
            for z in range(w):
                A[x][y][z]=a[x][y][z]
    RC=[0]*w
    for j in range(l+1):
        RC[(2**j)-1]=rc(j+7*i)
    for z in range(w):
        A[0][0][z]=xor(A[0][0][z],RC[z])
    return A
def rnd(a,i):
    return iota(chsy(pi(rho(theta(a)))),i)
def keccakp(b,n,s):
    w=b//25
    l=int(log2(w))
    for i in range(12+2*l-n,12+2*l):
        s=converta_s(int1(rnd(converts_a(s),i)))
    return s
def sponge(b,f,pad,r,s,d):
    P=s+pad(r,len(s))
    n=len(P)//(r)
    c=b-r
    S='0'*(b)
    for i in range(n):
        S=f(xors(S,P[i*r:r*(i+1)]+'0'*c))
    Z=''
    while len(Z)<d//4:
        Z=Z+b2h(S)[:r//4]
        S=f(S)
    return Z[:d//4]
def pad101(x,m):
    j=(-m-2)%x
    return '1'+j*'0'+'1'
g=lambda x:keccakp(1600,24,x)
def keccak(c,s,d):
    b=1600
    f=g
    pad=pad101
    r=b-c
    S=sponge(b,f,pad,r,s,d)
    return S
def binaire(m):
    l=list(m.encode())
    txt=''
    for e in l:
        s=len(bin(e))-2
        txt=txt+(8-s)*'0'+bin(e)[2:]
    return txt
def binhexa(s):
    s='0'*(-len(s)%4)+s
    b=''
    n=len(s)//4
    d={'0000':'0','0001':'1','0010':'2','0011':'3','0100':'4','0101':'5','0110':'6','0111':'7','1000':'8','1001':'9','1010':'a','1011':'b','1100':'c','1101':'d','1110':'e','1111':'f'}
    for i in range(n):
        b=b+d[s[4*i:4*i+4]]
    return b
def sha3_224(M):
    S=h2b(binhexa(binaire(M)))+'01'
    return keccak(448,S,224)
    
# Cryptage & Decryptage RSA    
    
def motnbr(texte):
    L=[]
    ch=""
    for i in range(len(texte)):
        indice=ord(texte[i])
        if 0<=indice<=9:
            ch=ch+"0"+"0"+str(indice)
        elif 10<=indice<=99:
            ch=ch+"0"+str(indice)
        else:
            ch=ch+str(indice)
    return ch    

def troisatrois(x):
    L=[]
    ch=""
    if len(x)%3==1:
        x="0"+"0"+x
    elif len(x)%3==2:
        x="0"+x
    for i in range(0,len(x)-2,3):
        for j in range(i,i+3):
            ch=ch+x[j]
        L.append(ch)
        ch=" "
    return L

def encrypter_rsa(pubkey,n, msg):    
    asci=motnbr(msg)
    L=troisatrois(asci)
    C=[]
    ch=""
    for i in range(len(L)):
        x=(int(L[i])**pubkey)%n                                                            
        if 0<=x<=9:
            ch=ch+"0"+"0"+str(x)
            C.append(ch)
            ch=""
        elif 10<=x<=99:
            ch=ch+"0"+str(x)
            C.append(ch)
            ch=""
        else:
            ch=ch+str(x)
            C.append(ch)
            ch=""
    return C
 
def dechiffrer(L):
    texte=""
    for i in range(len(L)):
        a=int(L[i])
        texte=texte+chr(a)
        a=0
    return texte

def decrypter_rsa(privkey, n, C):
    D=[]
    ch=""
    for i in range(len(C)):
        y=(int(C[i])**privkey)%n
        if 0<=y<=9:
            ch=ch+"0"+"0"+str(y)
            D.append(ch)
            ch=""
        elif 10<=y<=99:
            ch=ch+"0"+str(y)
            D.append(ch)
            ch=""
        else:
            ch=ch+str(y)
            D.append(ch)
            ch=""
    return dechiffrer(D)
    
    
    
# Le vote électronique
#connexion à la base de données déjà créée
def executer_requete_update(query):
    conn = sqlite3.connect('C:/Users/konan/Downloads/votesdb.db') #indiquer le chemin de la base de données votesdb.db
    cursor = conn.cursor()
    cursor.execute(query)
    conn.commit()
    conn.close()

def executer_requete_select(query):
    conn = sqlite3.connect('C:/Users/konan/Downloads/votesdb.db')#indiquer le chemin de la base de données votesdb.db
    cursor = conn.cursor()
    cursor.execute(query)
    resultats = cursor.fetchall()
    conn.close()
    return resultats
        

class Generateur : 

    def __init__(self):
        pass
        
    def reinitialiser_db(self):
        executer_requete_update("delete from electeurs")
        executer_requete_update("delete from empreintes")
        executer_requete_update("delete from votes")
        executer_requete_update("delete from Resultats")
        executer_requete_update("delete from codeN2crypte")
        
    def generer_codes(self):
        codes =  [['AHDETR', 'BHTDU'],['SAZERT', 'BNDAR'],['KILOPM', 'FREUI']]
        for c in codes :
            executer_requete_update("insert into electeurs values('"+c[0]+"','"+c[1]+"')")
    
    def generer_empreintes(self):
        # Remplir ( empreintes ) par la liste des codes N1 et les empreintes de N2
        codes =  [
            ['AHDETR', '3b3a467ded0cb3ad1db68e6e4a7be6858634838f81e005a283345c6e'],
            ['SAZERT', 'edf0da352cdc897b2a3831c73df66cacb5495319bfca6af5481e8cea'],
            ['KILOPM', '35760bee4b2fd1087364e1e55e74e358055b305da03f63ce11de73ed']
        ]
        for c in codes :
            executer_requete_update("insert into empreintes values('"+c[0]+"','"+c[1]+"',0)")
    
    def generer_partis(self):
        # Remplir ( partis ) par la liste des partis et le nbr de votes nul
        partis =  ['PAM','PJD','PI','UNL','NULL']
        
        for c in partis :
            executer_requete_update("insert into Resultats values('"+c+"',0)")        
            
class Commissaire : 
    
    def __init__(self):
        pass
    
    def verifier_authentification(self, N1):
        L = executer_requete_select("select * from empreintes where N1='"+N1+"' and VOTE=0")
        return len(L)!=0
        
    def verifier_electeur(self, N1):
        # Cette fonction vérifie si l'électeur n'a pas déja voté, 
        # si oui elle va marquer la ligne de ce code N1 comme déja voté et renvoie 
        # True si non elle renvoie False
        L = executer_requete_select("select * from empreintes where N1='"+N1+"' and VOTE=0")
        if len(L)==0:
            return False
        else :
            executer_requete_update("update empreintes set VOTE=1 where N1='"+N1+"'")
            return True
            
    def verifier_N2(self, N2):
        emp = sha3_224(N2)
        Empr = executer_requete_select("select * from empreintes where EMP='"+emp+"'")
        if len(Empr) == 0:
            return False
        else :
            return True
            


class Administrateur : 
    
    def __init__(self):
        pass
    
    def verifier_authentification(self,commissaire,  N1):
        return commissaire.verifier_authentification(N1)
        

        
        
        
class Anonymiseur : 
    
    def __init__(self):
        pass
    
    def traiter_vote(self, commissaire, vote, N2crypt, N1):
        # Cette fonction se charge de sauvegarder un vote s'il est valide  
        # cette fonction fait appel à la méthode "verifier_electeur" du commisaire si le résultat est True
        # elle sauvegarde le vote et renvoie True si non elle renvoie False
        c=Commissaire()
        l=c.verifier_electeur(N1)
        if l==True:
            for i in range (len(vote)):
                executer_requete_update("insert into votes values('"+str(vote[i])+"','"+N1+"')")
            for i in range (len(N2crypt)) :
                executer_requete_update("insert into codeN2crypte values('"+str(N2crypt[i])+"','"+N1+"')")
                
            
            return True
        else : 
            return False
        


class Decompteur : 
    
    def __init__(self):
        self.pub_key = 11
        self.__priv_key = 35
    
    def depouillement(self, administrateur, commissaire):
        # Retirer les bulletins de leur enveloppe (déchiffrer avec clé privée)

        L=executer_requete_select("select N1 from empreintes where VOTE=1")
        if (len(L)!=0):
            for i in range (len(L)):
                L1=executer_requete_select("select bulletin from votes where electeurs='"+L[i][0]+"'")
                L1bis = executer_requete_select("select N2crypt from codeN2crypte where N1 ='"+L[i][0]+"'")
                
                L2 = []
                L2bis = []
                chaine1=''
                chaine2=''
                for c in range(len(L1)):
                    sh=''
                    sh=str(L1[c])+sh
                    chaine1 = sh[2:len(sh)-3]
                    L2.append(float(chaine1))
                for i in range(len(L1bis)):
                    ch=''
                    ch=str(L1bis[i])+ch
                    chaine2 = ch[2:len(ch)-3]
                    L2bis.append(float(chaine2))
                    
                Parti = decrypter_rsa(self.__priv_key, 119, L2)
                N2 = decrypter_rsa(self.__priv_key, 119, L2bis)
                bool = commissaire.verifier_N2(N2)
                if bool == True :
                    
                    L3=executer_requete_select("select Nbr_de_votes from Resultats where Partis='"+Parti+"'")
                
                    conn = sqlite3.connect('C:/Users/konan/Downloads/votesdb.db')
                    cursor = conn.cursor()
                    cursor.execute("update Resultats set Nbr_de_votes=? where Partis='"+Parti+"'",(L3[0][0]+1,))
                    conn.commit()
                    conn.close()
        
                
        gagnant=executer_requete_select("select Partis from Resultats where Nbr_de_votes=(select max(Nbr_de_votes) from Resultats)")
        
        if(len(gagnant)>1):
            print("Egalite entre :")
            for i in range(len(gagnant)):
                print(gagnant[i][0]," ")
        else: 
            print("le gagnant est : ", gagnant[0][0])
        
        
class Electeur : 
    
    def __init__(self):
        self.N1= ''

    def authentification(self, administrateur, commissaire):
        # Cette fonction demande à l'utilisateur de saisir son code N1 et envoie une demande d'authentification
        # La fonction fait appel à la méthode "verifier_authentification" de l'administrateur
        # Cette fonction return True si l'authentification s'est bien passé, False si non
        self.N1 = input("Donner votre code Num 1 : ")
        verif = administrateur.verifier_authentification(commissaire, self.N1)
        return verif
        
    def voter(self, administrateur, decompteur, anonymiseur,commissaire) : 
        # Cette fonction demande la saisit du vote ( Un nombre par exemple) puis demande la saisie du code N2
        # puis elle construit un bulletin de vote masqué et demande à l'administrateur de le signer et récupère
        # le bulletin signé et masqué puis enlève le masque
        # la fonction ensuite chiffre avec la clé publique du decompteur le bulletin signé, 
        # et envoie le résultat avec le code N1 à l'anonymiseur
        # SI le vote est bien enregistré on retourne True , False si non
        
        self.vote = input("veuillez saisir votre vote : ")
        self.N2 = input("Veuillez saisir le code N2 : ")
        L=executer_requete_select("select Partis from Resultats where Partis='"+self.vote+"'")
        while(len(L)==0):
            print("Vote invalide, veuillez réessayer!\n")
            self.vote = input("veuillez saisir votre vote : ")
            L=executer_requete_select("select Partis from Resultats where Partis='"+self.vote+"'")        
        chiff_vote = encrypter_rsa(11, 119, self.vote)
        chiff_N2 = encrypter_rsa(11, 119, self.N2)         
        r = anonymiseur.traiter_vote(commissaire,chiff_vote,chiff_N2, self.N1)           
        return r
        
#======================= PROGRAMME PRINCIPAL =================================

generateur = Generateur()
administrateur = Administrateur()
electeur = Electeur()
commissaire = Commissaire()
decompteur = Decompteur()
anonymiseur = Anonymiseur()

# Préparer le vote
generateur.reinitialiser_db()
generateur.generer_codes()
generateur.generer_empreintes()
generateur.generer_partis()

nbr_electeurs=executer_requete_select("select count(N1) from empreintes")
i=0

while (i < nbr_electeurs[0][0]):

    # Authentification
    authorise = electeur.authentification(administrateur, commissaire)
    while authorise== False : 
        print("Code invalide\n")
        authorise = electeur.authentification(administrateur, commissaire)
    
    # Etape de vote
    print("###################################################\n")
    print("#    1- PAM                                       #\n")
    print("#    2- PJD                                       #\n")
    print("#    3- PI                                        #\n")
    print("#    4- UNL                                       #\n")
    print("#    5- NULL (Vote nul)                           #\n")
    print("###################################################\n")
    
    reponse = electeur.voter(administrateur, decompteur, anonymiseur,commissaire)
    while reponse==False: 
        print("Problème de sauvegarde de vote, veuillez réessayer!\n")
        reponse = electeur.voter(administrateur, decompteur, anonymiseur,commissaire)
    print("Votre vote a bien été sauvegardé, Merci.")
    
    i+=1
    
# Etape de dépouillement
decompteur.depouillement(administrateur, commissaire)