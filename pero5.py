#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Universidade Federal de Pelotas - Campus Capao do Leao
#
# Programa para avaliar Perovskitas
#
# Autor: Jonatas Favotto Dalmedico                                    data: Setembro 2022
# Orientador: Prf. Dr.  Maurício Jeomar Piotrowski                    Coorientador: Prof. Dr. Diego Guedes-Sobrinho


# In[2]:


# bibliotecas e pacotes
import numpy as np
import csv
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


# In[3]:


#funcoes

#busca o maximo e minimo de um vetor
def max_min_array(vetor_aux):
    
    
        #calculando o valor maximo do vetor

        #estipulando o valor maximo inicial como o primeiro elemento do vetor
        val_max = vetor_aux[0]

        #loop de procura pelo maior valor
        for i in range(len(vetor_aux)):

            #comparando cada elemento do vetor como numero maximo
            if (vetor_aux[i] > val_max):

                #se sim --> novo valor maximo abaixo
                val_max = vetor_aux[i]

        vetor_max = val_max

        #calculando o valor minimo do vetor

        #estipulando o valor minimo inicial como o primeiro elemento do vetor
        val_min = vetor_aux[0]

        #loop de procura pelo menor valor
        for i in range(len(vetor_aux)):

            #comparando cada elemento do vetor como numero minimo
            if (vetor_aux[i] < val_min):

                #se sim --> novo valor maximo abaixo
                val_min = vetor_aux[i]

        vetor_min = val_min
                
        return vetor_max, vetor_min

def direc_to_cart_coord(a_1, a_2, a_3, coord_direc):
    
    #1 parametros de rede
    
    #1.1 calculando as constantes de rede
    
    a = np.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
    b = np.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)
    c = np.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)
    
    #1.2 Calculando o ângulo alfa entre os vetores a2 e a3 (ou b e c);
    # o ângulo beta entre os vetores a1 e a3 (ou a e c); e o ângulo gama
    # entre os vetores a1 e a2 (ou a e b)

    alfa = np.arccos((a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2])/(b*c))
    beta = np.arccos((a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2])/(a*c))
    gama = np.arccos((a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2])/(a*b))
    
    #2 preparando as variaveis fracionais da matriz input "coord_direc"

    u = coord_direc[:,0]
    v = coord_direc[:,1]
    w = coord_direc[:,2]

    #3 calculnando o volume da celula
    
    omega = a*b*c*np.sqrt(1-np.cos(alfa)**2-np.cos(beta)**2-np.cos(gama)**2+2*np.cos(alfa)*np.cos(beta)*np.cos(gama))
    
    #4 Transformando as coordenadas frac para cart

    x = a*u + b*np.cos(gama)*v + c*np.cos(beta)*w
    y = b*np.sin(gama)*v + c*((np.cos(alfa)-np.cos(beta)*np.cos(gama))/(np.sin(gama)))*w
    z = ((omega)/(a*b*np.sin(gama)))*w

    #5 montando ma matriz de coordenadas finais

    coord_cart = np.array([x,y,z]).transpose()
    
    return coord_cart
    


# In[4]:


#identificando o numero de camadas
n_camada = 5


# In[5]:


# Abrindo arquivo de dados da perovskita para extrair a matriz de vetores da celula unitaria (lattice vecotors)

with open('CONTCAR','r') as inf: #errado! o certo seria CONTCAR ou o final do OUTCAR!
    #incsv = csv.reader(inf)
    #Lattice Constant
    latt_constant = inf.readlines()[1]
    
with open('CONTCAR','r') as inf: #errado! o certo seria CONTCAR ou o final do OUTCAR!
    #incsv = csv.reader(inf)
    #Lattice vector
    latt_vec = inf.readlines()[2:5]
    latt_vec = [line.split() for line in latt_vec] 

#print(latt_constant)    
#print(latt_vec)


# In[6]:


#trabalhando as coordenadas dos vetores a_1, a_2 e a_3
coord_latt_vec = np.zeros((3,3))
for i in range(3):
    coord_latt_vec[i][0] = latt_vec[i][0] #x
    coord_latt_vec[i][1] = latt_vec[i][1] #y
    coord_latt_vec[i][2] = latt_vec[i][2] #z

#print('Lattice Constant:',latt_constant)
#print('Lattice vectors:\na_1=',coord_latt_vec[0][:],'\na_2=',coord_latt_vec[1][:],'\na_3=',coord_latt_vec[2][:])


# In[7]:


# Abrindo arquivo de dados da perovskita

with open('CONTCAR','r') as inf: #errado! o certo seria CONTCAR ou o final do OUTCAR!
    #incsv = csv.reader(inf)
    reader = inf.readlines()[5:7]
    reader = [line.split() for line in reader] 
    
#print(reader)


# In[8]:


#preparando o calculo: separando os elementos quimicos
nome_atom = reader[0][:]
#print(nome_atom)

#preparando o calculo: separando o vetor de quantidade em funcao do elemento quimico
no_atom = reader[1][:]
#print(no_atom)

no_atom = [ int(x) for x in no_atom ]
##print(type(no_atom))

#calculando o no. total de atomos
no_atom_total=0
for i in range(5):
    no_atom_total = no_atom_total + no_atom[i]

#print(no_atom_total)


# In[9]:


#separando o cabecalho das coordenadas dos atomos
with open('CONTCAR','r') as inf:
    #incsv = csv.reader(inf)
    lista = inf.readlines()[9:]
    lista = [line.split()[:3] for line in lista]      
dados = lista
#dados = dados[:,-3]
##print(dados)


# In[10]:


#trabalhando as coordenadas ate o numero de atomos
coord_atom = np.zeros((3,no_atom_total))
coord_atom = np.zeros((no_atom_total,3))
for i in range(no_atom_total):
    coord_atom[i][0] = dados[i][0] #x
    coord_atom[i][1] = dados[i][1] #y
    coord_atom[i][2] = dados[i][2] #z
    
#print('Matriz de coordenadas dos atomos X e M em coordenadas "Direct".\n',coord_atom)

#estipulando o valor maximo nas direcoes x, y e z 

#em x
val_max = coord_atom[0][0]

#loop de procura pelo maior valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero maximo
    if (coord_atom[i][0] > val_max):
        
        #se sim --> novo valor maximo abaixo
        val_max = coord_atom[i][0]
    
dist_x_max = val_max

#print('\n x maximo: ',dist_x_max)

#estipulando o valor minimo 
val_min = coord_atom[0][0]

#loop de procura pelo menor valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero minimo
    if (coord_atom[i][0] < val_min):
        
        #se sim --> novo valor maximo abaixo
        val_min = coord_atom[i][0]
    
dist_x_min = val_min

#print('\n x minimo: ',dist_x_min)


#em y
val_max = coord_atom[0][1]

#loop de procura pelo maior valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero maximo
    if (coord_atom[i][1] > val_max):
        
        #se sim --> novo valor maximo abaixo
        val_max = coord_atom[i][1]
    
dist_y_max = val_max

#print('\n y maximo: ',dist_y_max)

#estipulando o valor minimo 
val_min = coord_atom[0][1]

#loop de procura pelo menor valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero minimo
    if (coord_atom[i][1] < val_min):
        
        #se sim --> novo valor maximo abaixo
        val_min = coord_atom[i][1]
    
dist_y_min = val_min

#print('\n y minimo: ',dist_y_min)

#em z
val_max = coord_atom[0][2]

#loop de procura pelo maior valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero maximo
    if (coord_atom[i][2] > val_max):
        
        #se sim --> novo valor maximo abaixo
        val_max = coord_atom[i][2]
    
dist_z_max = val_max

#print('\n z maximo: ',dist_z_max)

#estipulando o valor minimo 
val_min = coord_atom[0][2]

#loop de procura pelo menor valor
for i in range(no_atom_total):
    
    #comparando cada elemento do vetor como numero minimo
    if (coord_atom[i][2] < val_min):
        
        #se sim --> novo valor maximo abaixo
        val_min = coord_atom[i][2]
    
dist_z_min = val_min

#print('\n z minimo: ',dist_z_min)

#calculando o tamanho da celula

#em x
lx = 1

#em y
ly = 1

#em z
lz = 1

##print('\nLattice vectors')
##print('a=',lx,'\nb=',ly,'\nc=',lz)


# In[11]:


#separando coordenadas dos metais (considerando que esses sejam os dois atomos da lista 'reader')


#Metal cental M:
#print('\n Metal M:', reader[0][-2])
#de quatidade
no_M = int(reader[1][-2])
#print('\n Quatidade:', reader[1][-2])

#Metal de fronteira X:
#print('\n Metal X:', reader[0][-1])
#de quatidade
no_X = int(reader[1][-1])
#print('\n Quatidade:', reader[1][-1])

#no. total de metais M+X
no_metais = no_M + no_X
#print('\n Quatidade total M + X = ',no_metais)

##print(coord_atom[0][-10])

#separando coordenadas
#Metal cental M:
#coord_atom_M = np.zeros((3,no_M))
coord_atom_M = np.zeros((no_M,3))
for i in range(no_M):
    coord_atom_M[i][0] = coord_atom[-no_metais+i][0] #coord_atom[0][-no_metais+i] #x
    coord_atom_M[i][1] = coord_atom[-no_metais+i][1] #coord_atom[1][-no_metais+i] #y
    coord_atom_M[i][2] = coord_atom[-no_metais+i][2] #coord_atom[2][-no_metais+i] #z
    
#print('\n Coordenadas xyz de M:\n',coord_atom_M)

#Metal cental X:
coord_atom_X = np.zeros((no_X,3))
for i in range(no_X):
    coord_atom_X[i][0] = coord_atom[-no_metais+i+no_M][0] #x
    coord_atom_X[i][1] = coord_atom[-no_metais+i+no_M][1] #y
    coord_atom_X[i][2] = coord_atom[-no_metais+i+no_M][2] #z
    
#print('\n Coordenadas xyz de X:\n',coord_atom_X)


# In[12]:


#preparando dados de posicao de M

xm = np.transpose(np.delete(coord_atom_M,[1,2],1))
ym = np.transpose(np.delete(coord_atom_M,[0,2],1))
zm = np.transpose(np.delete(coord_atom_M,[0,1],1))

xmm = np.squeeze(xm) #xx
ymm = np.squeeze(ym) #yy
zmm = np.squeeze(zm) #zz

#preparando dados de posicao de X
xx = np.transpose(np.delete(coord_atom_X,[1,2],1))
yx = np.transpose(np.delete(coord_atom_X,[0,2],1))
zx = np.transpose(np.delete(coord_atom_X,[0,1],1))

xxx = np.squeeze(xx) #xx
yxx = np.squeeze(yx) #yy
zxx = np.squeeze(zx) #zz

#plt.title('Coordenadas')
#plt.plot(xxx,zxx,'o',color='blue',label='X')
#plt.plot(xmm,zmm,'o',color='black',label='M')
#plt.xlabel('x')
#plt.ylabel('z')
#plt.legend()
#plt.show()

#plt.title('Coordenadas')
#plt.plot(yxx,zxx,'o',color='blue',label='X')
#plt.plot(ymm,zmm,'o',color='black',label='M')
#plt.xlabel('y')
#plt.ylabel('z')
#plt.legend()
#plt.show()


#plt.title('Coordenadas')
#plt.plot(xxx,yxx,'o',color='blue',label='X')
#plt.plot(xmm,ymm,'o',color='black',label='M')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.legend()
#plt.show()

#conclusao

#direcao c: z
#plano ab: xy


# In[13]:


##caso 2: 2D (n=4)


#replicacao de supercelula 2x2x2

#Consiste em adicionar celulas unitarias nas fronteiras do plano ab (xy). 
#pensando em uma matriz, [0][0] seria a celula fundamental enquanto [0][1], [1][0] e [1][1] 
#as celulas de fronteira (pontos colaterais da rosa dos ventos: noroeste, nordeste, sudeste, 
#sudoeste). Para fazer isso basta adicionar essas celulas ao acrescentar +1 nas coordenadas 
#as quais correspondem aa direcao da fronteira.

#print('\n Em: \n [0][0] [0][1] \n [1][0] [1][1] \n ou: \n [Noroeste: NW] [Nordeste: NE] \n [Sudoeste: SW] [Sudeste: SE] \n')

#dados
aux_M = coord_atom_M
aux_X = coord_atom_X
#aux_M = coord_atom_M_cart
#aux_X = coord_atom_X_cart
coord_atom_M_NW = np.zeros((no_M,3))
coord_atom_X_NW = np.zeros((no_X,3))
coord_atom_M_NE = np.zeros((no_M,3))
coord_atom_X_NE = np.zeros((no_X,3))
coord_atom_M_SW = np.zeros((no_M,3))
coord_atom_X_SW = np.zeros((no_X,3))
coord_atom_M_SE = np.zeros((no_M,3))
coord_atom_X_SE = np.zeros((no_X,3))

#noroeste NO(ou NW) [0][0] (celula primordial)

coord_atom_M_NW = aux_M
coord_atom_X_NW = aux_X


#nordeste NE [0][1]

for i in range(no_M): #+a em x
    coord_atom_M_NE[i][0] = coord_atom_M_NE[i][0] + lx #a_uni
    
coord_atom_M_NE = coord_atom_M_NE + aux_M

for i in range(no_X): #+a em x
    coord_atom_X_NE[i][0] = coord_atom_X_NE[i][0] + lx #a_uni

coord_atom_X_NE = coord_atom_X_NE + aux_X


#sudoeste SO(ou SW) [1][0]

for i in range(no_M): #+ab em y
    coord_atom_M_SW[i][1] = coord_atom_M_SW[i][1] + ly #b_uni

coord_atom_M_SW = coord_atom_M_SW + aux_M 

for i in range(no_X): #+b em y
    coord_atom_X_SW[i][1] = coord_atom_X_SW[i][1] + ly #b_uni

coord_atom_X_SW = coord_atom_X_SW + aux_X

#sudeste SE [1][1]

for i in range(no_M): #+a e +b em x e y
    coord_atom_M_SE[i][0] = coord_atom_M_SE[i][0] + lx #a_uni
    coord_atom_M_SE[i][1] = coord_atom_M_SE[i][1] + ly #b_uni
    
coord_atom_M_SE = coord_atom_M_SE + aux_M

for i in range(no_X): #+b em x e y
    coord_atom_X_SE[i][0] = coord_atom_X_SE[i][0] + lx #a_uni
    coord_atom_X_SE[i][1] = coord_atom_X_SE[i][1] + ly #b_uni

coord_atom_X_SE = coord_atom_X_SE + aux_X

#variaveis finais
#print('\n--> celula primordial [0][0] \n')
#print('M noroeste final: \n',coord_atom_M_NW)
#print('X noroeste final: \n',coord_atom_X_NW)
#print('\n--> [0][1] \n')
#print('M nordeste final: \n',coord_atom_M_NE)
#print('X nordeste final: \n',coord_atom_X_NE)
#print('\n--> [1][0] \n')
#print('M sudoeste final: \n',coord_atom_M_SW)
#print('X sudoeste final: \n',coord_atom_X_SW)
#print('\n--> [1][1] \n')
#print('M sudeste final: \n',coord_atom_M_SE)
#print('X sudeste final: \n',coord_atom_X_SE)

#variaveis imutaveis
##print('\nVariaveis iniciais \n')
##print('aux de M \n',aux_M)
##print('aux de X \n',aux_X)
##print('coord_atom_M \n',coord_atom_M)
##print('coord_atom_X \n',coord_atom_X)


# In[14]:


#Passando a celula para cartesiana

#preparando os vetores 
a1 = coord_latt_vec[0][:]
a2 = coord_latt_vec[1][:]
a3 = coord_latt_vec[2][:]

B = coord_atom_M_NW
coord_atom_M_NW = direc_to_cart_coord(a1, a2, a3, B)
#print('\n--> celula primordial [0][0] \n')
#print('M noroeste final: \n',coord_atom_M_NW)

C = coord_atom_X_NW
coord_atom_X_NW = direc_to_cart_coord(a1, a2, a3, C)
#print('X noroeste final: \n',coord_atom_X_NW)

B = coord_atom_M_NE
coord_atom_M_NE = direc_to_cart_coord(a1, a2, a3, B)
#print('\n--> [0][1] \n')
#print('M nordeste final: \n',coord_atom_M_NE)

C = coord_atom_X_NE
coord_atom_X_NE = direc_to_cart_coord(a1, a2, a3, C)
#print('X nordeste final: \n',coord_atom_X_NE)

B = coord_atom_M_SW
coord_atom_M_SW = direc_to_cart_coord(a1, a2, a3, B)
#print('\n--> [1][0] \n')
#print('M sudoeste final: \n',coord_atom_M_SW)

C = coord_atom_X_SW
coord_atom_X_SW = direc_to_cart_coord(a1, a2, a3, C)
#print('X sudoeste final: \n',coord_atom_X_SW)

B = coord_atom_M_SE
coord_atom_M_SE = direc_to_cart_coord(a1, a2, a3, B)
#print('\n--> [1][1] \n')
#print('M sudeste final: \n',coord_atom_M_SE)

C = coord_atom_X_SE
coord_atom_X_SE = direc_to_cart_coord(a1, a2, a3, C)
#print('X sudeste final: \n',coord_atom_X_SE)


# In[15]:


#criando super celula para 2x2x1

#print('Parâmetros de rede da célula unitária:')
#parametros de rede da celula unitaria
a_uni = np.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)

b_uni = np.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)

c_uni = np.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)

'(a b c) = ({:^.2f} {:^.2f} {:^.2f})'.format(a_uni,b_uni,c_uni)


# In[16]:


#preparando dados para plotar o sistema

#preparando dados para os Ms
x_M_NW = np.transpose(np.delete(coord_atom_M_NW,[1,2],1))
x_M_NE = np.transpose(np.delete(coord_atom_M_NE,[1,2],1))
x_M_SW = np.transpose(np.delete(coord_atom_M_SW,[1,2],1))
x_M_SE = np.transpose(np.delete(coord_atom_M_SE,[1,2],1))
x_M = np.append(x_M_NW,x_M_NE)
x_M = np.append(x_M,x_M_SW)
x_M = np.append(x_M,x_M_SE)

y_M_NW = np.transpose(np.delete(coord_atom_M_NW,[0,2],1))
y_M_NE = np.transpose(np.delete(coord_atom_M_NE,[0,2],1))
y_M_SW = np.transpose(np.delete(coord_atom_M_SW,[0,2],1))
y_M_SE = np.transpose(np.delete(coord_atom_M_SE,[0,2],1))
y_M = np.append(y_M_NW,y_M_NE)
y_M = np.append(y_M,y_M_SW)
y_M = np.append(y_M,y_M_SE)

z_M_NW = np.transpose(np.delete(coord_atom_M_NW,[0,1],1))
z_M_NE = np.transpose(np.delete(coord_atom_M_NE,[0,1],1))
z_M_SW = np.transpose(np.delete(coord_atom_M_SW,[0,1],1))
z_M_SE = np.transpose(np.delete(coord_atom_M_SE,[0,1],1))
z_M = np.append(z_M_NW,z_M_NE)
z_M = np.append(z_M,z_M_SW)
z_M = np.append(z_M,z_M_SE)

xmm = np.squeeze(x_M) #xx #muda de formato [[...]] para []
ymm = np.squeeze(y_M) #yy
zmm = np.squeeze(z_M) #zz

#preparando dados dos Xs
x_X_NW = np.transpose(np.delete(coord_atom_X_NW,[1,2],1))
x_X_NE = np.transpose(np.delete(coord_atom_X_NE,[1,2],1))
x_X_SW = np.transpose(np.delete(coord_atom_X_SW,[1,2],1))
x_X_SE = np.transpose(np.delete(coord_atom_X_SE,[1,2],1))
x_X = np.append(x_X_NW,x_X_NE)
x_X = np.append(x_X,x_X_SW)
x_X = np.append(x_X,x_X_SE)

y_X_NW = np.transpose(np.delete(coord_atom_X_NW,[0,2],1))
y_X_NE = np.transpose(np.delete(coord_atom_X_NE,[0,2],1))
y_X_SW = np.transpose(np.delete(coord_atom_X_SW,[0,2],1))
y_X_SE = np.transpose(np.delete(coord_atom_X_SE,[0,2],1))
y_X = np.append(y_X_NW,y_X_NE)
y_X = np.append(y_X,y_X_SW)
y_X = np.append(y_X,y_X_SE)

z_X_NW = np.transpose(np.delete(coord_atom_X_NW,[0,1],1))
z_X_NE = np.transpose(np.delete(coord_atom_X_NE,[0,1],1))
z_X_SW = np.transpose(np.delete(coord_atom_X_SW,[0,1],1))
z_X_SE = np.transpose(np.delete(coord_atom_X_SE,[0,1],1))
z_X = np.append(z_X_NW,z_X_NE)
z_X = np.append(z_X,z_X_SW)
z_X = np.append(z_X,z_X_SE)
##print(z_X)

xxx = np.squeeze(x_X) #xx #muda de formato [[...]] para []
yxx = np.squeeze(y_X) #yy
zxx = np.squeeze(z_X) #zz

#plotando a super celula
#plt.title('Coordenadas')
#plt.plot(xxx,zxx,'o',color='blue',label='X')
#plt.plot(xmm,zmm,'o',color='black',label='M')
#plt.xlabel('x')
#plt.ylabel('z')
#plt.legend()
#plt.show()

#plt.title('Coordenadas')
#plt.plot(yxx,zxx,'o',color='blue',label='X')
#plt.plot(ymm,zmm,'o',color='black',label='M')
#plt.xlabel('y')
#plt.ylabel('z')
#plt.legend()
#plt.show()

#plt.title('Coordenadas')
#plt.plot(xxx,yxx,'o',color='blue',label='X')
#plt.plot(xmm,ymm,'o',color='black',label='M')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.legend()
#plt.show()


# In[17]:


#nova quantidade de M e X

no_M = 4*no_M
no_X = 4*no_X
##print(no_M)
#print(z_X)


# In[18]:


#separando coordenadas por camada

#para n=5
#separando Nde cima e de baixo

#dados e vetores para os calculos
#z_medio =
z_fifth = c_uni/n_camada #correspondente a metade da altura. A celula é dividida em tres camadas. n/3 Ms na camada 1 e n/2 na camada 2
no_n_M = int(no_M/n_camada)
##print(no_n_M)
M_camada_index = np.zeros((n_camada,no_n_M))

# supercelula 2x2x1
#Em M
x_M = np.squeeze(xmm)
y_M = np.squeeze(ymm)
z_M = np.squeeze(zmm)
##print(z_M )

#Em X
x_X = np.squeeze(xxx)
y_X = np.squeeze(yxx)
z_X = np.squeeze(zxx)
##print(xx)

#contadores
n1 = 0 
n2 = 0 
n3 = 0 
n4 = 0
n5 = 0

#iniciando a busca
for i in range(no_M):
    
    ##print('M_',i)
    #loop
    if (z_M[i] > 4*z_fifth) and (z_M[i] < 5*z_fifth):
        
        #M esta na camada n_camada 1
        M_camada_index[0][n1] = i
        n1 += 1
                    
    elif (z_M[i] > 3*z_fifth) and (z_M[i] < 4*z_fifth):
        
        #M esta na camada n_camada 2
        M_camada_index[1][n2] = i
        n2 += 1
        
    elif (z_M[i] > 2*z_fifth) and (z_M[i] < 3*z_fifth):
        
        #M esta na camada n_camada 3
        M_camada_index[2][n3] = i
        n3 += 1
        
    elif (z_M[i] > 1*z_fifth) and (z_M[i] < 2*z_fifth):
        
        #M esta na camada n_camada 4
        M_camada_index[3][n4] = i
        n4 += 1
        
    else: #(z_M[i] < 1*z_third):
        
        #M esta abaixo e na camada 5
        M_camada_index[4][n5] = i
        n5 += 1
        
#print('* Ms na primeira camada: ',M_camada_index[0][:])

#print('* Ms na segunda camada: ',M_camada_index[1][:])

#print('* Ms na terceira camada: ',M_camada_index[2][:])

#print('* Ms na quarta camada: ',M_camada_index[3][:])

#print('* Ms na quinta camada: ',M_camada_index[4][:])


# In[19]:


#pre calculo 

#dados
z_medio = 0
X_max_index = 0
X_min_index = 0
no_M_por_camada = len(M_camada_index[0][:])
##print(no_M_por_camada)
z_medio_por_camada = np.zeros((len(M_camada_index),2))
#nn = len(M_camada_index)
##print(nn)
##print(z_medio_por_camada)

# supercelula 2x2x1
#Em M
x_M = np.squeeze(xmm)
y_M = np.squeeze(ymm)
z_M = np.squeeze(zmm)
##print(xm)

#Em X
x_X = np.squeeze(xxx)
y_X = np.squeeze(yxx)
z_X = np.squeeze(zxx)
##print(xx)

#calculando a altura z media entre todos os M da mesma camada

#loop pelos M por camada 
for i in range(n_camada): 
    
    ##print('camada ',i)
    #zerando o valor medio
    z_medio = 0
    
    #loop em cada M do sistema
    for j in range(no_M_por_camada):
           
        #buscando as coordenadas de cada coluna (indice M) de cada linha (camada)
        ##print(' * M_',int(M_camada_index[i][j]))
        #calculando o z_medio nessa camada
        z_medio += (z_M[int(M_camada_index[i][j])])/(no_M_por_camada)
        
        ##print(z_medio,z_M[int(M_camada_index[i][j])])
    
    #armazenando (no_camada ; z_medio)
    z_medio_por_camada[i][0] = i
    z_medio_por_camada[i][1] = z_medio
    
##print(z_medio_por_camada)
    
#print('\n Coordenada medias dos M em z por camada (no_camada ; z_medio): \n',z_medio_por_camada)
##print('\n * Obs.: em coordenadas Direct!')


# In[20]:


#procurando o numero de X max e minimos em cada camada

## | camada |                vale da camada               |               topo da camada               |
##    n1:    (z_X > z_M_medio_n2) & (z_X < z_M_medio_n1);               (z_X > z_M_medio_n1)

##    n2:    (z_X > z_M_medio_n3) & (z_X < z_M_medio_n2) ;  (z_X > z_M_medio_n2) & (z_X < z_M_medio_n1)(=n1) 

##    n3:    (z_X > z_M_medio_n4) & (z_X < z_M_medio_n3) ;  (z_X > z_M_medio_n3) & (z_X < z_M_medio_n2)(=n2) 

##    n4:    (z_X < z_M_medio_n4) ;                         (z_X > z_M_medio_n4) & (z_X < z_M_medio_n3)(=n3) 

#dados e vetores importantes
colunas_M = len(M_camada_index[0][:])
##print(colunas_M)
linhas_M = n_camada+1
X_camadas_index = np.zeros((linhas_M,colunas_M))
b1 = 0
b2 = 0
b3 = 0
b4 = 0
b5 = 0
b6 = 0
##print(X_camadas_index)


    
#intervalo desconsiderado
z_max_n1 = z_medio_por_camada[0][1] + 1 #acima desse valor se encontram os X topo
z_min_n1 = z_medio_por_camada[0][1] - 1 #abaixo desse valor se encontram os X vale
z_max_n2 = z_medio_por_camada[1][1] + 1 #acima desse valor se encontram os X topo
z_min_n2 = z_medio_por_camada[1][1] - 1 #abaixo desse valor se encontram os X vale
z_max_n3 = z_medio_por_camada[2][1] + 1 #acima desse valor se encontram os X topo
z_min_n3 = z_medio_por_camada[2][1] - 1 #abaixo desse valor se encontram os X vale
z_max_n4 = z_medio_por_camada[3][1] + 1 #acima desse valor se encontram os X topo
z_min_n4 = z_medio_por_camada[3][1] - 1 #abaixo desse valor se encontram os X vale
z_max_n5 = z_medio_por_camada[4][1] + 1 #acima desse valor se encontram os X topo
z_min_n5 = z_medio_por_camada[4][1] - 1 #abaixo desse valor se encontram os X vale

##print("n1: [",z_max_n1,z_min_n1,"]")
##print("n2: [",z_max_n2,z_min_n2,"]")
##print("n3: [",z_max_n3,z_min_n3,"]")
##print("n4: [",z_max_n4,z_min_n4,"]")
##print("n5: [",z_max_n5,z_min_n5,"]")

##print('* X_',i,'z_X=',z_X[i],'> z_max',z_max_n1)
    
#Loop de busca:
for i in range(no_X):
    
    ##print('* X_',i)    
    #condicional para X topo da n1
    if (z_X[i] > z_max_n1):
        
        ##print('OK!')
        ##print('* X_',i,'z_X=',z_X[i],'> z_max',z_max_n1)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[0][b1] = i
            
        b1 += 1
            
    
    #condicional para X no vale de n1 e topo da n2 (entre-camada n1 e n2)
    if (z_X[i] > z_max_n2) & (z_X[i] < z_min_n1):
        
        ##print('* X_',i,'z_X=',z_X[i],'> z_max=',z_max_n2,'e < z_min=',z_min_n1)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[1][b2] = i
        
        b2 += 1
    
    #condicional para X vale da n2 e topo da n3 (entre-camada n2 e n3)
    if (z_X[i] > z_max_n3) & (z_X[i] < z_min_n2):
        
        ##print('* X_',i,'z_X=',z_X[i],'< z_min=',z_min_n2,'e > z_min=',z_min_n3)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[2][b3] = i
        
        b3 += 1
    
    #condicional para X vale da n3 e topo da n4 (entre-camada n3 e n4) 
    if (z_X[i] > z_max_n4) & (z_X[i] < z_min_n3):
        
        ##print('* X_',i,'z_X=',z_X[i],'< z_min=',z_min_n2)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[3][b4] = i
            
        b4 += 1
        
    #condicional para X vale da n4 e topo da n5 (entre-camada n4 e n5) 
    if (z_X[i] > z_max_n5) & (z_X[i] < z_min_n4):
        
        ##print('* X_',i,'z_X=',z_X[i],'< z_min=',z_min_n2)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[4][b5] = i
            
        b5 += 1
        
    #condicional para X vale da n5
    if (z_X[i] < z_min_n5):
        
        ##print('* X_',i,'z_X=',z_X[i],'< z_min=',z_min_n2)
        
        #guardando o indice desse X. O numero de topos e igual ao numero de Ms
        X_camadas_index[5][b6] = i
            
        b6 += 1

#print('Vetor dos X topo:               ',X_camadas_index[0][:])
#print('Vetor dos X entre camada n1/n2: ',X_camadas_index[1][:])
#print('Vetor dos X entre camada n2/n3: ',X_camadas_index[2][:])
#print('Vetor dos X entre camada n3/n4: ',X_camadas_index[3][:])
#print('Vetor dos X entre camada n4/n5: ',X_camadas_index[4][:])
#print('Vetor dos X vale:               ',X_camadas_index[5][:])


#print('\nPara entender: \n Posicao dos X de topo          --> X   X \n Posicao dos M de n=1           --> M   M ')
#print(' Posicao de X entre camada n1/n2 --> X   X \n Posicao dos M de n=2           --> M   M ')
#print(' Posicao de X entre camada n2/n3 --> X   X \n Posicao dos M de n=3           --> M   M ')
#print(' Posicao de X entre camada n3/n4 --> X   X \n Posicao dos M de n=4           --> M   M ')
#print(' Posicao de X entre camada n4/n5 --> X   X \n Posicao dos M de n=5           --> M   M ')
#print(' Posicao dos X de vale          --> X   X')


# In[21]:


#Transformando a matriz de posicao de X topos, vales e entre camada para em um unico vetor
#e organizando em ordem crescente
#print('Matrix de X topos, meio de camada e vales \n')
#separando os X proximos de M. Para isso, identifica-se os X pelos indices das matrizes X_max_index e X_min_index 
#e subtrai-se suas correspondentes nas matrizes de coordenadas x_X, y_X e z_X.

#dados
##print("\n Matrix de X no topo do sistema: ", X_max_index)
##print("\n Matrix de X no vale do sistema: ", X_min_index)
#print('Vetor dos X topo:               ',X_camadas_index[0][:])
#print('Vetor dos X entre camada n1/n2: ',X_camadas_index[1][:])
#print('Vetor dos X entre camada n2/n3: ',X_camadas_index[2][:])
#print('Vetor dos X entre camada n3/n4: ',X_camadas_index[3][:])
#print('Vetor dos X entre camada n4/n5: ',X_camadas_index[4][:])
#print('Vetor dos X vale:               ',X_camadas_index[5][:])

##print("\n -- Matrix de coordenadas X do sistema -- \n - [x]: \n", x_X,"\n - [y]: \n",y_X,"\n - [z]: \n",z_X)
i_max=0
i_entrec = 0
i_min=0
index_X_camadas = len(X_camadas_index[0][:])
#print(len(X_camadas_index))

#matriz dos X fora do plano equatorial ab
range_toral = (len(X_camadas_index))*index_X_camadas #len(X_camadas_index[0][:]) + len(X_camadas_index[1][:]) + len(X_camadas_index[2][:])
X_api_index = np.zeros(range_toral)

##print('range toral=',range_toral)

#transformando a matriz de 3 vetores em um unico vetor
vetor_X_camadas_index = X_camadas_index.flatten()
g = 0
##print(vetor_X_camadas_index)    

#organizando o vetor em ordem rescente
for i in range(range_toral):
    
    ##print('Indice: ',i)
    
    #passando por todos os outros X
    for j in range(i+1,range_toral):
        
        #condicional: qual serio o menor
        if (vetor_X_camadas_index[i] > vetor_X_camadas_index[j]):
                
            #novo indice
            index_min = vetor_X_camadas_index[i]
            vetor_X_camadas_index[i] = vetor_X_camadas_index[j]
            vetor_X_camadas_index[j] = index_min
                
X_api_index = vetor_X_camadas_index            

#print('\n Vetor de X na direcao apical (fora do plano equatorial): ',X_api_index)    


# In[22]:


#separando os X no plano equatorial 
##print('X no plano equatorial')
#isso é possivel ao subtrair os atomos da matriz de X calculada anteriormente
#contendo X topos, vales e entre camada

#variaveis para o loop:
#coord_X_plan_equa = np.zeros((len(x_X)-len(X_max_index)-len(X_min_index),3))
coord_X_plan_equa = 0
# X_equa_index: matriz calculada anteriormente
X_equa_x = x_X
X_equa_y = y_X
X_equa_z = z_X
tamanho = len(X_api_index)

##print("\n X_z inicial:",X_equa_z)

##print('no X: ',no_X,' ; no colunas em x_X: ',len(x_X))
#variaveis auxiliares
ii = 0

##print('\n Vetor de coord x de X equa inicial:\n ',X_equa_x)

#loop decrescente pelos vetores posicao correspondentes a todos os X(x,y,z)
for i in range(tamanho-1,-1,-1):
    
    ##print('\n --> valor de i nessa iteracao: ',i,'\n   ---> valor do indice nessa iteracao: ',int(X_api_index[i]))
    #retirando os X do topo e do vale
    ii = int(X_api_index[i])
    X_equa_x = np.delete(X_equa_x,ii,None)
    ##print('\n Vetor de coord x de X equa: \n',X_equa_x)
    X_equa_y = np.delete(X_equa_y,ii,None)
    ##print('\n Vetor de coord y de X equa: \n',X_equa_y)
    X_equa_z = np.delete(X_equa_z,ii,None)
    ##print('\n Vetor de coord z de X equa: \n',X_equa_z)
    
##print("\n X_z final:",X_equa_z)

##print('\n Coordenadas dos X pertencentes ao plano equatorial: \n --> [x]:\n',X_equa_x,'\n --> [y]:\n',X_equa_y,'\n --> [z]:\n',X_equa_z)       


# In[23]:


##Para uma camada (n=4)

#separando os X do plano equatorial por camada e M vizinho
#print('X do plano equatorial')
#Busca-se identificar os X ao redor de cada M a partir de uma distacia minima 
#e maxima como intervalo nas direcoes x e y.

#Todos os atomos X e M respeitam uma proporcao pois trata-se de uma celula

#variaveis e vetores utilizados
# X_equa_x (ou X_equa_y)
dist_max_x = 5
#print('d_max_x',dist_max_x)
dist_max_y = 5 
#print('d_max_y',dist_max_y)
dist_max_z = 4
#print('d_max_z',dist_max_z)

busca_M = len(x_M)
busca_X = len(X_equa_x)
##print(busca_M,busca_X)

intervalos_X_M = 0 #np.zeros((busca_X,4)) #4 colunas pois uma dela identifica o M
linha_auxiliar = 0

#contador para linhas da matriz intervalos_X_M
linha = 0

#loop que identifica os X_j por M_i na direcao x
for i in range(busca_M):
    
    for j in range(busca_X):
        
        #if (i == 0):
        #    a=X_equa_z[j]
        #    a_max=9
        #    a_min=5
            
        #    if (a < a_max) & (a > a_min):
                
        #        #print(i,"# M: (x,y,z)=",x_M[i],",",y_M[i],",",z_M[i])
        #        #print(j,"# X: (x,y,z)=",X_equa_x[j],",",X_equa_y[j],",",X_equa_z[j],"\n")

        #        #print(" -> (x_X) < (x_M max):",X_equa_x[j],"<",(x_M[i] + dist_max_x),"?")
        #        #print(" --> (y_X) < (y_M max):",X_equa_y[j],"<",(y_M[i] + dist_max_y),"?\n")
                
        #        #print(" -> (x_X) > (x_M min):",X_equa_x[j],">",(x_M[i] - dist_max_x),"?")
        #        #print(" --> (y_X) > (y_M min):",X_equa_y[j],">",(y_M[i] - dist_max_y),"?\n")

        #        #print(" ---> (z_X) < (z_M max):",X_equa_z[j],"<",(z_M[i] + dist_max_z),"?")
        #        #print(" ---> (z_X) > (z_M min):",X_equa_z[j],">",(z_M[i] - dist_max_z),"?\n")
        
        if (X_equa_x[j] < (x_M[i] + dist_max_x)) and (X_equa_y[j] < (y_M[i] + dist_max_y)):
            
            #se o x de X for maior que o intervalo minimo
            if (X_equa_x[j] > (x_M[i] - dist_max_x)) and (X_equa_y[j] > (y_M[i] - dist_max_y)):
                
                #se o z de X for menor que o limite maximo do intervalo
                #if (X_equa_z[j] < (z_M[i] + dist_max_z)) and (X_equa_z[j] < (z_M[i] + dist_max_z)):
                if (X_equa_z[j] < (z_M[i] + dist_max_z)):
                    
                    #se o z de X for maior que o intervalo minimo
                    #if (X_equa_z[j] > (z_M[i] - dist_max_z)) and (X_equa_z[j] > (z_M[i] - dist_max_z)):
                    if (X_equa_z[j] > (z_M[i] - dist_max_z)):
                
                        linha_auxiliar = np.array([[X_equa_x[j], X_equa_y[j], X_equa_z[j], i]])
                        ##print(linha_auxiliar)
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            #intervalos_X_M[0][:] = linha_auxiliar
                            intervalos_X_M = linha_auxiliar
                            ##print(intervalos_X_M)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            intervalos_X_M = np.append(intervalos_X_M,linha_auxiliar, axis=0)
                            ##print(intervalos_X_M)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1


##print('\n Matirix de M: \n',x_M,'\n',y_M,'\n',z_M) 
##print('\n Matirix de X: \n',X_equa_x,'\n',X_equa_y,'\n',X_equa_z)
##print('\n Matirix coordenadas de M (x,y,z,)(Direct): \n',coord_atom_M) 
##print(len(intervalos_X_M))
#print('\n Matirix de X separados por M (x,y,z,!numero de M!): \n',intervalos_X_M) 


# In[24]:


#calculo das distancias X-M dentro do plano equatorial
#print('Distancias no plano equatorial\n')
#vetores e variaveis essenciais para o loop
compr_Xs = len(intervalos_X_M)
##print(compr_Xs)
dist_X_M_equa = np.zeros((compr_Xs,2)) # distancia X-M no plano equatorial

#transformando as coordenadas x_M de Direct para cartesianas
#as coordenadas
coord_Ms = np.transpose(np.vstack([x_M,y_M,z_M]))

#chamando a matriz de coordenadas que foram transformadas de fracional para cartesiana
coord_M_cart = coord_Ms

#print('Matriz de coordenadas de posicao dos M em coordenadas cartesianas:\n',coord_M_cart)

#Loop que chama as coordenadas do atomo Xi (linha da matriz "Matirix de X separados por M" em "intervalos_X_M")
for i in range(compr_Xs):
    
    #chamando as coordenadas do atomo de M - corresponde a coordenada dentro da ultima coluna de cada linha
    #ponto final (X_i)
    x_f = intervalos_X_M[i][0]
    y_f = intervalos_X_M[i][1]
    z_f = intervalos_X_M[i][2]
    
    ##print('Coordenada do pto. final X_',i,': (',x_f,';',y_f,';',z_f,')')
        
    #identificando o atomo M_atom_M_j
    atom_M_j = int(intervalos_X_M[i][3])
    
    #ponto inicial (M_atom_M_j)
    x_i = coord_M_cart[atom_M_j][0] #x_i = x_M[atom_M_j]
    y_i = coord_M_cart[atom_M_j][1] #y_i = y_M[atom_M_j]
    z_i = coord_M_cart[atom_M_j][2] #z_i = z_M[atom_M_j]
    
    
    ##print('Coordenada do pto. inicial M_',atom_M_j,': (',x_i,';',y_i,';',z_i,') \n')

    
    #calculando a distancia (X_i <-> M_atom_M_j)
    
    #identificando o atomo M a qual a distancia seria atrelada
    dist_X_M_equa[i][0] = intervalos_X_M[i][3]
    
    #calculando a distancia
    dist_X_M_equa[i][1] = math.sqrt( (np.power(x_f-x_i,2)) + (np.power(y_f-y_i,2)) + (np.power(z_f-z_i,2)) )
    
    ##print('Distancia equatorial para X_',i,' e M_ ',int(dist_X_M_equa[i][0]),': ',dist_X_M_equa[i][1])
    
#print('\nMatriz das distancias X-M no plano equatorial (Cartesiano): \n', dist_X_M_equa)
#print('\n Onde a primeira coluna identifica o atomo M atrelado a distancia calculada na segunda coluna.')


# In[25]:


#Procurando o maior e o menor valor de X-M calculado (atribuido a camada e independente de qual for o atomo M)

#Separando a coluna de distancias e a transformando em vetor
dist_eq_aux = np.transpose(dist_X_M_equa)
dist_eq_aux = dist_eq_aux[1][:]
resultado_dist_eq_aux = max_min_array(dist_eq_aux)

#print('Vetor das distancias X-M no plano equatorial  (Cartesiano): \n',dist_eq_aux)

dist_eq_max = resultado_dist_eq_aux[0]

#print('\n Distancia maxima encontrada: ',dist_eq_max)

dist_eq_min = resultado_dist_eq_aux[1]

#print('\n Distancia minima encontrada: ',dist_eq_min)


# In[26]:


#separando coordenadas dos X na direcao apical
##print('Coordenadas dos X topo, meio de camada e vale (do vetor de indices crescentes)')
#variaveis para o loop:
# X_api_index: matriz calculada anteriormente
tamanho = len(X_api_index)
X_api_x = np.zeros((tamanho))
X_api_y = np.zeros((tamanho))
X_api_z = np.zeros((tamanho))

##print(X_api_x)

#variaveis auxiliares
ii = 0

##print('\n Vetor de indices de X apical:\n ',X_api_index)
##print('\n Vetor de x de X (todos):\n ',x_X)
##print('\n Vetor de y de X (todos):\n ',y_X)
##print('\n Vetor de z de X (todos):\n ',z_X)

#loop decrescente pelos vetores posicao correspondentes a todos os X(x,y,z)
for i in range(tamanho):
    
    ##print('\n --> valor de i nessa iteracao: ',i,'\n   ---> conteudo do vetor nessa posicao/indice: ',int(X_api_index[i]))
    
    #adicionar as coordenadas os X do topo e do vale
    ii = int(X_api_index[i])
    X_api_x[i] = x_X[ii]
    ##print('\n Vetor de coord x de X apical: \n',X_api_x[i])
    X_api_y[i] = y_X[ii]
    ##print('\n Vetor de coord y de X apical: \n',X_equa_y[i])
    X_api_z[i] = z_X[ii]
    ##print('\n Vetor de coord z de X apical: \n',X_equa_z[i])
    
##print('\n Coordenadas dos X topo e vale na direcao apical (Direct!): \n --> [x]:\n',X_api_x,'\n --> [y]:\n',X_api_y,'\n --> [z]:\n',X_api_z)       


# In[27]:


#calculo das distancias X-M na direcao apical
#print('Coordenadas dos X topo, meio de camada e vale por proximidade de M')
##Para uma camada (n=3)

#Todos os atomos X e M respeitam uma proporcao pois trata-se de uma celula
#calculando um intervalo medio X-M na direcao x e y

#vetores e variaveis essenciais para o loop
# X_api_index: seria o vetor de X pertencentes a direcao apical (apenas indices)
# X_api_x (ou X_api_y)

tolerancia_x = 1.5 #somado e subtraido da posicao x de M
tolerancia_y = 1.5 #somado e subtraido da posicao x de M
tolerancia_z = 1.5 #somado e subtraido ao z_medio de M da primeira e segunda fileira de Ms

z_fileira_1 = z_medio_por_camada[0][1] #z medio de Ms na fileira 1
z_fileira_2 = z_medio_por_camada[1][1] #z medio de Ms na fileira 2
z_fileira_3 = z_medio_por_camada[2][1] #z medio de Ms na fileira 3
z_fileira_4 = z_medio_por_camada[3][1] #z medio de Ms na fileira 4
z_fileira_5 = z_medio_por_camada[4][1] #z medio de Ms na fileira 4

#print(z_medio_por_camada)
#z_intervalo_M_M = z_fileira_1 - z_fileira_2 - z_fileira_3 #intervalo entre fileiras; serve como alcance
z_intervalo_M_M = z_fileira_5

#print('Intervalo',z_intervalo_M_M)

busca_M = len(x_M)
busca_X = len(X_api_x)
interv_api_X_M = 0 #np.zeros((busca_X,4)) #4 colunas pois uma dela identifica o M
linha_auxiliar = 0

#contador para linhas da matriz interv_api_X_M
linha = 0

#loop que identifica os X_j por M_i na direcao x
for i in range(busca_M):
    
    ##print('\n ////////////////////Iteracao para M_',i,'////////////////////\n ')
    
        
    ##print('==========Iteracao para X_',j,'==========')
    #identificando se o M da camada pertencente a fileira de cima ou de baixo 
    
    # atomos na camada n1
    if (z_M[i] > (z_fileira_1 - tolerancia_z)) and (z_M[i] < (z_fileira_1 + tolerancia_z)):
            
        #M pertence a fileira superior
        ##print('\n | Camada n1: M_',i)
        ##print('z=',z_M[i],'> max=',(z_fileira_1 - tolerancia_z))
        ##print('z=',z_M[i],'< min=',(z_fileira_1 + tolerancia_z))
            
        #loop que checa se x_X esta dentro do intervalo (x_M + dist_max ; x_M - dist_max)
        for j in range(busca_X):
            
            #Identificando se X é proximo de M na direcao x
            if (X_api_x[j] > (x_M[i] - tolerancia_x)) and (X_api_x[j] < (x_M[i] + tolerancia_x)):
                
                ##print('\n * X_',j)
                ##print('x=',X_api_x[j],'> min=',(x_M[i] - tolerancia_x))
                ##print('x=',X_api_x[j],'< max=',(x_M[i] + tolerancia_x))
                ##print('OK!')
                ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                
                ########um if em y aqui #############################################################
                if (X_api_y[j] > (y_M[i] - tolerancia_y)) and (X_api_y[j] < (y_M[i] + tolerancia_y)):
                
                    #Identificando se X esta contido na primeira fileira pela altura z
                    if (X_api_z[j] > (z_M[i]-z_intervalo_M_M)) and (X_api_z[j] < (z_M[i]+z_intervalo_M_M)):

                        ##print('\n * X_',j)
                        ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                        ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                        ##print('OK!')

                        linha_auxiliar = np.array([[X_api_x[j], X_api_y[j], X_api_z[j], i]])
                        ##print(linha_auxiliar)
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            interv_api_X_M = linha_auxiliar
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            interv_api_X_M = np.append(interv_api_X_M,linha_auxiliar, axis=0)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1
    # atomos na camada n2
    if (z_M[i] > (z_fileira_2 - tolerancia_z)) and (z_M[i] < (z_fileira_2 + tolerancia_z)):
        #M pertence a fileira inferior
        ##print('\n || Camada n2: M_',i)
        ##print('z=',z_M[i],'> max=',(z_fileira_2 - tolerancia_z))
        ##print('z=',z_M[i],'< min=',(z_fileira_2 + tolerancia_z))
            
        #loop que checa se x_X esta dentro do intervalo (x_M + dist_max ; x_M - dist_max)
        for j in range(busca_X):
            
            #Identificando se X é proximo de M na direcao x
            if (X_api_x[j] > (x_M[i] - tolerancia_x)) and (X_api_x[j] < (x_M[i] + tolerancia_x)):
                
                ##print('\n * X_',j)
                ##print('x=',X_api_x[j],'> min=',(x_M[i] - tolerancia_x))
                ##print('x=',X_api_x[j],'< max=',(x_M[i] + tolerancia_x))
                ##print('OK!')
                
                ########um if em y aqui #############################################################
                if (X_api_y[j] > (y_M[i] - tolerancia_y)) and (X_api_y[j] < (y_M[i] + tolerancia_y)): 
                
                    #Identificando se X esta contido na primeira fileira pela altura z
                    if (X_api_z[j] > (z_M[i] - z_intervalo_M_M)) and (X_api_z[j] < (z_M[i] + z_intervalo_M_M)):

                        ##print('\nM_',i)
                        ##print('\n * X_',j)
                        ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                        ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                        ##print('OK!')


                        linha_auxiliar = np.array([[X_api_x[j], X_api_y[j], X_api_z[j], i]])
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            interv_api_X_M = linha_auxiliar
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            interv_api_X_M = np.append(interv_api_X_M,linha_auxiliar, axis=0)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1
                        
    # atomos na camada n3                   
    if (z_M[i] > (z_fileira_3 - tolerancia_z)) and (z_M[i] < (z_fileira_3 + tolerancia_z)):
        #M pertence a fileira inferior
        ##print('\n ||| Camada n3: M_',i)
        ##print('z=',z_M[i],'> max=',(z_fileira_2 - tolerancia_z))
        ##print('z=',z_M[i],'< min=',(z_fileira_2 + tolerancia_z))
            
        #loop que checa se x_X esta dentro do intervalo (x_M + dist_max ; x_M - dist_max)
        for j in range(busca_X):
            
            #Identificando se X é proximo de M na direcao x
            if (X_api_x[j] > (x_M[i] - tolerancia_x)) and (X_api_x[j] < (x_M[i] + tolerancia_x)):
                
                ##print('\n * X_',j)
                ##print('x=',X_api_x[j],'> min=',(x_M[i] - tolerancia_x))
                ##print('x=',X_api_x[j],'< max=',(x_M[i] + tolerancia_x))
                ##print('OK!')
                
                ########um if em y aqui
                if (X_api_y[j] > (y_M[i] - tolerancia_y)) and (X_api_y[j] < (y_M[i] + tolerancia_y)):
                    
                    #Identificando se X esta contido na primeira fileira pela altura z
                    if (X_api_z[j] > (z_M[i] - z_intervalo_M_M)) and (X_api_z[j] < (z_M[i] + z_intervalo_M_M)):

                        ##print('\nM_',i)
                        ##print('\n * X_',j)
                        ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                        ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                        ##print('OK!')

                        linha_auxiliar = np.array([[X_api_x[j], X_api_y[j], X_api_z[j], i]])
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            interv_api_X_M = linha_auxiliar
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            interv_api_X_M = np.append(interv_api_X_M,linha_auxiliar, axis=0)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1  
    # atomos na camada n4                    
    if (z_M[i] > (z_fileira_4 - tolerancia_z)) and (z_M[i] < (z_fileira_4 + tolerancia_z)):
        #M pertence a fileira inferior
        ##print('\n ||| Camada n3: M_',i)
        ##print('z=',z_M[i],'> max=',(z_fileira_2 - tolerancia_z))
        ##print('z=',z_M[i],'< min=',(z_fileira_2 + tolerancia_z))
            
        #loop que checa se x_X esta dentro do intervalo (x_M + dist_max ; x_M - dist_max)
        for j in range(busca_X):
            
            #Identificando se X é proximo de M na direcao x
            if (X_api_x[j] > (x_M[i] - tolerancia_x)) and (X_api_x[j] < (x_M[i] + tolerancia_x)):
                
                ##print('\n * X_',j)
                ##print('x=',X_api_x[j],'> min=',(x_M[i] - tolerancia_x))
                ##print('x=',X_api_x[j],'< max=',(x_M[i] + tolerancia_x))
                ##print('OK!')
                
                
                if (X_api_y[j] > (y_M[i] - tolerancia_y)) and (X_api_y[j] < (y_M[i] + tolerancia_y)):
                    
                    #Identificando se X esta contido na primeira fileira pela altura z
                    if (X_api_z[j] > (z_M[i] - z_intervalo_M_M)) and (X_api_z[j] < (z_M[i] + z_intervalo_M_M)):

                        ##print('\nM_',i)
                        ##print('\n * X_',j)
                        ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                        ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                        ##print('OK!')

                        linha_auxiliar = np.array([[X_api_x[j], X_api_y[j], X_api_z[j], i]])
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            interv_api_X_M = linha_auxiliar
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            interv_api_X_M = np.append(interv_api_X_M,linha_auxiliar, axis=0)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1
    # atomos na camada n5                    
    if (z_M[i] > (z_fileira_5 - tolerancia_z)) and (z_M[i] < (z_fileira_5 + tolerancia_z)):
        #M pertence a fileira inferior
        ##print('\n ||| Camada n3: M_',i)
        ##print('z=',z_M[i],'> max=',(z_fileira_2 - tolerancia_z))
        ##print('z=',z_M[i],'< min=',(z_fileira_2 + tolerancia_z))
            
        #loop que checa se x_X esta dentro do intervalo (x_M + dist_max ; x_M - dist_max)
        for j in range(busca_X):
            
            #Identificando se X é proximo de M na direcao x
            if (X_api_x[j] > (x_M[i] - tolerancia_x)) and (X_api_x[j] < (x_M[i] + tolerancia_x)):
                
                ##print('\n * X_',j)
                ##print('x=',X_api_x[j],'> min=',(x_M[i] - tolerancia_x))
                ##print('x=',X_api_x[j],'< max=',(x_M[i] + tolerancia_x))
                ##print('OK!')
                
                
                if (X_api_y[j] > (y_M[i] - tolerancia_y)) and (X_api_y[j] < (y_M[i] + tolerancia_y)):
                    
                    #Identificando se X esta contido na primeira fileira pela altura z
                    if (X_api_z[j] > (z_M[i] - z_intervalo_M_M)) and (X_api_z[j] < (z_M[i] + z_intervalo_M_M)):

                        ##print('\nM_',i)
                        ##print('\n * X_',j)
                        ##print('z=',X_api_z[j],'> min=',(z_M[i]-z_intervalo_M_M))
                        ##print('z=',X_api_z[j],'< max=',(z_M[i]+z_intervalo_M_M))
                        ##print('OK!')

                        linha_auxiliar = np.array([[X_api_x[j], X_api_y[j], X_api_z[j], i]])
                        #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)

                        if linha == 0:

                            #declarando a linha "0" matriz
                            interv_api_X_M = linha_auxiliar
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha += 1

                        #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                        elif linha > 0:

                            interv_api_X_M = np.append(interv_api_X_M,linha_auxiliar, axis=0)
                            #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                            linha+=1
                            
#print('\n Matriz de X separados por M na apical (x,y,z,!numero de M!)(Direct!): \n',interv_api_X_M) 


# In[28]:


#calculo das distancias X-M na direcao apical
#print('Distancias na direcao apical')
#vetores e variaveis essenciais para o loop
comprim_Xs = len(interv_api_X_M)
##print(comprim_Xs)
dist_X_M_api = np.zeros((comprim_Xs,2)) # distancia X-M no plano equatorial
x_f = 0
y_f = 0
z_f = 0
x_i = 0
y_i = 0
z_i = 0
##print(x_M)

#print('Matriz de posicao X separados por M (cartesiana): \n',interv_api_X_M)
#print('\nMatriz de posicao M (cartesiana): \n',coord_M_cart)

#Loop que chama as coordenadas do atomo Xi (linha da matriz "Matirix de X separados por M" em "interv_api_X_M")
for i in range(comprim_Xs):
    
    #chamando as coordenadas do atomo de M - corresponde a coordenada dentro da ultima coluna de cada linha
    #ponto final (X_i)
    x_f = interv_api_X_M[i][0]
    y_f = interv_api_X_M[i][1]
    z_f = interv_api_X_M[i][2]
    
    ##print('Coordenada do pto. final X_',i,': (',x_f,';',y_f,';',z_f,')')
        
    #identificando o atomo M_atom_M_j
    atom_M_j = int(interv_api_X_M[i][3])
    
    #ponto inicial (M_atom_M_j)
    x_i = coord_M_cart[atom_M_j][0] #x_M[atom_M_j]
    y_i = coord_M_cart[atom_M_j][1] #y_M[atom_M_j]
    z_i = coord_M_cart[atom_M_j][2] #z_M[atom_M_j]
    
    ##print('Coordenada do pto. inicial M_',atom_M_j,': (',x_i,';',y_i,';',z_i,') \n')

    
    #calculando a distancia (X_i <-> M_atom_M_j)
    
    #identificando o atomo M a qual a distancia seria atrelada
    dist_X_M_api[i][0] = interv_api_X_M[i][3]
    
    #calculando a distancia
    dist_X_M_api[i][1] = math.sqrt( (np.power(x_f-x_i,2)) + (np.power(y_f-y_i,2)) + (np.power(z_f-z_i,2)) )
    
    ##print('Distancia equatorial para X_',i,' e M_ ',int(dist_X_M_api[i][0]),': ',dist_X_M_api[i][1])
    
##print('Matriz das distancias X-M na direcao apical (Direct!): \n', dist_X_M_api)
#print('\nMatriz das distancias X-M na direcao apical (Cartesiana): \n', dist_X_M_api)
#print('\n Onde a primeira coluna identifica o atomo M atrelado a distancia calculada na segunda coluna.')


# In[29]:


#Procurando o maior e o menor valor de X-M calculado (atribuido a camada e independente de qual for o atomo M)

#Separando a coluna de distancias e a transformando em vetor
dist_api_aux = np.transpose(dist_X_M_api)
dist_api_aux = dist_api_aux[1][:]
resultado_dist_api_aux = max_min_array(dist_api_aux)

#print('Vetor das distancias X-M na direcao apical (Cartesiana): \n',dist_api_aux)

dist_api_max = resultado_dist_api_aux[0]

#print('\n Distancia maxima encontrada: ',dist_api_max)

dist_api_min = resultado_dist_api_aux[1]

#print('\n Distancia minima encontrada: ',dist_api_min)


# In[30]:


#Calculo dos angulos M-X-M (plano equatorial)

##Caso para (n=2)

#Encontrar quais M tem o mesmo X em comum. Olhar na matriz "intervalos_X_M"
    #se os valores de (x,y,z) forem iguais: calcular o angulo pois temos ex. M_1 - X_3 - M_3
    #comparar partindo de i para os indices da ultima coluna. Ex.: para M_0: qual X em comum com M_1, M_2, etc.
    
    
#vetores e variaveis essenciais para o loop
# x_M, y_M e z_M: vetores posicao dos M
# intervalos_X_M: matriz de coordenadas de X separadas por vizinhos de M ( quarta coluna )
# busca_M: variavel com quantidade de M no sistema (<-- vetores posicao dos M)
busca_X_M = len(intervalos_X_M)
atomo_M_i = 0 #identifica o atomo M inicial para o arco M-X-M
atomo_M_f = 0 #identifica o atomo M final para o arco M-X-M
# dist_X_M_equa #distancias X-M
linha = 0
M_X_M_equa = 0
##print(dist_eq_aux)
#print('(Plano equatorial) => Matriz de coordenadas de X identificando o atomo M vizinho ( x, y, z, M_i): \n',intervalos_X_M,'\n')

#Procura por M_h vizinhos de um M_i (mesma camada) em um raio de alcance de 1 unidade de celula
for i in range(busca_X_M):
    
    ##print('Loop primario --> X_',i,' vizinho de M_',int(intervalos_X_M[i][3]))
    
    #identificando qual seria o atomo M_i
    atomo_M_i = int(intervalos_X_M[i][3])
    
    #quais outros X (ou outras linhas) possuem o mesmo M como vizinho?
    #procura-se entao por nas outras linhas (ignorando a linha i)
    for h in range(busca_X_M):
        
        #ignorando a mesma linha (M iguais)
        if (i != h):
            
            ##print('  >>Loop secundario --> X_',h,' vizinho de M_',int(intervalos_X_M[h][3]),'\n')
            
            #identificando qual seria o atomo M_h para o arco com M_i
            atomo_M_f = int(intervalos_X_M[h][3])
            
            #em M diferentes:
            if (atomo_M_i != atomo_M_f):
                
                ##print('M diferente localizado')
                
                #comparar atomos X de atomo_M_i com atomo_M_f. 
                #Se forem iguais (apenas uma possibilidade) entao calcula-se o angulo.
                #comparando vizinhos de atomo_M_i com atomo_M_f:
                
                #coordenadas de X de atomo_M_i
                x_vizinho_M_i = intervalos_X_M[i][0]
                y_vizinho_M_i = intervalos_X_M[i][1]
                z_vizinho_M_i = intervalos_X_M[i][2]
                
                #coordenadas de X de atomo_M_h
                x_vizinho_M_h = intervalos_X_M[h][0]
                y_vizinho_M_h = intervalos_X_M[h][1]
                z_vizinho_M_h = intervalos_X_M[h][2]
                
                #comparando a direcao x
                if (x_vizinho_M_i == x_vizinho_M_h):
                    
                    #comparando a direcao y
                    if (y_vizinho_M_i == y_vizinho_M_h):
                        
                        #comparando a direcao z
                        if (z_vizinho_M_i == z_vizinho_M_h):
                            
                            ##print('  --> M_',int(intervalos_X_M[i][3]),'- X_',i,'-M_',int(intervalos_X_M[h][3]))
                            ##print('  OU o mesmo que M_',int(intervalos_X_M[i][3]),'- X_',h,'-M_',int(intervalos_X_M[h][3]),'\n')
                            
                            #calculando o angulo
                            #dados para Mi
                            d_i = dist_X_M_equa[i][1]
                            
                            #dados para Mh
                            d_h = dist_X_M_equa[h][1]
                            
                            #distancia de M_i ate M_h
                            d_ih = math.sqrt( (np.power(d_i-d_h,2)) + (np.power(d_i-d_h,2)) + (np.power(d_i-d_h,2)) )
                            
                            #distancia de M_i ate M_h
                            delta_x_d_Mi_Mh = coord_M_cart[atomo_M_i][0]-coord_M_cart[atomo_M_f][0]
                            delta_y_d_Mi_Mh = coord_M_cart[atomo_M_i][1]-coord_M_cart[atomo_M_f][1]
                            delta_z_d_Mi_Mh = coord_M_cart[atomo_M_i][2]-coord_M_cart[atomo_M_f][2]
                            
                            d_ih = math.sqrt( delta_x_d_Mi_Mh**2 + delta_y_d_Mi_Mh**2 + delta_z_d_Mi_Mh**2 )
                            ##print("d_ih= ",d_ih)
                            ##print('d_',atomo_M_i,'=',d_i,'; d_',atomo_M_f,'=',d_h,'; d_{',atomo_M_i,',',atomo_M_f,'}=',d_ih)
                            
                            #calculando cosseno do angulo M-X-M
                            #cos_theta = ( np.power(d_i,2) + np.power(d_h,2) - np.power(d_ih,2) )/( 2*d_i*d_h ) 
                            aux_theta = (d_i**2 + d_h**2 - d_ih**2)/(2*d_i*d_h)
                            
                            theta = np.arccos(aux_theta)
                            ##print('Antes do arccos:',theta)
                            ##print('Cos(theta)=',cos_theta)
                            
                            #calculando arco (obtuso)
                            #theta = 180 - math.degrees(math.acos(cos_theta)) #corrigido com subtracao de 180 graus pois o "acos(theta)" calcula o angulo agudo.
                            theta = math.degrees(theta)
                            ##print('Em M_',int(intervalos_X_M[i][3]),'- X_',i,'-M_',int(intervalos_X_M[h][3]),': theta=',180-theta,'°')
                            
                                                       
                            #armazendo os thetas no formato ( M_i, X_i, X_h, Mh, theta)
                            linha_auxiliar = np.array([[intervalos_X_M[i][3], i, h, intervalos_X_M[h][3], theta]])
                            
                            ##print('Armazendo os thetas no formato ( M_i, X_i, X_h, M_h, theta):',linha_auxiliar)
                            
                            #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)
                            if linha == 0:

                                #declarando a linha "0" matriz
                                #interv__api_X_M[0][:] = linha_auxiliar
                                M_X_M_equa = linha_auxiliar
                                ##print(interv_api_X_M)

                                #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                                linha += 1

                            #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                            elif linha > 0:

                                M_X_M_equa = np.append(M_X_M_equa,linha_auxiliar, axis=0)
                                ##print(interv_api_X_M)

                                #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                                linha+=1
                            
    
#print('Angulos thetas no formato ( M_i, X_i, X_h, M_h, theta) (Cartesiana): \n',M_X_M_equa)
#print('Obs.: considerando que X_i e X_h sao o mesmo atomo;')
#print('Onde, na matriz de X vizinho de M, X_i vizinho de M_i e X_h vizinho de M_h.')


# In[31]:


#Procurando o maior e o menor valor de M-X-M calculado (plano equatorial)

#Separando a coluna de distancias e a transformando em vetor
ang_eq_aux = np.transpose(M_X_M_equa)
ang_eq_aux = ang_eq_aux[4][:]
resultado_ang_eq_aux = max_min_array(ang_eq_aux)

#print('Vetor das distancias X-M no plano equatorial  (Cartesiano): \n',ang_eq_aux)

ang_eq_max = resultado_ang_eq_aux[0]

#print('\n Distancia maxima encontrada: ',ang_eq_max)

ang_eq_min = resultado_ang_eq_aux[1]

#print('\n Distancia minima encontrada: ',ang_eq_min)


# In[32]:


#Calculo dos angulos M-X-M (direcao apical)

##Caso para (n=2)

#Encontrar quais M tem o mesmo X em comum. Olhar na matriz "intervalos_X_M"
    #se os valores de (x,y,z) forem iguais: calcular o angulo pois temos ex. M_1 - X_3 - M_3
    #comparar partindo de i para os indices da ultima coluna. Ex.: para M_0: qual X em comum com M_1, M_2, etc.
    
    
#vetores e variaveis essenciais para o loop
# x_M, y_M e z_M: vetores posicao dos M
# intervalos_X_M: matriz de coordenadas de X separadas por vizinhos de M ( quarta coluna )
# busca_M: variavel com quantidade de M no sistema (<-- vetores posicao dos M)
busca_X_M = len(interv_api_X_M)
atomo_M_i = 0 #identifica o atomo M inicial para o arco M-X-M
atomo_M_f = 0 #identifica o atomo M final para o arco M-X-M
# dist_X_M_equa #distancias X-M
linha = 0
M_X_M_api = 0
##print(dist_eq_aux)
#print('(Plano equatorial) => Matriz de coordenadas de X identificando o atomo M vizinho ( x, y, z, M_i): \n',interv_api_X_M,'\n')

#Procura por M_h vizinhos de um M_i (mesma camada) em um raio de alcance de 1 unidade de celula
for i in range(busca_X_M):
    
    ##print('Loop primario --> X_',i,' vizinho de M_',int(intervalos_X_M[i][3]))
    
    #identificando qual seria o atomo M_i
    atomo_M_i = int(interv_api_X_M[i][3])
    
    #quais outros X (ou outras linhas) possuem o mesmo M como vizinho?
    #procura-se entao por nas outras linhas (ignorando a linha i)
    for h in range(busca_X_M):
        
        #ignorando a mesma linha (M iguais)
        if (i != h):
            
            ##print('  >>Loop secundario --> X_',h,' vizinho de M_',int(intervalos_X_M[h][3]),'\n')
            
            #identificando qual seria o atomo M_h para o arco com M_i
            atomo_M_f = int(interv_api_X_M[h][3])
            
            #em M diferentes:
            if (atomo_M_i != atomo_M_f):
                
                ##print('M diferente localizado')
                
                #comparar atomos X de atomo_M_i com atomo_M_f. 
                #Se forem iguais (apenas uma possibilidade) entao calcula-se o angulo.
                #comparando vizinhos de atomo_M_i com atomo_M_f:
                
                #coordenadas de X de atomo_M_i
                x_vizinho_M_i = interv_api_X_M[i][0]
                y_vizinho_M_i = interv_api_X_M[i][1]
                z_vizinho_M_i = interv_api_X_M[i][2]
                
                #coordenadas de X de atomo_M_h
                x_vizinho_M_h = interv_api_X_M[h][0]
                y_vizinho_M_h = interv_api_X_M[h][1]
                z_vizinho_M_h = interv_api_X_M[h][2]
                
                #comparando a direcao x
                if (x_vizinho_M_i == x_vizinho_M_h):
                    
                    #comparando a direcao y
                    if (y_vizinho_M_i == y_vizinho_M_h):
                        
                        #comparando a direcao z
                        if (z_vizinho_M_i == z_vizinho_M_h):
                            
                            ##print('  --> M_',int(intervalos_X_M[i][3]),'- X_',i,'-M_',int(intervalos_X_M[h][3]))
                            ##print('  OU o mesmo que M_',int(intervalos_X_M[i][3]),'- X_',h,'-M_',int(intervalos_X_M[h][3]),'\n')
                            
                            #calculando o angulo
                            #dados para Mi
                            d_i = dist_X_M_api[i][1]
                            
                            #dados para Mh
                            d_h = dist_X_M_api[h][1]
                            
                            #distancia de M_i ate M_h
                            d_ih = math.sqrt( (np.power(d_i-d_h,2)) + (np.power(d_i-d_h,2)) + (np.power(d_i-d_h,2)) )
                            
                            #distancia de M_i ate M_h
                            delta_x_d_Mi_Mh = coord_M_cart[atomo_M_i][0]-coord_M_cart[atomo_M_f][0]
                            delta_y_d_Mi_Mh = coord_M_cart[atomo_M_i][1]-coord_M_cart[atomo_M_f][1]
                            delta_z_d_Mi_Mh = coord_M_cart[atomo_M_i][2]-coord_M_cart[atomo_M_f][2]
                            
                            d_ih = math.sqrt( delta_x_d_Mi_Mh**2 + delta_y_d_Mi_Mh**2 + delta_z_d_Mi_Mh**2 )
                            ##print("d_ih= ",d_ih)
                            ##print('d_',atomo_M_i,'=',d_i,'; d_',atomo_M_f,'=',d_h,'; d_{',atomo_M_i,',',atomo_M_f,'}=',d_ih)
                            
                            #calculando cosseno do angulo M-X-M
                            #cos_theta = ( np.power(d_i,2) + np.power(d_h,2) - np.power(d_ih,2) )/( 2*d_i*d_h ) 
                            aux_theta = (d_i**2 + d_h**2 - d_ih**2)/(2*d_i*d_h)
                            
                            theta = np.arccos(aux_theta)
                            ##print('Antes do arccos:',theta)
                            ##print('Cos(theta)=',cos_theta)
                            
                            #calculando arco (obtuso)
                            #theta = 180 - math.degrees(math.acos(cos_theta)) #corrigido com subtracao de 180 graus pois o "acos(theta)" calcula o angulo agudo.
                            theta = math.degrees(theta)
                            ##print('Em M_',int(intervalos_X_M[i][3]),'- X_',i,'-M_',int(intervalos_X_M[h][3]),': theta=',180-theta,'°')
                            
                                                       
                            #armazendo os thetas no formato ( M_i, X_i, X_h, Mh, theta)
                            linha_auxiliar = np.array([[interv_api_X_M[i][3], i, h, interv_api_X_M[h][3], theta]])
                            
                            ##print('Armazendo os thetas no formato ( M_i, X_i, X_h, M_h, theta):',linha_auxiliar)
                            
                            #quando nao ha nada no vetor inicialmente (vetor ainda nao criado)
                            if linha == 0:

                                #declarando a linha "0" matriz
                                #interv__api_X_M[0][:] = linha_auxiliar
                                M_X_M_api = linha_auxiliar
                                ##print(interv_api_X_M)

                                #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                                linha += 1

                            #quanto o vetor possui algum termo (indice maior que um dentro do vetor)
                            elif linha > 0:

                                M_X_M_api = np.append(M_X_M_api,linha_auxiliar, axis=0)
                                ##print(interv_api_X_M)

                                #passando para a proxima linha (referente ao proximo atomo X a ser encontrado)
                                linha+=1
                            
    
#print('Angulos thetas no formato ( M_i, X_i, X_h, M_h, theta) (Cartesiana): \n',M_X_M_api)
#print('Obs.: considerando que X_i e X_h sao o mesmo atomo;')
#print('Onde, na matriz de X vizinho de M, X_i vizinho de M_i e X_h vizinho de M_h.')


# In[33]:


#Procurando o maior e o menor valor de M-X-M calculado (direcao apical)

#Separando a coluna de distancias e a transformando em vetor
ang_api_aux = np.transpose(M_X_M_api)
ang_api_aux = ang_api_aux[4][:]
resultado_ang_api_aux = max_min_array(ang_api_aux)

#print('Vetor das distancias X-M na direcao apical (Cartesiana): \n',ang_api_aux)

ang_api_max = resultado_ang_api_aux[0]

#print('\n Distancia maxima encontrada: ',ang_api_max)

ang_api_min = resultado_ang_api_aux[1]

#print('\n Distancia minima encontrada: ',ang_api_min)


# In[34]:


#separando dador por camadas inter-face (intf) e entre-camada (entc)

#separando Ms por camada (identificacao)
camada_n1 = M_camada_index[0][:]
camada_n2 = M_camada_index[1][:]
camada_n3 = M_camada_index[2][:]
camada_n4 = M_camada_index[3][:]
camada_n5 = M_camada_index[4][:]

###X-M

#vetores e matrizes utilizados
# M_camada_index #matriz onde cada linha representa uma camada. Contem indices de M .
# dist_X_M_equa #matriz de X-M calculada no plano equatorial com uma coluna de indices de M vizinhos
# dist_X_M_api #matriz de X-M calculada na distancia apical com uma coluna de indices de M vizinhos

no_X_M_equa = len(dist_X_M_equa)
##print('quantidade de distancias calculadas no plano equatorial: ',no_X_M_equa)

no_X_M_api = len(dist_X_M_api)
##print('quantidade de distancias calculadas na direcao apical: ',no_X_M_api)
##print('Numero de Ms por camada: ',colunas_M) # quantidade de Ms por camada

c_intf_equa = 0 #contador para a quantidade de valores no plano equatorial
c_intf_api = 0 #contador para a quantidade de valores na direcao apical

c_entc_equa = 0 #contador para a quantidade de valores no plano equatorial
c_entc_api = 0 #contador para a quantidade de valores na direcao apical

#inter-face (intf) --> n1 e n5

#loop que passa pelas linhas 0 e 4
for i in range(n_camada):
    
    ##print('n',i+1)
    
    #pulando as camadas n2 (i=1), n3 (i=2) e n4 (i=3)
    if (i != 1) and (i != 2):
        
        if (i != 3):
        
            ##print('n',i+1)


            #loop que passa pelos Ms de cada camada
            for j in range(colunas_M):

                ##print('Camada n',i+1,'>> M_',int(M_camada_index[i][j]))

                #loop que passa pelas linhas dist_X_M_equa:
                for g in range(no_X_M_equa):

                    #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                    #aloca-se a coluna dist_X_M_equa[g][1]=dist
                    if (M_camada_index[i][j] == dist_X_M_equa[g][0]):


                        linha_auxiliar = dist_X_M_equa[g][1]

                        #criando vetor de distancias
                        if c_intf_equa == 0:

                            #declarando a linha "0" matriz
                            intf_X_M_equa = linha_auxiliar
                            c_intf_equa += 1

                        elif c_intf_equa > 0:

                            #adicionando valores
                            intf_X_M_equa = np.append(intf_X_M_equa,linha_auxiliar)
                            c_intf_equa += 1

                #loop que passa pelas linhas dist_X_M_api:
                for g in range(no_X_M_api):

                    #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                    #aloca-se a coluna dist_X_M_equa[g][1]=dist
                    if (M_camada_index[i][j] == dist_X_M_api[g][0]):


                        linha_auxiliar = dist_X_M_api[g][1]

                        #criando vetor de distancias
                        if c_intf_api == 0:

                            #declarando a linha "0" matriz
                            intf_X_M_api = linha_auxiliar
                            c_intf_api += 1

                        elif c_intf_api > 0:

                            #adicionando valores
                            intf_X_M_api = np.append(intf_X_M_api,linha_auxiliar)
                            c_intf_api += 1



#print('\nDistancias X-M de inter-face')
#print('\n* no plano equatorial: \n',intf_X_M_equa)
#print('\n* na direcao apical: \n',intf_X_M_api)
      
#entre-camada (entc) --> n2 e n3

#loop que passa pela linhas 1
for i in range(n_camada):
    
    ##print('n',i+1)
    
    #pulando a camada n2 (i=1)
    if (i != 0) and (i != 4):
        
        ##print('n',i+1)
        
        
        #loop que passa pelos Ms de cada camada
        for j in range(colunas_M):
            
            ##print('Camada n',i+1,'>> M_',int(M_camada_index[i][j]))
            
            #loop que passa pelas linhas dist_X_M_equa:
            for g in range(no_X_M_equa):
                
                #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                #aloca-se a coluna dist_X_M_equa[g][1]=dist
                if (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                    
                    
                    linha_auxiliar = dist_X_M_equa[g][1]

                    #criando vetor de distancias
                    if c_entc_equa == 0:

                        #declarando a linha "0" matriz
                        entc_X_M_equa = linha_auxiliar
                        c_entc_equa += 1

                    elif c_entc_equa > 0:

                        #adicionando valores
                        entc_X_M_equa = np.append(entc_X_M_equa,linha_auxiliar)
                        c_entc_equa += 1
                        
            #loop que passa pelas linhas dist_X_M_api:
            for g in range(no_X_M_api):
                
                #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                #aloca-se a coluna dist_X_M_equa[g][1]=dist
                if (M_camada_index[i][j] == dist_X_M_api[g][0]):
                    
                    
                    linha_auxiliar = dist_X_M_api[g][1]

                    #criando vetor de distancias
                    if c_entc_api == 0:

                        #declarando a linha "0" matriz
                        entc_X_M_api = linha_auxiliar
                        c_entc_api += 1

                    elif c_entc_api > 0:

                        #adicionando valores
                        entc_X_M_api = np.append(entc_X_M_api,linha_auxiliar)
                        c_entc_api += 1

#print('\nDistancias X-M de entre-camada')
#print('\n* no plano equatorial: \n',entc_X_M_equa)
#print('\n* na direcao apical: \n',entc_X_M_api)


# In[35]:


### X-M
#print('Distancias X-M')
#separando maximos e minimos

#inter-face (intf) --> n1 e n5 

#print('\nPlano equatorial na inter-face')
resultado_intf_X_M_equa = max_min_array(intf_X_M_equa)

intf_dist_eq_max = resultado_intf_X_M_equa[0]

#print('\n * max: ',intf_dist_eq_max)

intf_dist_eq_min = resultado_intf_X_M_equa[1]

#print('\n * min: ',intf_dist_eq_min)

#print('\nDirecao apical na inter-face')
resultado_intf_X_M_api = max_min_array(intf_X_M_api)

intf_dist_api_max = resultado_intf_X_M_api[0]

#print('\n * max: ',intf_dist_api_max)

intf_dist_api_min = resultado_intf_X_M_api[1]

#print('\n * min: ',intf_dist_api_min)


#entre-camada (entc) --> n2 até n4

#print('\nPlano equatorial na entre-camada')
resultado_entc_X_M_equa = max_min_array(entc_X_M_equa)

entc_dist_eq_max = resultado_entc_X_M_equa[0]

#print('\n * max: ',entc_dist_eq_max)

entc_dist_eq_min = resultado_entc_X_M_equa[1]

#print('\n * min: ',entc_dist_eq_min)

#print('\nDirecao apical na entre-camada')
resultado_entc_X_M_api = max_min_array(entc_X_M_api)

entc_dist_api_max = resultado_entc_X_M_api[0]

#print('\n * max: ',entc_dist_api_max)

entc_dist_api_min = resultado_entc_X_M_api[1]

#print('\n * min: ',entc_dist_api_min)


# In[36]:


###M-X-M

#vetores e matrizes utilizados
# M_camada_index #matriz onde cada linha representa uma camada. Contem indices de M .
# M_X_M_equa #matriz de M-X-M calculada no plano equatorial com uma coluna de indices de M vizinhos (4)
# M_X_M_api #matriz de M-X-M calculada na distancia apical com uma coluna de indices de M vizinhos (4)

no_M_X_M_equa = len(M_X_M_equa)
##print('quantidade de distancias calculadas no plano equatorial: ',no_X_M_equa)

no_M_X_M_api = len(M_X_M_api)
##print('quantidade de distancias calculadas na direcao apical: ',no_X_M_api)
##print('Numero de Ms por camada: ',colunas_M) # quantidade de Ms por camada

c_intf_equa = 0 #contador para a quantidade de valores no plano equatorial
c_intf_api = 0 #contador para a quantidade de valores na direcao apical

c_entc_equa = 0 #contador para a quantidade de valores no plano equatorial
c_entc_api = 0 #contador para a quantidade de valores na direcao apical

#inter-face (intf) --> n1 e n3 

#loop que passa pelas linhas 0 e 2
for i in range(n_camada):
    
    ##print('n',i+1)
    
    #pulando as camadas n2 (i=1), n3 (i=2) e n4 (i=3)
    if (i != 1) and (i != 2):
        
        if (i != 3):
        
            ##print('n',i+1)


            #loop que passa pelos Ms de cada camada
            for j in range(colunas_M):

                ##print('Camada n',i+1,'>> M_',int(M_camada_index[i][j]))

                #loop que passa pelas linhas dist_X_M_equa:
                for g in range(no_M_X_M_equa):

                    #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                    #aloca-se a coluna dist_X_M_equa[g][1]=dist
                    if (M_camada_index[i][j] == M_X_M_equa[g][0]):


                        linha_auxiliar = M_X_M_equa[g][4]

                        #criando vetor de distancias
                        if c_intf_equa == 0:

                            #declarando a linha "0" matriz
                            intf_M_X_M_equa = linha_auxiliar
                            c_intf_equa += 1

                        elif c_intf_equa > 0:

                            #adicionando valores
                            intf_M_X_M_equa = np.append(intf_M_X_M_equa,linha_auxiliar)
                            c_intf_equa += 1

                #loop que passa pelas linhas dist_X_M_api:
                for g in range(no_M_X_M_api):

                    #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                    #aloca-se a coluna dist_X_M_equa[g][1]=dist
                    if (M_camada_index[i][j] == M_X_M_api[g][0]):


                        linha_auxiliar = M_X_M_api[g][4]

                        #criando vetor de distancias
                        if c_intf_api == 0:

                            #declarando a linha "0" matriz
                            intf_M_X_M_api = linha_auxiliar
                            c_intf_api += 1

                        elif c_intf_api > 0:

                            #adicionando valores
                            intf_M_X_M_api = np.append(intf_M_X_M_api,linha_auxiliar)
                            c_intf_api += 1
        


#print('Ângulos M-X-M de inter-face')
#print('* no plano equatorial: \n',intf_M_X_M_equa)
#print('* na direcao apical: \n',intf_M_X_M_api)
      
#entre-camada (entc) --> n2 e n3

#loop que passa pela linhas 1
for i in range(n_camada):
    
    ##print('n',i+1)
    
    #pulando a camada n2 (i=1)
    if (i != 0) and (i != 4):
        
        ##print('n',i+1)
        
        
        #loop que passa pelos Ms de cada camada
        for j in range(colunas_M):
            
            ##print('Camada n',i+1,'>> M_',int(M_camada_index[i][j]))
            
            #loop que passa pelas linhas dist_X_M_equa:
            for g in range(no_M_X_M_equa):
                
                #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                #aloca-se a coluna dist_X_M_equa[g][1]=dist
                if (M_camada_index[i][j] == M_X_M_equa[g][0]):
                    
                    
                    linha_auxiliar = M_X_M_equa[g][4]

                    #criando vetor de distancias
                    if c_entc_equa == 0:

                        #declarando a linha "0" matriz
                        entc_M_X_M_equa = linha_auxiliar
                        c_entc_equa += 1

                    elif c_entc_equa > 0:

                        #adicionando valores
                        entc_M_X_M_equa = np.append(entc_M_X_M_equa,linha_auxiliar)
                        c_entc_equa += 1
                        
            #loop que passa pelas linhas dist_X_M_api:
            for g in range(no_M_X_M_api):
                
                #se os indices M forem iguais (M_camada_index[i][j] == dist_X_M_equa[g][0]):
                #aloca-se a coluna dist_X_M_equa[g][1]=dist
                if (M_camada_index[i][j] == M_X_M_api[g][0]):
                    
                    
                    linha_auxiliar = M_X_M_api[g][4]

                    #criando vetor de distancias
                    if c_entc_api == 0:

                        #declarando a linha "0" matriz
                        entc_M_X_M_api = linha_auxiliar
                        c_entc_api += 1

                    elif c_entc_api > 0:

                        #adicionando valores
                        entc_M_X_M_api = np.append(entc_M_X_M_api,linha_auxiliar)
                        c_entc_api += 1

##print('Distancias X-M de entre-camada')
#print('* no plano equatorial: \n',entc_M_X_M_equa)
#print('* na direcao apical: \n',entc_M_X_M_api)


# In[37]:


### M-X-M
#print('Ângulos M-X-M')
#separando maximos e minimos

#inter-face (intf) --> n1 e n3 

#print('\nPlano equatorial na inter-face')
resultado_intf_M_X_M_equa = max_min_array(intf_M_X_M_equa) #funcao busca max e min do vetor

intf_ang_eq_max = resultado_intf_M_X_M_equa[0]

#print('\n * max: ',intf_ang_eq_max)

intf_ang_eq_min = resultado_intf_M_X_M_equa[1]

#print('\n * min: ',intf_ang_eq_min)

#print('\nDirecao apical na inter-face')
resultado_intf_M_X_M_api = max_min_array(intf_M_X_M_api) #funcao busca max e min do vetor

intf_ang_api_max = resultado_intf_M_X_M_api[0]

#print('\n * max: ',intf_ang_api_max)

intf_ang_api_min = resultado_intf_M_X_M_api[1]

#print('\n * min: ',intf_ang_api_min)


#entre-camada (entc) --> n2

#print('\nPlano equatorial na entre-camada')
resultado_entc_M_X_M_equa = max_min_array(entc_M_X_M_equa)

entc_ang_eq_max = resultado_entc_M_X_M_equa[0]

#print('\n * max: ',entc_ang_eq_max)

entc_ang_eq_min = resultado_entc_M_X_M_equa[1]

#print('\n * min: ',entc_ang_eq_min)


#print('\nDirecao apical na entre-camada')
resultado_entc_M_X_M_api = max_min_array(entc_M_X_M_api)

entc_ang_api_max = resultado_entc_M_X_M_api[0]

#print('\n * max: ',entc_ang_api_max)

entc_ang_api_min = resultado_entc_M_X_M_api[1]

#print('\n * min: ',entc_ang_api_min)


# In[38]:


#salvando resultados

#n_camada = 5

#Inter-face
#vetor de dados de formato: [no_camadas; dist_max_equa; dist_min_equa; dist_max_api; dist_min_api]
result_interf = np.array([ n_camada, intf_dist_eq_min, intf_dist_eq_max, intf_dist_api_min, intf_dist_api_max ])
##print(result_interf)

#escrevendo arquivo de saida
data = np.squeeze(result_interf)
#data = list(data)
#print(data)

#with open("output_n1.dat", "w") as dat_file:
# Criando arquivo xyz com os resultado
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
#with open("output_XM_n4_interf.dat", "w") as f:
#with open("GeI_XM_n4_interf.dat", "w") as f:
#with open("GeBr_XB_n4_interf.dat", "w") as f:
with open("XM_n5_interf.dat", "w") as f:
#with open("Sn_XM_n4_interf.dat", "w") as f:
    writer = csv.writer(f, dialect='myDialect')
    writer.writerow(data)
    
#entre-camada
#vetor de dados de formato: [no_camadas; dist_max_equa; dist_min_equa; dist_max_api; dist_min_api]
result_interf = np.array([ n_camada, entc_dist_eq_min, entc_dist_eq_max, entc_dist_api_min, entc_dist_api_max ])
##print(result_interf)

#escrevendo arquivo de saida
data = np.squeeze(result_interf)
#data = list(data)
#print(data)

#with open("output_n1.dat", "w") as dat_file:
# Criando arquivo xyz com os resultado
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
#with open("output_XM_n4_entrec.dat", "w") as f:
#with open("Ge_XM_n4_entrec.dat", "w") as f:
#with open("GeI_XM_n4_entrec.dat", "w") as f:
#with open("GeBr_XB_n4_entrec.dat", "w") as f:
with open("XM_n5_entrec.dat", "w") as f:
#with open("Sn_XM_n4_entrec.dat", "w") as f:
    writer = csv.writer(f, dialect='myDialect')
    writer.writerow(data)


# In[39]:


#salvando resultados M-X-M

#n_camada = 5

#Inter-face
#vetor de dados de formato: [no_camadas; dist_max_equa; dist_min_equa; dist_max_api; dist_min_api]
result_interf = np.array([ n_camada, intf_ang_eq_min, intf_ang_eq_max , intf_ang_api_min, intf_ang_api_max ])
##print(result_interf)

#escrevendo arquivo de saida
data = np.squeeze(result_interf)
#data = list(data)
#print(data)

#with open("output_n1.dat", "w") as dat_file:
# Criando arquivo xyz com os resultado
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
#with open("output_MXM_n4_interf.dat", "w") as f:
#with open("GeI_MXM_n4_interf.dat", "w") as f:
#with open("GeBr_BXB_n4_interf.dat", "w") as f:
with open("MXM_n5_interf.dat", "w") as f:
#with open("Sn_MXM_n4_interf.dat", "w") as f:
    writer = csv.writer(f, dialect='myDialect')
    writer.writerow(data)
    
#entre-camada
#vetor de dados de formato: [no_camadas; dist_max_equa; dist_min_equa; dist_max_api; dist_min_api]
result_interf = np.array([ n_camada, entc_ang_eq_min, entc_ang_eq_max , entc_ang_api_min, entc_ang_api_max ])
##print(result_interf)

#escrevendo arquivo de saida
data = np.squeeze(result_interf)
#data = list(data)
#print(data)

#with open("output_n1.dat", "w") as dat_file:
# Criando arquivo xyz com os resultado
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
#with open("output_MXM_n4_entrec.dat", "w") as f:
#with open("Ge_MXM_n4_entrec.dat", "w") as f:
#with open("GeI_MXM_n4_entrec.dat", "w") as f:
#with open("GeBr_BXB_n4_entrec.dat", "w") as f:
with open("MXM_n5_entrec.dat", "w") as f:
#with open("Sn_MXM_n4_entrec.dat", "w") as f:
    writer = csv.writer(f, dialect='myDialect')
    writer.writerow(data)


# In[40]:


# Cálculo da média de distâncias d

somas_por_M=np.zeros((no_M,4))

# calculando o valor médio

for ocataedro_M in range(no_M): # separa distâncias r por átomo de metal M
    
    soma_dists0=0
    soma_dists1=0
    somas_por_M[ocataedro_M ][0]=ocataedro_M
    
    # Na formula são considerados 6 atomos de X no octaedro, 
    #porém nas celula unitária podemos encontrar menos de 6 atomos. Logo é necessário contar os X.
    Xs0=0 
    Xs1=0
    
    ##print(ocataedro_M)
    
    # somando as distâncias equatoriais por M
    
    for dists0 in range(len(dist_X_M_equa)): # varre linha por linha da matriz [octaedro, d_{X-M}]
        ##print(ocataedro_M,dists0)
        ##print(soma_dists0)
        if ocataedro_M==dist_X_M_equa[dists0][0]: # identificando a distancia pelo numero do octaedro

            soma_dists0=soma_dists0+dist_X_M_equa[dists0][1] # somatoria das distâncias equatoriais de cada octaedro
            ##print("Metal, soma de distâncias equa:",ocataedro_M,soma_dists0)
            ##print(" inicio",Xs0)
            #if ocataedro_M==0:
            #    #print("loop",dists0,";",ocataedro_M,"# octaedro; somatoria =",soma_dists0,"(+",dist_X_M_equa[dists0][1],")")
            Xs0=Xs0+1
            ##print(" fim",Xs0)
            
    # somando as distâncias apicais por M
    
    for dists1 in range(len(dist_X_M_api)): # vare linha por linha da matriz (r,M)
        
        if ocataedro_M==dist_X_M_api[dists1][0]:

            #soma_dists=soma_dists1+dist_X_M_api[dists1][1]
            soma_dists1=soma_dists1+dist_X_M_api[dists1][1]
            ##print("Metal, soma de distâncias api:",ocataedro_M,soma_dists1)
            ##print(" inicio",Xs1)
            #if ocataedro_M==0:
            #    #print("loop",dists1,";",ocataedro_M,"# octaedro; somatoria =",soma_dists1,"(+",dist_X_M_api[dists1][1],")")
            Xs1=Xs1+1
            ##print(" fim",Xs1)
    
    #if ocataedro_M==0:
    #        #print("\n soma_equa(",soma_dists0,")+soma_api(",soma_dists1,")=",soma_dists0+soma_dists1)
    #        ##print(Xs0,Xs1)
    
    ##print(Xs0,Xs1,Xs0+Xs1)
    if (Xs0+Xs1 != 0):
        somas_por_M[ocataedro_M][1]=(soma_dists0+soma_dists1)/(Xs0+Xs1)
    somas_por_M[ocataedro_M][2]=Xs0+Xs1
    
    #if ocataedro_M==0:
    #        #print("\n soma_equa(",soma_dists0,")+soma_api(",soma_dists1,")=",soma_dists0+soma_dists1)
            
    ##print(ocataedro_M,Xs0+Xs1)
    
##print(len(M_camada_index[0]))
##print(no_M)
# identificando em qual camada o octaedro se encontra.
for octa in range(no_M):
        
    for posi in range(len(M_camada_index[0])):
            
        ##print(posi)
        ##print(octa,posi)
        ##print(ixi)
        
        ##print(M_camada_index[0][octa],somas_por_M[posi][0])
        
        # identificando se está na camada 1
        if M_camada_index[4][posi] == somas_por_M[octa][0]:
            somas_por_M[octa][3]=1
            ##print(M_camada_index[0][posi],"(",posi,")",somas_por_M[octa][0],"(",octa,")")
        
        # identificando se está na camada 2
        if M_camada_index[3][posi] == somas_por_M[octa][0]:
            somas_por_M[octa][3]=2
            ##print(M_camada_index[0][octa],"(",octa,")",somas_por_M[posi][0],"(",posi,")")
            
        # identificando se está na camada 3
        if M_camada_index[2][posi] == somas_por_M[octa][0]:
            somas_por_M[octa][3]=3
            ##print(M_camada_index[0][posi],"(",posi,")",somas_por_M[octa][0],"(",octa,")")
        
        # identificando se está na camada 4
        if M_camada_index[1][posi] == somas_por_M[octa][0]:
            somas_por_M[octa][3]=4
            ##print(M_camada_index[0][octa],"(",octa,")",somas_por_M[posi][0],"(",posi,")")
            
        # identificando se está na camada 5
        if M_camada_index[0][posi] == somas_por_M[octa][0]:
            somas_por_M[octa][3]=5
            ##print(M_camada_index[0][octa],"(",octa,")",somas_por_M[posi][0],"(",posi,")")
            
#print("\n",somas_por_M)


# In[41]:


# calculando o desvio de distâncias \delta d

desvio_d_por_M=np.zeros((no_M,3))

for ocataedro_M in range(no_M): # separa distâncias r por átomo de metal M
    
    soma_dists0=0
    soma_dists1=0
    desvio_d_por_M[ocataedro_M ][0]=ocataedro_M 
    # Na forma são considerados 6 atomos de X no octaedro, 
    #porém nas celula unitária podemos encontrar menos de 6 atomos. Logo é necessário contar os X.
    Xs0=0 
    Xs1=0
    
    # somando as distâncias equatoriais por M
    
    for dists0 in range(len(dist_X_M_equa)): # vare linha por linha da matriz (r,M)
        
        if ocataedro_M==dist_X_M_equa[dists0][0]:

            soma_dists0=soma_dists0+((dist_X_M_equa[dists0][1]-somas_por_M[ocataedro_M][1])/somas_por_M[ocataedro_M][1])**2
            ##print("Metal, soma de distâncias equa:",ocataedro_M,soma_dists0)
            
            Xs0=Xs0+1
        
    # somando as distâncias apicais por M
    
    for dists1 in range(len(dist_X_M_api)): # vare linha por linha da matriz (r,M)
        
        if ocataedro_M==dist_X_M_api[dists1][0]:

            soma_dists1=soma_dists1+((dist_X_M_api[dists1][1]-somas_por_M[ocataedro_M][1])/somas_por_M[ocataedro_M][1])**2
            ##print("Metal, soma de distâncias api:",ocataedro_M,soma_dists1)
            #if ocataedro_M==0:
            #    #print("loop",dists1,";",ocataedro_M,"# octaedro; somatoria =",soma_dists1)
            #    #print("*somatoria =",soma_dists1,"=(",dist_X_M_api[dists1][1],"-",somas_por_M[ocataedro_M][1],")/",somas_por_M[ocataedro_M][1],"^2")
            Xs1=Xs1+1
    if (Xs0+Xs1 != 0):
        desvio_d_por_M[ocataedro_M][1]=(soma_dists0+soma_dists1)/(Xs0+Xs1)
    desvio_d_por_M[ocataedro_M][2]=somas_por_M[ocataedro_M][3]
    ##print("soma equa + soma api:",soma_dists0," + ",soma_dists1)
    
    #if ocataedro_M==0:
    #    #print("\n soma_equa(",soma_dists0,")+soma_api(",soma_dists1,")=",soma_dists0+soma_dists1)
    #    #print("\n Delta_d=",soma_dists0+soma_dists1,"/",Xs0+Xs1,"=",desvio_d_por_M[ocataedro_M][1])
#print("\n",desvio_d_por_M)


# In[42]:


# Separando apenas os octaedros completos (com 6 haletos)

completos=0

# contando a quantidade de octaedros completos
for octa in range(len(somas_por_M)):
    
    if (somas_por_M[octa][2] == 6.0):
        
        ##print(somas_por_M[octa][2],completos)
        
        completos=completos+1
        
# criando nova matriz de resultados
desvio_d_por_M_completos=np.zeros((completos,3))

##print(desvio_d_por_M_completos)

# reagrupando os resultados
cont_octa=0
for octa in range(len(desvio_d_por_M)):
    
    if (somas_por_M[octa][2] == 6.0):
    
        ##print(octa)

        desvio_d_por_M_completos[cont_octa][0]=somas_por_M[octa][0]
        desvio_d_por_M_completos[cont_octa][1]=desvio_d_por_M[octa][1]
        desvio_d_por_M_completos[cont_octa][2]=desvio_d_por_M[octa][2]
    
        cont_octa=cont_octa+1
    
#print(desvio_d_por_M_completos)


# In[43]:


# separando coordenadas por octaedro

desvio_theta_por_M=np.zeros((no_M,3))

for ocataedro_M in range(no_M): # buscar pelo índice de octaedros
    
    # Criando uma matriz auxiliar para armazenar cada octaedro
    matrix_aux=np.zeros((int(somas_por_M[ocataedro_M][2]),5)) # [r(x,y,z),distâmcia,número do octaedro]
    ##print(ocataedro_M,int(somas_por_M[ocataedro_M][2]),": \n",matrix_aux,"\n")
    i_aux=0
    conta=0
    
    for i in range(len(intervalos_X_M)): #loop pela matriz de coordenadas de X equatoriais

        if (ocataedro_M == intervalos_X_M[i][3]):
        ##print("Octaedro ",ocataedro_M,"; linha:",i)
        
            # preenchendo a matriz com coordenadas de X
            matrix_aux[i_aux][0]=intervalos_X_M[i][0]
            matrix_aux[i_aux][1]=intervalos_X_M[i][1]
            matrix_aux[i_aux][2]=intervalos_X_M[i][2]
            
            # preenchendo a distância equatorial
            matrix_aux[i_aux][3]=dist_X_M_equa[i][1]
            
            # preenchendo a quarta coluna com a identificação de cada X equatorial
            matrix_aux[i_aux][4]=i_aux
            
            ##print(i,matrix_aux,"\n")

            i_aux=i_aux+1

    for i in range(len(interv_api_X_M)): #loop pela matriz de coordenadas de X apicais

        if (ocataedro_M == interv_api_X_M[i][3]):
        ##print("Octaedro ",ocataedro_M,"; linha:",i)
            
            # preenchendo a matriz com coordenadas de X apical
            matrix_aux[i_aux][0]=interv_api_X_M[i][0]
            matrix_aux[i_aux][1]=interv_api_X_M[i][1]
            matrix_aux[i_aux][2]=interv_api_X_M[i][2]
            
            # preenchendo a distância apical
            matrix_aux[i_aux][3]=dist_X_M_api[i][1]
            
            # preenchendo a quarta coluna com a identificação de cada X
            matrix_aux[i_aux][4]=i_aux
            
            ##print(i,matrix_aux,"\n")

            i_aux=i_aux+1
            
    # Calculando os ângulos
    
    indice_thetas=0
    
    # Matriz auxiliar onde os ângulos são temporariamente conservados
    theta_temporario=np.zeros((int(somas_por_M[ocataedro_M][2]),int(somas_por_M[ocataedro_M][2])))
    
    # copiando as coordenadas dos X_i e X_h e calculando suas distancias
    # e preenchendo a matrix auxiliar
    for i_i in range(len(matrix_aux)): # vare linha por linha da matriz (x y z,X)
        #dados para X_i
        r_i = matrix_aux[i_i][0:3]
        d_i = matrix_aux[i_i][3]
        soma_theta=0
        i_theta=0
        
        
        for i_h in range(len(matrix_aux)):
            
            ##print("#",i_h," é diferente de ",i_i,"? \n")
            ##print("#",i_i," |<--->| ",i_h,": ")

            if i_i != i_h: # os átomos precisam ser diferentes para que o cálculo seja válido
                
               # #print("     *** i_i =",i_i," e i_h =",i_h," ***\n")
            
    
                #dados para X_h
                r_h = matrix_aux[i_h][0:3]
                d_h = matrix_aux[i_h][3]
                
                ##print("  R: sim! \n")
                ##print("  Então, \n   d_i =",d_i,"\n   d_h =",d_h,"\n")
                
                # calculando a distancia
                r_delta=r_i-r_h
                
                ##print("    d_i-d_h =",d_delta,"\n")
                d_ih = math.sqrt( r_delta[0]**2 + r_delta[1]**2 + r_delta[2]**2 )
                ##print("    d_ih =",d_ih,"\n")
                
                # calculando o angulo
                theta_cart= (d_i**2 + d_h**2 - d_ih**2)/(2*d_i*d_h)
                
                # O arcosseno compreende um intervalo de theta pertencente à [-1,1].
                # Logo, todos os valores fora desse intervalo serão descartados
                if (theta_cart <= 1):
                    
                       if (theta_cart >= -0.9): #(theta_cart >= -1): pois theta menor que -0.9 é maior que 150
                
                            theta = np.arccos(theta_cart)
                            theta = math.degrees(theta)

                            theta_temporario[i_i][i_h]=theta
                            
            #if (ocataedro_M == 0):
            #if (i_i == len(matrix_aux)-1):
            #        if (i_h == len(matrix_aux)-1):
            #            #print("\n",theta_temporario,"\n")
            #            ##print("     *** i_i =",i_i," e i_h =",i_h," ***\n")
                        
    # calculando a soma dos ângulos da matriz triangular
    sutien=0
    n_thetas=0
    
    ##print("\n",theta_temporario,"\n")
    
    for line in range(len(theta_temporario)):
        
        for column in range(conta,len(theta_temporario)):
            
            if (theta_temporario[line][column] != 0):
                
                sutien=sutien+((theta_temporario[line][column]-90)**2)
                n_thetas=n_thetas+1
                ##print(theta_temporario[line][column])
            ##print(sutien)
            ##print(line,column,conta)
            
        conta=conta+1
    
    ##print("Octaedro ",ocataedro_M)
        
    if n_thetas-1 == 0:
        desvio_angulo=sutien
    else:
        desvio_angulo=sutien/(n_thetas-1)
    ##print("\n",desvio_angulo,n_thetas-1,"\n")
    
    #
    desvio_theta_por_M[ocataedro_M][0]=ocataedro_M
    desvio_theta_por_M[ocataedro_M][1]=desvio_angulo


# In[44]:


# Separando apenas os octaedros completos (com 6 haletos)

# criando nova matriz de resultados
desvio_theta_por_M_completos=np.zeros((completos,3))

# reagrupando os resultados
cont_octa=0
for octa in range(len(desvio_d_por_M)):
    
    if (somas_por_M[octa][2] == 6.0):
    
        ##print(octa)

        desvio_theta_por_M_completos[cont_octa][0]=somas_por_M[octa][0] # octaedro
        desvio_theta_por_M_completos[cont_octa][1]=desvio_theta_por_M[octa][1] # theta
        #desvio_theta_por_M_completos[cont_octa][2]=desvio_d_por_M_completos[octa][2] # camada
        ##print(cont_octa,octa)
        ##print(desvio_d_por_M_completos)    
        cont_octa=cont_octa+1
        
for octa in range(len(desvio_theta_por_M_completos)):
    ##print(octa)
    desvio_theta_por_M_completos[octa][2]=desvio_d_por_M_completos[octa][2]

#print(desvio_theta_por_M_completos)


# In[45]:


#salvando resultados "bond length distortion - delta d (boleto)"

#n_camada = 4

#Inter-face e Entre-camada
n_Atomos = 1 #numero de atomos para o primeira linha
Cabec = [[n_Atomos],[]]
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
with open("interf_entrec_delta_d_n_5.dat", "w") as f:
    
    writer = csv.writer(f, dialect='myDialect',lineterminator = '\n')
            
    # escreve linhas no formato [ n_octaedro , delta d ]
    for row in desvio_d_por_M_completos:
        row = row.tolist()
        writer.writerow(row)


# In[46]:


#salvando resultados "bond angle variance - theta^2 (banva)"

#n_camada = 4

#Inter-face e Entre-camada
n_Atomos = 1 #numero de atomos para o primeira linha
Cabec = [[n_Atomos],[]]
csv.register_dialect('myDialect', delimiter = ' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
with open("interf_entrec_theta_n_5.dat", "w") as f:
    writer = csv.writer(f, dialect='myDialect')
    
    # escreve linhas no formato [ n_octaedro , theta^2 ]
    for row in desvio_theta_por_M_completos:
        row = row.tolist()
        writer.writerow(row)


# In[ ]:



