# -*- coding: utf-8 -*-
import matplotlib.pyplot as gph
 
import math

def solution(x):
    return (math.exp(x) + math.exp(-x) + 3.3*x*x - 3.3*x - 2)

def function(x,y):
    return y + 8.6 + 3.3 * x - 3.3*x*x

def Euler(N, M):
    #N = 20
    #M = 1# должно вводиться
    h = 1./N
    step = h
    #print step
    X = [0]
    Y = [M]
    Z = []
    #начальные условия, решение по методу Эйлера
    t_x = 0
    t_y = M
    t_z = t_y - 3.3
    Z.append(t_z)
    #print Z[0]
    #print(t_x)
    #print(t_y)
    #print(t_z)
    i = 1
    while(i <= N): # сделать цикл с условием на эпсилон
        t_z = t_z + h*function(t_x, t_y)
        Z.append(t_z)
        t_y = t_y + h * Z[i-1]
        Y.append(t_y)
        t_x += step
        X.append(t_x)
        i+=1
    #print "t_z"
    #print t_z
    #print "F(mu)"
    #print function(t_x,t_z)
    #gph.plot(X,Y)
    #print X
    #print Y

    return t_z
    #print function(5)
# конец метода Эйлера!

def vnesh_cycle(N): # для Эйлера
    eps = 0.000001
    M=[]
    M.append(0)
    M.append(1)
    i=1
    while(abs(M[i]-M[i-1]) > eps):     
        z1 = Euler(N, M[i-1])
        z2 = Euler(N, M[i])
        o1 = z1 - math.exp(1) + math.exp(-1) - 3.3
        o2 = z2 - math.exp(1) + math.exp(-1) - 3.3
        temp =  M[i] - (o2*(M[i]-M[i-1])/(o2-o1))
        #print temp
        M.append(temp)
        i+=1
    
    
    h = 1./N
    step = h
    t=M[i-1]
    #t = 0.1
    #print t
    #print step
    X = [0]
    Y = [0]
    Z = []
    #начальные условия, решение по методу Эйлера
    t_x = 0
    t_y = t
    t_z = t_y - 3.3
    Z.append(t_z)
    i = 1
    while(i <= N): # сделать цикл с условием на эпсилон
        t_z = t_z + h*function(t_x, t_y)
        Z.append(t_z)
        t_y = t_y + h * t_z
        Y.append(t_y)
        t_x += step
        X.append(t_x)
        i+=1        
    #print t_x
    #print t_y
    #print t_z
    print "Z:"
    #print Z
    print "Function equal: \t"
#Требуется вычислить значение в Методе Эйлера: f(Mu) = z(1,Mu) - e + 1/e - 3.3
    print function(t_x,1.30453) - math.exp(1) + math.exp(-1) - 3.3 #Если ищу значение в правильном месте и правильные аргументы подставляю в функцию? 
    print "при значении x = 1.0, mu = 1.30453" 
    gph.plot(X,Y)
    gph.show() 

def sol(x,N):
    h = 1./N
    step = h
    
    X=[0]
    Y=[]
    
    t = solution(0)
    #print t
    
    Y.append(t)
        
    i=1
    while(i <= N):
        x = x+step
        X.append(x)
        y = solution(x)
        Y.append(y)
        i += 1

    gph.plot(X,Y)
    #print X
    #print Y
    gph.show()

def progonka(N):
    h = 1./N
    a1= 0.
    A = [a1]
    k = 1
    while(k < N):
        A.append(1.)
        k += 1

    A.append(-1.)
    
    # test output
    #print A

    b1 = -(1 + h)
    B = [b1]
    k = 1
    while(k < N):
        temp = -(2 + h*h)
        B.append(temp)
        k +=1

    B.append(1.)
    
    # test output
    #print B
    c1 = 1.
    C = [c1]
    k = 1
    while(k < N):
        C.append(1.)
        k +=1

    C.append(0.)
    
    # test output 
    #print "Array C:"
    #print C
    d1 = -3.3*h
    D = [d1]
    k=1
    t_x = 0.
    while(k < N):
        t_x = t_x +h
        t_d = (8.6+3.3*t_x - 3.3*t_x*t_x)*h*h
        D.append(t_d)
        k+=1
        
    # test output
    #print D
    D.append((math.e - 1. / math.e + 3.3) * h)
    #print Lyambda[1] #
    
    # test output
    # The old factor of acceleration
    #  A = [0.]
    #  B = [-(1 + h)]
    #  C = [1.]
    #  D = [-3.3 * h]
    #
    #  x = 0.
    #  for i in range(N - 1):
    #      x += h
    #      A.append(1.)
    #      B.append(-(2 + h ** 2))
    #      C.append(1.)
    #      D.append((8.6+3.3*x - 3.3*x*x)*h*h)
    #
    #  A.append(-1.)
    #  B.append(1.)
    #  C.append(0.)
    #  D.append((math.e - 1. / math.e + 3.3) * h)

    Lyambda =[0.]
    Mu = [0.]
    
    # test output
    #lyambda = -C[0]/B[0]
    #mu=D[0]/B[0]
    #print "l"
    #print lyambda
    #print mu
    i=1
    while(i <= N+1):
        t_Lyamb = -C[i - 1] / (A[i - 1] * Lyambda[i - 1] + B[i - 1])
        t_Mu = (D[i - 1] - A[i - 1] * Mu[i - 1]) / (A[i - 1] * Lyambda[i - 1]+B[i - 1])
        Lyambda.append(t_Lyamb)
        Mu.append(t_Mu)
        i+=1
        
    # test output
    #print "Mu"
    #print Mu
    Y = [0]
    
    # test output
    #print Lyambda
    i = 1
    while(i <= N):
        Y.append(0)
        i+=1
# create array
    i = N
    Y[N]=Mu[N+1]
    while(i >=1):
        t_y = Lyambda[i] * Y[i] + Mu[i]
        Y[i-1] = t_y
        i-=1
    #print Y
    i = 1
    X=[0]
    t_x = 0
    while(i <= N):
        t_x += h
        X.append(t_x)
        i +=1
    ys = [solution(x) for x in X]
    gph.plot(X,Y, label="approx")
    gph.plot(X, ys, label="precise")
    gph.legend()
    gph.show()
def Main(N):

    # test output:
    #  sol(0, N)
    progonka(N)
    vnesh_cycle(N)
# To check:
#sol(0,40)
#Euler(10,0)
#progonka(40)
#vnesh_cycle(40)
#print solution(0)#math.exp(0)
Main(50)
