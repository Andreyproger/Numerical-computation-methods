#3_8
import math

print("input a")
a = input()

print("input b")
b = input()

print("input h")
h = input()

n = (b-a)/h

sum=0
temp = 0
i=1

while(i<n+1):
    sum += (math.sqrt(1+ pow(1 + (i-1)*h, 4))
    + 3*math.sqrt(1+ pow(1 + (i-1)*h + h/3, 4))
    + 3*math.sqrt(1+ pow(1 + i*h - h/3, 4)) +
    math.sqrt(1+ pow(1 + i*h, 4)))
    i+=1

def trapez_cycle(h):
	return h * sum/8

print (trapez_cycle(h))

#Gauss

print "input a (the left edge):"
a = math.ceil(input())
print a

print "input a (the right edge):"
b = math.ceil(input())
print b

 
print (math.sqrt( 1 + pow(1.5- math.sqrt(3)/6, 4)) + math.sqrt(1+ pow(1.5 + math.sqrt(3)/6,4)))/2


#trapezoid

print "Input step:"
h = input()
print h


print "input a (the left edge):"
a = math.ceil(input())
print a

print "input a (the right edge):"
b = math.ceil(input())
print b

m = (b-a)/h

"""
        sum = 0
        for i in range(1, m-1):
                sum += formula(1 + h*i)"""

def formula(x):
	return math.sqrt(1 + pow(x,4));

def sum(m,h):
        sum = 0
        i = 1
        while(i < m):
                sum += formula(1 + h*i)
                i += 1
        return sum

def trapez_cycle(a,b,m,h):
	return h * ((formula(a) + formula(b))/2 + sum(m,h))

print trapez_cycle(a,b,m,h)
