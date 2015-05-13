# -*- coding: utf-8 -*-

from itertools import product as cartesianProduct
from scipy.misc import comb
import numpy as np
from math import floor, sqrt

""" When set to true, more output is printed at run-time in order to facilitete
 debugging. """
verbose = True

#The maximum calue of n that we want to
Nmax = 10L

# What numeric value is used to represent that a value has not yet been
# computed
NA = float("nan")
def NAcheck(value):
    #returns True if "value" corresponds to NA
    return str(value) == str(NA)

#==============================================================================
# initialize a tables to store all values of T that have been computed so far.
#
#  t_*[i][j][k] stores t(n = i,l = j, L = k).
#  A ll values are initialized as nan indicates that these values have not been
#  computed
#==============================================================================
#t_table = [[[ nan for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
#table_t_p = np.array([[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)],ndmin=3)
table_t_p = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_np = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_s_p = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_s_np = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]

def validateInput(n,l,L):
    if type(n) != long or type(l) != long or type(L) != long:
        raise  Exception("Input passed to t(n,l,L) must have type 'long'.")
    if min(n,l,L+1) < 1:
        raise Exception("Inputs passed to t(n,l,L) too small.")

def divisors(n):
    factors = set()
    for x in range(1, int(sqrt(n)) + 1):
        if n % x == 0:
            factors.add(long(x))
            factors.add(long(n//x))
    return sorted(factors)

def binom(n,k):
    return long(comb(n,k,exact = True))

def lrange(start,stop):
    return map(long,range(start,stop))



def t_p(n,l,L):

    validateInput(n,l,L)

    if n>2 and n - l >= L and l>1:
        if verbose and NAcheck(table_t_s_p[n-1][l][L]):
            print "table_t_s_p[n-1][l][L] = NA at n,l,L = %i,%i,%i"%(n,l,L)
        if verbose and str(table_t_s_np[n-1][l-1][L]) == str(NA):
            print "table_t_s_np[n-1][l-1][L] = NA at n,l,L = %i,%i,%i"%(n,l,L)
        return table_t_s_p[n-1][l][L] + table_t_s_np[n-1][l-1][L]

    elif n==l and (n==1 or n==2) and L==0:
        return 1L

    else:
        return 0L

def t_s_p(n,l,L):

    validateInput(n,l,L)

    if n >= 2 and n - l + 1 >= L and L > 0:

        result = L*(table_t_p[n][l][L] + table_t_p[n][l][L-1])

        if NAcheck(result):
            errorString = "t_s_p(%i,%i,%i) evaluated to NA"%(n,l,L)
            raise Exception(errorString)

        return L*(table_t_p[n][l][L] + table_t_p[n][l][L-1])

    elif n==l and n==1 and L==0:
        return 1L

    else:
        return 0L

def t_s_np(n,l,L):

    validateInput(n,l,L)

    if l < 2:
        return 0L

    if n==(l+1) and n>2 and L==0:
        return 1L

    elif n>2 and n - l >= L and L>0:

        result = L*(table_t_np[n][l][L] + table_t_np[n][l][L-1])

        if NAcheck(result):
            errorString = "t_s_np(%i,%i,%i) evaluated to NA"%(n,l,L)
            raise Exception(errorString)

        return L*(table_t_np[n][l][L] + table_t_np[n][l][L-1])

    else:
        return 0L


def t_np(n,l,L):

    validateInput(n,l,L)

    if n < 3 or l < 2:
        return 0L

    elif n==l+1 and L==0:
        return 1L

    elif L+l>n-1:
        return 0L

    else:
        temp = 0L

        for m in range(1,n):
            for d in divisors(m):
                if d != (n-1) :

                    j = m//d
#                    if verbose and j == 0:
#                        print "(j,m,d)=(%i,%i,%i)"%(j,m,d)

                    ## case: planted, planted
                    for l2 in range(1, l//j + 1 +1):
                        l1 = l + 1 - j * (l2 - 1)
                        for L1 in range(0, min(L, n-m-l1) +1 ):
                            for L2 in range(L-L1,min(L, d-l2 +1)+1):
#                                temp += 1
                                temp += d * table_t_p[n-m][l1][L1] * table_t_s_p[d][l2][L2]
#                                temp += d * binom(L,L2) * table_t_p[n-m][l1][L1] * table_t_s_p[d][l2][L2]

                    ## case: planted, non-planted
                    for l2 in range(1, l//j+1):
                        l1 = l + 1 - j * l2
                        for L1 in range(0,min(L,n-m-l1)+1):
                            for L2 in range(L-L1,min(L,d-l2)+1):
#                                temp += 1
                                temp += d * table_t_p[n-m][l1][L1] * table_t_s_np[d][l2][L2]
#                                temp += d * binom(L,L2) * table_t_p[n-m][l1][L1] * table_t_s_np[d][l2][L2]

                    ## case: non-planted, planted
                    for l2 in range(1, (l-1)//j + 1 +1):
                        l1 = l - j * (l2 - 1)
                        for L1 in range(0,min(L,n-m-l1 -1)+1):
                            for L2 in range(L-L1,min(L,d-l2 +1)+1):
#                                temp += 1
                                temp += d * table_t_np[n-m][l1][L1] * table_t_s_p[d][l2][L2]
#                                temp += d * binom(L,L2) * table_t_np[n-m][l1][L1] * table_t_s_p[d][l2][L2]

                    ## case: non-planted, non-planted
                    for l2 in range(1, (l-1)//j +1):
                        l1 = l - j * l2
                        for L1 in range(0,min(L,n-m-l1-1)+1):
                            for L2 in range(L-L1,min(L, d-l2 )+1):
#                                temp += 1
                                temp += d * table_t_np[n-m][l1][L1] * table_t_s_np[d][l2][L2]
#                                temp += d * binom(L,L2) * table_t_np[n-m][l1][L1] * table_t_s_np[d][l2][L2]

        if verbose:
            control = "t_np(%i,%i,%i) = %f"%(n,l,L,float(temp)/(n-1))
            print control

        temp = temp//(n-1)
        return temp


#==============================================================================
#  BEGIN: fill out tables.
#==============================================================================

def generateTable():
    for n in lrange(1,Nmax+1):
        for L in lrange(0,Nmax+1):
#            table_t_p[n][0][L] = 0L
#            table_t_np[n][0][L] = 0L
            for l in lrange(1,Nmax+1):
#                print (n,l,L)
#                print (table_t_p[n][l][L],table_t_np[n][l][L])
                table_t_p[n][l][L] = t_p(n,l,L)
                table_t_np[n][l][L] = t_np(n,l,L)
#                print (table_t_p[n][l][L],table_t_np[n][l][L])
#                print ""
        for L in lrange(0,Nmax+1):
            for l in lrange(1,Nmax+1):
                table_t_s_p[n][l][L] = t_s_p(n,l,L)
                table_t_s_np[n][l][L] = t_s_np(n,l,L)

generateTable()
#==============================================================================
# BEGIN: setting initial values
#==============================================================================
#t_p[1][1][1] = 1L
##t_np[1][1][1] = 0L
#t_s_p[1][1][1] = 1L
##t_s_np[1][1][1] = 0L
#
### L = 0 implies t_*(n,l,L) = 0
#for n,l in cartesianProduct(range(1,Nmax+1),range(1,Nmax+1)):
#    t_p[n][l][0] = 0L
#    t_np[n][l][0] = 0L
#    t_s_p[n][l][0] = 0L
#    t_s_np[n][l][0] = 0L
#
### There are no non-planted trees with less than 3 nodes.
#for n,l,L in cartesianProduct(range(1,3),range(1,Nmax+1),range(1,Nmax+1)):
#    t_np[n][l][L] = 0L
#    t_s_np[n][l][L] = 0L
#
#for L in range(1,Nmax +1):
#    pass
#
#for n in range(1,Nmax+1):
#    t_p[n][1][1] = 1L

#==============================================================================
# END: setting initial values
#==============================================================================

#==============================================================================
# def t_p_recursion(n,l,L):
#
#     # Validate input
#     if type(n) != long or type(l) != long or type(L) != long:
#         raise  Exception("Input passed to t(n,l,L) must have type 'long'.")
#     if min(n,l,L+1) < 1:
#         raise Exception("Inputs passed to t(n,l,L) too small.")
#
#     #Handle base cases
#     if max(l,L) > n or n - l < L:
#         return 0L
#     if n==1:
#         if (l==1 and L==1):
#             return 1L
#         else:
#             return 0L
#
#==============================================================================