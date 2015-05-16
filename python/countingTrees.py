# -*- coding: utf-8 -*-

from scipy.misc import comb
import numpy as np
from math import sqrt


#==============================================================================
# BEGIN VARIABLES
#==============================================================================

# Print extra output when running
verbose = True # Print to terminal what is currently being compute
veryVerbose = False # print every non-trivial result of applying t_np(n,l,L).

#The maximum calue of n that we want to
Nmax = 50L

# What numeric value is used to represent that a value has not yet been
# computed.
#NA = float("nan")
NA = -1L

# initialize tables to store all values of recursions that have been computed
# so far.
#
#  t_*[i][j][k] stores t(n = i,l = j, L = k).
#
#  All values are initialized as NA indicating that these values are yet to be
#  computed
table_t_p = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_np = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_s_p = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]
table_t_s_np = [[[ NA for L in range(Nmax+1)] for l in range(Nmax+1)] for n in range(Nmax+1)]


#==============================================================================
# END VARIABLES
#
# BEGIN FUNCTIONS FOR PRINTING
#==============================================================================

#print all four tables to output for a fixed value of n
def printTables(n):
    print "print np.matrix(table_t_p[%i])"%n
    print np.matrix(table_t_p[n] , dtype=long )

    print "\nprint np.matrix(table_t_np[%i])"%n
    print np.matrix(table_t_np[n] , dtype=long )

    print "\nprint np.matrix(table_t_s_p[%i])"%n
    print np.matrix(table_t_s_p[n] , dtype=long )

    print "\nprint np.matrix(table_t_s_np[%i])"%n
    print np.matrix(table_t_s_np[n] , dtype=long )

# print a specific talbe (fixed n), formated as a latex-table.
# Supports printing to files.
def printTableLatex(n,planted=False,star=False,printToFile=False,path="./tables/LaTeX/"):
    if planted:
        if star:
            table = table_t_s_p
            rowLabel = lambda l : r"$t^{p \, \star}_"+"{(%i,%i,\\Lp)}$"%(n,l)
        else:
            table = table_t_p
            rowLabel = lambda l : r"$t^{p}_"+"{(%i,%i,\\Lp)}$"%(n,l)
    else:
        if star:
            table = table_t_s_np
            rowLabel = lambda l : r"$t^{\neg p \, \star}_"+"{(%i,%i,\\Lp)}$"%(n,l)
        else:
            table = table_t_np
            rowLabel = lambda l : r"$t^{\neg p}_"+"{(%i,%i,\\Lp)}$"%(n,l)

    output = "\\begin{tabular}{%s}\n  "%("r|"+"r"*(n+1),)
    for s in [ " & \\Lp=%i"%L for L in range(0,n+1)]:
        output += s
    output += " \\\\\n  \\hline\n  "
    for l in range(1,n+1):
        output += rowLabel(l)
        for s in [ " & $%i$"%table[n][l][L] for L in range(0,n+1)]:
            output += s
        output += " \\\\\n"
        if l < n : output += "  "
    output += "\\end{tabular}\n"

    if not printToFile:
        return output
    else:
        if planted:
            if star:
                filename = "table_t_s_p__n_is_%i.tex"%n
            else:
                filename = "table_t_p__n_is_%i.tex"%n
        else:
            if star:
                filename = "table_t_s_np__n_is_%i.tex"%n
            else:
                filename = "table_t_np__n_is_%i.tex"%n

        f1 = open(path+filename, "w")
        f1.write(output)
        f1.close()

# Print all tables to .tex files (stored in "./tables/LaTeX/")
def generateTablesLatex():
    for i in range(1,n+1):
        for planted in (True,False):
            for star in (True,False):
                printTableLatex(i,planted,star,printToFile=True,path = "./tables/LaTeX/")

# Print all tables as .csv-files (stored in "./tables/csv/")
def generateTablesCSV(path="./tables/csv/"):
    for i in range(1,Nmax+1):
        np.savetxt(path+"t_p__n_%i"%i, np.array(table_t_p[i]), delimiter=" , ")
        np.savetxt(path+"t_np__n_%i"%i, np.array(table_t_np[i]), delimiter=" , ")
        np.savetxt(path+"t_s_p__n_%i"%i, np.array(table_t_s_p[i]), delimiter=" , ")
        np.savetxt(path+"t_s_np__n_%i"%i, np.array(table_t_s_np[i]), delimiter=" , ")


#==============================================================================
# END PRINTING FUNCTIONS
#
# BEGIN AUXILIARY FUNCTIONS FOR RECURSIONS
#==============================================================================

#returns True if "value" corresponds to NA
def NAcheck(value):
    return str(value) == str(NA)

# Verifies that only admissible values are passed to enumeration-recursion.
def validateInput(n,l,L):
    if type(n) != long or type(l) != long or type(L) != long:
        raise  Exception("Input passed to t(n,l,L) must have type 'long'.")
    if min(n,l,L+1) < 1:
        raise Exception("Inputs passed to t(n,l,L) too small.")

# Return a list of all divisors of n in ascending order, stored as longs
def divisors(n):
    factors = set()
    for x in range(1, int(sqrt(n)) + 1):
        if n % x == 0:
            factors.add(long(x))
            factors.add(long(n//x))
    return sorted(factors)

#compute binomial coefficients (stored as longs)
def binom(n,k):
    return long(comb(n,k,exact = True))

#same as range(start,stop), but cast to long
def lrange(start,stop):
    return map(long,range(start,stop))

#==============================================================================
# END AUXILIARY FUNCTIONS
#
# BEGIN RECURSIONS
#==============================================================================

def t_p(n,l,L):

    validateInput(n,l,L)

    if n>2 and n - l >= L and l>1:
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

    if n < 3 or l < 2 or L+l>n-1:
        return 0L

    else:
        temp = 0L

        for m in range(1,n):
            for d in divisors(m):
                if d != (n-1) :

                    j = m//d

                    ## case: planted, planted
                    for l2 in range(1, l//j + int(d!=1) +1):
                        l1 = l + 1 - j * (l2 - int(d!=1))
                        for L1 in range(0, min(L, n-m-l1) +1 ):
                            for L2 in range(L-L1,min(L, d-l2 +int(d!=1))+1):
                                temp += d * binom(L,L1) * binom(L1,L1+L2-L) * table_t_p[n-m][l1][L1] * table_t_s_p[d][l2][L2]

                    ## case: planted, non-planted
                    for l2 in range(1, l//j +1):
                        l1 = l + 1 - j * l2
                        for L1 in range(0,min(L,n-m-l1)+1):
                            for L2 in range(L-L1,min(L,d-l2)+1):
                                temp += d * binom(L,L1) * binom(L1,L1+L2-L) * table_t_p[n-m][l1][L1] * table_t_s_np[d][l2][L2]

                    ## case: non-planted, planted
                    for l2 in range(1, (l-1)//j + int(d!=1) +1):
                        l1 = l - j * (l2 - int(d!=1))
                        for L1 in range(0,min(L,n-m-l1 -1)+1):
                            for L2 in range(L-L1,min(L,d-l2 +int(d!=1))+1):
                                temp += d * binom(L,L1) * binom(L1,L1+L2-L) * table_t_np[n-m][l1][L1] * table_t_s_p[d][l2][L2]

                    ## case: non-planted, non-planted
                    for l2 in range(1, (l-1)//j +1):
                        l1 = l - j * l2
                        for L1 in range(0,min(L,n-m-l1-1)+1):
                            for L2 in range(L-L1,min(L, d-l2 )+1):
                                temp += d * binom(L,L1) * binom(L1,L1+L2-L) * table_t_np[n-m][l1][L1] * table_t_s_np[d][l2][L2]

        if veryVerbose:
            control = "t_np(%i,%i,%i) = %f"%(n,l,L,float(temp)/(n-1))
            print control

        temp = temp//(n-1)
        return temp

def solveRecursionsIteratively():
    for n in lrange(1,Nmax+1):
        if verbose:
            print "Computing tables for n=%i ... "%n
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

        if verbose:
            print "n = %i"%n
            printTables(n)
            print "="*79+"\n"


#==============================================================================
# END RECURSIONS
#
# BEGIN ACTUALLY COMPUTING STUFF (AND STORING OUTPUT)
#==============================================================================

if verbose:
    print "Solving Recursions"
solveRecursionsIteratively()

if verbose:
    print "Generating CSV-tables"
generateTablesCSV()

if verbose:
    print "Generating LaTeX-tables"
generateTablesLatex()