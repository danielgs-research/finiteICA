# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:11:43 2018

@author: MateusMarcuzzo

Galois fields class and helpers

    News(02/11/2018)
    we can make polynomials with
    from numpy.polynomial import Polynomial as P

    p = P([GF(1),GF(2)])

    we can
    +, -, *, //, %, divmod, **, ==
    the polynomials

    Still, we couldn't evaluate using p(x), because it calls for
        float multiplication
    which does not make sense here.

    but with polyval:

    from numpy.polynomial.polynomial import polyval
    it is possible

    polyval is not evaluating correctly when doing composed polynomials in GF
num
TO-DO:

    (not urgent)Implementation of power of prime GFs.
        Needs knowledge of polynomials, suggestion: look the code books of
        1st semester of Master's graduation. (2017)

        maybe an extended class of this base one.

    (09/03/2019)
    Implementation of row-echelon-form mod P

    Do the try-catch on truediv
    Raise exception when not prime

"""
import numpy as np

if __name__=="__main__":
    from prime_numbers import is_prime
    
else:
    from prime_numbers import is_prime
    
    
    
class GaloisField():
    """
    Data Abstraction for Finite Field elements of GF(P^n)

    OBS:
        It just works for GF(P), not with power of primes (yet)(maybe never)
    """

    def __init__(self, number, prime = 2):
        if(is_prime(prime)):
            self.__prime = prime
        else:
            raise   ValueError("The second parameter must be a prime\n"+
                               "GF(number,prime)\n            ^")

        self.number = number % self.__prime



    def __int__(self):
        return int(self.number)

    def __eq__(self, other):
        if(isinstance(other, GaloisField)):
            return self.number == other.number

        return self.number == other


    def __add__(self, other):
        if(isinstance(other, GaloisField)):
            return GaloisField(self.number + other.number,\
                           self.__prime)
        elif(isinstance(other, int)):
            return GaloisField(self.number + other,\
                               self.__prime)


    def __radd__(self, other):
        return self.__add__(other)

    def __pow__(self, other):
        if (other < 0):
            return self.inverse() ** (-other)

        if(other > 0):
            return self**(other-1) * self
        else:
            return GaloisField(1,self.__prime)

    def __mul__(self, other):
        if(isinstance(other, GaloisField)):
            return GaloisField(self.number * other.number,\
                           self.__prime)
        elif (isinstance(other,int)):
            return GaloisField(self.number * other,
                               self.__prime)
        else:
            return other.__mul__(self)


    #since multiplication is commutative..
    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
        if(isinstance(other, GaloisField)):
            return GaloisField(self.number - other.number,\
                           self.__prime)

        elif(isinstance(other,int)):
            return GaloisField(self.number - other,\
                               self.__prime)
            
        else: 
            return other.__add__(self * -1)

    def __rsub__(self, other):
        return (-1)* self.__sub__(other)

    def __str__(self):
        return (str(self.number))

    def __repr__(self):
        return " GF(%s,%s) " % (self.number,
                           self.__prime)

    def __truediv__(self, other):
        if(isinstance(other, GaloisField)):
            return GaloisField(self.number, self.__prime) * other.inverse()

        elif(isinstance(other, int)):
            return self.__truediv__(GaloisField(other, self.__prime))
        else:
            print("it should happen?")
            raise TypeError("Strange division by invalid type")

    #retirado do wikipédia
    #https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    def inverse(self):
        t = 0
        newt = 1

        r = self.__prime
        newr = self.number


        if(newr == 0):
            raise   ValueError("It's not possible to invert zero")


        while (newr != 0):
            quotient = r // newr
            t, newt = newt, t - quotient * newt
            r, newr = newr, r - quotient * newr

        if (r > 1):
            raise   ValueError("It's not invertible. Check if prime is prime")

        if (t < 0):
            t = t + self.__prime

        return GaloisField(t,self.__prime)
    


GF = GaloisField


def getGFarray(the_array,prime):
    """
    getGFarray(the_array,prime/finite_field)
    
    Transforms element-wise an np.array to GF elements.
    the first parameter is which array, list you want to transform,
    the second parameter must be the finite field, default is 2
    """
    func = np.vectorize(GaloisField)
     
    return func(the_array,prime)

def getINTarray(the_array): 
    """
    Transforms element-wise in an np.array to integer_array
    """
    func = np.vectorize(int)
    
    return func(the_array)




def swap_row(array,i,j):
    """
    swap rows in an array
    used for row echelon form
    """
    if i==None or j==None:
        pass
    else:
        array[i], array[j]= array[j].copy(), array[i].copy()
        #it seems the swapping needs two copy() functions, so it does not bug
        # right after
    

def has_pivot(line,n_th):
    """
    On a given line, 
    checks if it has a non-zero element on n_th index position
    (It should be a pivot, but it's not on a standalone match)
    
    It's more useful as a support method for find_pivot_line_index
    """
    if line[n_th] != 0:
        return True
    else:
        return False

def find_pivot_line_index(array,n_th):
    """
    Given an array, (which must be 2d)
    searchs for the first line (above to bottom)
    which has the n_th pivot 
    """
    pivot_line = None
    for index,line in enumerate(array):
        if index >= n_th:
            if has_pivot(line,n_th):
                pivot_line = index
                break;
    return pivot_line

def leftmost_nonzero_index(line):
    """
    Finds the leftmost nonzero element and gives its index:
        
        example:
            
            [0,0,1,2,0,0]
                 ^
            would return 2.
            
            [1,0,0]
             ^
            would return 0.
    
    """
    for index, element in enumerate(line):
        if(element != 0):
            return index
        
    return None

def normalize_line_by_pivot_inv(line):
    """
    Given the leftmost_nonzero element, divide the entire line by its value.
    Must be a GF array to work, since int elements in python does not have
    inverse.
    
    Also, return this normalization term
    """
    index = leftmost_nonzero_index(line)
    norm_term = line[index]**(-1)
    line *= (norm_term)*1
    
    return norm_term


    

def row_echelon_form(GF_array,has_supp_matrix=False,supp_matrix=[]):
    """
    DOING! (09/03/2019)
    Re-doing (11/03/2019)
    
    (20/03/2019, correction, just row)
    Gives the row echelon form (forma escalonada), which is the
    final matrix when using gaussian elimination of a SQUARED matrix.
    
    It receives an array which has it's elements in GF form already
    If not a numpy array, transform into it.
    
    It may not work for float types or anything else.
    
    (22/03/2019)
     supp_matrix and its boolean checker are used in middle steps of
     reduced_row_echelon_form. so it eventually returns the invert of the 
     GF array.
    """
    # https://stackoverflow.com/questions/ \14933577/swap-slices-of-numpy-arrays/14933939#14933939
    # The above code helped to implement this function
    
    the_array = GF_array.copy()
    
    if(has_supp_matrix == False):
        tobe_inv = (GF_array).copy() * 0 + np.eye(the_array.shape[0],dtype=int)
    else:
        tobe_inv = supp_matrix
    
    assert(the_array.ndim == 2)
    
    #pseudo código
    # para cada linha da matriz, execute:
    #  verifique se a linha do n-ésimo pivô tem um pivô, senão, encontre o que tem o pivô, e troque a linha
    #  divida essa linha pelo inverso do pivô da linha.
    #  subtraia as linhas abaixo da linha do pivô pela linha do pivô até que "não tenha pivô" nas linhas abaixo
    #  continue esse processo até cada linha ter um pivô
    
    for index,line in enumerate(the_array):
        if(not has_pivot(line,index) ):
            pivot_line_index = find_pivot_line_index(the_array,index)
            if(pivot_line_index != None):
                swap_row(the_array,index,pivot_line_index)
                swap_row(tobe_inv,index,pivot_line_index)
            else:
                line *= 0
                
                continue
                 
        norm_term  = normalize_line_by_pivot_inv(line)
        tobe_inv[index] *= norm_term * 1
        
        # like the below, just multiplying by 1 to avoid some problems
    

            
        for index_,line_ in enumerate(the_array):
            if(index_ > index and has_pivot(line_,index)):
                scale = line_[index]*1
                line_ -= getINTarray(scale*getINTarray(line))
                
                tobe_inv[index_] -= getINTarray(scale*getINTarray(tobe_inv[index]))
                
                # observation:
                # doing the above with x[0]+x[0] gives None results.
                # Also: scale = line_[index] without the *1 gives wrong results
                # There's some reference problem going on here
                # so we did an adaption which does not show any problem
                
                # I think the problem was on the swap function, but we'll leave it
                # the way it is here.
                
                # DO NOT REMOVE the *1 S!!!
    return the_array, tobe_inv
        
   
def reduced_row_echelon_form(GF_array):
    """
    This makes a Squared matrix to be an identity matrix
    also, transforms an identity matrix, doing the same operations
    as on the GF array, so it will be the inverse of the GF_array
    
    it returns (premises right):
        row_echelon_form of the array, and the identity transformed
        given the same operations done in GF_array
    """
    
    the_array, the_inv = row_echelon_form(GF_array)
    
    the_array = np.fliplr(the_array)
    the_array = np.flipud(the_array)
    
    the_inv = np.fliplr(the_inv)
    the_inv = np.flipud(the_inv)
    
    
    the_array,the_inv = row_echelon_form(the_array,True,the_inv)
    
    the_array = np.flipud(the_array)
    the_array = np.fliplr(the_array)
    
    the_inv = np.flipud(the_inv)
    the_inv = np.fliplr(the_inv)
    
    
   # the_array  = row_echelon_form(the_array)
    
    return the_array, the_inv
        
def inv_mod(GF_array):
    """
    Returns the inverse, if it exists, of the GF_array
    
    """
    return reduced_row_echelon_form(GF_array)[1]       
              
def inv_mod_v2(int_array,prime):
    """
    Return the inverse matrix, if it exists, on a int array mod prime
    
    it is based on the following link algorithm:
        http://home.ie.cuhk.edu.hk/~wkshum/wordpress/?p=1574
        
    ## This is Matlab code
    #p = 11; % a prime number
    #A = floor(rand(4,4)*p)  % a randomly generated integer matrix
    #
    #determinant = round(mod(det(A),p));
    #
    #if determinant == 0
    #   disp('matrix is singular')
    #else
    #    [d r s] = gcd(determinant,p);
    
    OBS:
        The extended Euclidean algorithm is particularly useful when 
        a and b are coprime. With that provision, 
        x is the modular multiplicative
        inverse of a modulo b,
        
        ax + by = gcd(a,b)
        
        here, r is the 'x' of the above expression
    
    #    B =  mod(round(inv(A)*det(A)*r),p)  % B is the inverse of A mod p
    #    mod(A*B,p)  % check that B times A mod p is the identity matrix
    #end;
    
    
    (22/03/2019)
    well...it works, before exploding the int size, or some overflow.
    Octave shows the same behavior. Eventually, the algorithm explodes.
    
    Even if it is faster than reduced row echelon form.
    
    It raises overflow error in 
    prime = 7
    and
    n_sources = 13...
    """
    
    determinant = np.round(np.linalg.det(int_array)) % prime
    
    r = int(GF(determinant,prime)**(-1))
    
    the_inv = getINTarray(np.round(np.linalg.inv(int_array)*np.linalg.det(int_array)*r)) %  prime
    
    return the_inv
  
    
    

##If this module is executed as a __main__ module, it does tests.
if __name__ == "__main__":
    import unittest
    class test_GaloisField(unittest.TestCase):
        
        some_primes = [2,3,5,7,11,13,17,19,23,29,31]
        # use the rest if you want some workload for your computer
        # and exhaustive long, testing
        # some_primes = [37,41,34,47,53,59,61,67,71,73,79]

        
        def test_equality(self):        
            for prime in self.some_primes:
                for value in range(prime):
                    self.assertEqual(GF(value,prime) , value % prime)
    
        def test_inverse(self):
            for prime in self.some_primes:
                for value in range(prime):
                    if(value != 0):
                        self.assertEqual
                        (GF(value,prime) * GF(value,prime).inverse() , 1)
                        
        def test_power(self):
            powers = 100
            
            for p in self.some_primes:
                for v in range(p):
                    for pwr in range(powers):
                        self.assertEqual(GF(v,p)**pwr , (v**pwr) % p)
    
                        
        def test_array(self):
            array_test = np.arange(81).reshape((9,9))
            
            array_GF_creator = np.vectorize(GF)
            
            for p in self.some_primes:
                array_GF = array_GF_creator(np.arange(81),p).reshape((9,9))
                bool_array = (array_GF == (array_test % p))
                self.assertTrue( bool_array.all() )
                
                array_test_squared = array_test.dot(array_test) % p
                array_GF_squared = array_GF.dot(array_GF)
                bool_array = ( array_test_squared == array_GF_squared)
                
                self.assertTrue( bool_array.all() )
        
        def test_reduced_row_echelon_form(self):
            import BSS_functions as bss
            for prime in self.some_primes:
                for n_sources in np.arange(2,10):
                    a = getGFarray(bss.generateMixingMatrix(n_sources,prime),prime)
                    
                    idty = reduced_row_echelon_form(a)[0]
                    
                    bool_array = ( a == a.dot(idty) )
                    self.assertTrue( bool_array.all())
                    
                    bool_array = ( a == idty.dot(a))
                    self.assertTrue( bool_array.all() )
                    
        def test_inv_mod(self):
            import BSS_functions as bss
            for prime in self.some_primes:
                for n_sources in np.arange(2,10):
                    a = getGFarray(bss.generateMixingMatrix(n_sources,prime),prime)
                    
                    the_inv = inv_mod(a)
                    
                    bool_array = ( a.dot(the_inv) == np.eye(a.shape[0],dtype=int))
                
                    self.assertTrue( bool_array.all() )
          
        def test_inv_mod_v2_test(self):
            ## Not executing in the test, actually
            def test_inv_mod_v2():
                import BSS_functions as bss
                for prime in self.some_primes:
                    for n_sources in np.arange(2,15):
                        a = bss.generateMixingMatrix(n_sources,prime)
                        int_a = getINTarray(a)
                        
                        print("prime is: ",prime)
                        print("n_sources is: ",n_sources)
                        inv_a = inv_mod_v2(int_a,prime)
                        print(int_a.dot(inv_a) % prime)
 
    unittest.main()           
