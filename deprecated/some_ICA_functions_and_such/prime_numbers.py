# -*- coding: utf-8 -*-
from math import floor, sqrt

def is_prime(number):
    """
    Checks if number is a prime number.

    Naive implementation
    """
    for i in range(2, floor(sqrt(number)) + 1):
        if (number % i) == 0:
            return False

    if number < 2:
        return False

    return True


## Execute tests if module is executed as normal program
if __name__ == "__main__":
    from pocha import describe, it
    

    @describe('Prime module')
    def _():
        some_primes = [2, 3, 5, 7, 11, 13, 17, 23, 29, 31]
        some_non_primes = [1, 4, 8, 15, 21, 35]
    
        @it('is_prime should return True for primes', tags=['is_prime'])
        def check_primes():
            for prime in some_primes:
                assert is_prime(prime)
    
        @it('is_prime should return False for non-primes', tags=['is_prime'])
        def check_non_primes():
            for non_prime in some_non_primes:
                assert not is_prime(non_prime)
