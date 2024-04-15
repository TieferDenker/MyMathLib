###########################################################################################################
# Module name: Angad's Library of Mathematical Functions
# Version: 0
# Author: Angad Singh
# Date created: 06.05.2022
# Module description:
# This module contains the list of almost all the mathematical functions ever made.
###########################################################################################################

#Important Libraries
import cmath
import math
from math import ceil, cos, exp, factorial, floor, gcd, log, pi, sin, sqrt
import random

#Famous Mathematical Constants
pi =  3.1415926535897932384626433
tau =  6.28318530717958647692528676
e = 2.7182818284590452353602874
gamma = 0.5772156649015328606065120
phi = 1.6180339887498948482045868
sq_2 = 1.4142135623730950488016887
sq_3 = 1.7320508075688772935274463
sq_5 = 2.2360679774997896964091736
zeta_2 = 1.6449340668482264364724151
zeta_3 = 1.2020569031595942853997381

# Small primes for deterministic primality testing
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

#Function Definitions:

###########################################################################################################
# Function          : is_prime
# Description       : This function checks whether a given number 'n' is a prime or not.
# Input parameters  : A natural number
# Return value      : 1_p(n)
###########################################################################################################
def is_prime(n, k=10):
    if n in SMALL_PRIMES:
        return True
    if n <= 1 or n % 2 == 0 or n%3 == 0 or n%5 == 0 or n%7 == 0 or n%11 == 0 or n%13 == 0 or n%17 == 0 or n%19 == 0 or n%23 == 0 or n%29 == 0 or n%31 == 0 or n%37 == 0 or n%41 == 0 or n%43 == 0 or n%47 == 0:
        return False

    # Write n as d*2^r + 1
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop
    for _ in range(k):
        a = random.randint(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

###########################################################################################################
# Function          : d
# Description       : This function returns the total number of divisors of any natural number 'Num'.
# Input parameters  : A natural number
# Return value      : d(n)
###########################################################################################################
def d(Num):
    Divisors = 0
    for Index in range(1,Num):
        if ((Num % Index) == 0):
            Divisors = Divisors + 1

    return Divisors

###########################################################################################################
# Function          : sigma
# Description       : This function returns the sum of all the positive divisors of any natural number 'Num'.
# Input parameters  : A natural number
# Return value      : σ(n)
###########################################################################################################
def sigma(Num):
    Sum_Of_Divisors = 0
    for Index in range(1,Num):
        if ((Num % Index) == 0):
            Sum_Of_Divisors=Sum_Of_Divisors + Index

    return Sum_Of_Divisors

###########################################################################################################
# Function          : distinct_Prime_Factors
# Description       : This function returns the number of distinct prime factors of any natural number 'Num'.
# Input parameters  : A natural number
# Return value      : ω(n)
###########################################################################################################
def distinct_Prime_Factors(Num):
    Prime_Factors = 0
    if ((Num % 2) == 0):
        Prime_Factors = Prime_Factors + 1
        while((Num % 2) == 0):
            Num = Num / 2

    for(i=3;i*i<=Num;i=i+2)
        if ((Num % i) == 0):
            Prime_Factors = Prime_Factors + 1
            while ((Num % i) == 0):
                Num=Num/i

    if (Num != 1):
        Prime_Factors++

    return Prime_Factors

###########################################################################################################
# Function          : Total_Prime_Factors
# Description       : This function returns the number of prime factors of any natural number 'Num'.
# Input parameters  : A natural number
# Return value      : Ω(n)
###########################################################################################################
def Total_Prime_Factors(Num):
    Prime_factors = 0
    while (Num%2)==0:
        Prime_factors = Prime_factors + 1;
        Num = Num//2;
    for i in range(3,int(sqrt(Num))+1,2)
        if (Num%i==0):
            while (Num%i)==0:
                Prime_factors = Prime_factors + 1
                Num = Num//i
    if not(Num == 1):
        Prime_factors = Prime_factors + 1
    
    return Prime_factors

###########################################################################################################
# Function          : per_sq
# Description       : This function checks whether a number 'Num' is a perfect square or not.
# Input parameters  : A natural number
# Return value      : per_sq(n)
###########################################################################################################
def per_sq(Num)
    if (ceil(sqrt(Num)) == floor(sqrt(Num))):
        return True
    return false

###########################################################################################################
# Function          : per_cb
# Description       : This function checks whether a number 'Num' is a perfect cube or not.
# Input parameters  : A natural number
# Return value      : per_cb(n)
###########################################################################################################
def per_cb(Num)
    if (ceil(cbrt(Num)) == floor(cbrt(Num))):
        return True
    return false

###########################################################################################################
# Function          : btod
# Description       : This function converts a number 'Num' from binary to decimal.
# Input parameters  : A natural number
# Return value      : btod(n)
###########################################################################################################
def btod(Num)
    decimal = 0
    power = 1
    while not(Num == 0):
        digit = Num%10
        decimal = decimal + (power*digit)
        power = power<<1
        Num = Num//10
    return decimal

###########################################################################################################
# Function          : dtob
# Description       : This function converts a number 'Num' from decimal to binary.
# Input parameters  : A natural number
# Return value      : dtob(n)
###########################################################################################################
def dtob(Num)
    binary = 0
    power = 1
    while not(Num == 0):
        digit = Num%2
        binary = binary + (power*digit)
        power = power*10
        Num = Num//2
    return binary

# Start editing from here:
###########################################################################################################
# Function          : factorial
# Description       : This function converts a number 'Num' from decimal to binary.
# Input parameters  : A natural number
# Return value      : per_sq(n)
###########################################################################################################
int factorial(int n) //returns the product of all natural numbers from 1 to 'n'
{
    int i;
    if(n==0 || n==1)
    return 1;
    return n*factorial(n-1);
}

int num_of_digits(int n) //returns the sum of digits of any natural number 'n'
{
    int digit=0;
    while(n!=0)
    {
        digit++;
        n=n/10;
    }
    return digit;
}

int sum_of_digits(int n) //returns the sum of digits of any natural number 'n'
{
    int digit, sum=0;
    while(n!=0)
    {
        digit=n%10;
        sum=sum+digit;
        n=n/10;
    }
    return sum;
}

int prod_of_digits(int n) //returns the product of digits of any natural number 'n'
{
    int digit, pro=1;
    while(n!=0)
    {
        digit=n%10;
        pro=pro*digit;
        n=n/10;
    }
    return pro;
}

int reverse(int n) //returns the reverse of digits of any natural number 'n'
{
    int i,sum=0,digit,m=pow(10,num_of_digits(n)-1);
    while(n!=0)
    {
        digit=n%10;
        sum=sum+m*digit;
        n=n/10;
        m=m/10;
    }
    return sum;
}

int prime_counting_fx(int n) //returns the number of primes less than or equal to any natural number 'n'
{
    int prime_numbers=0;
    iloop(2,n)
    {
        if(prime(i))
        prime_numbers++;
    }
    return prime_numbers;
}

int nthprime(int n) //returns the 'n'th prime number where 'n' is any natural number
{
    int prime_numbers=0;
    if(n==1)
    return 2;
    iloop(2,n*n)
    {
        if(prime(i))
        prime_numbers++;
        if(prime_numbers==n)
        return i;
    }
}

int fibo(int n) //returns the 'n'th Fibonacci number where 'n' is any natural number (n<=46)
{
    int a=0,b=1,c; 
    if(n==0)
    return 0;
    else if(n==1 || n==2)
    return 1;
    else
    {
        while(n>=2)
        {
            c=a+b;
            a=b;
            b=c;
            n--;
        }
        return c;
    }
}
int gold_con(int n) //returns the number of ways to write wany natural number 'n' as the sum of two odd primes (order is important)
{
    int count=0;
    iloop(3,n-3)
    {
        if(prime(i) && prime(n-i))
        count++;
    }
    return count;
}

int sqr_2(int n) //returns the number of ways to write any natural number 'n' as the sum of two squares (order is important)
{
    int i,count=0;
    for(i=1;i*i<n;i++)
    {
        if(per_sq(n-i*i))
        count++;
    }
    return count;
}

int sqr_3(int n) //returns the number of ways to write any natural number 'n' as the sum of three squares (order is important)
{
    int i,count=0;
    for(i=1;i*i<n;i++)
    {
        count=count+sqr_2(n-i*i);
    }
    return count;
}

int sqr_4(int n) //returns the number of ways to write any natural number 'n' as the sum of four squares (order is important)
{
    int count=0;
    iloop(1,n-1)
    {
        count=count+(sqr_2(i))*(sqr_2(n-i));
    }
    return count;
}

int cbr_2(int n) //returns the number of ways to write any natural number 'n' as the sum of two cubes (order is important)
{
    int i,count=0;  
    for(i=1;i*i*i<n;i++)
    {
        if(per_cb(n-i*i*i))
        count++;
    }
    return count;
}

int cbr_3(int n) //returns the number of ways to write any natural number 'n' as the sum of three cubes (order is important)
{
    int i,count=0;  
    for(i=1;i*i*i<=n;i++)
    {
        count=count+cbr_2(n-i*i*i);
    }
    return count;
}

int cbr_4(int n) //returns the number of ways to write any natural number 'n' as the sum of four cubes (order is important)
{
    int count=0;    
    iloop(1,n-1)
    {
        count=count+(cbr_2(i))*(cbr_2(n-i));
    }
    return count;
}






def dig(n):
    s=0
    for i in str(n):
        s=s+1
    return s

def sumdig(n):
    s=0
    for i in str(n):
        s=s+int(i)
    return s

def prodig(n):
    s=1
    for i in str(n):
        if i!='0':
           s=s*int(i)
    return s

def pi(n):
    c=0
    for k in range(1,n+1):
        c=c+prime(k)
    return c

def nthpri(n):
    c=0
    q=0
    for k in range(1,int(n**1.5)):
        q=q+prime(k)
        if q==n:
           c=k+1
        elif q==n+1:
            break  
    return c
     
def rev(n):
    s=0
    i=int(math.log10(n))
    while i>=0: 
        m=n%10
        s=s+m*(10**i)
        i=i-1
        n=(n-m)/10
    return int(s)

def nthfibo(n):
    phi=(1+5**0.5)/2
    x=-1/phi
    m=int((phi**n-x**n)/(5**0.5))
    return m

def gcd(m,n):
    while n!=0:
        x=m
        m=n
        n=x%n
    return m

def powofpri(n,p):
    c=0
    while n%p==0:
          c=c+1
          n=n/p
    return c

def perpow(n):
    x=[]
    c=0
    for i in range(2,n+1):
        if prime(i)==1 and n%i==0:
           x.append(powofpri(n,i))
    x.sort()
    for i in range(0,len(x)):
        if x[0]>=2 and x[i]%x[0]==0:
           c=c+1
    if c==len(x):
       return 1
    else:
          return 0

def angpri(n):
    c=0
    for k in range(40000,n):
           p=900*k*k+600*k+101
           q=900*k*k+1200*k+401
           if prime(p)==1 and prime(p+2)==1 and prime(2*p+1)==1:
              print(p)
              c=c+1
              #if prime(q)==1 and prime(q+2)==1 and prime(2*q+1)==1:
                 #print(k)
                 #c=c+1
           elif prime(q)==1 and prime(q+2)==1 and prime(2*q+1)==1:
                print(q)
                c=c+1
    return c

def hcn(n):
    c=0
    for i in range(1,n):
            if d(n)>d(i):
                c=c+1
    if c==n-1:
       return 1
    else:
          return 0

def semip(n):
    for i in range(2,int(n**0.5)+1):
        if prime(i)==1 and n%i==0 and prime(n/i)==1:
           return 1
    return 0

def prigap(n):
    i=2
    while (True):
        if prime(i+n)==1 and prime(i)==1 and pi(i+n)-pi(i)==1 :
            return i
        i=i+1
    return 0

def omega(n):
    count=0
    for i in range(2,n+1):
        if prime(i)==1 and n%i==0:
            count=count+powofpri(n,i)
    return count

def angch(n):
    m=d(n)+1
    while prime(d(n)+1)!=1:
          m=d(n)+1
          if prime(m)==1:
             break
          else:
                n=2*d(m)
    return (d(n)+1)

def weirdp(n):
    if prime(n)==1 and semip(n-1)==1:
       return 1
    return 0


Number_Of_Terms = 5

###########################################################################################################
# Function          : SawTooth
# Description       : This function returns the value of the sawtooth function used in signal analysis
# Input parameters  : A real number
# Return value      : ((x))
###########################################################################################################
def SawTooth(x):
    if floor(x) == ceil(x):
        return 0
    else:
        return (x - floor(x) - 1/2)

###########################################################################################################
# Function          : D
# Description       : This function returns the value of the Dedekind sum used in analytic number theory
# Input parameters  : An integer
# Return value      : D(a,b;c)
###########################################################################################################
def DedekindSum(a,b,c):
    Sum = 0
    for u in range(1,c):
        Sum = Sum + SawTooth(a*u/c)*SawTooth(b*u/c)
    return Sum

###########################################################################################################
# Function          : A
# Description       : This is an auxillary function used in the computation of the value of p(n)
# Input parameters  : Both the arguements are natural numbers
# Return value      : D(a,b;c)
###########################################################################################################
def A(k,n):
    Trig_Sum = 0
    for m in range(0,k):
        if gcd(m,k) == 1:
            Arg = pi * (DedekindSum(1,m,k) - (2*n*m/k))
            Trig_Sum = Trig_Sum + cos(Arg)
        m = m + 1
    return Trig_Sum

###########################################################################################################
# Function          : Derivative
# Description       : This function calculates the derivative of a specific function
#                     NOTE: It's not the general derivative operator
# Input parameters  : Both the arguements are natural numbers
# Return value      : Derivative(k,n)
###########################################################################################################
def Derivative(k,n):
    Aux_Val1 = sqrt(n-1/24)
    Aux_Val2 = 1/sqrt(n-1/24)
    Aux_Val3 = (pi/k)*sqrt(2/3)*Aux_Val1
    Aux_Val4 = 0.5*(pi/k)*sqrt(2/3)*Aux_Val2
    Diff_wrt_n_at_n_and_k = (Aux_Val1*math.cosh(Aux_Val3)*Aux_Val4-math.sinh(Aux_Val3)*Aux_Val2*0.5)/pow(Aux_Val1,2)
    return Diff_wrt_n_at_n_and_k

###########################################################################################################
# Function          : p
# Description       : This function calculates the number of possible partitions of a non-negative integer n
# Input parameters  : A natural number
# Return value      : p(n)
###########################################################################################################
def p(n):
    if (n == 0) or (n == 1) or (n == 2):
        return 1
    else:
        k = 1
        Result = 0
        while k <= Number_Of_Terms:
            Result = Result + A(k,n)*sqrt(k)*Derivative(k,n)
            k = k + 1
        return round(Result / (pi * sqrt(2)))

###########################################################################################################
# Function          : Bernoulli
# Description       : This function calculates the kth Bernoulli nuber
# Input parameters  : A natural number
# Return value      : B_k
###########################################################################################################
def Bernoulli(k):
    if k%2 == 1:
        if k == 1:
            return 1/2
        else:
            return 0
    else:
        if k ==  0:
            return 1
        elif k == 2:
            return 1/6
        elif k == 4:
            return -1/30
        elif k == 6:
            return 1/42
        elif k == 8:
            return -1/30
        elif k == 10:
            return 5/66
        elif k == 12:
            return -691/2730
        elif k == 14:
            return 7/6

###########################################################################################################
# Function          : E2k
# Description       : This function calculates the value of the Eisenstein series of weight 2k and complex period tau
# Input parameters  : k-> A natural number, tau-> A complex number
# Return value      : E_{2k}
###########################################################################################################
def E2k(k,tau):
    if tau.imag > 0:
        N = 200
        Sum = 0
        q = cmath.exp(complex(0,2*pi)*tau)
        for n in range(1,N+1):
            Sum = Sum + pow(n,2*k-1)*pow(q,n)/(1-pow(q,n))
        return (1 - (4*k*Sum)/Bernoulli(2*k))
    else:
        return cmath.inf

###########################################################################################################
# Function          : j
# Description       : This function calculates the j-invariant of the complex number tau used in complex analysis
# Input parameters  : tau-> A complex number
# Return value      : j(tau)
###########################################################################################################
def j(tau):
    E4 = E2k(2,tau)
    E6 = E2k(3,tau)
    return (1728*(pow(E4,3)/(pow(E4,3)-pow(E6,2))))

def b(d):
    tau = complex(1/2,sqrt(d)/2)
    return cmath.sqrt(d*(1728-j(tau)))

def a(d):
    tau = complex(1/2,sqrt(d)/2)
    E2 = E2k(1,tau)
    E4 = E2k(2,tau)
    E6 = E2k(3,tau)
    return (b(d)/6)*(1-(E4/E6)*(E2-(6/(pi*sqrt(d)))))

def CalculatePi(d):
    N = 5
    Sum = 0
    tau = complex(1/2,sqrt(d)/2)
    A = a(d)
    B = b(d)
    J = j(tau)
    print(A, B, J)
    for n in range(0,N+1):
        F = factorial(6*n)/(factorial(3*n)*(pow(factorial(n),3)))
        G = (A+n*B)/pow(J,n)
        Sum = Sum + F*G
    return cmath.sqrt(-J)/Sum

HeegnerNum = [1, 2, 3, 7, 11, 19, 43, 67, 163]

for d in HeegnerNum:
    try:
        print(d)
        print(CalculatePi(d))
    except ZeroDivisionError:
        continue