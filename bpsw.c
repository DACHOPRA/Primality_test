//gcc -lm -o bpsw bpsw.c

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h>	
	
int main(int argc, char *argv[])
{
	long long int n; 
	int prime;			//0 = composite, 1 = prime
	for (n = 50001; n< 52500; n+=2)
	{	
	prime = BPSW(n);
	if (prime == 1)
	printf("number %lld, prime %lld\n",n,prime);
	}
}



int BPSW(long long int n)
{
int trialDiv(long long int n);	
long long int jacobi(long long int n);
int lucas (long long int d, long long int n);



long long int d = 0, temp;
int result, trialResult;

/* First eliminate all N < 3 and all even N. */

if(n < 2)return(0);
if(n == 2)return(1);
if(n%2 == 0)return(0);
temp = sqrt(n);
temp = temp*temp;
if (temp == n) return(0);


/* Check for small prime divisors p < 1000. */

trialResult = trialDiv(n);
if(trialResult==0) return(0);
if(trialResult==1) return(1);


/* Carry out the Miller-Rabin test with base 2. */

//if(MillerRabin(n, 2)==0) return(0);


/* Find d where the Jacobi value is -1 */

d = jacobi(n);

/* Do a Lucas probable prime test with d and n*/
return(lucas(d,n));
}


int trialDiv(long long int n)								//1 = composite, 0 = unknown
{
int i;
long long int sPrimes[168] = {						//First 168 primes
2,3,5,7,11,13,17,19,23,29, 						
31,37,41,43,47,53,59,61,67,71, 
73,79,83,89,97,101,103,107,109,113, 
127,131,137,139,149,151,157,163,167,173, 
179,181,191,193,197,199,211,223,227,229,
233,239,241,251,257,263,269,271,277,281, 
283,293,307,311,313,317,331,337,347,349, 
353,359,367,373,379,383,389,397,401,409, 
419,421,431,433,439,443,449,457,461,463, 
467,479,487,491,499,503,509,521,523,541, 
547,557,563,569,571,577,587,593,599,601, 
607,613,617,619,631,641,643,647,653,659, 
661,673,677,683,691,701,709,719,727,733, 
739,743,751,757,761,769,773,787,797,809, 
811,821,823,827,829,839,853,857,859,863, 
877,881,883,887,907,911,919,929,937,941, 
947,953,967,971,977,983,991,997};

for (i = 0; i < 167; i++)
{
	if (n == sPrimes[i]) return (1);
	if (n%sPrimes[i] == 0) return (0);			//composite
}

return(2);  /* No conclusion */
}


long long int jacobi(long long int n)
{
	int result, t;
	long long int d, dAbs,sign, temp, n1, d1;

	result = 0;
	dAbs = 5;
	sign = 1;

	while (result != -1)		//if result != -1, increment d and try again
	{
	d = dAbs*sign;	
	t = 1;
	n1 = n;				//reinitialize n1 to n
	d1 = d;				//reinitialize d1 to d
	d1 = d1 % n1;

	while (d1 != 0) 		
	{
	    while (d1 % 2 == 0)        //while d is even 
	    {
		d1 = d1 / 2;
		if (n1 % 8 == 3 || n1 % 8 == 5) t = -t;
	    }
	    temp = d1;
	    d1 = n1;
	    n1 = temp;    
	    if ((d1 % 4 == 3) &&  (n1 % 4 == 3)) t = -t;
	    d1 = d1 % n1;
	}
	if (n1 == 1) result = t;
	else result = 0;
	dAbs = dAbs + 2;
	sign = sign * -1;
	}
	return d;
}



int lucas (long long int d, long long int n)
{   
	int gcd (long long int a, long long int b);

	int p, i, length,k;
	long long int q, q2, u, u2, uold, v, v2, t;

	if (gcd(d,n) != 1)
		return (0);

	p = 1;
	q = (1-d)/4;
	u = 0;
	v = 2;
	u2 = 1;
	v2 = 1;
	q2 = 2*q;
	t = (n+1)/2;						//theta
	length = 64 - __builtin_clzll(t); //length of our number in bits. //clz(b00010010) = 3 	

	for (i = 0; i < length; i++)
	{	
		u2 = (u2 * v2) % n;
		v2 = (v2 * v2 - q2) % n;
		if (t & 1)				//mask = 1
		{
			uold = u;
			u = (u2 * v) + (u * v2);
			if (u % 2 == 1 )
			 u = u+n;
			u = (u / 2) % n;
			v = (v2 * v) + (u2 * uold * d);
			if (v % 2 == 1 ) v = v+n;  
			v = (v / 2) % n;
		}

		q = (q*q) % n;
		q2 = q + q;

		t = t >> 1;
	}
	return (u == 0);
}

int gcd (long long int a, long long int b)
{
  long long int c;
  while (a != 0) {
     c = a; 
	 a = b%a;  
	 b = c;
  }
  return b;
}
