/*
      2,3,5,7,11     13     17     19     23     29 
     31     37     41     43     47     53     59     61     67     71 
     73     79     83     89     97    101    103    107    109    113 
    127    131    137    139    149    151    157    163    167    173 
    179    181    191    193    197    199    211    223    227    229 
    233    239    241    251    257    263    269    271    277    281 
    283    293    307    311    313    317    331    337    347    349 
    353    359    367    373    379    383    389    397    401    409 
    419    421    431    433    439    443    449    457    461    463 
    467    479    487    491    499    503    509    521    523    541 
    547    557    563    569    571    577    587    593    599    601 
    607    613    617    619    631    641    643    647    653    659 
    661    673    677    683    691    701    709    719    727    733 
    739    743    751    757    761    769    773    787    797    809 
    811    821    823    827    829    839    853    857    859    863 
    877    881    883    887    907    911    919    929    937    941 
    947    953    967    971    977    983    991    997
	First primes < 1000	*/	
	
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>	
	
int main(int argc, char *argv[])
{
	int prime;			//0 = composite, 1 = prime
	unsigned long long int number = 97;
	prime = iBPSW(number);
	printf("number %lld, %d",number,prime);
}



int iBPSW(unsigned long long int n)
{
int iComp2, isPrime;
int trialResult;


/* First eliminate all N < 3 and all even N. */


if(n < 2)return(0);
if(n == 2)return(1);
if(n%2 == 0)return(0);

/* Check for small prime divisors p < 1000. */

trialResult=trialDiv(n);
if(trialResult==0)
	return(0);


/* Carry out the Miller-Rabin test with base 2. */

if(iMillerRabin(mpzN, 2)==0)return(0);

/* The rumored strategy of Mathematica could be imitated here by
 * performing additional Miller-Rabin tests. One could also
 * carry out one or more extra strong Lucas tests. See the
 * routine iPrP in trn.c for an implementation.
 *
 * Now N is a prime, or a base-2 strong pseudoprime with no prime
 * divisor < 1000. Apply the appropriate Lucas-Selfridge primality
 * test.
 */

return(iLucasSelfridge(mpzN));
}


int trialDiv(unsigned long long int n)								//1 = composite, 0 = unknown
{
int i;
int sPrimes[168] = {2,3,5,7,11,13,17,19,23,29, 						//First 168 primes
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

for (i = 0, i < 167, i++)
{
	if (n%sPrimes[i] == 0) return (0)			//composite
}

return(1);  /* No conclusion */
}