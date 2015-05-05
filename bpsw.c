#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h>	
	
int iMiller(uint64_t mpzN, long iB);
int iBPSW(uint64_t mpzN);
long mpz_scan1(uint64_t mpzNm1, int zero)
{
  
  long first_one = 0;
  while(mpzNm1%2 == zero)
  {
    mpzNm1 =  mpzNm1/2;
    first_one++;
  }
  return first_one;
}

uint64_t powm(long mpzB, uint64_t mpzd,  uint64_t mpzN)
{
  uint64_t expmod = 0, modulo=1;
  while(expmod++<mpzd)
  {
    modulo = (mpzB * modulo) % mpzN;
  }
  return modulo;
}

int main(int argc, char *argv[])
{
	int prime;			//0 = composite, 1 = prime
	uint64_t number =    32416166833;
	prime = iBPSW(number);
	printf("number %lld, %d\n",number,prime);
}



int iBPSW(uint64_t mpzN)
{
  printf("\n in bpsw\n");
int trialResult;
uint64_t mpzRem;

/* First eliminate all N < 3 and all even N. */


if(mpzN < 2)return(0);
if(mpzN == 2){ printf("\n1st check: prime");return(1);}
if(mpzN%2 == 0)return(0);

/* Check for small prime divisors p < 1000. */

trialResult=trialDiv(mpzN);
if(trialResult==0)
	return(0);

// printf("\nerror check\n");
/* Carry out the Miller-Rabin test with base 2. */

if(iMiller(mpzN, 2)==0)return(0);


/* Now N is a prime, or a base-2 strong pseudoprime with no prime
 * divisor < 1000. Apply the appropriate Lucas-Selfridge primality
 * test.
 */

// return(iLucasSelfridge(mpzN));
}


int trialDiv(uint64_t n)								//1 = composite, 0 = unknown
{
  // printf("\n in trialDiv\n");
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

for (i = 0; i < 167; i++)
{
	if (n%sPrimes[i] == 0 && n!=sPrimes[i]) return (0)	;		//composite
}

return(1);  /* No conclusion */
}

int iMiller(uint64_t mpzN, long iB)
{
  printf("\ninside miller\n");
/* Test N for primality using the Miller's strong probable prime
   test with base B. See Gary Miller's famous paper ("Riemann's
   hypothesis and tests for primality," Journal of Computation and
   System Science, 1976, Volume 13, pp 300-317).

   Returns 1 if N is a prime or a base-B strong probable prime.
   Returns 0 if N is definitely not a prime (composite or < 2).

   NOTE 1: Some will not regard this as a "pure" Miller's test with
   base B, since certain adjustments are made, prior to applying the
   algorithm, in order to work around invalid input values and
   improve performance:

   1) N < 3 and even N are screened out first.
   2) Multiples of the small primes (to qMax=# of binary digits in N)
      are returned as composite.
   3) Values of B < 2 are replaced by B=2.
   4) If N divides B exactly, B is replaced by B+1.

   If such adjustments are not made, a third return value (e.g., -1)
   must be allowed, indicating invalid input and an indeterminate result,
   and complicating the calling source code.

   NOTE 2: Not all authorities agree exactly on the conditions for
   Miller's test. Some allow B < 2 (Rosen, "Elementary number theory and
   its applications," 3rd ed., 1993, pp. 195-200), although there are good
   reasons to exclude such values. On the other hand, some require
   1 < B < N (Ribenboim, "The new book of prime number records,"
   3rd ed., 1996, pp. 113-116, 143-146). As far as I know, no one
   allows N to divide B, as this produces "false negatives"; e.g.,
   N=3 and B=6 fails Miller's test, thus indicating N=3 as composite.
   In practice, there appears to be no good reason to use any B < 2,
   and in fact its value is almost always chosen to be a small
   (positive) prime number. Furthermore, in my opinion, refusing to
   first screen out even values of N and N < 3 gratuitously adds
   unnecessary complication to the test.
*/

int mpzB = iB;
long s,j;
uint64_t mpzNm1, mpzd,mpzRem;

// if(mpz_divisible_p(mpzB, mpzN))mpz_add_ui(mpzB, mpzB, 1);

/* Now compute d and s, where d is odd and N - 1 = (2^s)*d. */

mpzNm1 =  mpzN -  1;
s=mpz_scan1(mpzNm1, 0); 
if(s==0)return(0);
uint64_t power = 1;
uint64_t product = 2;
printf("s = %lld\n", s);
while(power++<s)
{
  product = product*2;
}
mpzd = (uint64_t) mpzNm1/product;
  printf("\nd = %lld\n",mpzd);
/* Now proceed with the Miller's algorithm. First, if B^d is
   congruent to 1 mod N, N is a strong probable prime to base B. */

mpzRem =  powm(mpzB, mpzd, mpzN);
printf("Remainder = %lld\n", mpzRem);
if(mpzRem ==1){ printf("\n2nd check: prime");return(1);}


/* Now calculate B^((2^j)*d), for j=0,1,...,s-1 by successive
   squaring. If any of these is congruent to -1 mod N, N is a
   sprp base B. Start with j=0 and B^d, already computed.
   Miller's test uses repeated modular squaring in place of repeated
   modular exponentiation for speed (squaring is an order of
   magnitude faster). */

if(mpzRem ==  mpzNm1){ printf("\n3rd check: prime");return(1);} /* j=0 case */
uint64_t new_mod = 1;
for(j=1; j < s; j++)
  {
    printf("inside 4th check\n");
    while(new_mod++<mpzd+1)
{   mpzRem =  mpzB * mpzRem;
  mpzRem = (mpzRem)% mpzN;
  // printf("Remainder = %lld\n", mpzRem);
  if(mpzRem == mpzNm1 ){ printf("\n4th check: prime");return(1);}
}  new_mod = 1;

  }
  printf("error!\n");
return(0);
}

// int iLucasSelfridge(mpz_t mpzN)
// {
// /* Test mpzN for primality using Lucas's test with Selfridge's parameters.
//    Returns 1 if mpzN is prime or a Lucas-Selfridge pseudoprime. Returns
//    0 if mpzN is definitely composite. Note that a Lucas-Selfridge test
//    typically requires three to seven times as many bit operations as a
//    single Miller's test. The frequency of Lucas-Selfridge pseudoprimes
//    appears to be roughly four times that of base-2 strong pseudoprimes;
//    the Baillie-PSW test is based on the hope (verified by the author,
//    May, 2005, for all N < 10^13; and by Martin Fuller, January, 2007,
//    for all N < 10^15) that the two tests have no common pseudoprimes. */

// int iComp2, iP, iJ, iSign;
// long lDabs, lD, lQ;
// unsigned long ulMaxBits, ulNbits, ul, ulGCD;
// mpz_t mpzU, mpzV, mpzNplus1, mpzU2m, mpzV2m, mpzQm, mpz2Qm,
//       mpzT1, mpzT2, mpzT3, mpzT4, mpzD;

// #undef RETURN
// #define RETURN(n)           \
//   {                         \
//   mpz_clear(mpzU);          \
//   mpz_clear(mpzV);          \
//   mpz_clear(mpzNplus1);     \
//   mpz_clear(mpzU2m);        \
//   mpz_clear(mpzV2m);        \
//   mpz_clear(mpzQm);         \
//   mpz_clear(mpz2Qm);        \
//   mpz_clear(mpzT1);         \
//   mpz_clear(mpzT2);         \
//   mpz_clear(mpzT3);         \
//   mpz_clear(mpzT4);         \
//   mpz_clear(mpzD);          \
//   return(n);                \
//   }

// /* This implementation of the algorithm assumes N is an odd integer > 2,
//    so we first eliminate all N < 3 and all even N. As a practical matter,
//    we also need to filter out all perfect square values of N, such as
//    1093^2 (a base-2 strong pseudoprime); this is because we will later
//    require an integer D for which Jacobi(D,N) = -1, and no such integer
//    exists if N is a perfect square. The algorithm as written would
//    still eventually return zero in this case, but would require
//    nearly sqrt(N)/2 iterations. */

// iComp2=mpz_cmp_si(mpzN, 2);
// if(iComp2 < 0)return(0);
// if(iComp2==0)return(1);
// if(mpz_even_p(mpzN))return(0);
// if(mpz_perfect_square_p(mpzN))return(0);

// /* Allocate storage for the mpz_t variables. Most require twice
//    the storage of mpzN, since multiplications of order O(mpzN)*O(mpzN)
//    will be performed. */

// ulMaxBits=2*mpz_sizeinbase(mpzN, 2) + mp_bits_per_limb;
// mpz_init2(mpzU, ulMaxBits);
// mpz_init2(mpzV, ulMaxBits);
// mpz_init2(mpzNplus1, ulMaxBits);
// mpz_init2(mpzU2m, ulMaxBits);
// mpz_init2(mpzV2m, ulMaxBits);
// mpz_init2(mpzQm, ulMaxBits);
// mpz_init2(mpz2Qm, ulMaxBits);
// mpz_init2(mpzT1, ulMaxBits);
// mpz_init2(mpzT2, ulMaxBits);
// mpz_init2(mpzT3, ulMaxBits);
// mpz_init2(mpzT4, ulMaxBits);
// mpz_init(mpzD);

// /* Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
//    such that Jacobi(D,N) = -1 (Selfridge's algorithm). Although
//    D will nearly always be "small" (perfect square N's having
//    been eliminated), an overflow trap for D is present. */

// lDabs=5;
// iSign=1;
// while(1)
//   {
//   lD=iSign*lDabs;
//   iSign = -iSign;
//   ulGCD=mpz_gcd_ui(NULL, mpzN, lDabs);
//   /* if 1 < GCD < N then N is composite with factor lDabs, and
//      Jacobi(D,N) is technically undefined (but often returned
//      as zero). */
//   if((ulGCD > 1) && mpz_cmp_ui(mpzN, ulGCD) > 0)RETURN(0);
//   mpz_set_si(mpzD, lD);
//   iJ=mpz_jacobi(mpzD, mpzN);
//   if(iJ==-1)break;
//   lDabs += 2;
//   if(lDabs > ulDmax)ulDmax=lDabs;   tracks global max of |D| 
//   if(lDabs > INT32_MAX-2)
//     {
//     fprintf(stderr,
//       "\n ERROR: D overflows signed long in Lucas-Selfridge test.");
//     fprintf(stderr, "\n N=");
//     mpz_out_str(stderr, 10, mpzN);
//     fprintf(stderr, "\n |D|=%ld\n\n", lDabs);
//     exit(EXIT_FAILURE);
//     }
//   }

// iP=1;         /* Selfridge's choice */
// lQ=(1-lD)/4;  /* Required so D = P*P - 4*Q */

// /* NOTE: The conditions (a) N does not divide Q, and
//    (b) D is square-free or not a perfect square, are included by
//    some authors; e.g., "Prime numbers and computer methods for
//    factorization," Hans Riesel (2nd ed., 1994, Birkhauser, Boston),
//    p. 130. For this particular application of Lucas sequences,
//    these conditions were found to be immaterial. */

// mpz_add_ui(mpzNplus1, mpzN, 1); /* must compute U_(N - Jacobi(D,N)) */

// /* mpzNplus1 is always even, so the accumulated values U and V
//    are initialized to U_0 and V_0 (if the target index were odd,
//    U and V would be initialized to U_1=1 and V_1=P). In either case,
//    the values of U_2m and V_2m are initialized to U_1 and V_1;
//    the FOR loop calculates in succession U_2 and V_2, U_4 and
//    V_4, U_8 and V_8, etc. If the corresponding bits of N+1 are
//    on, these values are then combined with the previous totals
//    for U and V, using the composition formulas for addition
//    of indices. */

// mpz_set_ui(mpzU, 0);           /* U=U_0 */
// mpz_set_ui(mpzV, 2);           /* V=V_0 */
// mpz_set_ui(mpzU2m, 1);         /* U_1 */
// mpz_set_si(mpzV2m, iP);        /* V_1 */
// mpz_set_si(mpzQm, lQ);
// mpz_set_si(mpz2Qm, 2*lQ);

// ulNbits=mpz_sizeinbase(mpzNplus1, 2);
// for(ul=1; ul < ulNbits; ul++)  /* zero bit off, already accounted for */
//   {
// /* Formulas for doubling of indices (carried out mod N). Note that
//  * the indices denoted as "2m" are actually powers of 2, specifically
//  * 2^(ul-1) beginning each loop and 2^ul ending each loop.
//  *
//  * U_2m = U_m*V_m
//  * V_2m = V_m*V_m - 2*Q^m
//  */
//   mpz_mul(mpzU2m, mpzU2m, mpzV2m);
//   mpz_mod(mpzU2m, mpzU2m, mpzN);
//   mpz_mul(mpzV2m, mpzV2m, mpzV2m);
//   mpz_sub(mpzV2m, mpzV2m, mpz2Qm);
//   mpz_mod(mpzV2m, mpzV2m, mpzN);
//   if(mpz_tstbit(mpzNplus1, ul))
//     {
// /* Formulas for addition of indices (carried out mod N);
//  *
//  * U_(m+n) = (U_m*V_n + U_n*V_m)/2
//  * V_(m+n) = (V_m*V_n + D*U_m*U_n)/2
//  *
//  * Be careful with division by 2 (mod N)!
//  */
//     mpz_mul(mpzT1, mpzU2m, mpzV);
//     mpz_mul(mpzT2, mpzU, mpzV2m);
//     mpz_mul(mpzT3, mpzV2m, mpzV);
//     mpz_mul(mpzT4, mpzU2m, mpzU);
//     mpz_mul_si(mpzT4, mpzT4, lD);
//     mpz_add(mpzU, mpzT1, mpzT2);
//     if(mpz_odd_p(mpzU))mpz_add(mpzU, mpzU, mpzN);
//     mpz_fdiv_q_2exp(mpzU, mpzU, 1);
//     mpz_add(mpzV, mpzT3, mpzT4);
//     if(mpz_odd_p(mpzV))mpz_add(mpzV, mpzV, mpzN);
//     mpz_fdiv_q_2exp(mpzV, mpzV, 1);
//     mpz_mod(mpzU, mpzU, mpzN);
//     mpz_mod(mpzV, mpzV, mpzN);
//     }
// /* Calculate Q^m for next bit position, doubling the exponent.
//    The irrelevant final iteration is omitted. */
//   if(ul < ulNbits - 1)  /* Q^m not needed for MSB. */
//     {

//     mpz_mul(mpzQm, mpzQm, mpzQm);
//     mpz_mod(mpzQm, mpzQm, mpzN);  /* prevents overflow */
//     mpz_add(mpz2Qm, mpzQm, mpzQm);
//     }
//   }

// /* If U_(N - Jacobi(D,N)) is congruent to 0 mod N, then N is
//    a prime or a Lucas pseudoprime; otherwise it is definitely
//    composite. */

// if(mpz_sgn(mpzU)==0)RETURN(1);
// RETURN(0);
// }

