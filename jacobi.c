#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h> 

int main(int argc, char *argv[])
{
	int jacobi(unsigned long long int d, unsigned long long int n);
	int lucasPrime (unsigned long long int n, unsigned long long int d);

	unsigned long long int d, n, p, dAbs, sign;
	int result;
	
	

	n = 117;
	d = 5;
	dAbs = 5;
	sign = 1;
	result = 0;
/*
	while (result != -1)
	{
		d = dAbs*sign;	
		result = jacobi(d,n);
		dAbs = dAbs + 2;
		sign = sign * -1;
	}
*/
//	printf("d %lld, n %lld, result %d \n",d,n,result);
	result = lucasPrime(d,n);
	printf("d %lld, n %lld, result %d \n",d,n,result);
	
}


int jacobi(unsigned long long int d, unsigned long long int n)
{

int t = 1;
unsigned long long int temp;

d = d % n;
while (d != 0) 
{
    while (d % 2 == 0)        //while d is even 
    {
        d = d / 2;
        if (n % 8 == 3 || n % 8 == 5) t = -t;
    }
    temp = d;
    d = n;
    n = temp;    
    if ((d % 4 == 3) &&  (n % 4 == 3)) t = -t;
    d = d % n;
}
if (n == 1) return t;
return 0;
}

int lucasPrime (unsigned long long int d, unsigned long long int n)
{   
	int p, i, length,k;
	unsigned long long int q, q2, u, u2, uold, v, v2, t;

//	if (gcd(d,n) != 1)
//		return (0);

	p = 1;
	q = (1 - d)/4;
	u = 0;
	v = 2;
	u2 = 1;
	v2 = p;
	q2 = 2*q;

	t = (n+1)/2;						//theta
	length = sizeof(unsigned long long int)*8 - __builtin_clzll(t); //length of our number in bits. //clz(b00010010) = 3 	
//	printf("length = %d\n", length);	

	for (i = 0; i < length-1; i++)
	{	
		u2 = (u2 * v2) % n;
		v2 = (v2 * v2 - q2) % n;
		if (t & 1)				//mask = 1
		{
			uold = u;
			u = (u2 * v) + (u * v2);
			if (u % 2 == 1 ) u = u+n;
			u = (u / 2) % n;
			v = (v2 * v) + (u2 * uold * d);
			if (v % 2 == 1 ) v = v+n;  
			v = (v / 2) % n;
		}

		q = (q*q) % n;
		q2 = q + q;

		t = t >> 1;
		printf("u = %lld\n", u);	
	}
	return (u == 0);
}

int gcd (unsigned long long int a, unsigned long long int b)
{
  unsigned long long int c;
  while (a != 0) {
     c = a; 
	 a = b%a;  
	 b = c;
  }
  return b;
}
