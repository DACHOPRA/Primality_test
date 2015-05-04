#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h> 

int main(int argc, char *argv[])
{
  int jacobi(unsigned long long int d, unsigned long long int n);

  unsigned long long int d, n;
  int result;

  d = 111;
  n = 3337;
  result = jacobi(d,n);
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