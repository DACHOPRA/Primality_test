int lucasPrime (unsigned long long int n, unsigned long long int d)
{   
	int p, length;
	unsigned long long int q, u, u2, uold, v, v2, t;
	
	if (gcd(d, n) != 1)
		return (0);

    p = 1;
	q = (1 - d)/4;
    u = 0;
	v = 2;
	
	u2 = 1;
	v2 = 1;
	q2 = 2*q;
	
    t = (n + 1)/2;						//theta
	length = 64 - __builtin_clz(t); //length of our number in bits. //clz(b00010010) = 3 	
    h = 0
	
	for (i = 0; i < length; i++)
	{	
        u2 = (u2 * v2) % n;
        v2 = (v2 * v2 - q2) % n;
		if (t & 1)				//mask = 1
		{
			uold = u;
            u = u2 * v + u * v2;
			u = (u % 2 == 0) ? u : u+n;
            u = (u / 2) % n;
            v = (v2 * v) + (u2 * uold * d)
            v = (v % 2 == 0) ? v : v+n;  
            v = (v / 2) % n;
		}
		if (i < length - 1)
			q = (q*q) & n;
			q2 = 2*q;		
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