// Compile with the command
// c++ -o csidh csidh.cpp -lgmp
#include <gmp.h>
#include <stdio.h>
#include <math.h> //log exp and ceil
#include <string.h> // strcmp


/*  How to use:
 * 
 *  Alice runs the commands
 *  ./makekey bigp Akey
 *  ./csidh bigp Akey 0 Afile
 * 
 *  Bob runs the commands
 *  ./makekey bigp Bkey
 *  ./csidh bigp Bkey 0 Bfile
 * 
 *  Alice and Bob exchange Afile and Bfile
 * 
 *  Alice runs the command
 *  ./csidh bigp Akey Bfile Sfile
 * 
 *  Bob runs the command
 *  ./csidh bigp Bkey Afile Sfile2
 * 
 *  Sfile and Sfile2 will be identical
 */

const short MAX_PRIMES = 100; // maximum number of small prime factors of prime_p+1
mpz_t prime_p;                      // prime_p will be a global variable, set in the early stages of the program.
gmp_randstate_t random_seed;        // random number seed.

class elliptic_point {
public:
	mpz_t A, x, y;
	// Full Constructor (which ironically is never used)
	elliptic_point(const mpz_t A_value, const mpz_t x_value, const mpz_t y_value) {
		mpz_init(A);
		mpz_init(x);
		mpz_init(y);
		mpz_set(A, A_value);
		mpz_set(x, x_value);
		mpz_set(y, y_value);
	}
	~elliptic_point() {
		// free up the memory
		mpz_clear(A);
		mpz_clear(x);
		mpz_clear(y);
	}
	elliptic_point(const mpz_t A_value) {
		// if only A is given, use the point at infinity, which is represented by x and y = -1.
		mpz_init(A);
		mpz_init(x);
		mpz_init(y);
		mpz_set(A, A_value);
		mpz_set_si(x, -1);
		mpz_set_si(y, -1);
	}
	void cp(const elliptic_point& other){
		// copies one point into another.
		// does not allocate new memory
		mpz_set(A, other.A);
		mpz_set(x, other.x);
		mpz_set(y, other.y);
	}
	bool is_infinity(){
		// test if the point is infinity, which is the identity element.
		if(mpz_sgn(x) < 0) {
			return true;
		}
		return false;
	}
	void randomize() {
		// find a random point on the curve
		// uses Algorithm 1
		// (assumes prime_p mod 4 = 3)
		mpz_t a, b, t;
		mpz_init(a);
		mpz_init(b);
		mpz_init(t); // temporary variable
		bool point_found = false;
		while(point_found == false){
			// generates random integer from 0 to prime_p-1
			mpz_urandomm(x, random_seed, prime_p); 
			//a = x^3 + A*x^2 + x mod prime_p
			mpz_powm_ui(a,x,2,prime_p);    // a = x^2 mod prime_p
			mpz_mul(t,A,a);                // t = A*a
			mpz_mod(a,t,prime_p);          // a = t mod prime_p
			mpz_powm_ui(b,x,3,prime_p);    // b = x^3 mod prime_p
			mpz_add(t,a,b);                // t = a + b
			mpz_add(y,t,x);                // y = t + x
			mpz_mod(a,y,prime_p);          // a = y mod prime_p
			// b = a^((prime_p-1)/2) mod prime_p
			mpz_fdiv_q_ui(t,prime_p,2);    // t = (prime_p-1)/2
			mpz_powm(b,a,t,prime_p);       // b = a^t mod prime_p
			if(mpz_cmp_si(b,1)==0){        // if b == 1 (happens 50% of the time)
				// y = a^((prime_p+1)/4)
				mpz_cdiv_q_ui(t,prime_p,4);// t = (prime_p+1)/4
				point_found = true;
				mpz_powm(y,a,t,prime_p);   // y = a^t mod prime_p
			}
		}
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(t);
	}
	void printpt(){
		// prints out the point information
		printf("(");
		mpz_out_str(stdout,10,A);
		printf(") [");
		if(mpz_cmp_si(x, -1) == 0) {
			printf(" Infinity ");
		} else {
			mpz_out_str(stdout,10,x);
			printf(", ");
			mpz_out_str(stdout,10,y);
		}
		printf("]\n");
	}
	void doublept() {
		// doubles the value of the point.
		if(mpz_cmp_si(x, -1) == 0) {
			// The point is infinity, which is its own double
			return;
		}
        if(mpz_cmp_si(y, 0) == 0) {
			// The double will be the point at infinity
			mpz_set_si(x, -1);
			mpz_set_si(y, -1);
			return;
		}
		mpz_t a, b, c, d, t;
		mpz_init(a);
		mpz_init(b);
        mpz_init(c);
        mpz_init(d);
		mpz_init(t); // temporary variables
		// t = 1/(2y) mod prime_p
		mpz_mul_ui(a,y,2);          // a = 2*y
		mpz_invert(t,a,prime_p);    // t = 1/a mod prime_p
		// a = (((x+2A)x+6)x+2A)x+1 mod prime_p
		mpz_mul_ui(a,A,2);          // a = 2*A
		mpz_add(b,x,a);             // b = x+a
		mpz_mul(c,b,x);             // c = b*x
		mpz_mod(b,c,prime_p);       // b = c mod prime_p
		mpz_add_ui(c,b,6);          // c = b+6
		mpz_mul(b,c,x);             // b = c*x
		mpz_mod(c,b,prime_p);       // c = b mod prime_p
		mpz_add(b,c,a);             // b = c + a
		mpz_mul(c,b,x);             // c = b*x
		mpz_add_ui(b,c,1);          // b = c+1
		mpz_mod(a,b,prime_p);       // a = b mod prime_p
		// y = (x^2-1)*a*t^3 mod prime_p
		mpz_powm_ui(b,x,2,prime_p); // b = x^2 mod prime_p
		mpz_sub_ui(c,b,1);          // c = b-1
		mpz_mul(b,c,a);             // b = c*a
		mpz_mod(a,b,prime_p);       // a = b mod prime_p
		mpz_powm_ui(b,t,3,prime_p); // b = t^3 mod prime_p
		mpz_mul(d,a,b);             // d = a*b
		mpz_mod(y,d,prime_p);       // y = d mod prime_p
		// x = (c*t)^2 mod prime_p
		mpz_mul(a,c,t);             // a = c*t
		mpz_powm_ui(x,a,2,prime_p); // x = a^2 mod prime_p
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
        mpz_clear(t);
	}
	void addpt( const elliptic_point& P ) {
		// adds the point P to the current point, using the curve geometry.
		if(mpz_cmp(A, P.A) != 0) {
			printf("Points do not have the same equation\n");
			return;
		}
        if(mpz_cmp_si(P.x, -1) == 0) {
			// The point P is infinity, so leave the current point alone.
			return;
		}
		if(mpz_cmp_si(x, -1) == 0) {
			// The current point is infinity, so copy P to the current point.
			mpz_set(x, P.x);
			mpz_set(y, P.y);
			return;
		}
		if(mpz_cmp(x, P.x) == 0) {
			// if the x coordinates are the same, there are two posibilities
			// if the y coordinates are also the same, we double the point.
			if(mpz_cmp(y, P.y) == 0) {
				doublept();
				return;
			} else {
			// otherwise, the y coordinates must be negatives of each other, so return infinity
				mpz_set_si(x, -1);
				mpz_set_si(y, -1);
				return;
			}
		}
		mpz_t a, b, c, d, e, m;
		mpz_init(a);
		mpz_init(b);
        mpz_init(c);
        mpz_init(d);
		mpz_init(m); // temporary variables
		// m = (y - P.y)/(x - P.x) mod prime_p
		mpz_sub(a,x, P.x);          // a = x - P.x
		mpz_invert(b,a,prime_p);    // b = 1/a mod prime_p
		mpz_sub(a,y, P.y);          // a = y - P.y
		mpz_mul(c, a, b);           // c = a*b
		mpz_mod(m, c, prime_p);     // m = c mod prime_p
		// x = m^2 - A - x - P.x mod prime_p
		mpz_powm_ui(a,m,2,prime_p); // a = m^2 mod prime_p
		mpz_sub(b,a,A);             // b = a - A
		mpz_sub(a,b,x);             // a = b - x
		mpz_sub(b,a,P.x);           // b = a - P.x
		mpz_mod(x,b,prime_p);       // x = b mod prime_p
		// b = P.y - m P.x  mod prime_p
		mpz_mul(a,m,P.x);           // a = m*P.x
		mpz_sub(c,P.y, a);          // c = P.y - a
		mpz_mod(b,c,prime_p);       // b = c mod prime_p
		// y = -(m x + b) mod prime_p
		mpz_addmul(b,m,x);          // b = b + m x
		mpz_neg(a,b);               // a = -b
		mpz_mod(y,a,prime_p);       // y = a mod prime_p
		mpz_clear(a);
        mpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
        mpz_clear(m);
	}
	void multpt( const mpz_t count){
		if(mpz_sgn(count) <= 0) {
			printf("Argument of mult must be a positive integer\n");
			return;
		}
		mpz_t m;
		mpz_init(m);
		mpz_set(m, count); // copy argument, since it will be changed
		elliptic_point Q(A);
		mpz_set(Q.x, x);
		mpz_set(Q.y, y);   // Q starts as the original point
		mpz_set_si(x, -1);
		mpz_set_si(y, -1);  // clear the current point
		while(mpz_cmp_ui(m,1)>0) { //while m > 1
			if (mpz_even_p(m)==0) { // if m is odd
				addpt(Q);           // add Q to the total
			}
			Q.doublept();           // Q = 2*Q
			mpz_fdiv_q_2exp(m,m,1); // m = m/2, (right-shift 1 bit)
		}
		addpt(Q);	                // since m = 1, we add the final Q
		mpz_clear(m);
	}
};

bool is_supersingular( const short small_prime_list[], const mpz_t A) {
	// tests whether y^2 = x^3 + A x^2 + x is supersingular (mod prime_p).  
	// that is, there must be prime_p+1 points on the curve.
	mpz_t d, elp, rootp, t;
	mpz_init(d);
	mpz_init(elp);
	mpz_init(rootp);
	mpz_init(t);
	elliptic_point P(A);
	elliptic_point Q(A);
	mpz_sqrt(rootp, prime_p);  //rootp = sqrt(prime_p), truncated to nearest integer
	short i;
	bool infiniteQ;
	while(true) {  // Infinitesimal chance of more than one iteration
		P.randomize(); // pick a random point
		mpz_set_ui(d,1); // d = 1;
		i = 0;
		while(small_prime_list[i] > 0) {
			mpz_set_si(elp, (long)small_prime_list[i]); // convert short prime to mpz variable
			mpz_cdiv_q(t, prime_p, elp); // t = (prime_p+1)/elp, since it rounds up.
			Q.cp(P); // copies P into Q
			Q.multpt(t);  // Q = [t]Q = [(prime_p+1)/elp]Q
			infiniteQ = Q.is_infinity();  // test to see if Q is infinity at this point.
			Q.multpt(elp); // Q = [elp]Q
			if(Q.is_infinity() == false) {
				// Q is definitely not supersingular
				mpz_clear(d);
				mpz_clear(elp);
				mpz_clear(rootp);
				mpz_clear(t);
				return false;
			}
			if(infiniteQ == false) {
				// we found an element of order elp, so we can increase the lower bound on the order.
				mpz_mul(d,d,elp); // d = d*elp
				if(mpz_cmp(d, rootp) > 0) { // test if d is bigger that sqrt(prime_p)
					// Q is definitely supersingular
					mpz_clear(d);
				    mpz_clear(elp);
					mpz_clear(rootp);
					mpz_clear(t);
					return true;
				}
			}
			i++;
		}
	}
	return true; // never will get here, but all branches must return something.
}

void velu(mpz_t A, const short small_prime) {
	// Sets the value of A to phi_el(A), using Velu's algorithm
	// find an element of order small_prime
	mpz_t a, b, c, d, h1, hn1;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_init(d);
	mpz_init(h1);
	mpz_init(hn1);
	mpz_cdiv_q_ui(a, prime_p, (unsigned long)small_prime); // a = (prime_p+1)/small_prime, since it rounds up.
	elliptic_point P(A); // starts out as infinity
	while(P.is_infinity()) {
		P.randomize();  // random point
		P.multpt(a);    // [(prime_p+1)/small_prime] P
	}
	// P will now be a point of order small_prime.  
	elliptic_point R(A);
	mpz_set_ui(h1, 1);
	mpz_set_ui(hn1, 1);
	for(short n=1; n <= (small_prime-1)/2; n++) {
		R.addpt(P);
		mpz_neg(a, R.x);          // a = -(x coordinate)
		mpz_add_ui(b,a,1);        // b = a + 1
		mpz_mul(c, h1, b);        // c = h1*b
		mpz_mod(h1, c, prime_p);  // h1 = c mod prime_p
        mpz_sub_ui(b,a,1);        // b = a - 1
        mpz_mul(c, hn1, b);       // c = hn1*b
        mpz_mod(hn1, c, prime_p); // hn1 = c mod prime_p
	}
	// d = ((A-2)/(A+2))^small_prime (h1/hn1)^8 mod prime_p
	mpz_invert(a, hn1, prime_p);  // a = 1/hn1 mod prime_p
	mpz_mul(b,a,h1);              // b = h1*a
	mpz_mod(a,b,prime_p);         // a = b mod prime_p
	mpz_powm_ui(c,a,8,prime_p);   // c = a^8 mod prime_p
	mpz_add_ui(a, A, 2);          // a = A+2
	mpz_invert(b, a, prime_p);    // b = 1/a mod prime_p
	mpz_sub_ui(a, A, 2);          // a = A-2
	mpz_mul(d, a, b);             // d = a*b
	mpz_mod(a, d, prime_p);       // a = d mod prime_p
	mpz_powm_ui(b, a, (unsigned long)small_prime, prime_p);// b = a^small_prime mod prime_p
	mpz_mul(a, b, c);             // a = b*c
	mpz_mod(d, a, prime_p);       // d = a mod prime_p
	// B = 2(1+d)/(1-d) mod prime_p
	mpz_ui_sub(a, 1, d);          // a = 1-d
	mpz_invert(b, a, prime_p);    // b = 1/a mod prime_p
	mpz_add_ui(a, d, 1);          // a = d + 1
	mpz_mul(c,a,b);               // c = a*b
	mpz_mul_ui(a,c,2);            // a = c*2
	mpz_mod(A, a, prime_p);       // A = a mod prime_p   
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	mpz_clear(d);
	mpz_clear(h1);
	mpz_clear(hn1);
	return;
}

void apply_key(mpz_t A, const short small_prime_list[], const short ei[], const int m){
	// Uses algorithm 4 to apply a privite key to the elliptic curve y^2 = x^3 + A x^2 + x.
	// This is designed to take exactly the same time, doing the same steps, regardless of the key.
	// This prevents side-channel attacks. 
	mpz_t B, C;
	mpz_init(B); // used to catch the one value we want.
	mpz_init(C); // used to store a result that we will not need.
	short i = 0;
	short k = 0;
	while(small_prime_list[i] > 0) {
		mpz_set(B,A);  // sets B to A, in case small_prime_list[i] = 0.
		if(ei[i] < 0) {
			mpz_neg(C,A); // C = -A, Use the quadratic twist to go backwards.
			mpz_mod(A,C,prime_p); // A = C mod prime_p 
		} else {
			mpz_neg(C,A); // take the same time doing nothing important
			mpz_mod(C,C,prime_p);
		}
		for(k = 1; k <= m; k++ ){
			velu(A, small_prime_list[i]);
			if( k == abs(ei[i]) ) {
				mpz_set(B,A);  // this is the value we want, store it in B.
			} else {
				mpz_set(C,A);  // take the same time doing nothing important
			}
		}
		mpz_set(A,B); // get the value we caught earlier
		if(ei[i] < 0) {
			mpz_neg(C,A); // C = -A, Reverse the quadratic twist, if we did it earlier.
			mpz_mod(A,C,prime_p); // A = C mod prime_p
		} else {
			mpz_neg(C,A); // take the same time doing nothing important
			mpz_mod(C,C,prime_p);
		} 		
		fprintf(stdout, "."); // make a progress bar
		fflush(stdout); //flush the stream so that the . is immediately printed.
		i++;
	}
	printf("\n");
	mpz_clear(B);
	mpz_clear(C);
}
  
int main(int argc, char* argv[])  {

  short small_prime_list[MAX_PRIMES];
  short key[MAX_PRIMES];
  if (argc != 5 ) { 
	printf("Applies the key to the starting elliptic curve, to produce the output elliptic curve\n");
    printf("Command: csidh prime_filename key_filename start_filename output_filename\n");
    printf("start_filename can be 0 to start with y^2 = x^3 + x\n");
    printf("Examples: ./csidh bigp Akey 0 Afile\n");
    printf("          ./csidh bigp Bkey Afile Sfile\n");
    return -1;
  }
  FILE *primefile;
  primefile = fopen(argv[1], "r");
  if(primefile == NULL) {
	  printf("File %s does not exist, or cannot be read\n", argv[1]);
	  return -1;
  } 
  int output;
  bool eof = false;
  // read in all of the numbers in the primefile,
  // until we hit the EOF
  short i = 0;
  while(eof == false){
	if(fscanf(primefile, "%d", &output) == 1) { 
	    small_prime_list[i] = (short)output;
	    i++;
	} else {
		eof = true;
	}
  }
  fclose(primefile);
  short num_primes = i;
  small_prime_list[i] = 0; // so routines can find the last prime.
  /*
  // print out the small primes
  for(i = 0; i< num_primes; i++){
	  printf("%d ",small_prime_list[i]);
  }
  printf(": %d primes in all.\n",num_primes);
  */
  FILE *keyfile;
  // read in the keyfile
  keyfile = fopen(argv[2], "r");
  if(keyfile == NULL) {
	  printf("File %s does not exist, or cannot be read\n", argv[2]);
	  return -1;
  }
  for(i=0; i< num_primes; i++) {
	  if(fscanf(keyfile, "%d", &output) == 1) {
		  key[i] = (short)output;
	  } else {
		  printf("Problem reading file %s\n", argv[2]);
		  return -1;
	  }
  }
  fclose(keyfile);
  /*
  // printing the key
  printf("private key:\n");
  for(i=0; i < num_primes; i++) {
 	  printf("%d ", key[i]);
  }
  printf("\n");
  */
  //compute the prime prime_p, check that it is a prime, and prime_p mod 8 = 3.
  mpz_init_set_ui(prime_p, 4); // set prime_p = 4.
  for(i = 0; i< num_primes; i++){
	  mpz_mul_ui(prime_p,prime_p,small_prime_list[i]);  // multiply by l_i 
  }
  mpz_sub_ui(prime_p,prime_p,1); // subtract 1.
  printf("p = ");
  mpz_out_str(stdout,10,prime_p); // outputs in base 10
  printf("\n");
  if(mpz_probab_prime_p(prime_p,15)== 0){ // prime_p should now be prime.
	  printf("WARNING! prime_p is not prime!\n");
	  mpz_clear(prime_p);
	  return -1;
  }
  mpz_t r;
  mpz_init(r);
  long r2 = mpz_mod_ui(r,prime_p,8); // finds prime_p mod 8
  mpz_clear(r);
  if(r2 != 3){  //  prime_p mod 8 should be 3 for this to work.
	  printf("WARNING! prime_p mod 8 should be 3.\n");
	  mpz_clear(prime_p);
      return -1;
  }
  // read in starting value, if not 0
  mpz_t Astart;
  mpz_init(Astart);  // initialize the starting vaule to 0
  if(strcmp(argv[3], "0") != 0) {
	  // read in the starting value from a file.
	  FILE *startfile;
	  startfile = fopen(argv[3], "r");
	  if(startfile == NULL) {
		  printf("File %s does not exist, or cannot be read\n", argv[3]);
		  mpz_clear(Astart);
		  mpz_clear(prime_p);
		  return -1;
	  }
	  if(mpz_inp_str(Astart, startfile, 10)== 0) { // read in the number
		  printf("An error occured trying to read a number in file %s\n", argv[3]);
		  mpz_clear(Astart);
		  mpz_clear(prime_p);
		  fclose(startfile);
		  return -1;
	  }
	  fclose(startfile);
  }
  printf("Starting with the value\nA = ");
  mpz_out_str(stdout,10,Astart); // outputs in base 10
  printf("\n");
  gmp_randinit_default(random_seed); // initialize seed
  if(is_supersingular(small_prime_list, Astart)){
	  printf("This represents a supersingular curve y^2 = x^3 + A x^2 + x.\n");
  } else {
	  printf("Warning: Starting value does not represent a supersingular curve.\n");
	  mpz_clear(Astart);
	  mpz_clear(prime_p);
	  gmp_randclear(random_seed);
	  return -1;
  }
  // compute the value of max_key, the maximum number in the key
  double prd = log(4.0);
  for(i = 0; i< num_primes; i++){
	  prd = prd + log((double) small_prime_list[i]);
  }
  prd = prd/((double) num_primes)/2.0;
  prd = exp(prd);  // this will be sqrt(prime_p)^(1/n)
  prd = (prd - 1.0)/2.0;
  // we want the smallest m such that (2m+1)^n > sqrt(prime_p)
  prd = ceil(prd);
  int max_key = (int)prd;
  //printf("Maximum magnitude in key  = %d\n", max_key);
  apply_key(Astart, small_prime_list, key, max_key);
  printf("Output = ");
  mpz_out_str(stdout,10,Astart); // outputs in base 10
  printf("\n");
  FILE *outfile;
  // output the result to a file
  outfile = fopen(argv[4], "w");
  if(outfile == NULL) {
	  printf("Error in writing to file %s\n", argv[4]);
  } else {
	  mpz_out_str(outfile, 10, Astart);
	  fclose(outfile);
  }
  gmp_randclear(random_seed);
  mpz_clear(prime_p); // clean up the mpz_t handles or else they will leak memory.
}
