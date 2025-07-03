// Compile with the command
// c++ -o makekey makekey.cpp
#include <stdio.h>
#include <stdlib.h> // atoi 
#include <assert.h>
#include <math.h> //log exp and ceil
#include <time.h> //srand(time(NULL))

const short MAX_PRIMES = 100; // maximum number of small prime factors of prime_p+1

int main(int argc, char* argv[])  {

  short small_prime_list[MAX_PRIMES];
  short key[MAX_PRIMES];
  
  if (argc != 3 ) { 
	printf("Creates a random key for the prime set up in prime_filename\n");
    printf("Command: makekey prime_filename output_filename\n");
    printf("Example: ./makekey bigp Akey\n");
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
  // print out the small primes
  for(i = 0; i< num_primes; i++){
	  printf("%d ",small_prime_list[i]);
  }
  printf(": %d primes in all.\n",num_primes);
  double prd = log(4.0);
  for(i = 0; i< num_primes; i++){
	  prd = prd + log((double) small_prime_list[i]);
  }
  prd = prd/((double) num_primes)/2.0;
  prd = exp(prd);  // this will be sqrt(prime_p)^(1/n)
  prd = (prd - 1.0)/2.0;
  // we want the smallest m such that (2m+1)^n > sqrt(prime_p)
  prd = ceil(prd);
  int max_key = (int)prd; // this will be the largest magnitude number in the key
  printf("Maximum magnitude in key = %d\n", max_key);
  FILE *keyfile;
  keyfile = fopen(argv[2], "w");
  if(keyfile == NULL) {
	  printf("File %s cannot be created\n", argv[2]);
	  return -1;
  }
  srand(time(NULL));
  for(i = 0; i< num_primes; i++){
	  // give a random integer between -max_key and max_key, inclusive
	  int j = (rand()%(max_key*2+1)) - max_key; 
	  fprintf(keyfile, "%d ", j);
  }
  fclose(keyfile);
}
