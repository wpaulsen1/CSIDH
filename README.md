# CSIDH
C++ Source code for implementing the CSIDH cryptosystem

CSIDH stands for Commutative Supersingular Isogeny Diffie-Hellam.

CSIDH is likely to be quantum resistant, unlike other crptosystems such as RSA.

Files included:

bigp        -- list of the small primes that produce the large prime number used by the cryptosystem

makekey.cpp -- source code for creating a randomized private key.

csidh.cpp   -- source code for appling a private key to a supersingular curve.

Makes use of the Gnu Multiple Precision Library (GMP), which must first be installed.

To compile sourcecode:

    c++ -o makekey makekey.cpp

    c++ -o csidh csidh.cpp -lgmp

Usage:

In order for Alice and Bob to communicate securely, they need to agree on an encryption scheme key.

   Alice runs the commands
   
       ./makekey bigp Akey
   
       ./csidh bigp Akey 0 Afile
  
   Bob runs the commands
   
       ./makekey bigp Bkey
   
       ./csidh bigp Bkey 0 Bfile
  
   Alice and Bob exchange Afile and Bfile
  
   Alice runs the command
   
       ./csidh bigp Akey Bfile Sfile
  
   Bob runs the command
   
       ./csidh bigp Bkey Afile Sfile2
  
   Sfile and Sfile2 will be identical, so these can be used for an encryption scheme key.

   Note that an evedropper (Eve) cannot determine Sfile from the Afile and Bfile.
