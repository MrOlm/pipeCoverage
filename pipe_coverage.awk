#! /usr/bin/awk -f 

# Version 0.3
# 1/11/16
# Matt Olm

# VERSION 0.3 UPDATE
# When very few reads map, it was still determining the read length based on the RC called
# in the beginning. Now, I update RC at the end to become j (which will be RC if that value
# is reached). This effects the start of END


# VERSION 0.2 UPDATE
# NOTE When dealing with buckets, I'm subtracting 1 from the value before dividing. This is
# because the value 3000 belongs in bucket 0, not bucket 1. I'm essentially just making it
# a zero-based value, which is fine because we'll never get a length 0 scaffold or position.
# This effects lines 11 and 20. 


# BL is bucket length (window length), RC = read count (reads to count to figure out read length)
# variables can be changed be passing them into awk with -v (awk -v BL=5000)
BEGIN{
if (BL == "") BL = 3000;
if (RC == "") RC = 10000
j = 0
}

# initialize stl (scaffold to length) array, stc (scaffold to count) array, and stw (scaffold to windows) array 
/^@SQ/ {$2 = substr($2,4,length($2)); $3=substr($3,4,length($3)) ;stl[$2]=$3; stc[$2]=0; stw[$2]=int(($3-1)/BL);
# initialize stwc (scaffold to window count) arrays
for (i=0; i <= stw[$2]; i++) 
stwc[$2][i]= 0 }


# Count scaffolds and buckets
!/^@/ {stc[$3]++;stwc[$3][int(($4-1)/BL)]++;
# Update the average read length if needed
if (j < RC) {lengths += length($10) ;j++}}

END{
# If there are no reads, arbitrarily decide on a value (doesn't matter anyways)
if (j == 0) rd_ln = 150;
else rd_ln = (lengths/j);
print ("# window length: " BL);
print ("# read length: " rd_ln " (based on first " j " reads)");
print ("# scaffold \t coverage \t length")
for (s in stl) {
 # Print the whole scaffold values
 printf (s "\t" (stc[s]*rd_ln)/stl[s] "\t" stl[s] "\n")
 # print the bucket values
 for (w in stwc[s]) {
  # Set the bucket length length to BL- if it's the last bucket will be overwritten below
  l = BL; 
  # If it's the last bucket, set it to the remaining value after a modulo
  if (w + 1 == length(stwc[s])) l = stl[s] % BL; 
  # If the remaining value after the modulo is 0, set it back to BL
  if (l == 0) l = BL; # fixing division by 0 if scaffold length is BL divisor
  printf(">" s "_" w+1 "\t" ((stwc[s][w])*rd_ln)/l "\t" l "\n");
  }
 }
}