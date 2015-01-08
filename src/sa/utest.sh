#!/bin/bash

# testing each executable using input/hp1

./ad_sa1.x ../../input/hp1 2 123456 | grep Final | awk -v p="ad_sa1 p1 2" -v t=1 -f utest.awk
./ad_sa2.x ../../input/hp1 2 123456 | grep Final | awk -v p="ad_sa2 p1 2" -v t=1 -f utest.awk
./ad_rsa.x ../../input/hp1 2 123456 | grep Final | awk -v p="ad_rsa p1 2" -v t=1 -f utest.awk

# testing each executable using input/hp9

./ad_sa1.x ../../input/hp9 2 123456 | grep Final | awk -v p="ad_sa1 p9 2" -v t=40 -f utest.awk
./ad_sa1.x ../../input/hp9 3 123456 | grep Final | awk -v p="ad_sa1 p9 3" -v t=56 -f utest.awk
./ad_sa2.x ../../input/hp9 2 123456 | grep Final | awk -v p="ad_sa2 p9 2" -v t=87 -f utest.awk
./ad_sa2.x ../../input/hp9 3 123456 | grep Final | awk -v p="ad_sa2 p9 3" -v t=119 -f utest.awk
./ad_rsa.x ../../input/hp9 2 123456 | grep Final | awk -v p="ad_rsa p9 2" -v t=35 -f utest.awk

# EOF


