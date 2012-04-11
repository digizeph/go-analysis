#!/usr/bin/perl -w

system("rm -rf prep/* reports/* && perl extract.pl && perl func-ana-100.pl")
