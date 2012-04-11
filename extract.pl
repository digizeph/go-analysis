#!/usr/bin/perl
use strict;
use warnings;
use Tie::File;

sub isnan { ! defined( $_[0] <=> 9**9**9 ) }



open (ENTREZ, "data/MAP/entrez.txt") or die "Cannot open entrez.txt";
open (SYMBOL, "data/MAP/gene.txt") or die "Cannot open gene.txt";

# Look gene database, find every enterz->symbol mapping.
my (%ent2sym,%loc2ent);     # define entrez2symbol and local2entrez variable
while(<SYMBOL>){
   my @in = split("\t",$_);
   $ent2sym{$in[10]} = $in[1] unless $in[10] eq "";
}

# map every local column to entrezID
while(<ENTREZ>){
    my @in = split(" ",$_);
    $loc2ent{$in[1]+1} = $in[0];    # localID start from 0, have to add one
}

close ENTREZ;
close SYMBOL;

# check hash
my ($wc, $miss);
$wc = `wc -l data/MAP/entrez.txt`;
$miss = 0;
foreach my $key (sort {$a<=>$b}  keys %loc2ent){
    unless (exists $ent2sym{$loc2ent{$key}}){
        $miss++;
    #    print "symbol missing : $key -> $loc2ent{$key}\n";
    }
    else{
        #print "$key -> $ent2sym{$loc2ent{$key}}\n";
    }
}
print "\n$miss/$wc\n";


# Get sorted gene-symbol list
open (TARGET,"data/TARGETS/tumor.txt") or die "Cannot open tumor.txt";

tie my @db, 'Tie::File', "data/CIPHERDB/dcipher.txt" or die "Cannot tie :$!";

unless(-d "prep"){
    mkdir "prep" or die;
}

my $count=1;
while(<TARGET>){
    last if /^\s\n/;
    my %sorthash;
    my @tuple = split(" ",$_);
    print "$count : $tuple[0] $tuple[1]\n";
    $count++;
    open (OUT100, ">prep/$tuple[0]-$tuple[1].symbol-100");
    open (OUT, ">prep/$tuple[0]-$tuple[1].symbol");
    # Read into sorthash.
    my $id = 0;
    my @values = split(" ",$db[$tuple[1]-1]); # line of values;
    foreach my $value (@values)
    {   
        $id++;
        if(! isnan($value))
        {   
            $sorthash{$id}=$value;
        }
    }

    # Sort the hash! Largest first.
    my $round = 0;
    foreach my $col(sort{$sorthash{$b} <=> $sorthash{$a}} keys %sorthash)
    {   

        print OUT (" $ent2sym{$loc2ent{$col}}") if exists($ent2sym{$loc2ent{$col}});
        print OUT ("\n");

        $round++;
#        last if ($round==$LIMIT);
        if($round<=100){
            print OUT100 (" $ent2sym{$loc2ent{$col}}") if exists($ent2sym{$loc2ent{$col}});
            print OUT100 ("\n");
        }
    }

    close OUT;
    close OUT100;

}

close TARGET;
