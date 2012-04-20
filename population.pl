open (ENTREZ, "data/MAP/entrez.txt") or die "Cannot open entrez.txt";
open (SYMBOL, "data/MAP/gene.txt") or die "Cannot open gene.txt";
open (OUT, ">data/MAP/population.txt");

# Look gene database, find every enterz->symbol mapping.
my (%ent2sym,%loc2ent);     # define entrez2symbol and local2entrez variable
while(<SYMBOL>){
   my @in = split("\t",$_);
   $ent2sym{$in[10]} = $in[1] unless $in[10] eq "";
}

# map every local column to entrezID
while(<ENTREZ>){
    my @in = split(" ",$_);
    print OUT "$ent2sym{$in[0]}\n" unless  $ent2sym{$in[0]}eq "";
#    $loc2ent{$in[1]+1} = $in[0];    # localID start from 0, have to add one
}

close ENTREZ;
close SYMBOL;
