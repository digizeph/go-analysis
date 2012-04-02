#!/usr/bin/perl


use strict;
use warnings;
use diagnostics;

use File::Find;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;
use GO::TermFinderReport::Text;
use GO::Utils::File qw (GenesFromFile);

my $ontologyFile = 'data/GODB/gene_ontology.obo';
my $aspect = 'F';
my $annotationFile = 'data/GODB/gene_association.goa_human';

my $totalNum = 18410;
my $ontology   = GO::OntologyProvider::OboParser->new(ontologyFile => $ontologyFile,
                              aspect       => $aspect);
my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider => $annotation,
                     ontologyProvider   => $ontology,
                     totalNumGenes      => $totalNum,
                     aspect             => $aspect);


my @files;
my $indir="prep";
my $outdir="reports2";

mkdir $outdir if !(-d $outdir);
# add input files into array
#sub wanted{
#	push (@files, $_) if $_=~/.*\.symbol-100/;
#}
#find (\&wanted, $indir);

push (@files,'TA-drugable');
push (@files,'TP-drugable');

foreach my $file (sort @files){
	
    my $outFile = $file.'.txt';


    
    my @genes = GenesFromFile($file);

    my @pvalues    = $termFinder->findTerms(genes        => \@genes,
                        calculateFDR => 1);

    # now just print the info back to the client

    my $report = GO::TermFinderReport::Text->new();

    my $cutoff = 0.05;


    open (OUT, ">$outdir/$outFile") or die "open outfile failed!";

    my $numHypotheses = $report->print(pvalues  => \@pvalues,
                       numGenes => scalar(@genes),
                       totalNum => $totalNum,
                       cutoff   => $cutoff,
                       fh       => \*OUT);

    # if they had no significant P-values

    if ($numHypotheses == 0){
        
        print "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
        
    }
    close OUT;
}
