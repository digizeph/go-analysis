#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;

use GO::TermFinderReport::Text;

use GO::Utils::File qw (GenesFromFile);

print "Enter the fully qualified name of your obo file:\n";

chomp(my $ontologyFile = <STDIN>);

print "What is the aspect of this ontology (F, P, or C)?\n";

chomp(my $aspect = uc(<STDIN>));

print "Enter the fully qualified name of your associations file:\n";

chomp(my $annotationFile = <STDIN>); 

print "Enter a the fully qualified name of your file with a list of genes for which to find term:\n";

chomp(my $genesFile = <STDIN>);

print "How many genes (roughly) exist within the organism?\n";

chomp(my $totalNum = <STDIN>);

print "Finding terms...\n";

my $ontology   = GO::OntologyProvider::OboParser->new(ontologyFile => $ontologyFile,
						      aspect       => $aspect);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider => $annotation,
				     ontologyProvider   => $ontology,
				     totalNumGenes      => $totalNum,
				     aspect             => $aspect);

my @genes = GenesFromFile($genesFile);

my @pvalues    = $termFinder->findTerms(genes        => \@genes,
					calculateFDR => 1);

# now just print the info back to the client

my $report = GO::TermFinderReport::Text->new();

my $cutoff = 0.05;

my $numHypotheses = $report->print(pvalues  => \@pvalues,
				   numGenes => scalar(@genes),
				   totalNum => $totalNum,
				   cutoff   => $cutoff);

# if they had no significant P-values

if ($numHypotheses == 0){
    
    print "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
    
}

