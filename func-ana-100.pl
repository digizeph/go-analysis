#!/usr/bin/perl


use strict;
use warnings;
use diagnostics;

use File::Find;

use CGI qw/:all :html3/;

use IO::File;
use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;
use GO::View;
use GO::TermFinderReport::Html;
use GO::TermFinderReport::Text;
use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

my $ontologyFile = 'data/GODB/gene_ontology.obo';
my $aspect = 'F';
my $annotationFile = 'data/GODB/gene_association-large.goa_human';
#my $annotationFile = 'data/GODB/gene_association.goa_human';

my $totalNum = 48410;



# Create GO objects.

my $ontology   = GO::OntologyProvider::OboParser->new(ontologyFile => $ontologyFile,
                              aspect       => $aspect);
my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider => $annotation,
                     ontologyProvider   => $ontology,
                     totalNumGenes      => $totalNum,
                     aspect             => $aspect);

# HTMLs
my $confFile = "GoView.conf";
my $conf = &ReadConfFile($confFile);

$conf->{'totalNumGenes'} ||= $annotation->numAnnotatedGenes;
#push(@additionalArgs, ('totalNumGenes', $conf->{'totalNumGenes'}));
my $report = GO::TermFinderReport::Html->new();

&GenerateFrameset;
my $htmlFile = $conf->{'outDir'}.'batchGOViewList.html';
my $listFh = IO::File->new($htmlFile, q{>} )|| die "Cannot make $htmlFile : $!";




my @files;
my $indir="prep";
my $outdir="reports";

mkdir $outdir if !(-d $outdir);
# add input files into array
sub wanted{
	push (@files, $_) if $_=~/.*\.symbol-100/;
}
find (\&wanted, $indir);




########
########
# go through each list
########
########


foreach my $file (sort @files){
    my $outFile = $file.'.txt';

    my @genes = GenesFromFile("$indir/$file");

    my @pvalues    = $termFinder->findTerms(genes        => \@genes,
                        calculateFDR => 1);

    # now just print the info back to the client

    my $report = GO::TermFinderReport::Text->new();

    my $cutoff = 0.05;

    $file =~ /(.*)\.symbol-100/;
    my $name = $1;
    mkdir "$outdir/$name" unless -d "$outdir/$name";
    open (OUT, ">$outdir/$name/$outFile") or die "open outfile failed!";

    my $numHypotheses = $report->print(pvalues  => \@pvalues,
                       numGenes => scalar(@genes),
                       totalNum => $totalNum,
                       cutoff   => $cutoff,
                       fh       => \*OUT);
     $conf->{'outDir'}="$outdir/$name/";
                       
	my $goView = GO::View->new(-ontologyProvider   => $ontology,
			       -annotationProvider => $annotation,
			       -termFinder         => \@pvalues,
			       -aspect             => $conf->{'aspect'},
			       -configFile         => $confFile,
			       -imageDir           => $conf->{'outDir'},
			       -imageLabel         => "Batch GO::View",
			       -nodeUrl            => $conf->{'goidUrl'},
			       -geneUrl            => $conf->{'geneUrl'},
			       -pvalueCutOff       => $conf->{'pvalueCutOff'},
			       -maxTopNodeToShow   => 4);   # 4 is specifically suitable
	my $imageFile;
	if ($goView->graph) {
	
	$imageFile = $goView->showGraph;
	
    }
    
    my $htmlFile = &GenerateHTMLFile($file, $goView->imageMap, \@pvalues,
				     scalar($termFinder->genesDatabaseIds), "Terms for $file"); 

    print $listFh a({-href   => $htmlFile,
		     -target => 'result'}, $htmlFile), br;

    # if they had no significant P-values

    if ($numHypotheses == 0){
        
        print "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
        
    }
    close OUT;

}

$listFh->close;

########
########
# Missce Functions
########
########


sub GenerateHTMLFile{

    my ($file, $map, $pvaluesRef, $numGenes, $title) = @_;

    # work out name of html file
    
    my $htmlFile = $file;

    # delete anything up to and including the last slash

    $htmlFile =~ s/.*\///;

    # delete anything following the last period

    $htmlFile =~ s/\..*//;

    # now add an html suffix

    $htmlFile .= ".html";

    my $fullHtmlFile = $conf->{'outDir'}.$htmlFile;

    my $htmlFh = IO::File->new($fullHtmlFile, q{>} )|| die "Cannot make $fullHtmlFile : $!";

    print $htmlFh start_html(-title=>$title);

    print $htmlFh center(h2($title)), hr;

    print $htmlFh $map if defined $map;

    my $numRows = $report->print(pvalues      => $pvaluesRef,
				 aspect       => $conf->{'aspect'},
				 numGenes     => $numGenes,
				 totalNum     => $conf->{'totalNumGenes'},
				 fh           => $htmlFh,
				 pvalueCutOff => $conf->{'pvalueCutOff'},
				 geneUrl      => $conf->{'geneUrl'},
				 goidUrl      => $conf->{'goidUrl'});

    if ($numRows == 0){

	print $htmlFh h4(font({-color=>'red'}),
			 center("There were no GO nodes exceeding the p-value cutoff of $conf->{'pvalueCutOff'} for the genes in $file."));

    }

    print $htmlFh end_html;

    $htmlFh->close;

    return ($htmlFile);

}

sub ReadConfFile{

    my $confFile = shift;

    my %conf;

    my $confFh = IO::File->new($confFile, q{<} )|| die "cannot open $confFile : $!";

    while (<$confFh>){

	next if /^\#/; # skip comment lines

	chomp;

	next if /^\s*$/; # skip blank lines, or those without content

	next unless /(.+) = (.+)/;

	my ($param, $value) = ($1, $2);

	$value =~ s/\s$//;

	$conf{$param} = $value;

    }

    $confFh->close;

    if (!exists $conf{'annotationFile'} || !defined $conf{'annotationFile'}){

	die "Your conf file must specify an annotation file entry.";

    }elsif (!exists $conf{'ontologyFile'} || !defined $conf{'ontologyFile'}){

	die "Your conf file must specify an ontology file entry.";

    }elsif (!exists $conf{'aspect'} || !defined $conf{'aspect'}){

	die "Your conf file must specify an aspect entry.";

    }

    if (!exists $conf{'totalNumGenes'} || !defined $conf{'totalNumGenes'}){

	$conf{'totalNumGenes'} = ""; # simply make it the empty string for now

    }

    if (!exists $conf{'outDir'} || !defined $conf{'outDir'}){

	$conf{'outDir'} = ""; # set to empty string for now

    }

    $conf{'geneUrl'} ||= "";
    $conf{'goidUrl'} ||= "";

    $conf{'pvalueCutOff'} ||= 1;

    $conf{'calculateFDR'} ||= 0;

    # now make sure that file paths are treated relative to the conf file

    my $confDir = "./"; # default

    if ($confFile =~ /(.+)\//){

	$confDir = $1."/"; # adjust if necessary

    }

    foreach my $file ($conf{'annotationFile'}, $conf{'ontologyFile'}, $conf{'outDir'}){

	# $file is an alias for the hash entry

	if ($file !~ /^\//){ # if it's not an absolute path

	    $file = $confDir.$file; # add the confDir on the front

	}

    }

    # return a reference to the hash

    return \%conf;

}

sub GenerateFrameset{

# start an index file that a user can use to browse the output data,
# using frames

    my $framesFile = $conf->{'outDir'}."batchGOView.html";

    my $framesFh = IO::File->new($framesFile, q{>} )|| die "Cannot create $framesFile : $!";

    print $framesFh frameset({-cols         => "100, *",
			      -marginheight => '0',
			      -marginwidth  => '0',
			      -frameborder  => '1',
			      -border       => '1'},
			  
			     frame({'-name'       => "list",
				    -src          => "batchGOViewList.html",
				    -marginwidth  => 0,
				    -marginheight => 0,
				    -border       => 1}),
		   
			     frame({'-name'       =>'result',
				    -marginwidth  => 0,
				    -marginheight => 0,
				    -border       => 1}));

    $framesFh->close;

    return;

}

