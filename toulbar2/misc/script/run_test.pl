#!/usr/bin/perl

=pod

=head1 run_test.pl 

=over 4

=back
=head2 USAGE :

 run_test.pl -w // 

=head2 DESCRIPTION
	run_test.pl is a perl script to execute test with toulbar2 program. 
	depending parameters the script allows to :
	
	i) exécute toulbar2  
	2) teste l'expression régulière
	3) génère un fichier de sortie contenant les résultats
=over 2 
=item exécute toulbar2:
=item teste l'expression régulière:
=item génère un fichier de sortie contenant les résultats:



=head2 OPTION LISTE:
		
	-wcsp				=> wcsp file 
	-uai				=> uai file 
    -help				=> print the current message
    -option				=> option du programme (string)
	-ub	 				=> ub the current directory
	-verbose			=> verbose    # flag ==> verbose mode
	-rank				=> date en second ajoutée au nom du fichier sortie 
	-timeout			=> temps maximum pour trouver la solution
	-help				=> help,     # flag
	-ub					=> ub,  #flag ==> ub the current directory
	-verbose			=> $verbose   # flag*
	-path_project		=> Chemin absolu du projet
________________________________________________________



=cut

use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use strict;

# dependance
# apt-get install libdatetime-per

# declaration des variables globales;
{    # debut main
    my $wcsp;
    my $uai;
    my $ub;
    my $default;
    my $regexp;
    my $rank;
    my $help;
    my $option;
    my $verbose;
    my $timeout;
	my $path_project="./";
    my $instance;
    

    # hash necessaire pour creer une ligne de command getopt
    # cf perldoc -f getopt
   
    my $result = GetOptions(
        "wcsp=s"    	=> \$wcsp,        # string
        "uai=s"    	=> \$uai,        # string
        "help"      	=> \$help,        # flag
        "option=s"  	=> \$option,      # option du programme (string )
        "ub=i" 			=> \$ub,  		  # flag ==> ub the current directory
        "default"		=> \$default,	  # flag ==> default the current directory
        "verbose"  	 	=> \$verbose,     # flag
        "rank=i"    	=> \$rank,        # integer
        "timeout=i" 	=> \$timeout,     # timeout du programme (string)
        "regexp"		=> \$regexp,	  # flag ==> default, ub or enum
        "path_project=s"=> \$path_project # path of project
    );

    my %options = (
        "wcsp"    		=> $wcsp,          # string
        "uai"    		=> $uai,          # string
        "help"    		=> $help,          # flag
        "option"  		=> \$option,       # option du programme (string )
        "ub"      		=> $ub,            # flag ==> ub the current directory
        "default"		=> $default,	   # flag ==> default the current directory
        "verbose" 		=> $verbose,       # flag
        "rank"    		=> \$rank,         # integer
        "timeout" 		=> \$timeout,      # timeout ( parametre de time out du test)
        "regexp=s"		=> \$regexp,	   # flag ==> default, ub or enum (expression reguliere du programme)
        "path_project"	=> $path_project   # string
    );

	my $test_directory = "$path_project/Test_done";
	my $evid;

############################################
    # creation du dossier de repertoire des tests si il n existe pas
    if ( !( -e $test_directory ) ) { `mkdir -p $test_directory` }
   
#############################################
    # affichage de l aide du programme
#############################################

    if ( defined $help ) {

        # partie execute sur la machine de lancement
        exists $options{file} || &usage('missing required -w option');
        &usage('Help:') if ($help);
        -r $wcsp || &usage("Could not open wcsp file ==> add -w option!");

        die("Could not open wcsp file ==> add -w option!");
    }
    
################
    my $toulbar2_bin = "$path_project/bin/Linux/toulbar2";

    # test existance du fichier exec toulbar2
    if ( !( -e $toulbar2_bin ) ) {
        print "binary toulbar2 not found";
        exit(-1);
    }
    my $problem;
    if(defined $wcsp ) {
    $problem = basename( $wcsp, ".wcsp" );
	$instance = $wcsp;
	}
	if (defined $uai) {
     $problem = basename( $uai, ".uai" );
	$instance = $uai;
	my $f = $uai;
	if( -e $uai.".evid" ) { $evid = $uai.".evid" }
	}

#################

    # definition d'un rank = date en second (timestamp);
    my $output = "$test_directory/$problem.out";

    my $timestamp=&MyTimeStamp;

    if ( defined $rank ) { 
    $rank.= $timestamp;
    $output.= $rank;
   } else {

    $output.= $timestamp;

  } 

    if ( defined $ub ) { $option.= "-ub=$ub" }
    if( $option=~m/TOULBAR2_OPTION/ ) { $option=~s/TOULBAR2_OPTION//;};
	
	open( FIC, ">$output" ) || die(" impossible d'ouvrir le fichier sorti ");
	print FIC "command_line = $option";
	print FIC "\n";
	print FIC "timeout = $timeout\n"; 
	close FIC;

	my $command = "$toulbar2_bin $instance $evid $option >> $output";
    print " affiche command = $command\n";
    my $command_out = `$command`;
#############################################
    # recupération de l'expression reguliere
#############################################
	
	
			my ( $etat, $res_simple ) = &simple_test_regex( $option, $output, $timeout );
			print "\n************\n$res_simple\n***************\n";
			print "$etat";



}    #fin du programme ( main )
#############################################################
sub usage ( $ ) {
    my ($msg) = @_;
    print STDERR "$msg\n";
    system("/usr/bin/pod2text $0 ");
    exit(1);
}

############
sub simple_test_regex {
    my $option     = shift;
    my $output     = shift;
    my $timeout    = shift;
    my $optimum    = "-";
    my $backtracks = "-";
    my $nodes      = "-";
    my $seconds;
    my $resultat;
    my $etat_test = "failed";    # msg
    open FIC, "$output" or warn "Problème pour ouvrir $output: $!\n";

# pour chaque ligne du fichier on test si elle est equivalente a l'expression reguliere
#"Optimum: 48 in 16 backtracks and 18 nodes and 0.05 seconds.";

    while (<FIC>) {

        if ( $_ =~ m/(Optimum:)\s(\d+)/ ) {
            $optimum = $2;

            if ( $_ =~ m/(\d+)\s+(backtracks)/ ) {
                $backtracks = $1;
            }
            if ( $_ =~ m/(\d+)\s+(nodes)/ ) {
                $nodes = $1;
            }
            if ( $_ =~ m/(\d+.\d*|\d+)\s+(seconds.)/ ) {
                $seconds = $1;
            }

        }

    }
    if ( $optimum ne "-" ) { $etat_test = "test ok"; }
    if ( not defined $seconds ) { $seconds = $timeout; }
    $resultat ="$output : optimum  $optimum , backtracks  $backtracks , nodes  $nodes , seconds  $seconds , option  $option";

    close FIC;
    return ( $etat_test, $resultat );
}

sub MyTimeStamp {
my @T=localtime(time);
my $timesamp=$T[0]+$T[1]*60+$T[2]*3600+$T[3]*3600*24+$T[5]*3600*24*365+$T[7]*3600*24;
return $timesamp;
}

