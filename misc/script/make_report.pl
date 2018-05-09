#!/usr/bin/perl

=pod

=head1 make_report.pl 

=over 4

=back

=head2 DESCRIPTION
	make_report.pl is a perl script that generate report  . 
	
	1) the script will scan -outpath directory and search for *.out* file (suposed to be solver output for a given instance) 
	2) it will extract information for each output found, information are defined by default in MatchRegex.txt file 
	3) it will print extracted information inside a report file (-name can change report name) and  -rapdir can change default location of this report

=head2 LISTE DES OPTIONS:

	
________________________________________________________

        -help   		=>  curent help
        -verbose	      	=> verbose  flag

        -name 			=> report file basename (default rapport_)
        -outpath		=> location of test output  (default ./Test_done )
	-rapdir			=> location of the resulting report (./Rapport)
	-MatchRegex		=> file default regular expression (./MatchRegexp.txt)

	-score			=> scoring flag trigger score calculation ( uai score require -baseline file ) (default off)
        -baseline 		=> filename contening  reference result for score calculation ( file structure = instance_name;energy)

        -fregexp		=> regular expresion for test output selection ( regexp on  filename) exp
        -timeout 		=> timeout  in second ( ex: -timeout 30 ) used in average calculation or for replacing execussion time for unterminated test.
________________________________________________________


=head2 USAGE : make_report.pl -rapdir MY_report_directory -name My_report_basename  -outpath My_toulbar2_instance_output   

=head2 Default command : make_report.pl 
	=+> 

=cut


use File::Basename;
use List::Util ;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use strict;


# declaration des variables globales;
{    # debut main

    my $path_project;
    my $report_directory;
    my $file_report;
    my $rep;
    my $help;
    my $fic;
    my @data;
	my $fregexp ; # expression reguliere de filtrage de fichier input
	my $MatchRegexp="./MatchRegexp.txt" ;
	my $report;
	my $head;
	my $cmd;
	my %structure;
	my  $solved=0;
	my $nsf=0;
	my $format_average;
	my $maximum;
	my $verbose;
	my $baseline;
	my %Hash_score_ref;
	my @score;
        my $timeout;
        my @Optimum;
	my $With_score;
	
    

    my $result = GetOptions(
        "help"      	=> \$help,        		# flag
        "verbose"      	=> \$verbose,        		# flag
        "score"      	=> \$With_score,        		# flag
        "outpath=s"		=> \$path_project, 		# path of project
		"rapdir=s"		=> \$report_directory,
		"name=s"		=> \$report,
		"baseline=s"		=> \$baseline,
		"MatchRegexp=s"	=> \$MatchRegexp,		# file default regular expression
		"timeout=i"	=> \$timeout,		# file default regular expression
        "fregexp=s"		=> \$fregexp 			# option regular expression
        );

    my %option = (
        "help"   		=> $help,          			# flag
        "verbose"      	=> \$verbose,        		# flag
        "score"      	=> $With_score,        		# flag scoring calculation
        "outpath"		=> $path_project,   		# string
		"rapdir=s"		=> $report_directory,  		# path of Report produce
		"name=s"		=> \$report,
		"MatchRegexp=s"	=> $MatchRegexp,			# file default regular expression
        "fregexp" 		=> $fregexp 				# file regular expression
        );
#########################
# base line
#################
if ( defined $baseline) {

		   open REF, "$baseline" or warn "Problème pour ouvrir $baseline: $!\n";
	while ( <REF> ) {
	chomp($_);
	my ($bench,$log)=split(";",$_);
	$bench=~s/\.\d+$//;
	$Hash_score_ref{$bench}=$log;
	if( defined $verbose ) { print "$bench ==> $log \n"; }

	}
	close REF;
}

#############################################
    # affichage de l aide du programme
#############################################
if(not defined $report ) {$report = "report";}

  if ( defined $help ) {

        # partie execute sur la machine de lancement
        exists $option{help} || &usage(' $0 option');
        &usage('Help:') if ($help);
    }

#############################################
# Chemin répertoire Test_done
   if( not defined $path_project) { $rep = "./Test_done"; }
   else { $rep = "$path_project" ; }
  print "out put location $rep \n";
   
# Chemin répertoire Rapports
   if( not defined $report_directory) { $report_directory 	= "./Rapports"; }
   
# Creation du repertoire s'il n'existe pas   
   if( !-e $report_directory ) {  `mkdir -p $report_directory`;}
   
#############################################

# ouverture du repertoire Test_done
   opendir( REP, $rep ) or die "Problème pour ouvrir $rep: $!\n";

   my $inputregexp= '.*\.out*';
   if( defined $fregexp ) {  $inputregexp= $fregexp;}
   print "-------------------------\n";
   print " file filetring -rgexp => $inputregexp \n";
   print "-------------------------\n";
   while ( defined( my $fic = readdir REP ) ) {
	   my $etat_test  = "failed";    # msg
	   my $resultats;
	   my $f = "${rep}/$fic";
	   my $option     = "-";

	   if ( $fic =~m/$inputregexp/ ) {

# recuperation du time out
		   open DATA, "$f" or warn "Problème pour ouvrir $fic: $!\n";
		   while (<DATA>) {
			   if ( $_ =~ m/^timeout =/ ) {
				   $timeout = $_;
				   $timeout =~ s/timeout = //;
				   chomp $timeout;
			   }

			   if ( $_ =~ m/^command_line/ ) {
				   $option = $_;
				   $option =~ s/command_line = //;
				   chomp $option;
				  $cmd=$option;
			   }
		   }

		   close DATA;


		   my $name = $fic  ; 
		   $name =~s/.out././;
		$name=~s/\.\d+$//;
		   
		# if Ref energy


# pour chaque ligne du fichier on test si elle est equivalente a l'expression reguliere
#"Optimum: 48 in 16 backtracks and 18 nodes and 0.05 seconds.";		   

		   open FIC, "$f" or warn "Problème pour ouvrir $fic: $!\n";

########################################
# lecture des expressions regulieres
#########################################
		my %hash_reg; # contient les regexp
		my %hash_rank; # contient les rank des données recupéré $1 $2...
		my %hash_mandatory;# si madatory = 1 ==> status = test ok sinon test failed;
		my %hash_results;

		   open FREG, "$MatchRegexp" or warn "regexp file no found \n";
		   while (<FREG>) {
	
			   if(($_ ne '\n') and !(/^(\s)*$/) and !(/^#/)) 
			   {
				   chomp $_;
				   my ($key,$rank,$reg,$mandatory) = split(";", $_);
				   $hash_reg{$key}=$reg;
				   $hash_rank{$key}=$rank;
				   $hash_mandatory{$key}=$mandatory;
				   if($key eq 'seconds' ) {
				   $hash_results{$key}=$timeout;
				} else { 
				   $hash_results{$key}='-';
				}
			   } else { if($verbose) {print "Match line skip :\n $_";} }
		   }
		   close FREG;

		   while (<FIC>) {
			   my $line = $_;
			   foreach my $r (keys %hash_reg)
			   {
				   my @data1;
				   my $rank = $hash_rank{$r}; 					#recuperation de la variable de rank $rank
					   (@data1)=($line=~m/$hash_reg{$r}/);
				   if(defined $data1[$rank-1] ) {
					   $hash_results{$r}=$data1[$rank-1];
					   

					   if ($hash_mandatory{$r} >0 ) { $etat_test= "test ok" ; $solved++ }
				   }
			   }
			   if($line=~m/^No solution in/) {
				   $hash_results{"Optimum"}="nsf";
				   $etat_test= "test ok" ; 
				   $solved++;
				   $nsf++;


			   }
		   }

		   $resultats =  "$name; $option";
		   $head="instance;option";
		   my $position=2;
		   while(my($key,$value) = each %hash_results){
#~ print "$key ==> $value";
			   $resultats .= ";$value";
			   $head.=";$key";
			   $structure{"$key"}=$position;
			   $position++;

		   }
		   $head.=";status\n";
		   $resultats.=";$etat_test\n";

		   print "$resultats";
		   push(@data, $resultats);

		   close FIC;
	   }
   }

   {
# Stocker les données dans un tableau
	   my @datas;
	   my @instance;
	   my @optimum;
	   my @backtracks;
	   my @nodes;
	   my @pretimes;
	   my @times;
	   my $aver_time;
	   my $aver_backtracks;
	   my $aver_nodes;
	   my $average;

#instance;option;seconds;backtracks;Pretime;nodes;Optimum;status
#~ shift(@data);
	   for (@data) { @datas = split(';', $_ );
		   print "@datas";

# Données utiles pour les calculs
#~ push( @instance,   $datas[1] );
#~ push( @optimum,    $datas[2] );
#	Optimum backtracks nodes seconds Pretime

		 if($verbose) {
		   print "seconds =".$structure{"seconds"}." ==> ".$datas[$structure{"seconds"}]."\n";
		   print "Pretime =".$structure{"Pretime"}." ==> ".$datas[$structure{"Pretime"}]."\n";
		   print "bactrack=".$structure{"backtracks"}." ==> ".$datas[$structure{"backtracks"}]."\n";
		}
		   my $name = $datas[0];
		
			$name=~s/\.\d+$//;
		  if($verbose) {
		   print "nom trouvé = $name \n";
		   }
		   push( @times,      $datas[$structure{"seconds"}] );
		   if(  $datas[$structure{"Optimum"}] ne '-' ) {
		   push( @Optimum,      $datas[$structure{"Optimum"}] );
		}
		   push( @backtracks,$datas[$structure{"backtracks"}] );
		   push( @nodes,       $datas[$structure{"nodes"}] );
		   push( @pretimes,   $datas[$structure{"Pretime"}] );
			if(defined $With_score) {
			   if(  $datas[$structure{"logLike"}] ne '-' ) {
				   push( @score,  score( $datas[$structure{"logLike"}], $Hash_score_ref{$name}));
			   }
			}
	   }


#	calcul moyenne et max du temps, des nodes, des backtracks 
	   print "aver size of pretime = ".scalar(@pretimes)."\n";

	   my $timemax = 0;			
	   $timemax = &max(@times);
	   $aver_time = &average(\@times);

	   my  $aver_pretime ;
	   $aver_pretime = &average(\@pretimes);

	   my  $max_pretime=0;
	   $max_pretime=&max(@pretimes);

	   my $backtrackmax = 0;
	   $backtrackmax = &max(@backtracks);
	   $aver_backtracks = &average(\@backtracks);

	   my $nodemax = 0;
	   $nodemax = &max(@nodes);
	   $aver_nodes = &average(\@nodes);

	   my $aver_score;
	   my $score_display;
	if(defined $With_score) {
	    print "Aver score =";
	   $aver_score=&average(\@score);
	}
	   $average = "aver_time;$aver_time;aver_backtracks;$aver_backtracks;aver_node;$aver_nodes;aver_pretime;$aver_pretime"."; solved; ".scalar(@Optimum). " \n";
	   $maximum = "max_time;$timemax;max_backtracks;$backtrackmax;max_node;$nodemax;max_pretime; $max_pretime;\n";
	if(defined $With_score) {
	   $score_display="Score total = ".&sum(\@score)." ; average score = $aver_score ; nb score =". scalar(@score) ."; nb optimum =".scalar(@Optimum). ";sovled;$solved;no solution;$nsf";
	   print "$score_display";
	}
	   print "$average";
	   print "$maximum";

	   $format_average = sprintf("aver_time;%5.3f",$aver_time);
	   $format_average .= sprintf(";aver_backtracks;%5.3f",$aver_backtracks);
	   $format_average .= sprintf(";aver_node;%5.2f",$aver_nodes);
	   $format_average .= sprintf(";aver_pretime;%5.2f",$aver_pretime);
	   $format_average .= ";sovled;$solved;no solution;$nsf";
	   $format_average .=";$cmd\n";
	if(defined $With_score) {
	   $format_average .="$score_display";
	}
   }
print "$format_average\n";

# générer un numéro pour le prochain rapport
opendir (REP, $report_directory) or die(" impossible d'ouvrir le dossier $report_directory\n");
my @contenu = grep { !/^\.\.?\z/ } readdir REP;

closedir REP;

my @num;
my $report1;
foreach (@contenu){
	$report1 = basename($_, ".txt");
	my @value = split ('_', $report1);
	push (@num, $value[1]);
}

my @num_triee = sort(tri_numerique @num);
my $num_max = @num_triee;
my $num_report = $num_max + 1;

$file_report = $report_directory."/".$report."_".$num_report.".txt";


#générer le fichier rapport
open( REPORT, ">>$file_report" )
|| die(" impossible de creer le fichier $file_report: $!\n");

#ecriture dans le fichier rapport
my @Data = sort mycriteria @data;
unshift(@Data, $head);
push (@Data, $format_average);
push (@Data, $maximum);

foreach my $v (@Data) { print REPORT $v ; }


#fermeture du fichier rapport
close REPORT;

print "------------------\nnom du rapport : $file_report\n";	

}    #fin du programme ( main )
#############################################################
sub usage ( $ ) {
	my ($msg) = @_;
	print STDERR "$msg\n";
	system("/usr/bin/pod2text $0 ");
	exit(1);
}

sub tri_numerique {
	if($a < $b)
	{ return -1; }
	elsif($a == $b)
	{ return 0; }
	else
	{ return 1; }
}

sub mycriteria {
	my ($aa) = $a =~ /_\d*_\d*_\d+_(\d+)/;
	my ($bb) = $b =~ /_\d*_\d*_\d+_(\d+)/;
#	my (@i) = split("_",$aa);
#	my (@j) = split("_",$aa);
#	$i[1] <=> $j[1];
#	print "my i $i[2] size = ".scalar(@i)."\n" ;
	$aa <=> $bb;
#	print "my a $aa my bb = $bb \n";
}

sub max {
	my $max = shift;
	$_ > $max and $max = $_ for @_;
	return $max
}


sub count_solve{
	my($data) = @_;
	my $n=0;
	if (not @$data) {
		die("Empty array\n");
	}
	my $total = 0;
	foreach (@$data) {
		if( /\d+/){
			$n++;
		}
	}
	return $n;
}
sub average{
	my $total = 0;
	my $n=0;
	my $average;
	my($data) = @_;
	if (not @$data) {
		die("Empty array\n");
	}
	foreach (@$data) {

	if( /\d+/ ){
		$total += $_;
		$n++;

	}
}
print "average calcul=> total = $total => over $n value \n";
if( $n > 0 ) {
	$average = $total / $n;
	return $average;
} else { exit "averge error  n= 0\n";}
}
sub stdev{
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}
sub min {
	my($a,$b) = @_;     # $a and $b are private variables
# and get values from the array @_

		if ($a < $b) {
			return $a;           # don't need parentheses
		} else {
			return($b);          # but you can use parentheses
		}
}

sub score{
	my $Ex=shift;
	my $refscore=shift;
	my $score;
	if(($refscore ne "-" ) and  ( $refscore != 0   )){
		$score = ($Ex-$refscore) / abs($refscore);
	} else { $score=0 };
	return($score);
}

sub sum{
	my($data) = @_;
	my $total=0;
	foreach  (@$data) { 
		$total += $_;
	}

	return($total);
}
