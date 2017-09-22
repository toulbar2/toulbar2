#!/usr/bin/perl

=pod

=head1 exp_opt.pl

=over 4

=back

=head2 DESCRIPTION
        the script allows to mesure efficienty of different command line option declare in a given texte file
        
        i) create option file ==> each line will correspond to a experimentation
	 on a instances declared in a template file
		exemple :
		-l=2 -L=200 
		-l=4 -L=200
		....

        2) template file : by default named  CTestTestfile.cmake.template
	include a list of the benchmark that you want to iteratively used foreach condition
	if you want generate easly this file you just have to launch
	 ccmake ..	
	set the benchmark location ( benchdir)
	set test option with the string ("TOULBAR2_OPTION")


=head2 LISTE DES OPTIONS:

        
________________________________________________________

	    

        -help               =>  print curent help
        -verbose            => verbose  flag

        -opt_file           => filename containing list of command line option  (Ex:  -opt_file My_list_of_option  )  (format 1 line -> 1 option

				

        -name               => report file basename (default rap_opt_r$i) where $i = line rank in the command line list ( declared in -opt_file)

        -data_root          => root of tree directory containing test output (default ./All_res_tuned_option)

        -name               => report file basename (default rap_opt_r$i) where $i = line rank in the command line list ( declared in -opt_file)
        -rapdir             => location of the resulting report ( Default ./Rapport_all_opt/ )

        -MatchRegex         => file default regular expression
        -template 	    => template file containing test list in CTEST format (usually generate by cmake .. command ) default name ( CTestTestfile.cmake.template)
        -fregexp            => regular expresion for test output selection ( regexp on  filename) (Default *.out.*)

        -timeout            => timeout  in second ( ex: -timeout 30 )
        -score              => score flag ( require baseline file )   (default OFF)
        -baseline           => tyref result for score calculation ( structure = filename;energy)
        -exec               =>  execussion flag ( exec ctest command for each element defined in the command line option ( -opt_file filename ) (default exec is OFF) 
        -njob               =>  thread number usd by ctest for test execusion ( ex: -njob 4 ==> 4 jobs are launched simultantely)
________________________________________________________


=head2 USAGE : exp_opt.pl -opt_file My_option_list  -data_root ALL_TEST_ROOT_directory -template CTestTestfile.cmake.templeate



=cut


use File::Basename;
use Getopt::Long;
use Pod::Usage;
use strict;

my $opt_file;

my $template="CTestTestfile.cmake.template";
my $size=0;
my $exec; 
my @optionl; #list of option
my $data_root="All_res_tuned_option";
my $nbjob=4;
my $help        ;  # flag
my $verbose     ;  # flag
my $rapdir = "Rapport_all_opt"      ; 
my $name        ; 
my $timeout=20     ; # file default regular expression
my $fregexp     ; # option regular expression
my $report="rap_opt_r";
my $baseline="ref_best_restart.txt";
$baseline="resultat_ref.txt"; # default release 
my $score;

	    my $result = GetOptions(
        			"help"       => \$help,                # current help
				"verbose"       => \$verbose,             # verbose flag on
				"exec"       => \$exec,                   # exec ctest
				"data_root=s"     => \$data_root,         # path of project
				"opt_file=s"     => \$opt_file,           # file containing list option
				"rapdir=s"      => \$rapdir,
				"name=s"        => \$report,			#report basename
				"baseline=s"    => \$baseline,
				"template=s" => \$template,               # file default regular expression
				"timeout=i"     => \$timeout,         		# file default regular expression
				"job=i"     => \$nbjob,          		# thread number used during ctest execution
				"fregexp=s"     => \$fregexp                    # option regular expression
			    );

	    my %option = (
			"help"        => $help,                               # flag
			"verbose"     => $verbose,                   # flag
			"exec"        => $exec,                   # execustion flag ==> ctest execution
			"data_root"   => $data_root,               # string
			"rapdir"      => $rapdir,           # path of Report produce
			"name"        => $report,
			"template" => $template,                        # file default regular expression
			"job"      => 	$nbjob,          		# thread number used during ctest execution
			"fregexp"     => $fregexp                             # file regular expression
			 );


#####################################################################
  if ( defined $help ) {

        # partie execute sur la machine de lancement
        exists $option{help} || &usage(' $0 option');
        &usage('Help:') if ($help);
    }

#####################################################################

	    if(defined $opt_file ) {

		    open (FOPT, $opt_file) or die(" $0 enable to open $opt_file\n");
		    while (<FOPT>) {
			    if( $_=~m/\s*-\w/) 
			    {
				    chomp $_;
				    push(@optionl,$_);

			    }


		    }
	    } else {

       &usage();


	   } 

close FOPT;

foreach my $opt (@optionl) {
	open (TMP, $template) or die(" impossible d'ouvrir le dossier $template\n");
	print "building option $opt => rank =$size \n";
	open REPORT, "> CTestTestfile.cmake" or die $! ;

	while (<TMP>) {

		my $line = $_;
		$line=~s/"TOULBAR2_OPTION"/\"$opt\"/;

		print REPORT $line;

	}
	close REPORT;


	if( !(-e "$data_root" ) ) { `mkdir $data_root` };

	if( defined $exec ) {
		my $out =`ctest -j$nbjob`;
		`mv Test_done  $data_root/Test_done_r$size`;
	}


	my $ref;
	if (defined $score)  {
     	$baseline = "-score -baseline $baseline";
	}
	

		my $o = "./make_report.pl  $baseline -timeout $timeout -name ".$report."_$size -rapdir $rapdir -outpath $data_root/Test_done_opt_$size";
	print $o."\n";
	my	@s = `$o`;
	print @s;


	$size++;
}


#############################################################
sub usage ( $ ) {
        my ($msg) = @_;
        print STDERR "$msg\n";
        system("/usr/bin/pod2text $0 ");
        exit(1);
}


