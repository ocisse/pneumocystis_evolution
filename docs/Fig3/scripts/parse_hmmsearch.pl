#!/Users/cisseoh/anaconda3/envs/snakemake/bin/perl -w

use Data::Dumper;
use feature 'say';
use IO::All;
use Carp;
use Getopt::Long;

my $cutoff = 1e-5;
my $verbose;

GetOptions ("cutoff=s" => \$cutoff,
			"verbose" => \$verbose)
or die ("Error in command line arguments\n");

my %prot2domains = (); # key = "protein|domain" value = domain
#my @doms;
my %description = ();
my $f = io($ARGV[0]);
$f->autoclose(0);
while(my $l = $f->getline || $f->getline){
	chomp $l;
		next if $l =~m/^#/;
		my @data = split /\s+/, $l;
		my ($pid,$des,$domain,$e) = ($data[0],$data[3],$data[4],$data[6]);
		unless ($e > $cutoff){
			$prot2domains{"$pid|$domain"} = $domain;
			#(@doms,$domain);
		}
		$description{$domain} = $des;
}
my $d = "";
my %doms = ();
foreach $d (values %prot2domains){
	$doms{$d}++;
}

my $d1 = "";
foreach $d1 (keys %doms){
	say "$d1,$doms{$d1},$description{$d1}";
}