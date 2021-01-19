#!/Users/cisseoh/anaconda3/envs/snakemake/bin/perl -w 

use Data::Dumper;
use IO::All;
use Carp;
use feature 'say';


my $f = io($ARGV[0]);
$f->autoclose(0);
while(my $l = $f->getline || $f->getline){
	chomp $l;
		if ($l =~ m/^>/){
				my @data = split /\_/, $l;
				my $id = join("_",@data[0..1]);
				say $id;
			} else {
				say $l;
			}
}
