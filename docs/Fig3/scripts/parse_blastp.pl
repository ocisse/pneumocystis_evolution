#!/Users/cisseoh/anaconda3/envs/snakemake/bin/perl -w

use Data::Dumper;
use feature 'say';
use IO::All;
use Carp;

my %prot = ();
# load blastp report
my $f = io($ARGV[0]);
$f->autoclose(0);
while(my $l = $f->getline || $f->getline){
	chomp $l;
		next if $l =~m/^#/;
		#say "TEST\t$l";
		my @data = split/\t/, $l;		
		my ($q,$s,$e) = ($data[0],$data[1],$data[10]);
		#say "TEST2\t$q,$s,$e";
		if($prot{$q}){
			my @old = @{$prot{$q}};
			#say "TEST3\t$q\t@old";
			my ($s1,$e1) = @old[0..1];
			if ($e < $e1){
				my @updated = ();
				push(@updated,$1,$e1);
				@{$prot{$q}} = @updated;
			} else {
				@{$prot{$q}} = @old;
			}
		} else {
			my @old = ();
			push(@old, $s,$e);
			@{$prot{$q}} = @old;
		}
}
# extract the annotation from tcdb header
#say Dumper \%prot;

extract_info(%prot);

# subs
sub extract_info {
	my (%h) = @_;

	my $p = "";
	foreach $p (keys %h){
		my @data1 = @{$h{$p}};
		#say "TEST4\t$p\t@data1";
		my @data2 = split/\|/, $data1[0];
		#say "TEST5\t$p\t@data2";
		say "$p\t$data2[3]";
	}
}
