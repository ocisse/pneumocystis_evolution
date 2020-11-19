#!/Users/cisseoh/anaconda3/envs/snakemake/bin/perl -w

use Data::Dumper;
use feature 'say';
use IO::All;
use Carp;



my %pj = load_data($ARGV[0]);
my %pmac = load_data($ARGV[1]);
my %po = load_data($ARGV[2]);
my %pcanCk1 = load_data($ARGV[3]);
my %pcanCk2 = load_data($ARGV[4]);
my %pcanA = load_data($ARGV[5]);
my %pc = load_data($ARGV[6]);
my %pm = load_data($ARGV[7]);
my %pw = load_data($ARGV[8]);

# merge All
my %all = (%pj, %pmac, %po, %pcanCk1, %pcanCk2, %pcanCk2, %pcanA, %pc, %pm, %pw);

# load tcdb annotation
my %ann = load_ann($ARGV[9]);
# itirate
say "ID,pj,pmac,pory,pcanCk1,pcanCk2,pcanA,pc,pm,pw,ANN";
my $t = "";
foreach $t (keys %all){
	unless ($pj{$t}){$pj{$t} = 0};
	unless ($pmac{$t}){$pmac{$t} = 0};
	unless ($po{$t}){$po{$t} = 0};
	unless ($pcanCk1{$t}){$pcanCk1{$t} = 0};
	unless ($pcanCk2{$t}){$pcanCk2{$t} = 0};
	unless ($pcanA{$t}){$pcanA{$t} = 0};
	unless ($pc{$t}){$pc{$t} = 0};
	unless ($pm{$t}){$pm{$t} = 0};
	unless ($pw{$t}){$pw{$t} = 0};
	
	say "$t,$pj{$t},$pmac{$t},$po{$t},$pcanCk1{$t},$pcanCk2{$t},$pcanA{$t},$pc{$t},$pm{$t},$pw{$t},$ann{$t}";
}
#
sub load_data {
	my @tcdb = ();
	my $f = io(@_);
	$f->autoclose(0);
	while(my $l = $f->getline || $f->getline){
		chomp $l;
			my @data = split /\t/, $l;
			push(@tcdb,$data[1]);
	}
	my $c = "";
	my %count = ();
	foreach $c (@tcdb){
		$count{$c}++;
	}
	return(%count);
}
sub load_ann {
	my %h1 = ();
	my $f1 = io(@_);
	$f1->autoclose(0);
	while (my $l1 = $f1->getline || $f1->getline) {
		chomp $l1;
			my @data1 = split /\,/, $l1;
			my ($acc,$des) = @data1[0..1];
			$h{$acc} = $des;
	}
	return(%h1);
}