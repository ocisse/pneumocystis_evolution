#!/Users/cisseoh/anaconda3/envs/snakemake/bin/perl -w

use Data::Dumper;
use feature 'say';
use IO::All;
use Carp;

my ($pj,$pmac,$po,$pck1,$pck2,$pA,$pc,$pm,$pw) = @ARGV;

my %pj = load_data($pj);
my %pmac = load_data($pmac);
my %po = load_data($po);
my %pck1 = load_data($pck1);
my %pck2 = load_data($pck2);
my %pA = load_data($pA);
my %pc = load_data($pc);
my %pm = load_data($pm);
my %pw = load_data($pw);

# merge
my %all = (%pj, %pmac, %po, %pck1, %pck2, %pA, %pc, %pm, %pw);

say "DOM,Pj,PjDes,Pmac,PmacDes,Pory,PoryDes,PcanCk1,PcanCk1Des,PcanCk2,PcanCk2Des,PcanA,PcanADes,Pcar,PcarDes,Pmur,PmurDes,Pwk,PwkDes";
my $d = "";
foreach $d (keys %all){
	unless ($pj{$d}){ $pj{$d} = '0,0'}
	unless ($pmac{$d}){ $pmac{$d} = '0,0'}
	unless ($po{$d}){ $po{$d} = '0,0'}
	unless ($pck1{$d}){ $pck1{$d} = '0,0'}
	unless ($pck2{$d}){ $pck2{$d} = '0,0'}
	unless ($pA{$d}){ $pA{$d} = '0,0'}
	unless ($pc{$d}){ $pc{$d} = '0,0'}
	unless ($pm{$d}){ $pm{$d} = '0,0'}
	unless ($pw{$d}){ $pw{$d} = '0,0'}
	say "$d,$pj{$d},$pmac{$d},$po{$d},$pck1{$d},$pck2{$d},$pA{$d},$pc{$d},$pm{$d},$pw{$d}";

}

# sub
sub load_data {
	my %h = ();
	my $f = io(@_);
	$f->autoclose(0);
	while (my $l = $f->getline || $f->getline) {
		chomp $l;
			my @data = split /\,/, $l;
			$h{$data[0]} = join(",",@data[1..2]);
	}
	return(%h);
}