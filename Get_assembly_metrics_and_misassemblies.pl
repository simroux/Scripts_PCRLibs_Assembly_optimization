#!/usr/bin/env perl
use strict;
use autodie;
use File::Spec::Functions;
use Getopt::Long;
my $h='';
my $wdir='';
GetOptions ('help' => \$h, 'h' => \$h, 'd=s'=>\$wdir);
if ($h==1 || $wdir eq ""){ # If asked for help or did not set up any argument
	print "# Script to parse assembly
# Arguments : 
# -d: working directory
# Must include the following files:
# A fasta file of contigs ≥ 1kb (Contigs_1k.fna)
# A fasta file of reference contigs ≥ 1kb (Ref_Contigs_1k.fna) ## Note: if no reference file, there will be no estimation of error rate (through QUAST)
## REUIQREMENTS
## bbtools (v37+)
## quast (v5+)
";
	die "\n";
}


my $bb_path=&run_cmd('which stats.sh');
if ($bb_path eq ""){die("pblm, we did not find the stats.sh path\n");}
my $spade_path=&run_cmd('which quast');
if ($spade_path eq ""){die("pblm, we did not find the quast path\n");}


my $fna_file=catfile($wdir, 'Contigs_1k.fna');
my $fna_file_ref=catfile($wdir, 'Ref_Contigs_1k.fna');

my $out_file_stat=catfile($wdir, 'Contigs_1k_stats.csv');
if (-e $out_file_stat){print "$out_file_stat already here\n";}
else{
	my %count;
	open my $fa,"<",$fna_file;
	print "\t\tReading $fna_file ... \n";
	my $seq_c="";
	while(<$fa>){
		chomp($_);
		if ($_=~/^>.*/){
			$count{"total_contigs"}++;
			$count{"total"}+=length($seq_c);
			if (length($seq_c)>=10000){
				$count{"10k"}+=length($seq_c);
				$count{"total_contigs_10k"}++;
			}
			$seq_c="";
		}
		else{
			$seq_c.=$_;
		}
	}
	close $fa;
	if (length($seq_c)>=10000){
		$count{"10k"}+=length($seq_c);
		$count{"total_contigs_10k"}++;
	}
	my $stat_file=catfile($wdir,"stats.txt");
	if (-e $stat_file){}
	else{&run_cmd("stats.sh $fna_file format=5 > $stat_file");}
	open my $txt,"<",$stat_file;
	while(<$txt>){
		chomp($_);
		my @tab=split("\t",$_);
		if ($tab[0] eq "n_contigs"){next;}
		$count{"L50"}=$tab[4];
		$count{"L90"}=$tab[6];
		$count{"max"}=$tab[7];
	}
	close $txt;
	open my $s1,">",$out_file_stat;
	print $s1 "# Total_contigs,Total_contigs_bp,Total_contigs_10k_n,Total_contigs_10k_bp,L50_bp,L90_bp,Max_contig_bp\n";
	print $s1 $count{"total_contigs"}.",".$count{"total"}.",".$count{"total_contigs_10k"}.",".$count{"10k"}.",".$count{"L50"}.",".$count{"L90"}.",".$count{"max"}."\n";
	close $s1;
}


if (!(-e $fna_file_ref)){
	die("No file $fna_file_ref, we stop here\n");
}
### If we have a reference file, we run QUAST
my $quast_out_dir=catdir($wdir,"QUAST/");
if (-d $quast_out_dir){print "QUAST was already run\n";}
else{
	&run_cmd("quast --fast -L -o $quast_out_dir -r $fna_file_ref $fna_file");
}

my $out_parsed_QUAST=catfile($wdir,"Misassemblies_QUAST.tsv");

my %store;
my %col_name;
my $base_file=catfile($quast_out_dir,"report.tsv");
if (!(-e $base_file)){die("We didn't find the expected QUAST output ($base_file), something went wrong\n");}
open my $tsv,"<",$base_file;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Assembly"){
		for (my $i=1;$i<=$#tab;$i++){
			$col_name{$i}=$tab[$i];
		}
	}
	elsif($tab[0] eq "Total length (>= 5000 bp)"){
		for (my $i=1;$i<=$#tab;$i++){
			$store{$col_name{$i}}{"size_5kb"}=$tab[$i];
		}
	}
	elsif($tab[0] eq "Total length (>= 1000 bp)"){
		for (my $i=1;$i<=$#tab;$i++){
			$store{$col_name{$i}}{"size_1kb"}=$tab[$i];
		}
	}
}
close $tsv;
%col_name=();
my $mis_file=catfile($quast_out_dir,"contigs_reports/misassemblies_report.tsv");
if (!(-e $mis_file)){die("We didn't find the expected QUAST output ($mis_file), something went wrong\n");}
open my $tsv,"<",$mis_file;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
# 				print "===$tab[0]===\n";
	if ($tab[0] eq "Assembly"){
		for (my $i=1;$i<=$#tab;$i++){
			$col_name{$i}=$tab[$i];
		}
	}
	elsif($tab[0] eq "    # c. relocations"){
		for (my $i=1;$i<=$#tab;$i++){
			$store{$col_name{$i}}{"relocation"}=$tab[$i];
			$store{$col_name{$i}}{"total_mis"}+=$tab[$i];
		}
	}
	elsif($tab[0] eq "    # c. inversions"){
		for (my $i=1;$i<=$#tab;$i++){
			$store{$col_name{$i}}{"inversion"}=$tab[$i];
			$store{$col_name{$i}}{"total_mis"}+=$tab[$i];
		}
	}
}
close $tsv;


open my $s1,">",$out_parsed_QUAST;
foreach my $asb (keys %store){
	print $s1 "$asb Total Size\t".$store{$asb}{"size_1kb"}."\n";
	print $s1 "$asb Total Miassemblies\t".$store{$asb}{"total_mis"}."\n";
	print $s1 "$asb Total Relocation\t".$store{$asb}{"relocation"}."\n";
	print $s1 "$asb Total Inversion\t".$store{$asb}{"inversion"}."\n";
	print $s1 "### NOTE : TRANSLOCATION ARE NOT CONSIDERED HERE, BECAUSE THE REFERENCE IS NOT A REAL REFERENCE BUT INSTEAD ANOTHER DRAFT ASSEMBLY\n";
}
close $s1;





sub run_cmd{
	my $cmd=$_[0];
	if ($_[1] ne "veryquiet"){print "$cmd\n";}
	my $out=`$cmd`;
	if ($_[1] ne "quiet" && $_[1] ne "veryquiet"){
		if ($_[1] eq "stderr"){print STDERR "$out\n";}
		else{print "$out\n";}
	}
	return($out);
}
