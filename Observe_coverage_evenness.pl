#!/usr/bin/env perl
use strict;
use autodie;
use File::Spec::Functions;
use Getopt::Long;
my $h='';
my $wdir='';
GetOptions ('help' => \$h, 'h' => \$h, 'd=s'=>\$wdir);
if ($h==1 || $wdir eq ""){ # If asked for help or did not set up any argument
	print "# Script to parse read mapping 
# Arguments : 
# -d: working directory
# Must include the following files:
# A sorted bam file (Reads_to_ref-1k_sorted.bam) with corresponding index (Reads_to_ref-1k_sorted.bam.bai), obtained by mapping reads to contigs ≥ 1kb, then sort and index using samtools
# A fasta file of contigs ≥ 1kb (Contigs_1k.fna)
";
	die "\n";
}
my $fna_file=catfile($wdir, 'Contigs_1k.fna');
my $bam_file=catfile($wdir, 'Reads_to_ref-1k_sorted.bam');


### Test that samtools is available
my $test=`which samtools`;
if ($test eq ""){die("pblm, we did not find the samtools path\n");}



### First, get the contig length, average depth of coverage, median depth of coverage, and associated standard deviation and coefficient of variation (for each contig)
my $out_file=catfile($wdir, 'Samtools_coverage_with_variation.csv');
if (-e $out_file){print "$out_file already here\n";}
else{
	print "Getting the basic metrics on contig coverage\n";
	my %seq_len_selected=&get_seq_len($bam_file);
	open my $s1,">",$out_file;
	print $s1 "Contig,Length,Average,Median,Stdev,Coeff\n";
	foreach my $contig (sort keys %seq_len_selected){
		my $metric=&get_coverage_metrics($bam_file,$contig);
		if ($metric==0){}
		else{
			print $s1 $contig.",".$seq_len_selected{$contig}.",".$metric."\n";
		}
	}
	close $s1;
}

### Then, look for duplicated reads (based on mapping) and get the corresponding insert size and GC% of these duplicated regions
my $dup_file=catfile($wdir,"Coverage-vs-insert_gc.csv");
if (-e $dup_file){print "$dup_file already processed\n";}
else{
	print "Getting the metrics on specific regions with (seemingly) duplicated reads\n";
	my $th=&get_len_th($bam_file);
	my $count=0;
	my %store_windows;
	my %store_contigs;
	my $c_c="";
	my $seq_c="";
	print "Reading $fna_file\n";
	open my $fa,"<",$fna_file;
	while(<$fa>){
		chomp($_);
		if ($_=~/^>(\S+)/){
			my $id_temp=$1;
			if ($seq_c ne ""){
				if (length($seq_c)>$th){
					$store_contigs{$c_c}=$seq_c;
				}
			}
			$c_c=$id_temp;
			$seq_c="";
		}
		else{
			$seq_c.=$_;
		}
	}
	close $fa;
	if ($seq_c ne ""){
		if (length($seq_c)>$th){
			$store_contigs{$c_c}=$seq_c;
		}
	}
	my %store_tmp;
	$c_c="";
	my $max=0;
	open my $s1,">",$dup_file;
	print $s1 "Contig,Contig_length,Position,Fragment_length,GC,Duplication,Coverage,Relative_coverage\n";
	foreach my $contig (sort {length($store_contigs{$b}) <=> length($store_contigs{$a})} keys %store_contigs){
		print "## $contig \t".length($store_contigs{$contig})."\n";
		my %depth;
		my @t_file=split("\n",&run_cmd("samtools view $bam_file $contig","quiet"));
		%store_tmp=();
		$max=&load_coverage($contig,$bam_file,\%store_tmp);
		my $c_code="";
		my $dup=0;
		foreach my $line (@t_file){
			chomp($line);
			my @tab=split("\t",$line);
			if ($tab[8]<=0){next;} 
			## If 0 we don't care, no length estimated
			## If < 0 that's the second time we see this fragment, and we don't want to count it twice
			my $code=$tab[2]."|".$tab[3]."|".$tab[8];
			if ($code ne $c_code){
				if ($c_code ne ""){
					my @t=split(/\|/,$c_code);
					if ($t[2]>0 && $t[2]<1000){ ## We look at "reasonable" insert sizes here
						## We look at the coverage and GC content of this window
						my $cov=&get_coverage($t[1],$t[2],\%store_tmp);
						my $gc=&get_gc_fragment($store_contigs{$contig},$t[1],$t[2]);
						print $s1 $contig.",".length($store_contigs{$contig}).",".$t[1].",".$t[2].",".$gc.",".$dup.",".$cov.",".$cov/$max."\n";
					}
				}
				$c_code=$code;
				$dup=1;
			}
			else{
				$dup++;
			}
		}
	}
	close $s1;
}





sub get_seq_len{
	my $in=$_[0];
	my %seq_len;
	my %total;
	print "\tSelecting contigs\n";
	open my $txt,"samtools view -H $in |";
	while(<$txt>){
		chomp($_);
		if ($_=~/^\@SQ/){
			if ($_=~/SN:(\S+)\tLN:(\d+)/){
				$seq_len{$1}=$2;
				if ($2>=10000){$total{"10k"}+=$2}
				if ($2>=2000){$total{"2k"}+=$2}
			}
			else{
				print "?!?!?! PBLM WITH $_\n";
				<STDIN>;
			}
		}
	}
	close $txt;
	print "\t".$total{"10k"}." bp in 10k+, ".$total{"2k"}." bp in 2k+\n";
	my %filtered_len;
	if ($total{"10k"}>=50000){
		# We take 10k contigs only
		print "\t\t we select 10k+\n";
		foreach my $seq (sort keys %seq_len){
			if ($seq_len{$seq}>=10000){
				$filtered_len{$seq}=$seq_len{$seq};
# 				print "\t\t\t $seq comes with us\n";
			}
		}
	}
	else {
		# Or only 2k if not enough in 10k+
		print "\t\t we select 2k+\n";
		foreach my $seq (sort keys %seq_len){
			if ($seq_len{$seq}>=2000){
				$filtered_len{$seq}=$seq_len{$seq};
# 				print "\t\t\t $seq comes with us\n";
			}
		}
	}
	return %filtered_len;
}

sub get_len_th {
	my $in=$_[0];
	my %total;
	open my $txt,"samtools view -H $in |";
	while(<$txt>){
		chomp($_);
		if ($_=~/^\@SQ/){
			if ($_=~/SN:(\S+)\tLN:(\d+)/){
				if ($2>=10000){$total{"10k"}+=$2}
				if ($2>=2000){$total{"2k"}+=$2}
			}
			else{
				print "?!?!?! PBLM WITH $_\n";
				<STDIN>;
			}
		}
	}
	close $txt;
	my $th=-1;
	if ($total{"10k"}>=50000){$th=10000;}
	else{$th=2000;}
	return($th);
}


sub get_coverage_metrics{
	my $in=$_[0];
	my $contig=$_[1];
	print "\tParsing contig $contig\n";
	open my $tsv,"samtools depth -a -r $contig $in |";
	my @tab_cover;
	while(<$tsv>){
		chomp($_);
		my @tab=split("\t",$_);
		if ($tab[0] ne $contig){
			print "?!?! pblm with $contig - $tab[0] ?\n";
		}
		else{
			push(@tab_cover,$tab[2]);
		}
	}
	close $tsv;
	if ($#tab_cover<=0){
		print "$contig is totally empty ??\n";
		return 0;
	}
	else{
		my $avg=&average(\@tab_cover);
		my $median=&median(\@tab_cover);
		my $stdev=&stdev(\@tab_cover);
		my $coeff=$stdev/$avg;
		return $avg.",".$median.",".$stdev.",".$coeff;
	}
}



sub load_coverage{
	my $c_c=$_[0];
	my $mapping_file=$_[1];
	my $hash=$_[2];
	my $max=0;
	my @t=split("\n",&run_cmd("samtools depth -r $c_c $mapping_file","quiet"));
	foreach my $line (@t){
		chomp($line);
		my @tab=split("\t",$line);
		$$hash{$tab[1]}=$tab[2];
		if ($tab[2]>$max){$max=$tab[2];}
	}
	return($max);
}

sub get_coverage{
	my $start=$_[0];
	my $length=$_[1];
	my $hash=$_[2];
	my $total=0;
	for (my $i=$start;$i<$start+$length;$i++){$total+=$$hash{$i};}
	return $total/$length;
}

sub get_gc_fragment{
	my $contig=$_[0];
	my $start=$_[1];
	my $length=$_[2];
	my $seq_c=substr($contig,$start,$length);
	$seq_c=~tr/atcgn/ATCGN/;
	$seq_c=~s/N//g;
	$seq_c=~tr/C/G/;
	my $gc_pcent=($seq_c =~ tr/G//)/$length;
	## 
# 	my $seq=substr($contig,$start,$length);
# 	print $seq."\n".$seq_c."\n".$length."\t".$gc_pcent."\n";
# 	<STDIN>;
	##
	return $gc_pcent;
}


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


sub average{ ## Stolen from https://edwards.sdsu.edu/research/calculating-the-average-and-standard-deviation/
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub median{ ## Stolen from https://edwards.sdsu.edu/research/calculating-the-average-and-standard-deviation/
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $med="NA";
        @$data=sort {$a <=> $b} (@$data);
	if (scalar(@$data) % 2 ==0 ){
		$med=(@{$data}[scalar(@$data)/2]+@{$data}[scalar(@$data)/2-1])/2;
	}
	else{
		$med=@{$data}[int(scalar(@$data)/2)];
	}
        return $med;
}



sub stdev{ ## Stolen from https://edwards.sdsu.edu/research/calculating-the-average-and-standard-deviation/
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
