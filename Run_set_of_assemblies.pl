#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h='';
my $lib='';
my $mem=120;
my $n_cpu=8;
my $lib='';
my $num=-1;
my $fastq_filtered="";
GetOptions ('help' => \$h, 'h' => \$h, 'c=s'=> \$lib, 'p=s'=>\$fastq_filtered, 't=s'=>\$n_cpu, 'm=s'=>\$mem);
if ($h==1 || $fastq_filtered eq "" || $lib eq ""){ # If asked for help or did not set up any argument
	print "# Script to get ref contigs for a specific spid
# Arguments : 
# -c : name of the library
# -p : path to the input reads (interleaved)
# optionals:
# -m: memory limit (Gb) - default 120Gb
# -t: number of threads - default 8
";
	die "\n";
}
my $dir_assembly="Assemblies/";


my $bfc_path=&run_cmd('which bfc');
if ($bfc_path eq ""){die("pblm, we did not find the bfc path\n");}
my $spade_path=&run_cmd('which spades.py');
if ($spade_path eq ""){die("pblm, we did not find the spades path\n");}
my $bbmap_path=&run_cmd('which bbmap.sh');
if ($bbmap_path eq ""){die("pblm, we did not find the bbmap path\n");}

if (!(-e $fastq_filtered)){
	die("\n\n\nwe stop here, the file that we need ".$fastq_filtered." is not available\n");
}

## Create the two other reads set if needed
my $dir_out=$dir_assembly.$lib;
if (!(-d $dir_out)){&run_cmd("mkdir $dir_out");}
$dir_out.="/";
my $log_out=$dir_out."Processing_".$lib.".log";

# Corrected
my $bfc_reads=$dir_out.$lib."_bfc_corrected.fastq.gz";
if (!(-e $bfc_reads)){
	&run_cmd("bfc -1 -s 10g -k 21 -t $n_cpu $fastq_filtered | seqtk dropse - | gzip -c > $bfc_reads");
}

my $tadpole_reads=$dir_out.$lib."_tadpole_corrected.fastq.gz";
if (!(-e $tadpole_reads)){
	&run_cmd("tadpole.sh mode=correct ecc=t prefilter=2 out=$tadpole_reads in=$fastq_filtered -threads=$n_cpu pigz unpigz ow=t");
}

# Normalized
my $bbnorm_reads=$dir_out.$lib."_bbnormed.fastq.gz";
if (!(-e $bbnorm_reads)){
	&run_cmd("bbnorm.sh in=$fastq_filtered out=$bbnorm_reads bits=32 min=2 target=100 pigz unpigz ow=t >> $log_out 2>&1");
}

# .. then corrected
my $bbnorm_bfc_reads=$dir_out.$lib."_bbnormed_bfc_corrected.fastq.gz";
if (!(-e $bbnorm_bfc_reads)){
	&run_cmd("bfc -1 -s 10g -k 21 -t $n_cpu $bbnorm_reads | seqtk dropse - | gzip -c > $bbnorm_bfc_reads");
}

my $bbnorm_tadpole_reads=$dir_out.$lib."_bbnormed_tadpole_corrected.fastq.gz";
if (!(-e $bbnorm_tadpole_reads)){
	&run_cmd("tadpole.sh mode=correct ecc=t prefilter=2 out=$bbnorm_tadpole_reads in=$bbnorm_reads -threads=$n_cpu pigz unpigz ow=t");
}

# Clumpified
my $dedupe_reads=$dir_out.$lib."_clump_deduped.fastq.gz";
if (!(-e $dedupe_reads)){
	&run_cmd("clumpify.sh in=$fastq_filtered out=$dedupe_reads dedupe subs=0 passes=2");
}

# .. then corrected
my $dedupe_bfc_reads=$dir_out.$lib."_clump_deduped_bfc_corrected.fastq.gz";
if (!(-e $dedupe_bfc_reads)){
	&run_cmd("bfc -1 -s 10g -k 21 -t $n_cpu $dedupe_reads | seqtk dropse - | gzip -c > $dedupe_bfc_reads");
}

my $dedupe_tadpole_reads=$dir_out.$lib."_clump_deduped_tadpole_corrected.fastq.gz";
if (!(-e $dedupe_tadpole_reads)){
	&run_cmd("tadpole.sh mode=correct ecc=t prefilter=2 out=$dedupe_tadpole_reads in=$dedupe_reads -threads=$n_cpu pigz unpigz ow=t");
}

# Corrected then clumpified
my $tadpole_dedupe_reads=$dir_out.$lib."_tadpole_corrected_clump_deduped.fastq.gz";
if (!(-e $tadpole_dedupe_reads)){
	&run_cmd("clumpify.sh in=$tadpole_reads out=$tadpole_dedupe_reads dedupe subs=0 passes=2");
}

my $bfc_dedupe_reads=$dir_out.$lib."_bfc_corrected_clump_deduped.fastq.gz";
if (!(-e $bfc_dedupe_reads)){
	&run_cmd("clumpify.sh in=$bfc_reads out=$bfc_dedupe_reads dedupe subs=0 passes=2");
}

# Corrected then normalized
my $tadpole_bbnorm_reads=$dir_out.$lib."_tadpole_corrected_bbnormed.fastq.gz";
if (!(-e $tadpole_bbnorm_reads)){
	&run_cmd("bbnorm.sh in=$tadpole_reads out=$tadpole_bbnorm_reads bits=32 min=2 target=100 pigz unpigz ow=t >> $log_out 2>&1");
}

my $bfc_bbnorm_reads=$dir_out.$lib."_bfc_corrected_bbnormed.fastq.gz";
if (!(-e $bfc_bbnorm_reads)){
	&run_cmd("bbnorm.sh in=$bfc_reads out=$bfc_bbnorm_reads bits=32 min=2 target=100 pigz unpigz ow=t >> $log_out 2>&1");
}


## Now doing the assemblies

## Normalized and corrected reads
my $dir_spades=$dir_out."SPAdes_meta_normalized_corrected";
my $asbly_meta_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $bbnorm_bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_normalized_corrected";
my $asbly_sc_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $bbnorm_bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}

## Normalized and corrected reads - Tadpole
my $dir_spades=$dir_out."SPAdes_meta_normalized_tadpole_corrected";
my $asbly_meta_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $bbnorm_tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_normalized_tadpole_corrected";
my $asbly_sc_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $bbnorm_tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}



## Dedupe and corrected reads
my $dir_spades=$dir_out."SPAdes_meta_dedupe_corrected";
my $asbly_meta_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $dedupe_bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_dedupe_corrected";
my $asbly_sc_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $dedupe_bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}

## Dedupe and corrected reads - Tadpole
my $dir_spades=$dir_out."SPAdes_meta_dedupe_tadpole_corrected";
my $asbly_meta_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $dedupe_tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_dedupe_tadpole_corrected";
my $asbly_sc_normcorr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_normcorr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $dedupe_tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}



## Corrected reads
my $dir_spades=$dir_out."SPAdes_meta_corrected";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_corrected";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $bfc_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


## Tadpole Corrected reads
my $dir_spades=$dir_out."SPAdes_meta_tadpole_corrected";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_tadpole_corrected";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}



### Corrected then clumpified
# meta_corrected_dedupe
# meta_tadpole_corrected_dedupe
# sc_corrected_dedupe
# sc_tadpole_corrected_dedupe

## Tadpole corrected dedupe reads
my $dir_spades=$dir_out."SPAdes_meta_tadpole_corrected_dedupe";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_dedupe_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_tadpole_corrected_dedupe";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_dedupe_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


## bfc corrected dedupe reads
my $dir_spades=$dir_out."SPAdes_meta_corrected_dedupe";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $bfc_dedupe_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_corrected_dedupe";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $bfc_dedupe_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}



### Corrected then normalized
# meta_corrected_normalized
# meta_tadpole_corrected_normalized
# sc_corrected_normalized
# sc_tadpole_corrected_normalized


## Tadpole corrected normalized reads
my $dir_spades=$dir_out."SPAdes_meta_tadpole_corrected_normalized";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_bbnorm_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_tadpole_corrected_normalized";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $tadpole_bbnorm_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


## bfc corrected normalized reads
my $dir_spades=$dir_out."SPAdes_meta_corrected_normalized";
my $asbly_meta_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_meta_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --meta --only-assembler -k 21,33,55,77,99,127 --12 $bfc_bbnorm_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
}


my $dir_spades=$dir_out."SPAdes_sc_corrected_normalized";
my $asbly_sc_corr=$dir_spades."/contigs.fasta";
if (!(-e $asbly_sc_corr)){
	if (-d $dir_spades){
		my $last_step=&get_last_step_SPAdes($dir_spades);
		&run_cmd("spades.py --restart-from $last_step -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	else{
		&run_cmd("spades.py --sc --only-assembler -k 21,33,55,77,99,127 --12 $bfc_bbnorm_reads -t $n_cpu -m $mem -o $dir_spades >> $log_out 2>&1");
	}
	&clean_dir($dir_spades);
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


sub clean_dir{
	my $out_dir=$_[0];
	my $contig_result=$out_dir."/contigs.fasta";
	## Clean the output directory
	if (-e $contig_result){
		my $contig_result_filtered=$out_dir."/contigs_1k.fasta";
		&filter_on_size($contig_result,$contig_result_filtered,1000);
		$contig_result_filtered=$out_dir."/contigs_10k.fasta";
		&filter_on_size($contig_result,$contig_result_filtered,10000);
		&run_cmd("rm -rf $out_dir/K* $out_dir/before* $out_dir/first* $out_dir/mismatch* $out_dir/tmp $out_dir/misc $out_dir/*.fastg $out_dir/*.paths $out_dir/*.info $out_dir/*.yaml $out_dir/*.txt $out_dir/split_input $out_dir/assembly_graph_with_scaffolds.gfa");
	}
	else{
		print "No contig generated in $out_dir \n";
	}
}
