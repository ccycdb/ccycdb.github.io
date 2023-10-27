#!/usr/bin/perl
use strict;
use List::Util qw(shuffle sum);
use Getopt::Long;

my ($situation, $workdir, $method, $filetype, $seqtype, $identity, $evalue,$tpm, $norm, $randomsample, $thread, $outdir,$help);
GetOptions(
  "situation=s" => \$situation,
  "wd=s" => \$workdir, 
  "m=s" => \$method,
  "f=s" => \$filetype, 
  "s=s" => \$seqtype,
  "id=s" => \$identity,
  "e=s" => \$evalue,
  "tpm=s" => \$tpm,
  "norm=s" => \$norm,
  "rs=s" =>  \$randomsample,
  "thread=s" => \$thread,
  "od=s" => \$outdir,
  "h=s" => \$help
);
if (   !defined $situation
	|| !defined $workdir
	|| !defined $method
	|| !defined $filetype
	|| !defined $seqtype
	|| !defined $seqtype
	|| $situation   !~ /^read-based|assembly-based|tabular$/
	|| $method   !~ /^diamond|usearch|blast$/
	|| $filetype !~ /^fastq|fastq.gz|fasta|fasta.gz|fq|fq.gz|fa|fa.gz$/
	|| $seqtype !~ /^nucl|prot$/
	|| (defined $norm && $norm !~ /0|1/)
	|| (defined $randomsample && $norm == 0)
	)
{ 
  print "Incorrect parameters!\n";
  &PrintHelp();
  die;
}

#software
my $blast = "$workdir/data/Softwares/ncbi-blast-2.13.0+/bin";
my $seqkit = "$workdir/data/Softwares/seqkit";
my $diamond = "$workdir/data/Softwares/diamond";
my $usearch = "$workdir/data/Softwares/usearch11";
my $csvtk = "$workdir/data/Softwares/csvtk";

#option
if (-e $outdir){
	print "Attation!!  \"$outdir\" already exists, be cautious as $outdir may be overwritten!\n";
};

if (defined $outdir && !-e $outdir){
	system ("mkdir $outdir");
}elsif(!defined $outdir && !-e $outdir){
	system ("mkdir $workdir/output")
}

my $dthread;
if (defined $thread){
	$dthread =$thread;
}else{
	$dthread=4;
}

my $didentity;
if (defined $identity){
	$didentity =$identity;
}else{
	$didentity=30;
}

my $devalue;
if (defined $identity){
	$devalue =$evalue;
}else{
	$devalue=1e-5;
}

my $diamond_parameters = "-k 1 -e $devalue --id $didentity -p $dthread --sensitive";
my $uidentity = $didentity / 10; 
my $usearch_parameters = "-id $uidentity";
my $blast_parameters = "-evalue $devalue -num_threads $dthread -max_target_seqs 1 -outfmt 6";

system ("mkdir $outdir/ORF2GENE") if (!-e "$outdir/ORF2GENE" && $situation eq "assembly-based"); #s2
system ("mkdir $outdir/SEQ2GENE") if (!-e "$outdir/SEQ2GENE" && $situation ne "assembly-based");
system ("mkdir $outdir/tmp") if (!-e "$outdir/tmp" && $situation eq "assembly-based" && $tpm == 1);


#### Referseq2Gene--id2genemap ####
my %id2gene;
my (@genes,%ha);
open(FILE, "$workdir/data/id2gene.CCycDB.map") || die "#1 cannot open id2gene.map\n";
print "Reading id2genemap ...\n";
while (<FILE>) {
  chomp;
  my @items = split("\t", $_);
  $id2gene{$items[0]} = $items[1];
  push(@genes,$items[1]);
}
close FILE;
print "Done!\n";
my @genelist=grep{++$ha{$_}<2}@genes;

### Annotation ###
print "Reading annotation ...\n";
my %annotation;
open(ANNO, "$workdir/data/annotations.txt") || die "# cannot open annotation.txt";
while (<ANNO>) {
  chomp;
  my @items = split("\t", $_);
  $annotation{$items[0]}{"Annotation"}     = $items[1];
  $annotation{$items[0]}{"KOID"}           = $items[2];
  $annotation{$items[0]}{"Brenda"}         = $items[3];
  $annotation{$items[0]}{"CAZY"}            = $items[4];
  $annotation{$items[0]}{"CAZY_group"}            = $items[5];
  $annotation{$items[0]}{"COG"} = $items[6];
=pod
  my @path = split(",", $items[1]);
  foreach my $path (@path) {
	$path=~s/"//g;
    $pathway{$path}{$items[0]} = 1;
  }
=cut
}
close ANNO;


#### Calculate the sample info ####
print "Calculating the number of sequences ... \n";
system ("rm $workdir/sampleinfo.txt") if -e "$workdir/sampleinfo.txt";
system ("$seqkit -j $dthread stat $workdir/*$filetype >>$workdir/sampleinfo.txt");


my (%size,@sizes);
open( FILE, "$workdir/sampleinfo.txt" ) || die "#1\n";
while (<FILE>) {
  chomp;
  $_=~ s/,//g;
  my @items = split( " ", $_);
  $items[0] =~ /$workdir\/(.*?).$filetype/;
  my $s = $1;
  $size{$s} = $items[3] if $items[3] =~ /^\d/;
  push(@sizes,$items[3]) if $items[3] =~ /^\d/;
}
close FILE;
print "Done!\n";


#### Alignment #####
print "Mapping ... \t";
my @files;
if ($situation ne "tabular" && $method eq "diamond") {
  @files = glob("$workdir/*$filetype");
  my $diamond_db = "$diamond makedb --in $workdir/data/CCycDB.database --db $workdir/data/CCycDB.database" if (!-e "$workdir/data/CCycDB.database.dmnd");
  system("$diamond_db") if (!-e "$workdir/data/CCycDB.database.dmnd");
  
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/$filetype/diamond/;
	print "Mapping:$file\n";
	if(!-e $out){
	print "$diamond blastx $diamond_parameters -d $workdir/data/CCycDB.database -q $file -o $out\n";
	system("$diamond blastx $diamond_parameters -d $workdir/data/CCycDB.database -q $file -o $out") if $seqtype eq "nucl";
	system("$diamond blastp $diamond_parameters -d $workdir/data/CCycDB.database -q $file -o $out") if $seqtype eq "prot";
	}
  }
  
}elsif ($situation ne "tabular" && $method eq "usearch") {
  die "Please specify the location of usearch!" if !-e $usearch;
	if ($filetype =~ /gz/) {
	  die "Only fastq and fasta files are supported by usearch!\n";
	}
	@files = glob("$workdir/*$filetype");
	foreach my $file (@files) {
	my $out = $file;
	$out =~ s/$filetype/usearch/;
	if(!-e $out){
	system("$usearch -usearch_global $file -db $workdir/data/CCycDB.database $usearch_parameters -blast6out $out");
	}
  }
  
}elsif ($situation ne "tabular" && $method eq "blast") {
  die "Please specify the location of blast and/or formatdb!" if !-e $blast;
  if ($filetype ne /fasta|fa/) {
  #if ($filetype =~ /gz|fastq/) {
    die "Only fasta files(filetype is \"fasta|fa\") are supported by blast program!";
  }
  @files = glob("$workdir/*$filetype");
  system("$blast/makeblastdb -in $workdir/data/CCycDB.database -dbtype prot -out $workdir/data/CCycDB.blast");
  
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/$filetype/blast/;
	if(!-e $out){
    system("$blast/blastp -db $workdir/data/CCycDB.blast -query $file -out $out $blast_parameters") if $seqtype eq "prot";
    system("$blast/blastx -db $workdir/data/CCycDB.blast -query $file -out $out $blast_parameters") if $seqtype eq "nucl";
	}
  }
}
print "done!\n";


my (%identity,%abundance,%samples);
my @sfiles = glob("$workdir/*diamond") if $method eq "diamond";
@sfiles = glob("$workdir/*usearch") if $method eq "usearch";
@sfiles = glob("$workdir/*blast")   if $method eq "blast";
die "No diamond/usearch/blast files were detected!\n" if $#sfiles == -1;

if($situation ne "assembly-based"){
  foreach my $file (@sfiles) {
	my (%diamond,%hit);
	$file =~ /$workdir\/(.*?)\.$method/;
	my $sample = $1;
	$samples{$sample} = 1;

	open(DIA, "$file") || die "#cannot open $file\n";
	open(SEQ2G, ">$outdir/SEQ2GENE/$sample.SEQ2G.txt") || die "#cannot open $outdir/SEQ2GENE/$sample.SEQ2GENE.txt\n";
	while (<DIA>) {
	  chomp;
	  my @items = split("\t", $_);
	  my $gene  = $id2gene{$items[1]};
	  $diamond{$items[0]} = $items[1];
	  #if (!$hit{$items[0]}) {
		$abundance{$sample}{$gene}++         if $gene;
		push(@{$identity{$gene}}, $items[2]) if $gene;
		#$hit{$items[0]} = 1;
	  #}
	  print SEQ2G "$items[0]\t$id2gene{$diamond{$items[0]}}\n" if $id2gene{$diamond{$items[0]}};
	}
	close DIA;
	close SEQ2G;
  }
}elsif ($situation eq "assembly-based"){
  foreach my $file (@sfiles) {
  my (%diamond,%hit);
  $file =~ /$workdir\/(.*?)\.$method/;
  my $sample = $1;
  $samples{$sample} = 1;

  open( DIA, "$file" ) || die "#cannot open ORF2GENE.txt\n\n";
  open(ORF2G, ">$outdir/ORF2GENE/$sample.ORF2GENE.txt") || die "#cannot open $outdir/ORF2GENE/$sample.ORF2GENE.txt\n";
  while (<DIA>) {
	chomp;
	my @items = split( "\t", $_ );
	my $gene  = $id2gene{$items[1]};
	$diamond{$items[0]} = $items[1];
	if (!$hit{$items[0]}) {
	$abundance{$sample}{$gene}++         if $gene;
	push(@{$identity{$gene}}, $items[2]) if $gene;
	$hit{$items[0]} = 1;
	}
	print ORF2G "$items[0]\t$id2gene{$diamond{$items[0]}}\n" if $id2gene{$diamond{$items[0]}};
  }
  close DIA;
  close ORF2G;
  
  if ($tpm == 1){
	  print "Please make sure the file of $workdir/$sample.tpm exists...\n";
	my (@tgene,%ta);
	open( VTPM, "$workdir/$sample.tpm" ) || die "#$workdir/$sample.tpm needed if -tpm 1\n";
	
	while (<VTPM>) {
	chomp;
	my @items = split( "\t", $_ );
	  if ($diamond{$items[0]}&&$id2gene{$diamond{$items[0]}}){
		  $abundance{$id2gene{$diamond{$items[0]}}} += $items[3];
		  push(@tgene,$id2gene{$diamond{$items[0]}});
	  }
	}
	
	my @uniqtg=grep{++$ta{$_}<2}@tgene;
	open( TPM, ">$outdir/tmp/$sample\.tpm.txt" ) || die "#1\n";
	print TPM "Gene\t$sample\n";
	foreach my $key (sort @uniqtg){
	  print TPM "$key\t$abundance{$key}\n" if $abundance{$key};
	  print TPM "$key\t0\n" if !$abundance{$key};
	}
	close TPM;
	close VTPM;
  }
 }
 system ("$csvtk join -t -f 1 -O --na 0 $outdir/tmp/*txt.tmp >$outdir/ORF2GENE.tpm") if -e "$outdir/tmp";
 #system ("rm $outdir/tmp/*tpm.tmp");
}

@sizes = sort { $a <=> $b } @sizes;
my $rs= $sizes[0] if (!defined $randomsample && $norm == 1);
$rs = $randomsample if (defined $randomsample && $norm == 1);

my @samplelist = sort keys %samples;

#### Not randomsampling ####
if ($norm == 0){
 print "not random sampling ...\t";
 open(FOUT, ">$outdir/FunProfile_$situation\_$method\_norandom.txt") || die "#5 cannot write profile.txt\n";
 print FOUT "#Not random sampling\n";
 print FOUT "Gene\tMean Identity\tAnnotation\tKO\tBrenda\tCAZY\tCAZY_group\tCOG\t",join("\t", @samplelist), "\n";
 foreach my $gkey (sort keys %identity){
  print FOUT "$gkey";
  if ($identity{$gkey}){
	  my $identity = sprintf("%.2f",sum(@{$identity{$gkey}}) / scalar(@{$identity{$gkey}}));
	  print FOUT "\t", sum(@{$identity{$gkey}}) / scalar(@{$identity{$gkey}});
  }else{
	  print FOUT "\t0";
  }
  #annotation
  print FOUT "\t$annotation{$gkey}{'Annotation'}\t$annotation{$gkey}{'KOID'}\t$annotation{$gkey}{'Brenda'}\t$annotation{$gkey}{'CAZY'}\t$annotation{$gkey}{'CAZY_group'}\t$annotation{$gkey}{'COG'}";
  foreach my $skey (@samplelist){
	print FOUT "\t$abundance{$skey}{$gkey}" if $abundance{$skey}{$gkey};
	print FOUT "\t0" if !$abundance{$skey}{$gkey};
  }
  print FOUT "\n";
 }
 close FOUT;
 
 #### randomsampling ####
}elsif(!defined $norm or $norm == 1){
 print "random sampling....\t";
 my %abundance_rs = &RandomSampling(\%abundance, \%size, $rs);
 open(FOUT, ">$outdir\/FunProfile_$situation\_$method\_randomTO$rs.txt") || die "#4 cannot write outfile\n";
 print FOUT "#Random sampling: $rs\n";
 print FOUT "Gene\tMean Identity\t",join("\t", @samplelist), "\n";
 foreach my $gkey (sort keys %abundance_rs) {
  print FOUT "$gkey\t";
  my $identity = sprintf("%.2f",sum(@{$identity{$gkey}}) / scalar(@{$identity{$gkey}}));
  print FOUT sum(@{$identity{$gkey}}) / scalar(@{$identity{$gkey}});
  print FOUT "\t$annotation{$gkey}{'Annotation'}\t$annotation{$gkey}{'KOID'}\t$annotation{$gkey}{'Brenda'}\t$annotation{$gkey}{'CAZY'}\t$annotation{$gkey}{'CAZY_group'}\t$annotation{$gkey}{'COG'}";
  foreach my $skey (@samplelist) {
    print FOUT "\t$abundance_rs{$gkey}{$skey}" if $abundance_rs{$gkey}{$skey};
    print FOUT "\t0" if !$abundance_rs{$gkey}{$skey};
  }
  print FOUT "\n";
 }
 close FOUT;
}
print "done!\n";

#### Reading categories... ####
print "Reading categories...\t";
my %pathway;
open(SUBIII, "$workdir/data/subcategoryII.txt") || die "# cannot open subcategoryIII.txt";
while (<SUBIII>) {
  chomp;
  my @items = split("\t", $_);
  push (@{$pathway{$items[1]}{$items[2]}{$items[3]}},$items[0]);
}
close SUBIII;
print "done.\n";

#### calculating the abundance of pathways ####
print "calculating the abundance of pathways ...\n";
open(PATH, ">$outdir/Fun_gene2pathway.txt") || die "#Line356 cannot write Fun_gene2pathway.txt\n";
print PATH "#The number of genes in pathways\n";
print PATH "Category\tSubCategoryI\tSubCategoryII\tGenes\t",join("\t", @samplelist), "\n";

open(PATH2A, ">$outdir/Fun_abundance2pathway.txt") || die "#Line360 cannot write Fun_abundance2pathway.txt\n";
print PATH2A "#The abundance of pathway\n";
print PATH2A "Category\tSubCategoryI\tSubCategoryII\t",join("\t", @samplelist), "\n";
foreach my $path1 (sort keys %pathway){
  #print PATH "$gkey";
  foreach my $path2 (sort keys %{$pathway{$path1}}){
	foreach my $path3 (sort keys %{$pathway{$path1}{$path2}}){
	  
	  print PATH "$path1\t$path2\t$path3\t", scalar(keys @{$pathway{$path1}{$path2}{$path3}});
	  print PATH2A "$path1\t$path2\t$path3";
	  my @genes = @{$pathway{$path1}{$path2}{$path3}};
	  
	  foreach my $skey (@samplelist) {
		my ($genesum, %abundance2path);
		foreach my $gene (@genes){
		  if($gene && $abundance{$skey}{$gene}){
			$genesum++;
			$abundance2path{$skey} += $abundance{$skey}{$gene}
		  }
		}
		print PATH "\t$genesum";
		print PATH "\t0" if !$genesum;
		print PATH2A "\t$abundance2path{$skey}" if $abundance2path{$skey};
		print PATH2A "\t0" if !$abundance2path{$skey};
	  }
	  print PATH "\n";
	  print PATH2A "\n";
	  
	}
  }
}
close PATH;
close PATH2A;
print "Done.\n";


sub RandomSampling() {
  my ($abundance, $size, $randomsample) = @_;
  my %abundance = %$abundance;
  my %size      = %$size;
  my %sum;
  foreach my $sample (keys %abundance) {
    foreach my $gene (keys %{$abundance{$sample}}) {
      $sum{$sample} += $abundance{$sample}{$gene};
    }
  }
  my %abundance_rs;
  foreach my $sample (keys %size) {
    my @array = shuffle(1 .. $size{$sample});
    @array = @array[0 .. $randomsample - 1];
    @array = grep { $_ <= $sum{$sample} } @array;
    my $i = 1;
    foreach my $gene (keys %{$abundance{$sample}}) {
      my @tmp
        = grep { $_ >= $i && $_ <= ($abundance{$sample}{$gene} + $i - 1) }
        @array;
      $abundance_rs{$gene}{$sample} = @tmp;
      $i += $abundance{$sample}{$gene};
    }
  }
  return %abundance_rs;
}

sub PrintHelp() {
  print "Manunal:\n";
  print
    "perl CCycdb.pl [-situation situation] [-wd work_directory] [-m diamond|usearch|blast] [-f filetype] [-s seqtype] [-tpm] [-norm] [-rs random_sampling_size] [-thread] [-od out_directory] [-h help]\n\n";
  print "E.g. perl CCycdb.pl -situation read-based -wd ./ -m diamond -f fasta -s nucl -norm 0 -thread 20 -od ./\n\n";
  
  print "E.g. perl CCycdb.pl -situation assembly-based -wd ./ -m diamond -f fasta -s nucl -tpm 1 -thread 60 -od ./\n\n";
  
  print "-situation\tthe situation for input files(read-based|assembly-based|tabular)\n";
  print "-wd\twork directory or current directory. <data> directory from github and your fasta/fastq files should be included in this directory\n";
  print "-od\toutfile directory\n";
  print "-m\tdatabase searching program you plan to use(diamond|usearch|blast)\n
  \t(diamond: -k 1 -e 1e-5)\n
  \t(usearch: -id 0.3)\n
  \t(blast: -evalue 1e-5)\n";
  print "-f\tspecify the extensions of your sequence files(fastq, fastq.gz, fasta, fasta.gz, fq, fq.gz, fa, fa.gz) or (faa, fna).\n\tPlease make sure that the filetype is correct according to the tool selected by -m.\n\tE.g. if -m usearch , filetype is support the (fastq|fasta);\n\tif -m blast, filetype is support fasta\n";
  print "-s\tsequence type (nucl|prot)\n";
  print "-tpm\t(0|1)(default: 0). \"1\" need \$sample.tpm exist in the work directory\n";
  print "-norm\t(0|1).0: don`t need random sampling; 1: need random sampling\n
		\tNote: [-situation assembly-based] is a prerequisite for this parameter.";
  print "-rs\tthe number of sequences for random subsampling. (default: the lowest number of sequences).\n\tNote: [-norm 1] is a prerequisite for this parameter.\n";
  print "-thread\tnumber of CPUs (default: 4)\n";
  print "-h\thelp documentation\n";
}
