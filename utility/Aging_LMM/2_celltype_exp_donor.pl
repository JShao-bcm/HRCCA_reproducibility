#!/usr/bin/perl -w
use strict;

my $sample_list= "snRNA_donor_majorclass_num_atlas"; #$ARGV[0]; #snRNA_donor_majorclass_num_atlas
my $ct_list = "ct_list_majorclass"; $ARGV[1]; #ct_list_majorclass
my $indir= "genexp_donor_raw_atlas"; #$ARGV[2]; #genexp_donor_raw_atlas
my $outdir= "2_genexp_donor_count_matrix", #$ARGV[3]; #genexp_donor_raw_atlas_2

my %sample;
open(INPUT,$sample_list);
my $header = <INPUT>;
chomp $header;
my @sample_list1 = split(/\,/,$header);
my @sample_list = @sample_list1; #@sample_list1[0..30];
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\,/,$line);
for(my $i=1; $i<=$#info; $i++){

if(($sample_list[$i] =~ /mix/) || ($sample_list[$i] =~ /\,/)){
next;
}
if($info[$i]>=20){

$sample{$info[0]}->{$sample_list[$i]}=$info[$i];
}
}
}

#my @ct_name = ("AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte","Microglia"); #,"RPE");
#my @ct_name = ("Oligodendrocyte","Astrocyte","Fibroblast","Rod","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","Melanocyte","HC","Cone","RGC","Schwann_cell");
my %ct_name;
open(INPUTct, $ct_list);
while(my $line = <INPUTct>){
	chomp $line;
	$ct_name{$line}=1;
	#	push(@ct_name,$line);
}
my %people;
for(my $s=1;$s<=$#sample_list;$s++){


my $sam = $sample_list[$s];
if($sam =~ /mix/){
next;
}

#my $file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas/$sam".".txt.gz";
my $file = $indir."/$sam".".txt.gz";

open(INPUT,"gunzip -c $file|");
my $line = <INPUT>;
chomp $line;
my @ct = split(/\s+/,$line);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i=1; $i<=$#info;$i++){
$people{$ct[$i]}->{$info[0]}->{$sam} = $info[$i];
}
}
}
#`mkdir -p /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_celltype_atlas_sample`;
`mkdir -p $outdir`;

for my $key (keys %people){
if(($key eq "unassigned")|| (!(exists $ct_name{$key} ))){
next;
}
#my $output = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_celltype_atlas_sample/exp_"."$key";
my $key1 = $key;
$key1 =~ s/\//_/g;
#my $output = $outdir."/exp_"."$key";
my $output = $outdir."/exp_"."$key1";

open(OUTPUT,">$output");
my @sample_ct = sort (keys %{$sample{$key}});
my $line=join("\n", @sample_ct);
my $header = join("\t",@sample_ct);
print OUTPUT "$header\n";
for my $gene (sort keys %{$people{$key}}){
my $zero_gene=0;
    for(my $i=0; $i<=$#sample_ct; $i++){
        if(!(defined $people{$key}->{$gene}->{$sample_ct[$i]})){
                  $people{$key}->{$gene}->{$sample_ct[$i]}=0;
        $zero_gene++;                  
        }
    }

print OUTPUT "$gene\t";

    for(my $i=0; $i<=$#sample_ct-1; $i++){
        if(!(defined $people{$key}->{$gene}->{$sample_ct[$i]})){
                  $people{$key}->{$gene}->{$sample_ct[$i]}=0;
        }
     print OUTPUT "$people{$key}->{$gene}->{$sample_ct[$i]}\t";
    }


  my $tmp = $#sample_ct;
  if(!(defined $sample_ct[$#sample_ct])){
  print "$key\t$tmp\n";
  exit;
  }
  if(!(defined $people{$key}->{$gene}->{$sample_ct[$#sample_ct]})){
                  $people{$key}->{$gene}->{$sample_ct[$#sample_ct]}=0;
   }
     print OUTPUT "$people{$key}->{$gene}->{$sample_ct[$#sample_ct]}\n";
}
}
