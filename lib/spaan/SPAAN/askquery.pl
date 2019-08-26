#!/usr/bin/perl
use File::Copy;
use File::Path;
use File::Basename;
use Cwd;

$SPAAN_dir=dirname(__FILE__);

$tmp_dir="$SPAAN_dir/".int(rand(1000000));
$cur_dir=cwd();
$query_file = $ARGV[0];
$output_file = $ARGV[1];

print ">moving files...\n";

mkdir("$tmp_dir");
copy("$query_file", "$tmp_dir/input.fasta");
chdir("$tmp_dir");

print ">standardization going on...\n";
system("$SPAAN_dir/standard.o input.fasta");
move("input.std", "query.ext");

print(">filteration going on...\n");
system("$SPAAN_dir/filter.o query.ext");

print(">annotation going on...\n");
system("$SPAAN_dir/annotate.o query.flt 1gi\\|");

print(">recognition going on...\n");

print(">AAcompo recognition...\n");
system("$SPAAN_dir/AAcompo/AAcompo.o query.ant query.A.NN1in");
system("$SPAAN_dir/AAcompo/recognize.o $SPAAN_dir/AAcompo/spec.reco");

print(">charge recognition...\n");
system("$SPAAN_dir/charge/charge.o query.ant query.c.NN1in");
system("$SPAAN_dir/charge/recognize.o $SPAAN_dir/charge/spec.reco");

print(">hdr recognition...\n");
system("$SPAAN_dir/hdr/hdr.o query.ant query.h.NN1in");
system("$SPAAN_dir/hdr/recognize.o $SPAAN_dir/hdr/spec.reco");

print(">multiplets recognition...\n");
system("$SPAAN_dir/multiplets/multiplets.o query.ant query.m.NN1in");
system("$SPAAN_dir/multiplets/recognize.o $SPAAN_dir/multiplets/spec.reco");

print(">dipep recognition...\n");
chdir("$SPAAN_dir/dipep/");
system("$SPAAN_dir/dipep/querydipep.o $tmp_dir/query.ant $tmp_dir/query.d.NN1in");
#system("$SPAAN_dir/dipep/querydipep.o query.ant query.d.NN1in");

chdir("$tmp_dir");
system("$SPAAN_dir/dipep/recognize.o $SPAAN_dir/dipep/spec.reco");
#system("$SPAAN_dir/dipep/recognize.o $SPAAN_dir/dipep/spec.reco");

print(">final output...\n");
system("$SPAAN_dir/finalp1.o query.A.NN2in query.c.NN2in query.h.NN2in query.m.NN2in query.d.NN2in query.flt query.out");

move("query.out", "$output_file");
chdir("$cur_dir");
rmtree("$tmp_dir");

#unlink("query.std query.ext");
#unlink("query.ant", "query.flt");
#unlink("query.A.
