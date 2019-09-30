#!/usr/bin/perl
use File::Copy;
use File::Path;
use File::Basename;
use Cwd;

$SPAAN_dir=dirname(__FILE__);

$cur_dir=cwd();
$query_file = $ARGV[0];
$output_file = $ARGV[1];

print ">moving files...\n";

chdir("$SPAAN_dir");
copy("$query_file","query.fasta");

print ">standardization going on...\n";
system("./standard.o query.fasta");
move("query.std", "query.ext");

print(">filteration going on...\n");
system("./filter.o query.ext");

print(">annotation going on...\n");
system("./annotate.o query.flt 1gi\\|");

print(">recognition going on...\n");

print(">AAcompo recognition...\n");
chdir("AAcompo");
system("./AAcompo.o ../query.ant query.NN1in");
system("./recognize.o spec.reco");

print(">charge recognition...\n");
chdir("../charge");
system("./charge.o ../query.ant query.NN1in");
system("./recognize.o spec.reco");

print(">hdr recognition...\n");
chdir("../hdr");
system("./hdr.o ../query.ant query.NN1in");
system("./recognize.o spec.reco");

print(">multiplets recognition...\n");
chdir("../multiplets");
system("./multiplets.o ../query.ant query.NN1in");
system("./recognize.o spec.reco");

print(">dipep recognition...\n");
chdir("../dipep");
system("./querydipep.o ../query.ant query.NN1in");
system("./recognize.o spec.reco");

print(">final output...\n");
chdir("..");
system("./finalp1.o AAcompo/query.NN2in charge/query.NN2in hdr/query.NN2in multiplets/query.NN2in dipep/query.NN2in query.flt query.out");

copy("query.out", "$output_file");
unlink("query.ant", "query.flt");
chdir("$cur_dir");
