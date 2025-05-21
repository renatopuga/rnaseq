#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;

# config groups
my %groups = (
    '0'  => [qw(S01 S02 S03 S04)],
    '10' => [qw(S01 S02 S03 S04)],
    '20' => [qw(S01 S02 S03 S04)],
    '50' => [qw(S01 S02 S03 S04)],
);

my $reads_per_sample = 7500000;
my $replicates = 5;
my $seed_base = 42;

my $concat_dir = "concat";
my $output_dir = "fake_pools";
my $log_dir    = "logs";

make_path($concat_dir, $output_dir, $log_dir);

# 1. Concat L001 + L002 by sample
print ">> Concat lanes...\n";
for my $group (keys %groups) {
    for my $sample (@{$groups{$group}}) {
        my $merged = "$concat_dir/${sample}-${group}.fastq.gz";
        next if -e $merged;

        print "  + $sample-$group\n";
        my @files = glob("${sample}-${group}_S*_L001_R1_001.fastq.gz ${sample}-${group}_S*_L002_R1_001.fastq.gz");

        unless (@files) {
            warn "  ! File not found $sample-$group\n";
            next;
        }

        system("zcat @files | gzip > $merged") == 0 or die "Error concat $sample-$group com zcat\n";
	#system("zcat @files > $merged") == 0 or die "Error concat $sample-$group\n";
    }
}

# 2. creating fake pools
print ">> Creating fake pools...\n";
for my $group (keys %groups) {
    for my $rep (1..$replicates) {
        my $pool_dir = "$output_dir/group$group/rep$rep";
        make_path($pool_dir);
        my $pool_fastq = "$pool_dir/fake_pool_group${group}_rep${rep}.fastq";
        open my $fh_out, "|-", "gzip > $pool_fastq.gz" or die "Error  $pool_fastq.gz: $!\n";

        print "  > Grupo $group - Réplica $rep\n";

        for my $sample (@{$groups{$group}}) {
            my $input_fastq = "$concat_dir/${sample}-${group}.fastq.gz";
            my $seed = $seed_base + $rep;

            print "    + Subsamples $sample-$group with seed $seed\n";
            my $cmd = "seqtk sample -s$seed $input_fastq $reads_per_sample";
            open my $fh_in, "-|", $cmd or die "Error $cmd: $!\n";
            while (<$fh_in>) {
                print $fh_out $_;
            }
            close $fh_in;
        }

        close $fh_out;
    }
}

print "Fake pools done!\n";
