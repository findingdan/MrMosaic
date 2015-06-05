#!/usr/bin/env perl

# v5: For 4000 probands, includes opt for pos_path
# Min_depth 7 and min_het_baf 0.05 were selected in powerpoint 29_4_2014
# NOTE: These values (10) were picked because it appears that prepared_BLd_runs.txt file, 
# the file used to generate the average depths, used this.

# Load Modules
use warnings;
use strict;
use Data::Dumper; 
use Getopt::Long;
use Vcf;
use Iterator::Simple qw( iter );
use List::MoreUtils  qw( any uniq all );
use Bio::DB::Sam;
use chr2num;

# Define Global Variables	
my $pileup_depth_min_threshold = 7;
my $min_het_baf                = 0.06;
my $min_base_qual              = 10;
my $min_mapping_qual           = 10;


# STRICTLY REQUIRE DEFINITION
my %opt             = process_options();
my $bam_path        = $opt{bam_path} || die;
my $targets_and_l2r = $opt{l2r_path} || die;
my $position_list   = $opt{pos_path} || die;

# MAIN
{
	my $pileup_iter = generate_pileup_iter_from_BAM_at_polymorphic_pos($bam_path, $position_list);
	my $ETR_iter    = generate_ETRs_iter($targets_and_l2r);
	my $pileup_inst = $pileup_iter->next;
	my $ETR_inst    = $ETR_iter->next;

	print "Name\tChr\tPosition\tLog.R.Ratio\tB.Allele.Freq\tGType\n";

	while ($pileup_inst && $ETR_inst) {
		last if ($pileup_inst->{chr} =~ /[XYM]/); 

		# Workflow:
		# for each high-coverage polymorphic position,
		# with sufficient depth,
		# get pileup, 
		# use # of alt alleles to calculate BAF and genotype

		my $depth = ($pileup_inst->{numalt} +  $pileup_inst->{numref});
		if ($depth < $pileup_depth_min_threshold) {	# tis depth points at the number of all hq bases, regardless of which bas
			$pileup_inst = $pileup_iter -> next;
			next;
		}

		# get the L2R from the ETR region holding this position
		(my $l2r, $ETR_inst, $ETR_iter) = get_L2R( $ETR_iter, $ETR_inst, $pileup_inst );

		my $chr = format_chr_for_gada($pileup_inst->{chr});
		my $pos = $pileup_inst->{pos};
		my $baf = sprintf ("%.3f" , $pileup_inst->{numalt} / $depth );	
		my $genotype = convert_baf_to_genotype($baf);
		my $snp_name = "e" . "." . $chr . "." . $pos;

		print join ("\t", $snp_name, $chr, $pos, $l2r, $baf, $genotype), "\n";

		$pileup_inst = $pileup_iter->next;
	}
}

# SUBS
sub format_chr_for_gada {
	my $chr = shift;
	$chr =~ s/23/X/;
	$chr =~ s/24/Y/;
	return $chr;
}

sub process_options {
	my %opt;
	my @options = qw( bam_path=s l2r_path=s pos_path=s );
    GetOptions(\%opt, @options);
	return ( %opt );
}

sub generate_pileup_iter_from_BAM_at_polymorphic_pos {
    my ($bam_path, $position_list_file) = @_;
	Bio::DB::Sam->max_pileup_cnt(250);
    my $bam_obj = Bio::DB::Sam->new( -bam => $bam_path);
	open (my $pos_fh, "<", $position_list_file) or die "$!: Problem, can't open position list\n";
	<$pos_fh>; # skip header
	return iter sub {
		my $pileup_pos_needed = <$pos_fh>;
		return if (! defined $pileup_pos_needed);
		my ($chr, $pos, $ref, $alt) = split /\t/, $pileup_pos_needed;
		my $x = retrieve_pileup_counts($bam_obj, $chr, $pos, $ref, $alt);
		my $numalt = $x->{alt};
		my $numref = $x->{ref};
		if (! defined $numalt) { $numalt = 0 };
		if (! defined $numref) { $numref = 0 };
		return {
			chr         => chr2num($chr),
			pos         => $pos,
			numalt      => $numalt,
			numref		=> $numref
		}
	}
}

sub retrieve_pileup_counts {
    my ($bam_obj, $needed_chr, $needed_pos, $ref, $alt) = @_;
    my %count;
    my $cb = sub {
        my ($seqid, $pos, $pileup_aref, $sam) = @_;
        if ($pos == $needed_pos) {
            for my $pileup (@{$pileup_aref}){
                my $al = $pileup->alignment;
				
				next unless $al->proper_pair == 1; 				# Returns true if mate and pair are both mapped 

                my $baseQual = $al->qscore->[$pileup->qpos];    # Phred qual of base (in numeric format)
                next unless $baseQual >= $min_base_qual;

				my $alignment_mapping_Qual = $al->qual;
				next unless $alignment_mapping_Qual >= $min_mapping_qual;  	# "white" or "yellow" in tview

				(my $cigar_str_letters = $al->cigar_str) =~ s/\d//g;
				next unless $cigar_str_letters =~ m/^M$/;		# Only use alignments with cigar strings of only matches (no soft clipped bases or indels)

                my $qBase = substr($al->qseq, $pileup->qpos, 1);# The base at the needed position
                $count{$qBase}++;
            }
        }
    };
    $bam_obj->fast_pileup("$needed_chr:$needed_pos-$needed_pos", $cb);
	return { alt => $count{$alt}, ref => $count{$ref} } ;
}

sub convert_baf_to_genotype {
	my $baf = shift;
	
	if ($baf < $min_het_baf) {
		return "AA"
	}
	elsif ($baf > (1-$min_het_baf)) {
		return "BB"
	}
	else {
		return "AB"
	}
}

sub get_L2R {
	my ($ETR_iter, $ETR_inst, $BAF_inst) = @_;
	# I think i can use perl's "state" function here instead of returning the insts
	while ($ETR_inst && $BAF_inst) {
		my ($chr_cmp, $pos_region_cmp) = cmp_region_and_pos ( $BAF_inst, $ETR_inst );
		if ($chr_cmp != 0 ) {
			# Need to catch up the ETR 
			$ETR_inst = $ETR_iter->next;
			next;
		}
		if ( $pos_region_cmp == 2) {
			$ETR_inst = $ETR_iter->next;
			next;
		}
		elsif ( $pos_region_cmp >= -1 && $pos_region_cmp <= 1) {
			return ($ETR_inst->{l2r}, $ETR_inst, $ETR_iter);
			# Use the L2R from this region
		}
		else {
			print "there shouldn't be positions before any ETR regions\n";
			die Dumper $ETR_inst, $BAF_inst;
		}
	}
}

sub cmp_region_and_pos {
	my ($position_inst, $region_inst) = @_;
	my $pos_chr = $position_inst->{chr};
	my $pos_pos = $position_inst->{pos};
	
	# Calculate chromosome comparison
	my $region_chr = $region_inst->{chr};
	my $chr_cmp = $pos_chr <=> $region_chr;

	# Calculate position comparison
	my $start = $region_inst->{start};
	my $end = $region_inst->{end};

	my $pos_to_start = $pos_pos <=> $start;
	my $pos_to_end   = $pos_pos <=> $end;
	my $total = $pos_to_start + $pos_to_end;
	
	return ($chr_cmp, $total);
}

sub convert_genotype {
	my $genotype = shift;
	return "NA" if ($genotype eq "NA");
	if ($genotype !~ m/^[01][01]$/) { die "I see $genotype genotype!\n" };
	$genotype =~ s/00/AA/;
	$genotype =~ s/11/BB/;
	$genotype =~ s/01|10/AB/;
	return $genotype;
}

sub generate_ETRs_iter {	
	my $ETR_file = shift;
	open (my $ETR_fh, "<", $ETR_file) or die "$!: Can't open ETR file\n";
	return iter sub {
		my $region_line = <$ETR_fh>;
		return if (! defined $region_line);
		chomp $region_line;
		my ($chr, $start, $end, $l2r_corrected) =  split /\t/, $region_line ;
		$l2r_corrected = sprintf ("%.3f", $l2r_corrected);
		return {
			chr => chr2num( $chr ) ,
			start => $start,
			end => $end,
			l2r => $l2r_corrected
		}
	}	
}

sub cmp_pos {
    my ( $x, $y ) = @_;
	# If the chromosomes are the same, then compare the positions.  
	# <=> comparison operator; returns 0 if x == y, 1 if x > y, -1 if x < y
	# ( ex: 14 <=> 14 = 0.  The left side is 0, so the || then compares the right side )
	# this is a Ray Miller coding genius tip
		#	if ($x->{chr} - $y->{chr} > 1) { die "Files not sorted correctly\n" };
		#	if ($y->{chr} - $x->{chr} > 1) { die "Files not sorted correctly\n" };
		#   Files must be sorted correctly!  Bam should be 1 through 22 in order.
    return $x->{chr} <=> $y->{chr} || $x->{pos} <=> $y->{pos};
}

__END__

=head1 NAME

 BAF_extractor.pl

=head1 VERSION

 Version 0.5

=head1 SYNOPSIS

 perl exome_LRRs_and_BAFs_generator.v5.pl --bam_path <file> --l2r_path <file> --pos_path <file>
 
=head1 DESCRIPTION
 
 This script is designed to extract the fraction of non-reference alleles from a BAM files at positions in target regions.

=head2 Arguments

=over 20

=item C<bam_path> 

a bam file

=item C<l2r_path>

a file containing target regions

=item C<pos_path> 

a file containing positions, reference allele, alt allele

=back

=head1 EXAMPLE

 perl ./BAF_extractor.pl --bam_path my_bam.bam --l2r_path my_l2r_file.txt --pos_path my_pos_file.txt

=head1 AUTHORS

 Dan King < dk6@sanger.ac.uk >
 Project oversight by Matt Hurles

=head1 INSTITUTION

 Wellcome Trust Sanger Institute
 Wellcome Trust Genome Campus,
 Hinxton, Cambridge, CB10 1SA, UK
 
=head1 LICENSE

 Copyright (c) 2014 Genome Research Ltd.
 Author: Dan King <dk6 [at] sanger.ac.uk>
 
 Any redistribution or derivation in whole or in part including any
 substantial portion of this code must include this copyright and
 permission notice.
 
 THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Founl2rion, either version 3 of the License, or
 (at your option) any later version.
 
=cut
