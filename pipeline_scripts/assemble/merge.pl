#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Parallel::ForkManager; 
no warnings 'experimental::smartmatch';

# DNA sequences of reference, AA sequences of reference, input of exonerate_best.pl result, sample list, output dir of merged contigs
my ($query, $queryp, $exonerateout, $samplelist, $outdir, $help);
my $error_trsd = 0.05; # maximum error rate allowed to consider two sequences aligned 
my $olp_lower_trsd = 10; # minimum overlap between contig
my $olp_upper_trsd = 24; # maximum overlaps found in exonerate_best.pl
my $simi_trsd = 60; # minimum similarity required between query and contig
my $cov_trsd = 80; # minimum coverage required for contig
my $score_trsd = 0; # minimum score required for contig
my $nflen_trsd = 0; # minimum non-flanking sequence required for contig
my $amb_trsd = 1; # maximum percentage of ambiguous nucleotide allowed
my $st_trsd = 0; # maximum stop codon allowed in contig
my $abortnum = 500; # maximum path number allowed to feed into subroutine get best
my $process = 1; # number of cpu used

my $getbesttemp = "getbesttemp"; # run dir for subroutine getbest

my $opt = GetOptions( 'queryn:s', \$query,
					  'queryp:s', \$queryp,
                      'filtered:s', \$exonerateout,
                      'samplelist:s', \$samplelist,
                      'merged:s', \$outdir,
                      'error_rate:s', \$error_trsd,
                      'min_overlap:s', \$olp_lower_trsd,
                      'similarity:s', \$simi_trsd,
                      'len_percentage:s', \$cov_trsd,
                      'raw_score:s', \$score_trsd,
                      'nonflank_len:s', \$nflen_trsd,
                      'ambiguous_percentage:s', \$amb_trsd,
                      'stop_codons:s', \$st_trsd,
                      'cpu:s', \$process,
                      'help|h!', \$help) or usage();
 
# display help message if one of the following option is not specified
if (!($opt && $query && $queryp && $exonerateout && $samplelist && $outdir) || $help) {
usage();
} 

# check the existence of DNA sequences of reference
if (!(-e $query)) {
die ("ERROR: Cannot open $query ($!)");
}

# check the existence of AA sequences of reference
if (!(-e $queryp)) {
die ("ERROR: Cannot open $queryp ($!)");
}

# check the existence of output of exonerate_best.pl
if (!(-e $exonerateout)) {
die ("ERROR: Cannot open $exonerateout ($!)"); 
}

# get the name of sample from sample list
my %sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
$sample{$samplename} = "";
}
close SAMPLELIST;

# create output dir
mkdir $outdir; # outdir
mkdir "$outdir/nf"; # non-flanking dir
mkdir "$outdir/f"; # flanking dir
mkdir "$outdir/p"; # amino acid dir
mkdir "$getbesttemp"; # run dir for subroutine getbest

# input dir 
my $filtered = "$exonerateout/filtered"; # dir containing filtered contigs
my $exonerate = "$exonerateout/exonerate"; # dir containing exonerate result
my $overlap = "$exonerateout/overlap"; # dir containing overlapping information

# only analysis sample exist in sample list
my @exoneratesub; 
opendir FILTEREDDIR, $filtered;
while (my $exoneratesub = readdir FILTEREDDIR) {
push @exoneratesub, $exoneratesub if (exists $sample{$exoneratesub});
}
closedir FILTEREDDIR;

# save reference sequence into hash $query{loci name} = sequence
my %query;
open QUERY, $query;
while (my $line = <QUERY>) {
	if ($line =~ />(\S+)/) {
	chomp(my $seq = <QUERY>);
    $query{$1} = $seq;
	}
}
close QUERY;

# create dir of amino acid sequences of reference and construct a hash $queryp{loci name} = sequence
my %queryp;
if (-e "queryp") { # save hash only if queryp exists
	open QUERYP, $queryp;
	while (my $line = <QUERYP>) {
		if ($line =~ />(\S+)/) {
		chomp(my $pseq = <QUERYP>);
		$queryp{$1} = $pseq;
		}
	}
	close QUERYP;
} else { # save hash and create dir only if queryp not exists
	mkdir "queryp";
	open QUERYP, $queryp;
	while (my $line = <QUERYP>) {
		if ($line =~ />(\S+)/) {
		chomp(my $pseq = <QUERYP>);
		open QUERYPOUT, ">queryp/$1.fas";
		$queryp{$1} = $pseq;
		print QUERYPOUT ">$1\n$pseq\n";
		close QUERYPOUT;
		}
	}
	close QUERYP;
}

# split reads into $split set
my $split = $process;

my @queryp = sort keys %queryp;
my $splited_array = split_array(\@queryp, $split);

function("assembling contigs further and retrieving best contigs of each gene");

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($process/2));

# merge contigs foreach chunk in parallel
DATA_LOOP: while (my $array = shift @$splited_array) {
$pm->start and next DATA_LOOP;
	
	# merge contigs for each gene in each chunk
	while (my $gene = shift @$array) {
	
		my %result;
		foreach my $exoneratesub (@exoneratesub) {
			if (open INFILE, "$filtered/$exoneratesub/$gene.fas") {

				# get the information of each contig
				my %sub;
				while (my $line = <INFILE>) {
					if ($line =~ />(\S+)\|(\S+)\|(\S+)\|(\S+)\|(\S+)\|(\S+)/) {
					$sub{$1}->{stq} = $2; # start point on query
					$sub{$1}->{edq} = $3; # end point on query
					$sub{$1}->{stsm} = $4; # start point on contig
					$sub{$1}->{edsm} = $5; # end point on contig
					$sub{$1}->{score} = $6; # alignment score
					chomp(my $seq = <INFILE>);
					$sub{$1}->{seqf} = $seq; # contig
					}
				}
				close INFILE;
				
				# read the overlapping information
				my %overlap;
				if (open OVERLAP, "$overlap/$exoneratesub/$gene.txt") {
					while (my $line = <OVERLAP>) {
					$line =~ /(\S+)\|(\S+)\s(\S+)/; # contig1|contig2 length of overlap
					$overlap{"$1\|$2"} = $3;
					}
					close OVERLAP;
					
					# if found overlap in this gene, construct overlapping graph and find all possible path in the graph
					if (scalar keys %overlap > 0) {
						open EXONERATE, "$exonerate/$exoneratesub/$gene.txt";
						while (my $line = <EXONERATE>) {
							if ($line =~ /sugar:\s\S+\s(\d+)\s(\d+)\s\S+\s(\S+)\s\d+\s\d+\s[+|-]\s\d+\s\S+/) {
								my $stq = $1;
								my $edq = $2;
								$edq--;
								my $subname = $3;
								if ((exists $sub{$subname})&&($sub{$subname}->{stq} == $stq)&&($sub{$subname}->{edq} == $edq)) { # ascertain the contig in overlapping information by its name start and end point on query
								# save the information of per transition score (alignment score for each nucleotide pair)
								chomp(my $pt = <EXONERATE>);
								my @score = $pt =~ /\S+\_\S+\_\-*\d+/g;
								pop @score;
								shift @score;
								$sub{$subname}->{pt} = \@score;
								}
							}
						}
						close EXONERATE;
						
						# merge sequences base on overlap grpah
						my $graph = graph(\%overlap, \%sub, $gene); # construct overlap graph
						my ($reducedgraph) = graphreduction($graph, \%sub, $gene); # reduce the transitive edges	
						my $allpath = findallpath($reducedgraph, $gene, $exoneratesub); # work through all path in the resulting graph
						merge($allpath, $reducedgraph, \%sub); # merge each path into longer contig	                
					}
				}
						
				my ($bestsub, $simi, $bestseqnf, $bestprotseq) = getbest($gene, \%sub); # get the best contig
				
				# save best result into %result
				if ($bestsub) {
				$result{$exoneratesub}->{nf} = uc($bestseqnf);
				$result{$exoneratesub}->{f} = uc($sub{$bestsub}->{seqf});
				$result{$exoneratesub}->{p} = $bestprotseq;
				$result{$exoneratesub}->{simi} = $simi;
				}
			}
		}
		
		# write output of non-flanking, flanking and amino acid sequences
		if (scalar keys %result > 0) {
		open NF, ">$outdir/nf/$gene.fas";
		open F, ">$outdir/f/$gene.fas";
		open P, ">$outdir/p/$gene.fas";
		print NF ">$gene\n$query{$gene}\n";
		print F ">$gene\n$query{$gene}\n";
		print P ">$gene\n$queryp{$gene}\n";
			map {
			print NF ">$_\n$result{$_}->{nf}\n";
			print F ">$_\n$result{$_}->{f}\n";
			print P ">$_\n$result{$_}->{p}\n";
			} sort keys %result;
		close NF;
		close F;
		close P;
		}
	}
$pm->finish;
}
$pm->wait_all_children;

print "All coding sequences has been written to $outdir/nf/ \n";
print "All flanking sequences has been written to $outdir/f/ \n";
print "All amino acid sequences has been written to $outdir/p/ \n";

# `rm -rf query`;
# `rm -rf queryp`;
`rm -rf $getbesttemp`;

##################################################
# subroutine
##################################################

# subroutine to translate nucleotide into amino acid sequence
sub translate {
my $seq = shift @_; # DNA sequence

my $start_codon = shift @_;
$start_codon --; 
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X");
$pro=$pro->seq; 
return($pro);
}

# subroutine to report the start of a step of program
sub function {
my $function = shift;

if ($process > 1){
print "Start $function in parallel\n";
} else {
print "Start $function in serial\n";
}
}

# subroutine to construct overlap graph
sub graph {
my $overlap = shift; # overlapping information
my $sub = shift; # hash of contigs and their information
my $gene = shift; # loci name

# construct overlapping graph

# save the next node of each node and overlapping information in graph
my (%graph, %startemp, %endtemp); # hash of overlapping graph, contigs at the start of graph, contigs at the end of graph
foreach my $overlapinfo (sort keys %$overlap) {
my ($cur, $nex) = $overlapinfo =~ /(\S+)\|(\S+)/; # get the overlap pair
my ($assication, $content) = edgeinfo($sub, $overlap->{$overlapinfo}, $cur, $nex, $gene); # find exact overlap between the pair
	
	# if they are qualified overlap save the information of assication(overlap) and content(how many overlapped nucleotide) between one node and its next node
	if (($assication)&&($content)) {
	$graph{$cur}->{nex}->{$nex}->{assication} = $assication; 
	$graph{$cur}->{nex}->{$nex}->{content} = $content;
	$startemp{$cur} += 1;
	$endtemp{$nex} += 1;
	} else {
	delete $overlap->{$overlapinfo};
	}
}

# if a contig always on the right, then it must the end of a graph		
foreach (sort keys %endtemp) {
$graph{$_}->{nex}->{end} = "" if (!(exists $startemp{$_}));
}

# save the previous node of each node in graph

# if a contig always on the left, then it must the start of a graph		
foreach (sort keys %startemp) {
$graph{$_}->{pre}->{start} = "" if (!(exists $endtemp{$_}));
}

# save the previous node of internal node
foreach my $overlapinfo (sort keys %$overlap) {
my ($pre, $cur) = $overlapinfo =~ /(\S+)\|(\S+)/;
$graph{$cur}->{pre}->{$pre} = "";
}	

# return hash of graph			
return(\%graph);
}

# subroutine to get the exact overlap information between to contigs
sub edgeinfo {
my $sub = shift; # hash of contigs and their information
my $overlapinfo = shift; # overlapping information
my $cur = shift; # contig name of current contig
my $nex = shift; # contig name of next contig
my $gene = shift; # loci name

# get the information of the sequence
my $seq1qst = $sub->{$cur}->{stq}; # start point on reference of seq1
my $seq1qed = $sub->{$cur}->{edq}; # end point on reference of seq1
my $seq1sst = $sub->{$cur}->{stsm}; # start point on contig of seq1
my $seq1sed = $sub->{$cur}->{edsm}; # end point on contig of seq1
my $seq1lf = $seq1sst; # left unaligned sequence length of seq1
my $seq1rf = length($sub->{$cur}->{seqf})-$seq1sed-1; # right unaligned sequence length of seq1
my $seq1 = $sub->{$cur}->{seqf}; # sequence 1
my @score1 = @{$sub->{$cur}->{pt}}; # per transition score of seq1
# same for seq2
my $seq2qst = $sub->{$nex}->{stq}; 
my $seq2qed = $sub->{$nex}->{edq};
my $seq2sst = $sub->{$nex}->{stsm};
my $seq2sed = $sub->{$nex}->{edsm};
my $seq2lf = $seq2sst;
my $seq2rf = length($sub->{$nex}->{seqf})-$seq2sed-1;
my $seq2 = $sub->{$nex}->{seqf};
my @score2 = @{$sub->{$nex}->{pt}};
	
	# scalar for how two contigs connect ("o" only, means overlap), scalar for how many overlap between two contigs
	my ($assication, $content);
	if ($overlapinfo ne "na") { # if this overlap is detected by sga
		my $olp = $overlapinfo;
		if (fp($seq2sst, $olp) == fp($seq1sst, length($seq1))) { # if these two sequence is also overlapped on codon position
		# calculate score of overlapped portion of 2 contigs in alignment with reference
		my @rscore1 = reverse @score1;
		my $leftscore = olpscore($olp-$seq1rf, \@rscore1); # score of current contig
		my $rightscore = olpscore($olp-$seq2lf, \@score2); # score of next contig
		($assication, $content) = ("o", "$olp $leftscore $rightscore");     
		}
	} else { # if this overlap is detected by exonerate_best.pl
	($assication, $content) = olp($seq1, $seq2, $seq1sst, $seq2sst, \@score1, \@score2, $seq1rf, $seq2lf); # find exact overlap
	}

# return assication and content of overlap
return ($assication, $content);
}

# subroutine to find exact overlap between to contigs
sub olp {
my $seq1 = shift; # current contig
my $seq2 = shift; # next contig
my $seq1stsm = shift; # start position on current contig
my $seq2stsm = shift; # start position on next contig
my $score1 = shift; # pt score on current contig 
my $score2 = shift; # pt score on next contig
my $seq1rf = shift; # length of unaligned sequence in current contig 
my $seq2lf = shift; # length of unaligned sequence in next contig 
	
# find the overlap by compare the right of current contig and the left of next contig
for (my $i = $olp_upper_trsd; $i >= $olp_lower_trsd; $i--) {
	# get assumed overlap
	my $firolp = substr $seq1, -$i; 
	my $nexolp = substr $seq2, 0, $i;
	if ((length($seq1) > $i)&&(length($seq2) > $i)) { # if sequences are at least longer than overlap
		if (diff($firolp, $nexolp) <= int($error_trsd*length($firolp))) { # if number mismatch is lower than threshold
			if (fp($seq2stsm, length($nexolp)) == fp($seq1stsm, length($seq1))) { # codon position is also overlapped
			# calculate score of overlapped portion of 2 contigs in alignment with reference
			my @rscore1 = reverse @$score1;
			my $leftscore = olpscore((length($firolp)-$seq1rf), \@rscore1);
			my $rightscore = olpscore((length($nexolp)-$seq2lf), $score2);
			my $len = length($firolp);
			return("o", "$len $leftscore $rightscore");     
			last;
			}
		}
	}
}
}

# subroutine to find the codon position of provided base
sub fp {
my $frst = shift; # position of first codon
my $rp = shift; # postulate the codon positon of this nucleotide

my $dis = $rp-$frst; # distance between first codon and query postion

# find the position
my $pos = do {
	if ($dis > 0) { # if query postion is after position of first codon
	$dis%3;
	} else { # if query postion is before position of first codon
	(3-abs($dis)%3)%3;
	}
};

# return codon position
return $pos;
}

# subroutine to find the number of different nucleotide in two sequences
sub diff {
my $seq1 = shift; # seq1
my $seq2 = shift; # seq2

# compare each column of nucleotide
my $diff = 0;
for (my $i=0; $i<length($seq1); $i++) {
my $nucl1 = substr $seq1, $i, 1;
my $nucl2 = substr $seq2, $i, 1;
$diff++ if ($nucl1 ne $nucl2);
}

return $diff;
}

# subroutine to calculate the alignment score of overlapping sequences to query
sub olpscore {
my $olplen = shift; # length of overlap in nucleotide
my $score = shift; # pt score

$olplen = int($olplen/3); # length of overlap in AA

# get score
my $olpscore = do { # if overlap length > 0
	if ($olplen > 0) { 
		my ($count, $all);
		# add up score of each positon covered in overlap region
		foreach (@$score) {
		my ($nucleo, $score) = $_ =~ /\S+\_(\S+)\_(\-*\d+)/;
		$count++ if ($nucleo =~ /[A-Z]{3}/);
		$all += $score;
		last if ($count == $olplen);
		} 
		$all;
	} else {
	0;
	}
};

# return overlap score
return($olpscore);
}

# subroutine to reduce furcation in graph
sub graphreduction {
my $graph = shift; # input graph
my $sub = shift; # hash of contig
my $gene = shift; # loci name
	
my $lag_allfurcation = 0;
while () {
	# get all nodes connect with forward furcation 
	my %furcation;
	foreach my $node (sort keys %$graph) {
	$furcation{$node} = "" if ((scalar keys %{$graph->{$node}->{nex}}) > 1);
	}
	
	# get all nodes connect with reverse furcation 
	my %reversefurcation;
	foreach my $node (sort keys %$graph) {
	$reversefurcation{$node} = "" if ((scalar keys %{$graph->{$node}->{pre}}) > 1);
	}
	
	# stop reduction if there's no furcation or it cannot be reduced further
	my $allfurcation = (scalar keys %furcation) + (scalar keys %reversefurcation);
	if (($allfurcation == 0)||($allfurcation == $lag_allfurcation)) {
	return ($graph);
	last;
	}

	$lag_allfurcation = $allfurcation;
	
	# reduce forward furcation 
	foreach my $furcation (sort keys %furcation) {
	$graph = simple_furcation($furcation, $graph, $sub);
	}
	
	# reduce reverse furcation
	foreach my $reversefurcation (sort keys %reversefurcation) {
	$graph = simple_reversefurcation($reversefurcation, $graph, $sub);
	}	
}
}

# subroutine to remove forward tips and bubbles in overlap graph
sub simple_furcation {
my $furcation = shift; # nodes with furcations
my $graph = shift; # graph
my $sub = shift; # hash of contigs
    
# find all paths of forward tips and bubble
my (%endchainhash, %convergechainhash);
foreach my $node (sort keys %{$graph->{$furcation}->{nex}}) {
	# start from node which connect with forward furcation
	my @chain;
	push @chain, $furcation;        
	while () {
		push @chain, $node;
		my @nexnode = sort keys %{$graph->{$node}->{nex}}; # get the next node
		my @prenode = sort keys %{$graph->{$node}->{pre}}; # get the previous node
		if ((@nexnode == 1) && (@nexnode == @prenode)&&(!(exists $graph->{$node}->{nex}->{end}))) { # if there's only one pre-node -> one current-node -> one next-node, then this node is in the path of a furcation
		$node = $nexnode[0];
		next;
		}
		if ((exists $graph->{$node}->{nex}->{end})&&(@prenode == 1)) { # if next node is a end node and only one pre-node, then this chain must a tip
		push @{$endchainhash{$node}}, \@chain;
		last;
		}
		if (@prenode > 1) {
		push @{$convergechainhash{$node}}, \@chain; # find the node which pre-node is more than 1, it probably form a bubble
		last;
		}
		if (@nexnode > 1) { # stop walking if find a furcation
		last;
		}
	}		
}

# find the bubble
foreach my $converge (sort keys %convergechainhash) {
	# if more than one node converge at node $converge, then it must have form a bubble
	if (@{$convergechainhash{$converge}} > 1) {
		my $max = 0;
		my $maxchain;
		foreach my $chain (@{$convergechainhash{$converge}}) {
		my $chainscore = total_score($chain, $graph, $sub); # get the total score of the chain
			if ($chainscore > $max) { # save only the best chain
			$max = $chainscore;
			$maxchain = $chain;
			}
		}
		
		# delete all the chain except the best chain
		my @deletechain = grep {$_->[1] ne $maxchain->[1]} @{$convergechainhash{$converge}};
		foreach my $deletechain (@deletechain) {
			delete $graph->{$deletechain->[0]}->{nex}->{$deletechain->[1]};
			delete $graph->{$deletechain->[-1]}->{pre}->{$deletechain->[-2]};
			for (1..(@$deletechain-2)) {
			delete $graph->{$deletechain->[$_]};
			}
		}
	}
}

# get the forward tips
my @end = grep {@{$endchainhash{$_}}  == 1} sort keys %endchainhash;
my $max = 0;
my $maxchain;
foreach my $endnode (@end) {
	my $chainscore = total_score($endchainhash{$endnode}->[0], $graph, $sub); # get score of all chain
	if ($chainscore > $max) { # get best chain
	$max = $chainscore;
	$maxchain = $endchainhash{$endnode}->[0];
	}
}
# delete all tips except the best one
my @deletendnode = grep {$endchainhash{$_}->[0]->[1] ne $maxchain->[1]} @end;
foreach my $deletendnode (@deletendnode) {
	my @deletearray = @{$endchainhash{$deletendnode}->[0]};
	delete $graph->{$deletearray[0]}->{nex}->{$deletearray[1]};
	for (1..(@deletearray-1)) {
	delete $graph->{$deletearray[$_]};
	}
}
 
# return simplified graph   
return $graph;
} 

# subroutine to remove reverse tips in overlap graph
sub simple_reversefurcation {
my $furcation = shift; # node with reverse furcations
my $graph = shift; # graph
my $sub = shift; # hash of contigs
    
# get all reverse tips
my %reversechainhash;
foreach my $node (sort keys %{$graph->{$furcation}->{pre}}) {
	my @chain;
	push @chain, $furcation;        
	while () {
		push @chain, $node;
		my @prenode = sort keys %{$graph->{$node}->{pre}};
		my @nexnode = sort keys %{$graph->{$node}->{nex}};
		if ((@prenode == 1) && (@prenode == @nexnode)&&(!(exists $graph->{$node}->{pre}->{start}))) { # if there's only one pre-node -> one current-node -> one next-node, then this node is in the path of a reverse furcation
		$node = $prenode[0];
		next;
		}
		if ((exists $graph->{$node}->{pre}->{start})&&(@nexnode == 1)) { # end if there's only one next node and pre node is start node, then this chain must a tip
		push @{$reversechainhash{$node}}, \@chain;
		last;
		}
		if ((@prenode > 1)||(@nexnode > 1)) { # do not work if we find a reverse or forward furcation
		last;
		}
	}		
}

# get the reverse tips
my @end = grep {@{$reversechainhash{$_}} == 1} sort keys %reversechainhash;
my $max = 0;
my $maxchain;
foreach my $startnode (@end) {
	my @reversechain = reverse @{$reversechainhash{$startnode}->[0]}; 
	my $chainscore = total_score(\@reversechain, $graph, $sub); # get the score of chain
	if ($chainscore > $max) { # save only best chain
	$max = $chainscore;
	$maxchain = $reversechainhash{$startnode}->[0];
	}
}
# delete all chain except the best one
my @deletendnode = grep {$reversechainhash{$_}->[0]->[1] ne $maxchain->[1]} @end;
foreach my $deletendnode (@deletendnode) {
	my @deletearray = @{$reversechainhash{$deletendnode}->[0]};
	delete $graph->{$deletearray[0]}->{pre}->{$deletearray[1]};
	for (1..(@deletearray-1)) {
	delete $graph->{$deletearray[$_]};
	}
}

# return graph
return $graph;
}

# subroutine to get the score of provided chain
sub total_score {
my $chain = shift; # one-way connected node
my $graph = shift; # graph
my $sub = shift; # hash of contig
	
# intergrate over the score of each contig in the chain
my $count = 0;
my $totalscore = $sub->{$chain->[$count]}->{score};
while ($count <= (@$chain-2)) {
	$totalscore+=$sub->{$chain->[$count+1]}->{score}; # add the score of two contig together
	my $assication = $graph->{$chain->[$count]}->{nex}->{$chain->[$count+1]}->{assication}; # get the assocation(overlap only) 
	my $content = $graph->{$chain->[$count]}->{nex}->{$chain->[$count+1]}->{content}; # get the amount of overlap
	if ($assication eq "o") {
	my ($leftscore, $rightscore) = $content =~ /\S+\s(\S+)\s(\S+)/;
	($leftscore, $rightscore) = ($rightscore, $leftscore) if ($leftscore < $rightscore);
	$totalscore-=$rightscore; # minus total by lower score
	}
	$count++;
}

# return score of chain
return $totalscore;
} 

# subroutine to find all path in a graph
sub findallpath {
my $graph = shift; # simplified graph
my $gene = shift; # loci name
my $sub = shift; # hash of contigs

# get the start contig
my @start = grep {exists $graph->{$_}->{pre}->{start}} sort keys %$graph;

# find all path from each start postion
my @allpath;
foreach my $start (@start) {
	my %chain;
	push @{$chain{0}}, $start; # push start contig into chain first
	while () {
		foreach my $subchainnum (sort keys %chain) {
			my $node = $chain{$subchainnum}->[-1];
			while () {
				my @nex = sort keys %{$graph->{$node}->{nex}};
				if (@nex == 1) { # if there's only one next contig
					if ($nex[0] ne "end") { # if next node is not end node
						if ($nex[0] ~~ @{$chain{$subchainnum}}) { # if a node was found in a previously bulit chain, it may form a loop.
						push @allpath, \@{$chain{$subchainnum}}; # save the chain
						delete $chain{$subchainnum};
						last;
						} else {
						push @{$chain{$subchainnum}}, $nex[0]; # continue to walk in the path
						$node = $nex[0];
						}	
					} else { # if it go to the end node
					push @allpath, \@{$chain{$subchainnum}}; # save the chain
					delete $chain{$subchainnum}; 
					last;
					}
				} else { # if there are more than one next contig
					for (1..(@nex-1)) { # duplicate the chain and add first node of the furcation to the end of each duplicated chain
					my @num = sort {$a<=>$b} keys %chain;
					my $last = $num[-1]+1;
					push @{$chain{$last}}, @{$chain{$subchainnum}};
					push @{$chain{$last}}, $nex[$_];
					}
					push @{$chain{$subchainnum}}, $nex[0];
					last;	
				}	
			}
		}
		my $chainnum = scalar keys %chain;
		last if (!(%chain)); # last if all path has been found
		if ($chainnum > $abortnum) { # abort if there's too many paths
		print STDERR "Too much contigs($chainnum) start with $start in $gene\n";
		last;
		}
	}
}

# return all paths
return \@allpath;				
}

# subroutine to merge contigs in a path
sub merge {
my $allpath = shift; # all path
my $graph = shift; # simplified graph
my $sub = shift; # hash of contigs
    
# merge for each path
foreach my $path (@$allpath) {
	my @path = @$path;
	
	my $merge = $sub->{$path[0]}->{seqf}; # add the first contig
	for my $nodenum (0..(@path-2)) { # merge from second contig to the end
		my $curnode = $path[$nodenum];
		my $nexnode = $path[$nodenum+1];
		my $assication = $graph->{$curnode}->{nex}->{$nexnode}->{assication}; # get the association (overlap only)
		my $content = $graph->{$curnode}->{nex}->{$nexnode}->{content}; # get the exact amount of overlap
		if ($assication eq "o") {
			my ($olp, $score1, $score2) = $content =~ /(\S+)\s(\S+)\s(\S+)/;
			$merge = do{ # keep the overlap which got higher alignment score
				if ($score1>=$score2) {
				my $seq2 = substr $sub->{$nexnode}->{seqf}, $olp;
				$merge.$seq2;
				} else {
				my $seq1 = substr $merge, 0, (length($merge)-$olp);
				$seq1.$sub->{$nexnode}->{seqf};
				}
			};
		} 
	}
			
	my $subname = join "\|", @path;
	$sub->{$subname}->{seqf} = $merge; # join all contig name together as the name of the path
}
}

# subroutine to split array in to several sub-array
sub split_array {
my $array = shift; # array of loci name
my $split = shift; # parts to split
	
$split++ if ($split == 0);
my $block_num = int(@$array/$split);
$block_num++ if ($block_num == 0);

my @splited_array;
my $size_cnt = 0;
my $block_cnt = 0;

while (my $file = shift @$array) {
	if ($size_cnt >= $block_num) {
	$size_cnt = 0;
	$block_cnt ++;
	}
push @{$splited_array[$block_cnt]}, $file;
$size_cnt ++;
}

# return splited dataset
return(\@splited_array);
}

# subroutine to get best contigs of each gene
sub getbest {
my $query = shift; # reference sequence
my $sub = shift; # hash of merged assemblies

# write the contigs to temporary file 
open TEMP, ">$getbesttemp/$query.fas";
map {print TEMP ">$_\n$sub->{$_}->{seqf}\n"} sort keys %$sub;
close TEMP;

# run exonerate to align query against contigs
my $exonerate = `exonerate -m protein2dna queryp/$query.fas $getbesttemp/$query.fas --subopt FALSE --showvulgar FALSE --showalignment FALSE --score 0 --ryo ">%ti %s %ps %qal %ql %td\n%tas"`;
my %bestseq;
while ($exonerate =~ />\S+\s\d+\s\d+\.\d+\s\d+\s\d+\s*\[*[a-z]*\]*\n[A-Z\*\n]+/gs) {
next if ($& =~ /revcomp/);
# get contig name, score, similarity, alignment length, query length, sequence
my ($id, $score, $similarity, $qalnlen, $qlen, $seq) = $& =~ />(\S+)\s(\d+)\s(\d+\.\d+)\s(\d+)\s(\d+)\s*\[*[a-z]*\]*\n([A-Z\*\n]+)/s;
my $coverage = sprintf("%.2f", $qalnlen/$qlen*100); # compute coverage alignment length/query length
next if (($similarity < $simi_trsd)||($coverage < $cov_trsd)||($score < $score_trsd)); # skip unqualified contigs
	
	# keep the best contig with best hit if a contig has multiple hit
	if (exists $bestseq{$id}) {
	next if ($score < $bestseq{$id}->{score});
	}

$seq =~ s/\n//g;
my ($protseq) = translate($seq, 1) =~ /(\S+\w+)\**/; # translate nucleotide into amino acid
my $nfseq = substr $seq, 0, length($protseq)*3; # get coding sequence
my $nflen = length($nfseq); # length of non-flanking sequence
my @ambiguous = $sub->{$id}->{seqf} =~ /[^ATCG]/g; # get number of ambiguous bases in sequence
my $flen = length($sub->{$id}->{seqf});  # length of flanking sequence
my $ambiguous = sprintf("%.2f", @ambiguous/$flen*100); # percentage of ambiguous bases
my @stcodon = $protseq =~ /\*/g; # number of stop codon 

next if ((@stcodon > $st_trsd)||($nflen < $nflen_trsd)||($ambiguous > $amb_trsd)); # skip unqualified contigs

# save the qualified sequences
$bestseq{$id}->{seq} = $nfseq;
$bestseq{$id}->{protseq} = $protseq;
$bestseq{$id}->{wholelen} = length($sub->{$id}->{seqf});
$bestseq{$id}->{score} = $score;
$bestseq{$id}->{simi} = $similarity;
$bestseq{$id}->{cov} = $coverage;
}

# if there are seq left after filter
if (%bestseq) { 
my @score = sort {$bestseq{$a}->{score} <=> $bestseq{$b}->{score}} keys %bestseq;
@score = reverse @score;
	
	# find the seq with highest score, if score is the same, then keep the longer one
	my $bestsub = $score[0];
	for (1..(@score-1)) {
		if ($bestseq{$bestsub}->{score} == $bestseq{$score[$_]}->{score}) {
			if ($bestseq{$bestsub}->{wholelen} < $bestseq{$score[$_]}->{wholelen}) {
			$bestsub = $score[$_];
			}
		} else {
		last;
		}
	}

# return contig name, similarity, DNA and AA sequence	
return ($bestsub, $bestseq{$bestsub}->{simi}, $bestseq{$bestsub}->{seq}, $bestseq{$bestsub}->{protseq});
}
}

sub usage {
print STDERR "
Script name: merge.pl

This is a script to assemble contigs further and retrieve best contigs for each locus

Dependencies:
(1) Exonerate v2.2.0 or higher
(2) Perl Module:
	1. Bio::Seq (Included in Bioperl)
	2. Parallel::ForkManager

Example usage:

	perl merge.pl --queryn query.dna.fas --queryp query.aa.fas --filtered exonerateout --merged mergeout --samplelist samplelist.txt

Input files:
(1) query.dna.fas
(2) query.aa.fas
(3) exonerateout
(4) samplelist.txt

Output files:
(1) mergeout

Options:
--queryn
  Full coding nucleotide sequences of reference in fasta format
--queryp
  Amino acid sequences of reference in fasta format 
--filtered
  Directory containing filtered and potentially overlapped contigs
--samplelist
  A list of sample name. Only samples in the list will be assembled
--merged
  Directory containing further assembled nucleotide coding sequences(nf), whole sequences with flanking region(f), and amino acid sequences(p)
--error_rate
  Maximum error rate allowed to consider two contigs aligned, $error_trsd in default
--min_overlap
  Minimum overlap length between contigs required for further assemble, $olp_lower_trsd in default
--similarity
  Minimum similarity required between contigs and query, $simi_trsd in default
--len_percentage
  Keep contigs if length of aligned query in alignment is more than X percentage of whole query length, $cov_trsd in default
--raw_score
  Minimum raw score required for query-contig alignment, $score_trsd in default
--nonflank_len
  Minimum length of non-flanking sequences required, $nflen_trsd in default
--ambiguous_percentage
  Maximum allowed percentage of ambiguous nucleotides in contigs, $amb_trsd in default
--stop_codons
  Maximum allowed number of stop codon in a contig, $st_trsd in default
--cpu
  Limit the number of CPUs, $process in default.
--help , -h
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}
