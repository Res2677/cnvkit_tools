use Data::Dumper;
open IN,'gencode.v22.annotation.gtf',or die $!;
sub min { 
my $currentMinCnt = shift @_;
foreach ( @_ ){
	if ( $_ < $currentMinCnt )
	{
		$currentMinCnt = $_;
	}
}
return $currentMinCnt;
}

while(<IN>)
{
	chomp;
	next if (/##/);
	@a = split/\t/,$_;
	if ($a[8]=~/gene_name "(\S+)";/){$gene = $1;}
	if ($a[8]=~/transcript_name "(\S+)";/){$transcript = $1;}
	#if ($a[2]=~/exon/){push @{$hash{$gene}{$a[3]}},$a[4];}
	if ($a[2]=~/transcript/){$transsta{$gene}{$transcript} = $a[3];}
	if ($a[2]=~/transcript/){$transend{$gene}{$transcript} = $a[4];}
	if ($a[2]=~/exon/){$hash{$gene}{$transcript}{$a[3]} = $a[4];}
}

foreach $gene (keys %hash)
{
	foreach $tr (keys %{$hash{$gene}})
	{
		foreach $ee (keys %{$hash{$gene}{$tr}})
		{
			$num{$gene}{$tr} ++;
		}
	}
}

$tr_n =2;


foreach $gene (keys %hash)
{
	foreach $cc (sort {$num{$gene}{$a} <=> $num{$gene}{$b}} keys %{$num{$gene}})
	{
		$ttt = $cc;
	}
	$gt{$gene} = $ttt;
}
#print Dumper \%gt;

foreach $gene (keys %hash)
{
	$ks = $transsta{$gene}{$gt{$gene}};
	$ke = $transend{$gene}{$gt{$gene}};
	foreach $start (sort keys %{$hash{$gene}{$gt{$gene}}})
	{
		#print "$start\n";
		if ($start >$ks){$end = $start-1;print "$gene\t$kk\t$end\tintron\n"};
		print "$gene\t$start\t$hash{$gene}{$gt{$gene}}{$start}\texon\n";
                $kk = $hash{$gene}{$gt{$gene}}{$start}+1;
		#print "$gene\t$kk\t$end\tintron\n";
		#$kk = $end+1;
	}
	if ($hash{$gene}{$start} > $ke){print "$gene\t$kk\t$ke\tintron\n";}
}
