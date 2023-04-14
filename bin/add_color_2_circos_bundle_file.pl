#!/usr/bin/perl
use strict;
use warnings;

=head1 DESCRIPTION

Adds colours to a CIRCOS bundle file.

=head1 AUTHOR

Original: Ross Crowhurst L<mailto:ross.crowhurst@plantandfood.co.nz>
Modified: Usman Rashid L<mailto:usman.rashid@plantandfood.co.nz>

=cut

my $low = 0;

my %bundleColorsRGB = (
        3000 => "128,0,0,0.5",
        2000 => "229,0,10,0.5",
        1500 => "229,19,9,0.5",
        1000 => "216,38,8,0.5",
        500  => "210,57,7,0.5",
        250  => "204,76,6,0.5",
        100  => "198,95,5,0.5",
        50   => "192,114,4,0.5",
        25   => "186,113,3,0.5",
        10   => "180,152,2,0.5",
        5    => "174,171,1,0.5",
        0    => "168,191,0,0.5"
);

my %bundleColorsRGBLow = (
        55 => "128,0,0,0.5",
        50 => "229,0,10,0.5",
        45 => "229,19,9,0.5",
        40 => "216,38,8,0.5",
        35  => "210,57,7,0.5",
        30  => "204,76,6,0.5",
        25  => "198,95,5,0.5",
        20   => "192,114,4,0.5",
        15 => "186,113,3,0.5",
        10 => "180,152,2,0.5",
        5    => "174,171,1,0.5",
        0    => "168,191,0,0.5"
);

sub usage {
	print "USAGE: $0 -i=bundle_file_in -o=colored_bundle_file_out [-low]\n";
	print "To get colors:\n\n";
	print "       $0 -colorsRGB [or -colorsRGBAsHTMLTable] [-low]\n";
	print "or\n";
	print "       $0 -colorsHex [-low]\n";
	print "or\n";
	print "       $0 -colorsHexAsHTMLKeyTable [-low]\n";
	exit(0);
}

sub exportRGB {
	if ($low)
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGBLow)
		{
			print "$threshold\t$bundleColorsRGBLow{$threshold}\n";
		}
	}
	else
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGB)
		{
			print "$threshold\t$bundleColorsRGB{$threshold}\n";
		}
	}
	exit(0);
}

sub exportRGBHTMLTable {
	print "<table border=1>\n";
	print "<tr><th>Bundled Links</th><th>RGB</th></tr>\n";
	if ($low)
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGBLow)
		{
			my $cellBgColor = rgbToHex($bundleColorsRGBLow{$threshold});
			print qq{<tr><td>$threshold</td><td bgcolor="$cellBgColor">$bundleColorsRGBLow{$threshold}</td></tr>\n};
		}
	}
	else
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGB)
        	{
			my $cellBgColor = rgbToHex($bundleColorsRGB{$threshold});
                	print qq{<tr><td>$threshold</td><td bgcolor="$cellBgColor">$bundleColorsRGB{$threshold}</td></tr>\n};
        	}
	}
	print "</table>\n";
	exit(0);
}

sub exportAsHTMLKeyTable {
        print "<table border=1>\n";
        print "<tr><th>Bundled Links</th></tr>\n";
        if ($low)
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGBLow)
		{
			my $cellBgColor = rgbToHex($bundleColorsRGBLow{$threshold});
			print qq{<tr><td bgcolor="$cellBgColor">&nbsp;<span style="color:white">$threshold</span></td></tr>\n};
		}
	}
	else
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGB)
        	{
                	my $cellBgColor = rgbToHex($bundleColorsRGB{$threshold});
                	print qq{<tr><td bgcolor="$cellBgColor">&nbsp;<span style="color:white">$threshold</span></td></tr>\n};
        	}
	}
        print "</table>\n";
	exit(0);
}
sub rgbToHex {
	my ($r, $g, $b) = split/,/, $_[0];
	return sprintf ("#%2.2X%2.2X%2.2X", $r, $g, $b);
}

my $bundleFileIn = "";
my $bundleFileOut = "";

(@ARGV) or usage();
foreach my $arg (@ARGV)
{
	($arg =~ m/^-(h|help)$/) and usage();
	($arg =~ m/^-low$/) and $low = 1;
	($arg =~ m/^-colorsRGB$/) and exportRGB();
	($arg =~ m/^-colorsRGBAsHTMLTable$/) and exportRGBHTMLTable();
	($arg =~ m/^-colorsHexAsHTMLKeyTable$/) and exportAsHTMLKeyTable();
	($arg =~ m/^-i=(.+)$/) and $bundleFileIn = $1;
	($arg =~ m/^-o=(.+)$/) and $bundleFileOut = $1;
}

open(OUT, ">$bundleFileOut") or die "ERROR: can not open bundle out file $bundleFileOut $!\n";
open(IN, "<$bundleFileIn") or die "ERROR: can not open bundle in file $bundleFileIn $!\n";
while (my $line = <IN>)
{
#ASB_LG19 13470754 14218750 Ss262 2177839 2976275 nlinks=672,bsize1=150447,bsize2=150419,bidentity1=0.201133,bidentity2=0.188392,depth1=0,depth2=0,
#ASB_LG19 14250080 15061508 Ss262 1303606 2191377 nlinks=1076,bsize1=279892,bsize2=278553,bidentity1=0.344937,bidentity2=0.313766,depth1=0,depth2=0,
#ASB_LG19 14314359 14314420 Ss262 7198136 7198167 nlinks=9,bsize1=62,bsize2=32,bidentity1=1.000000,bidentity2=1.000000,depth1=1,depth2=1,
#ASB_LG19 15064224 15625360 Ss262 672993 1254783 nlinks=881,bsize1=305520,bsize2=304727,bidentity1=0.544466,bidentity2=0.523774,depth1=0,depth2=0,
#ASB_LG19 15650721 16282135 Ss262 8995 672359 nlinks=786,bsize1=199405,bsize2=198505,bidentity1=0.315807,bidentity2=0.299239,depth1=0,depth2=0,
#ASB_LG19 17026943 17042421 Ss262 965 7848 nlinks=35,bsize1=7610,bsize2=4363,bidentity1=0.491634,bidentity2=0.633788,depth1=0,depth2=0,
	chomp $line;
	my @data = split/\s+/, $line;
	my @bundleFields = split/,/, $data[6];
	my ($label, $count) = split/=/, $bundleFields[0];
	my $colorText = "color=(168,191,0)";
	if ($low)
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGBLow)
		{
			if ($count > $threshold)
			{
				$colorText = "color=($bundleColorsRGBLow{$threshold})";
			}
		}
	}
	else
	{
		foreach my $threshold (sort {$a <=> $b} keys %bundleColorsRGB)
		{
			if ($count > $threshold)
			{
				$colorText = "color=($bundleColorsRGB{$threshold})";
			}
		}
	}
	my $newline = join(" ", $data[0], $data[1], $data[2], $data[3], $data[4], $data[5], $colorText, $data[6]);
	select OUT; print OUT "$newline\n";
}
close(OUT);
exit(0);
