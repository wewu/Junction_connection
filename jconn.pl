#!perl
                              
$fname = $ARGV[0];
$maxexon=$ARGV[1];

check_param();

my @farr = read_file($fname);                                       
my @farr1 = step1_pre_filter(@farr);

my @farr1 = read_file("$fname.step1");
my @farr2 = step2_jconn(@farr1);



## ============== step functions ================

sub comp_coverage
{
    my ($chr, $beg, $end, @carr) = @_;

    my $clen = @carr;
    my $covlen = 0; 

    my $cbegline = get_cov_beg_line($beg, $end, $clen, @carr);
    for (my $k=$cbegline; $k<$clen; $k++)
    {
	my ($cchr, $cbeg, $cend, $cnum) = read_cov_line($carr[$k]);

      if ($cchr eq $chr) 
      {
	  my $ext = 0;
	  if (($cbeg <= $beg) && ($cend >= $beg) && ($cend <= $end))    { $covlen = $covlen + $cend - $beg ;  $ext = 1;}
	  elsif (($cbeg >= $beg) && ($cbeg <= $end) && ($cend >= $end)) { $covlen = $covlen + $end - $cbeg ;  $ext = 1;}
	  elsif (($cbeg <= $beg) && ($cend >= $end))                    { $covlen = $covlen + $end - $beg ;   $ext = 1;}
	  elsif (($cbeg >= $beg) && ($cend <= $end))                    { $covlen = $covlen + $cend - $cbeg ; $ext = 1;}

         if ($ext == 1) 
         {
	     my $tbeg = $beg+50;
	     my $tend = $end - 50;
            #show("$tbeg,$tend, case=($case1,$case2,$case3,$case4), $cbeg, $cend, $covlen");
         }
      }
      if ($cbeg > $end) 
      { 
	  $k = $clen + 1; 
      }
    }
    my $totlen = $end - $beg + 1;
   #show("totlen=$totlen");
    my $perc = round100($covlen / $totlen);
    $perc;
}

sub step2_jconn
{
    my @farr1 = @_;
    my $flen1= @farr1;

    my @garr, $glen = 0;
    my @uarr=0, $ulen=$flen1;
    my @rarr, $rlen=0;
    my $add, $beg;
    my $hgaplen = 50, $tgaplen = 50;

    my @covarr = read_file("$fname.cov");
    my $covlen = @covarr;
    my $covcutoff = 0.95;

    for (my $u=0; $u<$ulen; $u++)
    {
	$uarr[$u] = 0;
    }

    $beg = get_next_begin(@uarr);   
   while ($beg >= 0)
   {
       $garr[0] = $beg;
       $glen = 1 ;
       $uarr[$beg] = 1;

       for (my $g=0; $g<$glen; $g++)
       {
	   my $cidx = get_last_idx($garr[$g]);
	   my ($cchr, $cbeg, $cend, @tmparr) = read_bed_line($farr1[$cidx]);
	   $cbeg = $cbeg + $hgaplen;
	   $cend = $cend - $tgaplen;
	   $add = 0;

         ## find all next lines
	   for (my $j=$cidx+1; $j<$flen1; $j++)
	   {
	       my ($nchr, $nbeg, $nend, @tmparr) = read_bed_line($farr1[$j]);
	       $nbeg = $nbeg + $hgaplen;
	       $nend = $nend - $tgaplen;

	       my $gaplen = $nbeg - $cend;
            if (($gaplen <= $maxexon) && ($gaplen > 0))
            {
		my $hasone = has_between($cidx, $j, $hgaplen+$tgaplen);
               if ($hasone == 0)
               {
		   my $covper = comp_coverage($nchr, $cend, $nbeg, @covarr);
                  #show("covper=$covper, $g, $j,  $cname, $nchr, $cend, $nbeg");
                  if ($covper > $covcutoff)
                  {
		      $garr[$glen] = $garr[$g].",$j";
		      $glen++;
		      $uarr[$j] = 1;
		      $add = 1;
                  }
               }
            }
            elsif ($gaplen > $maxexon)
            {
		$j=$flen1;
            }
	   }

	   if ($add == 0) { $rarr[$rlen] = $garr[$g]; $rlen++; }
       }

       $beg = get_next_begin(@uarr);
      #show("beg=$beg, rlen=$rlen");
   }

   #write_arr_to_file("$fname.step1.idx", @rarr);
   #print "Output file: $fname.step1.idx\n";

    my @retarr, $k=0;
    for (my $r = 0; $r<$rlen; $r++)
    {
	my @iarr = split(/\,/, $rarr[$r]);
	my $ilen = @iarr;

	my ($schr, $sbeg, $send, $sname, $sscore, $sstrnd, $stbeg, $stend, $srgb, $sbkcnt, $sbksz, $sbkend) = read_bed_line($farr1[$iarr[0]]);
	my ($echr, $ebeg, $eend, @tmparr) = read_bed_line($farr1[$iarr[$ilen-1]]);

	my $wchr = $schr;
	my $wbeg = $sbeg;
	my $wend = $eend;
	my $num=$r+1;
	my $wname= "j".$num;
	my $wscore=$sscore;
	my $wstrnd=$sstrnd;
	my $wbeg = $sbeg;
	my $wend = $eend;
	my $wrgb = $srgb;
	my $wbkcnt = $ilen + 1;

	my $wbksz = "";

	my $bklen = $send - $sbeg - $hgaplen;
	my $wbkend = "0,".$bklen;
	for (my $i=1; $i<$ilen; $i++)
	{
	    my ($nchr, $nbeg, $nend, @tmparr) = read_bed_line($farr1[$iarr[$i]]);
	    my $bksz = $nbeg - $send + $hgaplen + $tgaplen;
	    $wbksz = $wbksz.",".$bksz;

	    $bklen = $nend - $sbeg - $hgaplen;
	    $wbkend = $wbkend.",".$bklen;

	    $send = $nend;
	}

	$wbksz = $hgaplen.$wbksz.",".$tgaplen;

	$retarr[$k] = comb_col($wchr, $wbeg, $wend, $wname, $wscore, $wstrnd, $wbeg, $wend, $wrgb, $wbkcnt, $wbksz, $wbkend);
	$k++;
    }

    write_arr_to_file("$fname.step2", @retarr);
    print "Junction connection output file: $fname.step2\n";

    @retarr;
}

sub step1_pre_filter
{
    my (@farr) = @_;
    my $flen = @farr;

   # Codification of initial filtering to remove low scoring stuff that does not seem
   # like it should be part of a transcript:
   #
   # Score <= 2, unknown, overlaps with something with score >= 20.
   # unless it is a perfect alternate splice form without being ridiculously long

    my %junction_num;
    my %kept_junctions;
    my $scorefilter = 3;
    my $scorefilter_max = 15;
    my %junctions_file;
    my %starts;
    my %ends;

    for (my $i=0; $i<$flen; $i++)
    {
	$line = $farr[$i];
	chomp($line);
	my @a = split(/\t/,$line);
	$junction_num{$a[0]}=$junction_num{$a[0]}+0;
	$junctions_file{$a[0]}[$junction_num{$a[0]}][0] = $a[1]+50;
	$junctions_file{$a[0]}[$junction_num{$a[0]}][1] = $a[2]-50;
	$junctions_file{$a[0]}[$junction_num{$a[0]}][2] = $a[4];
	$junctions_file{$a[0]}[$junction_num{$a[0]}][3] = $line;
      if($a[4] >= $scorefilter) 
      {
	  $starts{$a[0]}{$a[1]+50}=1;
	  $ends{$a[0]}{$a[2]-50}=1;
      }
	$junction_num{$a[0]}++;
    }

    my @cnter_arr, $cnt_idx = 0;
   foreach my $chr (keys %junction_num) 
   {
       my $start = 0;
       my $kept_counter=0;
       for(my $i=0; $i<$junction_num{$chr}; $i++) 
       {
	   if($junctions_file{$chr}[$i][2] > $scorefilter) 
	   {
            # keep becuase it exceeds the minimum score cutoff for unqualified inclusion
	       $kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
	       $kept_counter++;
	   } 
	   elsif($starts{$chr}{$junctions_file{$chr}[$i][0]}+0==1 && $ends{$chr}{$junctions_file{$chr}[$i][1]}+0==1) 
	   {
            # keep because it shares start and end with junctions that exceed the minimum
            # score cutoff for unqualified inclusion
	       $kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
	       $kept_counter++;
	   } 
	   elsif($junctions_file{$chr}[$i][3] =~ /24,116,205/ || $junctions_file{$chr}[$i][3] =~ /16,78,139/) 
	   {
            # keep regardless of score because it's a known junction
	       $kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
	       $kept_counter++;     
	   } 
         else 
         {
            # Going to see if it overlaps with a junction of score above the min, if not then
            # we will keep it.
	     my $flag = 0;
	     for(my $j=$start; $j<$junction_num{$chr}; $j++) 
	     {
		 if($junctions_file{$chr}[$j][2] >= $scorefilter_max 
		    && $junctions_file{$chr}[$j][0] <= $junctions_file{$chr}[$i][1] 
		    && $junctions_file{$chr}[$j][1] >= $junctions_file{$chr}[$i][0]) 
		 {
		     $flag = 1;
		 }
		 if($junctions_file{$chr}[$i][1] < $junctions_file{$chr}[$j][0]) 
		 {
                  # we've checked everything that can overlap, set $j to jump out of loop
		     $j = $junction_num{$chr};
		 }
	     }
            if($flag == 0) 
            {
		$kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
		$kept_counter++;
            }
	     while($junctions_file{$chr}[$start][1] < $junctions_file{$chr}[$i][0] && $start < $junction_num{$chr}) 
	     {
		 $start++;
	     }        
         }
       }
       $cnt_arr[$cnt_idx] = $kept_counter;
       $cnt_idx ++;
   }
   
    my $idx = 0;
    my @retarr, $k=0;
   foreach my $chr (keys %kept_junctions) 
   {
       for (my $i=0; $i<$cnt_arr[$idx]; $i++)
       {
	   $retarr[$k] = $kept_junctions{$chr}[$i];
	   $k++;
       }
       $idx ++;
   }
   #write_arr_to_file("$fname.step1", @retarr);
   #print "Preprocessing output file $fname.step1\n";

    return @retarr;
}




sub get_cov_beg_line
{
    my ($beg, $end, $clen, @carr) = @_;

    my $first = 0;
    my $last = $clen;
    my $half = int(($first + $last)/2);
    my ($cchr, $cbeg, $cend, $cnum) = read_cov_line($carr[$half]);

    for (my $loop = 0; $loop<=10; $loop++)
    {
      if ($cend < $beg)
      {
	  $first = $half; 
      }
      elsif ($cbeg > $end)
      {
	  $last = $half; 
      }
      else
      {
	  $loop = 100;
      }
      $half = int(($first + $last)/2);
      ($cchr, $cbeg, $cend, $cnum) = read_cov_line($carr[$half]);
  }
    $first;
}


## ============== small functions ================
sub check_param
{
   if (length($fname) == 0)
   {
       print "Usage: input_file(in BED format)   maxexon_length(default 2000)\n";
       exit;
   }
                                
   if (length($maxexon) == 0)
   {
       $maxexon=2000;
   }
}


sub get_next_begin
{
    my (@uarr) = @_;
    my $ulen = @uarr;

    my $ret = -1;
    for (my $u=0; $u<$ulen; $u++)
    {
      if ($uarr[$u] == 0)
      {
	  $ret = $u;
	  $u = $ulen;
      }
  }
    $ret;
}

sub has_between
{
    my ($beg, $end, $gaplen) = @_;

    my $found = 0;
    my ($bchr, $bbeg, $bend, @tmparr) = read_bed_line($farr[$beg]);
    my ($echr, $ebeg, $eend, @tmparr) = read_bed_line($farr[$end]);

    for (my $t = $beg; $t<$end; $t++)
    {
	my ($mchr, $mbeg, $mend, @tmparr) = read_bed_line($farr[$t]);

      if  ((($bend-$mbeg) < $gaplen ) && (($mend - $ebeg)<$gaplen))
      {
	  $found = 1;
	  $t = $end;
      }
    }
    $found;

}

sub get_last_idx
{
    my ($path) = @_;

    my @parr = split(/\,/,$path);
    my $plen = @parr;
    $parr[$plen-1];
}

sub read_bed_line
{
    my ($ln) = @_;

    my $lncur = $ln;
    chomp($lncur);
    $lncur =~  s/\r//g;

    my @larr = split(/\t/, $lncur);
    my $chr = $larr[0];
    my $beg = $larr[1];
    my $end = $larr[2];
    my $name =$larr[3];
    my $score=$larr[4];
    my $strnd=$larr[5];
    my $tbeg =$larr[6];
    my $tend =$larr[7];
    my $rgb  =$larr[8];
    my $bkcnt=$larr[9];
    my $bksz =$larr[10];
    my $bkend=$larr[11];
    ($chr, $beg, $end, $name, $score, $strnd, $tbeg, $tend, $rgb, $bkcnt, $bksz, $bkend);
}

sub read_cov_line
{
    my ($ln) = @_;

    my $lncur = $ln;
    chomp($lncur);
    $lncur =~  s/\r//g;

    my @larr = split(/\t/, $lncur);
    my $chr = $larr[0];
    my $beg = $larr[1];
    my $end = $larr[2];
    my $cov = $larr[3];
    ($chr, $beg, $end, $cov);
}

sub read_file
{
    my ($fname) = @_;
   
    open (Fmsg, "< $fname" ) || print_and_append_log("\t\tCan't open $fname file!\n");
    my @fmsg = <Fmsg>;
    my $len = @fmsg;

    for (my $i=0; $i<$len; $i++)
    {
	chomp($fmsg[$i]);
    }
    close Fmsg;
    @fmsg; 
}


sub create_empty_file
{
    my ($fname) = @_;
    open (Fmsg, "> $fname" ) || print_and_append_log("Can't create empty $fname file!\n");
    close Fmsg;
}

sub write_file
{
    my ($fname, $msg) = @_;
    open (Fmsg, "> $fname" ) || print_and_append_log("Can't write $fname file!\n");
   if (length($msg)>0)
   {   
       print Fmsg "$msg\n";
   }
    close Fmsg;
}

sub write_arr_to_file
{
    my ($fname, @cont) = @_;
    my $lencont = @cont;
    my $str="";   
   
    for (my $i=0; $i<$lencont-1; $i++)
    {
	$str = $str.$cont[$i]."\n";
    }
   if ($lencont > 0)
   {
       $str = $str.$cont[$lencont-1];
   }
    &write_file($fname, $str);
}

sub append_file
{
    my ($fname, $msg) = @_;
    open (Fmsg, ">> $fname" ) || print_and_append_log("Can't append $fname file1!\n");
   if (length($msg)>0)
   {   
       print Fmsg "$msg\n";
   }
    close Fmsg;
}

sub append_arr_to_file
{
    my ($fname, @cont) = @_;
    my $lencont = @cont;
    my $str="";   
   
    for (my $i=0; $i<$lencont-1; $i++)
    {
	$str = $str.$cont[$i]."\n";
    }
   if ($lencont > 0)
   {
       $str = $str.$cont[$lencont-1];
   }
    &append_file($fname, $str);
}

sub create_empty_file
{
    my ($fname) = @_;
    open (Fmsg, "> $fname" ) || print_and_append_log("Can't create empty $fname file!\n");
    close Fmsg;
}


sub show
{
    my ($str) = @_;
    print "\n$str\n";
}

sub show_arr
{
    my ($name, @arr) = @_;
    my $len = @arr;
   
    print "$name\t";
    for (my $k=0; $k<$len; $k++)
    {
	print "$k->$arr[$k]; ";
    }
    print "\n";
}

sub sort_array_str
{
    my (@arr1) = @_;
   
    my @arr2 = sort {lc $a cmp lc $b} @arr1; 
    @arr2;
}

sub sort_array_num
{
    my (@arr1) = @_;
   
    my @arr2 = sort {$a <=> $b} @arr1; 
    @arr2;
}

sub comb_arr
{ 
    my ($sep, @arr) = @_;
    my $str = "";
    my $len = @arr;

   if ($len > 0)
   {
       $str = $arr[0];
   }
    for (my $y=1; $y<$len; $y++)
    {
	$str = $str.$sep.$arr[$y];
    }
    $str;
}


sub comb_col
{ 
    my (@arr) = @_;
    my $str = comb_arr("\t", @arr);
    $str;
}

sub round
{
    my ($val) = @_;
   if ($val > 0)
   {
       $val = int($val + 0.5);
   }
   else
   {
       $val = int($val - 0.5);
   }
    $val;
}



sub round100
{
    my ($val) = @_;
    $val = &round($val*100)/100;

    my ($h, $t) = split(/\./, $val);
    $t = &pad0($t, 2, "r");   # t return with ".t" if it has no . before
   if ($t =~ /\./)
   {
       $val = $h.$t;   
   }
   else
   {
       $val = $h.".".$t;   
   }
    $val;
}

sub pad0
{
    my ($val, $len, $dir) = @_;

   if (($dir eq "r") && (!($val =~ /\./)) && (($len - length($val)) > 1))
   {
       $val = $val.".";
   }
    my $curlen = length($val);
   while ($curlen < $len)
   {
      if ($dir eq "l")
      {
	  $val = "0".$val;
      }
      elsif ($dir eq "r")
      {
	  $val = $val."0";
      }
      $curlen++;
  }
    $val;
}

sub min
{
    my ($v1, $v2) = @_;

    my $ret;
   if ($v1 < $v2)
   {
       $ret = $v1;
   }
   else
   {
       $ret = $v2;
   }
    $ret;
}

sub max
{
    my ($v1, $v2) = @_;
    my $ret;
   if ($v1 > $v2)
   {
       $ret = $v1;
   }
   else
   {
       $ret = $v2;
   }
    $ret;
}


