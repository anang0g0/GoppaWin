$str = shift || ' 8x^3+12x^2-2x-3 ';        #<----- 

$str =~ s/\s//g;
print factor($str)."\n";
sub factor{
    my @a = arrange(shift);     # 式の整理 次数でまとめる
    if((my $i = arrayGcd(@a)) > 1){
        for(@a){ $_ /= $i;}
        return ed2(ed($i).factor(array2expr(@a)));
    }
    return factor_div(@a) if($#a > 2);  #３次以上の時は因数定理？除算の定理？で求める
    return ed2("x(".ed($a[2])."x".ed($a[1]).")") unless($a[0]);
    my @a0 = measure($a[0],1);  # 約数を求める a0 = (b1,b2)...
    my @a2 = measure($a[2],0);  #   "          a2 = (c1,c2)...
    for (@a0){                  # (a0,a1,a2)=(b1,c1)(b2,c2)
        my ($b1,$b2) = split ",";   # a1=b1*c2+b2*c1となるb1,b2,c1,c2を求める
        for(@a2){
            my ($c1,$c2) = split ",";
            return ed2("(".ed($c1)."x".ed($b1).")(".ed($c2)."x".ed($b2).")")
                    if($a[1] == $b1*$c2+$b2*$c1);
            return ed2("(".ed($c2)."x".ed($b1).")(".ed($c1)."x".ed($b2).")")
                    if($a[1] == $b1*$c1+$b2*$c2);
        }
    }
    return "(".array2expr(@a).")";  # 数値を式に変換 -> a2 x^2 + a1 x + a0 
}
sub factor_div{             # 因数定理？　除算の定理？
    my @a = @_;
    my ($tmp,$tmp1,$k);
    for my $z (1..abs($a[$#a])){
        for(my $i=1;$i<=abs($a[0]);$i++){
            $tmp = $tmp1 = 0;
            $k = 0 - $i;
            for(my $j=0;$j<=$#a;$j++){
                my $z2 = ($j) ? $z : 1;
                $tmp += $a[$j]*($i**$j)/($z2**$j);
                $tmp1 += $a[$j]*($k**$j)/($z2**$j);
            }                                                      # 誤差調整の為、sprintfを追加
            return ed2("(".ed($z)."x".ed(0-$i).")".factor(array2expr(div($z,0-$i,@a)))) 
                                                                   unless(sprintf("%f",$tmp)+0);
            return ed2("(".ed($z)."x".ed(0-$k).")".factor(array2expr(div($z,0-$k,@a)))) 
                                                                   unless(sprintf("%f",$tmp1)+0);
        }
    }
    if($a[0]==0){
        shift @a;
        return ed2("x".factor(array2expr(@a)));
    }
    if(@a==5 and $a[1]==0 and $a[2]==0 and $a[3]==0 and
           $a[4] >= 0 and $a[4] == (int(sqrt($a[4])))**2){
        my @b = (@a[0],0,1);
        my $n=sqrt(@a[4]);
        my $ans = factor(array2expr(@b));
        if($n == 1){$n="";};
        if($ans =~ s/x([+-])/${n}x\^2$1/g){
            return $ans;
        }else{
            return "(".array2expr(@a).")";
        }
    }elsif(@a==5 and $a[1]==0 and $a[3]==0){
        my @b = @a[0,2,4];
        my $ans = factor(array2expr(@b));
        if($ans !~ /(\(.+\)){2,}/){            # X^2をXと置いてa^2-b^2=(a+b)(a-b)の公式にて ※1
            if($b[0] >= 0 and $b[2] == 1 and $b[1]+1 == sqrt($b[0])*2){ 
                return "(" . array2expr(sqrt($b[0]),1,1) . ")(" . array2expr(sqrt($b[0]),-1,1) . ")" ;
            }
            return "(" . array2expr(@a) . ")";
        }
        $ans =~ s/x/x^2/g;
        return $ans;
    }else{
        return "(".array2expr(@a).")";
    }
}
sub div{
    my ($x,$z,@a) = @_;
    my @b;
    for (my $i = $#a;$i>0;$i--){
            $b[$i-1]=$a[$i]/$x;
            $a[$i-1]-= $z*$b[$i-1];
    }
    return @b;
}
sub arrange{                # 式の整理
    my $str = shift;
    my @a = (0,0,0);
    my ($b,$c);
    while($str =~ s/([+-]*[\d\.]*)x\^([\d\.]*)//){
        $b=$2;
        $c=$1;
        $a[$b] += ($c =~ /^[+-]*$/) ? "$c" . "1" : $c;
    }
    while($str =~ s/([+-]*[\d\.]*)x//){
        $c = $1;
        $a[1] += ($c =~ /^[+-]*$/) ? "$c" . "1" : $c;
    }
    while($str =~ s/([+-]*[\d\.]{1,})//){
        $a[0] += $1;
    }
    return @a;
}
sub measure{               # 約数を求める
    my ($a,$flag) = @_;
    my $n = abs($a);
    my @ans = ();
    for (my $j = int(sqrt($n)); $j > 0;$j--){
        unless ($n % $j){
            if($a<0){
                push(@ans,join(",",$j,$n/$j*(-1)));
                push(@ans,join(",",$j*(-1),$n/$j));
            }else{
                push(@ans,join(",",$j,$n/$j));
                push(@ans,join(",",$j*(-1),$n/$j*(-1))) if($flag);
            }
        }
    }
    return @ans;
}
sub array2expr{            # 数値を式に変換
    my @a = @_;
    my $new = "";
    for (my $i = $#a ; $i>=0;$i--){
        next if($a[$i] ==0);
        $new .= ed($a[$i]);
        if($i == 1){
            $new .= "x";
        }elsif($i > 1){
            $new .= "x^".$i;
        }
    }
    return ed2($new);
}
sub arrayGcd{
#    my @a = @_;
    my $ans = 0;
    for(@_){ $ans = gcd($ans,$_);}
    return abs($ans);
}
sub gcd{
    my ($x,$y) = @_;
    return $y ? gcd($y,$x % $y) : $x ;
}

sub ed{
    my $n = shift;
    return ($n >= 0) ? "+" . $n : $n;
}
sub ed2{
    my $str = shift;
    $str =~ s/([+-])1(x)/$1$2/g;
    $str =~ s/^\+//;
    $str =~ s/\(\+/\(/g;
    my $sub = 1;
    if ($str =~ s/xx(\^(\d+))/XX/){
        $sub = $2;
    }
    while($str =~ s/xx/X/i){$sub++;}
    $str =~ s/X/x^$sub/ ;
    return $str;
}
