#!/usr/bin/perl 
#use strict;
use Math::Trig;


#### Simulation parameters
$tstep = 10;  # in Myr
$tmax = 500;

# collisional cascade 
$RadLargeRem = 0.7937;  # largest remnant: Bottke 0.5M -> 0.79 in size
$RadLargeFrag = 0.7937; # largest fragment.
$RadSecLarge = $RadLargeFrag; #deprecated terminology
$FamSFDslope = -3.5;



# simulation resolution
$npt = 1000;
$SmallLimit = 0.5;  # km don't go below this in the loop below
$MaxNew = 1999;
$MinSize = $SmallLimit;
$MaxSize = 80;


#family and asteroid properties
$pV = 0.045;
$ainit = 2.4189; #Polana is at 2.4189

#Yarko drift properties
$dadtbennu_measured = 19.0e-4; #au/Myr; Drift rate from Spotto + Bolin's TI adjustment
$dadtbennu = $dadtbennu_measured; 
$bennudiam = 0.492;
$bennua    = 1.12;

#counters and internals - code check
$toadd=0;
$random_number = rand();
print "Test Random: $random_number \n";


#File opening.
open(OUT,">out.out") || die("Could not open file out.out \n");

#time zero to output file
print OUT "0 $npt \n";

## Loop over initial particles to initialize key values.
for ($i=0;$i<$npt;$i++)  {

    #diameter and albedo
    $D[$i]=getD();  #km
    $H[$i] =  2.5*(6.244 - log($pV) - 2*log($D[$i]));

    #select the pole. Ymult is deprecated.
    $polgood[$i]=acos(1-2*rand());
    $polsin[$i]=cos($polgood[$i]);
    $Ymult[$i]=1.0;
    
    #set counters. last collision. when formed. which generation.
    $col[$i] = 0;
    $form[$i] = 0;
    $gen[$i] = 1;
    $restime[$i] = 0;
    $a[$i] = $ainit;
    $J72[$i] = 0;
    print OUT "$i $a[$i] $H[$i] 0 0 $Ymult[$i] $D[$i] $polsin[$i] $J72[$i]\n";
}

######### BIG LOOP    ########
for ($t=0;$t<$tmax;$t=$t+$tstep) {
    print OUT "$t $npt \n";
    open(TEMP,">temp.temp") || die("Could not open file temp.temp \n");
    $npt = $npt + $toadd;
    $toadd = 0;
    print "TIME ",$t," ",$npt,"\n";
    for ($i=0;$i<$npt;$i++) {
	
	if ($a[$i] > 0){

	    $adot = (sqrt($bennua)/sqrt($a[$i]))*(($bennudiam/$D[$i]))*$dadtbennu;
	    $da =  $tstep*$adot*$polsin[$i]*$Ymult[$i];
      	    $a[$i]=$a[$i]+($da);
		
	    #check if object is destroyed
	    $destroy = getColl_Ballouz($tstep,$D[$i]);
	    if ($destroy>0) {
		#track destroyed objects index
		$origindex=$i;
		#how many new things to get added.
		$newnumb=getNumbNew($D[$i]);
		#print "NewNumb: ",$newnumb,"\n";
		for($j=($npt+$toadd);$j<($npt+$toadd+$newnumb);$j++){

		    getNewD($j,$D[$i],($j-($npt+$toadd)),$origindex);
		    print "Adding ",$j," ",$D[$i]," ",$D[$j]," ",($j-$npt)," ",$a[$i]," ",$form[$j]," ",$gen[$j],"\n";

		}

		$col[$i] = $t;
		$a[$i]=-8;
		$restime[$i] = $t;

		# If we change npt right here... that will be a pain.
		# So we save that number for later.
		$toadd = $toadd + $newnumb;
		
	    }
	}
	if ($a[$i] < 2.15 && $a[$i] > 0) {
	    $a[$i]=-6;
	    $restime[$i]=$t;
	}
	if ($a[$i] < 2.2569 && $a[$i] > 0 && $J72[$i] < 1) {
	    #one time only to roll the dice on J72
	    $J72[$i] = 1;
	    #Now, roll the dice. Bottke 2015 says its a 1/3rd chance
	    if (rand() < 0.33333333) {
		$a[$i]=-7;
		$restime[$i]=$t;
	    }
	}
	if ($a[$i] > 2.49 && $a[$i] > 0) {
	    $a[$i]=-3;
	    $restime[$i]=$t;
	}
	
       
	print OUT "$i $a[$i] $H[$i] $da $adot $Ymult[$i] $D[$i] $restime[$i] $col[$i] $form[$i] $polsin[$i] $gen[$i]\n";
	print TEMP  "$i $a[$i] $H[$i] $da $adot $Ymult[$i] $D[$i] $restime[$i] $col[$i] $form[$i] $polsin[$i] $gen[$i]\n";;
	
    }


    $outps=sprintf("%05d.ps",$t);
    $outpng=sprintf("%05d.png",$t);
    $outjpg=sprintf("%05d.jpg",$t);

    close TEMP;

}
close OUT;
    print "FINISH ",$npt,"\n";


#### seeing QstarD_Bennu_Equations_MKS in BallouzBoulder folder.
sub getColl_BallouzMKS {
    my $ts = $_[0]; #timestep in year?
    my $di = $_[1]; #dimaeter in km
    my $mu_g=0.33;

    my $U = 5000.0; # m/s MBA impact speed
    my $rho_target = 2000.0; #kg/m**3
    my $rho_imp = 2000.0;
    my $Rw = 80.0; #meters - Otohime.
    
    my $rt = $di*1000.0/2.0; # into MKS

#    my $first = 
    $QstarD = (3.45e2)*$rt**(-0.47) + (2.72e-1)*$rt**(0.99);

}

sub getColl_Ballouz {
    my $ts = $_[0]; #timestep
    my $di = $_[1]; #dimaeter in km
    my $rtarg = $di*100000.0/2.0; #cm and diameter
    my $rtarg_m = $rtarg/100.0;
    my $colverb = 0;
    
    #print "Entering getColl_Ballouz ",$ts," ",$di," ",$rtarg," ",$rtarg_m,"\n";
    
    ###Target Density - need to document in paper
    my $rho_g=1.19; #g/cm3
    my $rho_s=2.2; #g/cm3

    my $mu_g=0.33;

    ### mu_s, we find that it is 0.47 +/- 0.02.
    my $mu_s=0.47;


    #### Q1_func
    my	$G=6.67e-8;
    my	$K1=0.12;
    my	$Kr=1.2;
    my	$mu=$mu_g;
    my  $Rw=8e3; #Otohime radius
    my $ratio_g = 1.0;
    my $delta = 3.0; #gcc
    my $U_mba= 2.5e5; #cm/ss
    ### Q0_func
    my $b=3.0*$mu_g;
    my $a= (-1.0)*$mu_s;

    ### QstarRecon

    

    my $Q1 = 0.5*(($Kr*$K1**(1.0/3.0)**(3.0*(2.0+$mu)/2.0)))*(($G)**(3.0*$mu/2.0))*((4.0*3.14159/3.0)**($mu-1.0))*($rho_g**(5.0*$mu/2.0))*($delta**(-$mu))*($ratio_g**(3.0*(2.0+$mu)/2.0));
    my $Q0 = -$Q1*$b*($Rw**($b-$a))/($a*$U_mba**(3*$a+$b));

    my $QstarD = ($Q0*($rtarg**$a)*($U_mba**(2.0-(3.0*$mu_s))))+($Q1*($rtarg**$b)*($U_mba**(2.0-(3.0*$mu_g))));

    my $impactor_radius_mba = (($QstarD*($rho_s/$delta)*(2.0/($U_mba)**2))**(1.0/3.0)*$rtarg)/100.0;

    if ($colverb > 1) {
	print "Rtarg Q1 ",$rtarg," ",$Q1,"\n";
	print "Rtarg Q0 ",$rtarg," ",$Q0,"\n";
	print "Rtarg QstarD ",$rtarg," ",$QstarD,"\n";
	print "Rtarg imp_rad_mba ",$rtarg," ",$impactor_radius_mba,"\n";
    }

    
    my $alpha = -2.564813222553626;
    my $C1=13998520270481.287; ##m2.8
    my $Dip1=$impactor_radius_mba*(10.0**0.1);
    
    my $Pi_mba = 2.9e-24; #per m2 per year

    my $N_mba_law = (($C1*$impactor_radius_mba**(-2.8))-($C1*$Dip1**(-2.8)));
    
    my $np_mba = ($N_mba_law*$Pi_mba*(($impactor_radius_mba)+($rtarg/100.0))**2)/2.0;
    
    #lifetime_mba from simplelifetime
    $lifetime = (1.0/$np_mba/1.0e6);
    #####
    
    if ($colverb > 1) {
	print "N_mba_law ",$impactor_radius_mba," ",$N_mba_law,"\n";
	print "Np_mba ",$rtarg," ",$np_mba,"\n";
    }
    if ($colverb > 0){
	$farinella = 16.8*($rtarg_m)**(1.0/2.0); #Myr and meters
	print "Farinella v Ballouz: ",$rtarg_m," ",$farinella," ",$lifetime,"\n";
	#print "Ron Lifetime ",$rtarg," ",$lifetime,"\n";
    }
    
    if ($lifetime < $tstep) { # where we need to do fancy stats
	# HACK
	return 1;
    } elsif ($lifetime > $tstep) {
	$prob=($ts/$lifetime); # prob for destruction during this interval
	$f=rand();
	if ($prob > $f) {
	    print "remove ",$i," ",$di," ",$lifetime," ",$prob," ",$f," ",$t,"\n"; 
		return 1;
	} else {
	    return 0;
	}
    } else {
	return 0;
    }

}

#farinella based lifetime estimates - for comparison.
sub getColl {
    my $ts = $_[0];
    my $di = $_[1];

    # lifetime = 16.8 Myr * sqrt(R)
    $lifetime = 16.8*($di*500.0)**0.5;
    #print "getColl ",$di," ",$lifetime," ",$ts/$lifetime,"\n";
    if ($lifetime < $tstep) { # where we need to do fancy stats
	return 1;
    } elsif ($lifetime > $tstep) {
	$prob=($ts/$lifetime); # prob for destruction during this interval
	$f=rand();
	if ($prob > $f) {
	    print "remove ",$i," ",$di," ",$lifetime," ",$prob," ",$f," ",$t,"\n"; 
		return 1;
	} else {
	    return 0;
	}
    } else {
	return 0;
    }    
}


##################################################
# Here we want to be able to designate a SFD
#Durda(pg5): 10000 of 2km  1 20km  Super Catastrophic SLOPE = 
# x=0.30103, y=4   x=1.30103 y=0  SLOPE = -4  CUMULATIVE
#Durda(pg5):
### km? YES
#This function gets some variables from the top global settings.
#MinSize
#MaxSize
#It also gets slope from FamSFDslope

sub getD { 
#    $r = rand()*18 + 1;  # sizes from 1-20 km
    
    $small=$MinSize;
    $big=$MaxSize;
    $q=$FamSFDslope;   ### THIS IS THE SLOPE!!!!
    $a=$small**($q);
    $b=$big**($q);
    $f=rand();
    $rd = ((1-$f)*$a + $f*$b)**(1/$q);
#    print "$rd $f $q $a $b \n";
    if ($rd <= 0) {print "Diam = 0"; exit;}
    return $rd;
}


sub getNumbNew {
    my $initD = $_[0];
    my $k=0;
    my $freshD = 0;
    my $count = 0;
    my $massconserve = 0;
    #Question is size of 2nd largest.
    # Answered below. I use the RadSecLarge* the bro$RadLargeRem*$initD;

    #Bootstrap to start with Largest Remnant
    $freshD=$RadLargeRem*$initD;
    $massconserve = ($massconserve**(3.0) + $freshD**(3.0))**(1.0/3.0);
    if ($freshD < $SmallLimit) {
	print "FreshD (RadLargeRem) is below limit - maxed out new guys at ",$count," ",$initD," ",$freshD," ",$massconserve,"\n";
	return $count;
    }

    ## If the first remnant is large enough.. we keep it and start counting.
    $count++;

    for($k=1;$k<$MaxNew;$k++) {

	$freshD = 10**((log(($k+1))/log(10))/($FamSFDslope) + log($RadLargeFrag*$RadLargeRem*$initD)/log(10));
	if ($freshD > $SmallLimit) {
	    $count++;
	    $massconserve = ($massconserve**(3) + $freshD**(3))**(1.0/3.0);
	} else {
	    print "FreshD (RadLargeFrag) is below limit - maxed out new guys at ",$count," ",$initD," ",$freshD," ",$massconserve,"\n";
	    return $count;
	}

	if ($massconserve > $initD) {
	    print "massconserve limit met getNumbNew ",$initD," ",$freshD," ",$count," ",$massconserve,"\n";
	    return $count;
	}
    }
    print "getNumbNew ",$initD," ",$freshD," ",$count,"\n";

    return $count;
}

    
sub getNewD {
    my $index = $_[0];
    my $initD = $_[1];
    my $cindex = $_[2];
    my $oi = $_[3];

    # Bottke2002 style
    # Largest Remnant first.
    # Then make the largest fragment and teh slope off of it.
    if ($cindex == 0) {
	$D[$index] = $RadLargeRem*$initD;
    } else {
	$D[$index] = 10**((log(($cindex+1))/log(10))/($FamSFDslope) + log($RadLargeFrag*$RadLargeRem*$initD)/log(10));
    }
	
    $H[$index] =  2.5*(6.244 - log($pV)/log(10) - 2*log($D[$index])/log(10));
    print "getNewD ",$initD," ",$D[$index]," ",$index," ",$cindex,"\n";


    $polgood[$index]=acos(1-2*rand());
    $polsin[$index]=cos($polgood[$index]);
    $Ymult[$index]=1.0;

    
    $col[$index] = 0;
    $form[$index] = $t;
    $restime[$index] = 0;
    $a[$index] = $a[$oi];

    $gen[$index]=$gen[$oi]+1;
    #print OUT "$i $a[$i] $H[$i] 0 0 $Ymult[$i] $D[$i] $polsin[$i]\n";

}

# From DCR ssic.c
# q = p->surf_den_exp + 2;
# a = pow(p->r_inner,q);
# b = pow(p->r_outer,q);
#   f = ran(0);
# sma = pow((1 - f)*a + f*b,1/q);

############################################################

sub parse {
    my(@parsed);
    @parsed = split(/\s+/,$_[0]);
    chomp(@parsed);
  #remove whitespace at beginning of line, if any
    while ($parsed[0] eq "") {shift(@parsed)}
    return(@parsed);
}
