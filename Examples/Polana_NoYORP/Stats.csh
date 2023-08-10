#!/bin/csh

#set num = $argv[1]
set num = `awk 'NF ==2 {print $2}' out.out | tail -n1`
set initn = `awk 'NF ==2 {print $2}' out.out | head -n1`
set initnn = `echo " " | awk '{print '$initn' -1 }' `
echo $num total final particles

tail -n $num out.out > final.out

awk '$2 == -8 {i=i+1} END {print "Collision: ",i}'  final.out
awk '$2 == -6 || $2 == -7 {i=i+1} END {print "Inward: ",i}'  final.out
awk '$2 == -6 {i=i+1} END {print "Inward nu6: ",i}'  final.out
awk '$2 == -7 {i=i+1} END {print "Inward J7:2: ",i}'  final.out
awk '$2 == -3 {i=i+1} END {print "Outward: ",i}'  final.out
awk '$2 > 2  {i=i+1} END {print "Alive: ",i}'  final.out
awk '$2 > 2 && $10 == 0 {i=i+1} END {print "Alive0: ",i}'  final.out
awk '$2 > 2 && $10 > 0 {i=i+1} END {print "AliveColl: ",i}'  final.out
awk '$2 == -6 && $10 == 0 {i=i+1} END {print "Inward0: ",i}'  final.out
awk '$2 == -6 && $10 > 0 {i=i+1} END {print "InwardColl: ",i}'  final.out
awk '$2 == -7 && $10 == 0 {i=i+1} END {print "Inward J7:2: ",i}'  final.out
awk '$2 == -7 && $10 > 0 {i=i+1} END {print "InwardColl J7:2: ",i}'  final.out




awk '$2 == -6 {print $0}' final.out > final.inward.out
awk '$2 == -7 {print $0}' final.out > final.inwardJ72.out
awk '$2 > 0 {print $0}' final.out > final.alive.out
awk '$2 == -3 {print $0}' final.out > final.outward.out
awk '$2 == -8 {print $0}' final.out > final.collision.out

awk '$10 == 0 {print $0}' final.inward.out  > final.inward.orig.out
awk '$10 == 0 {print $0}' final.inwardJ72.out  > final.inwardJ72.orig.out

awk '$10 > 0 {print $0}' final.inward.out  > final.inward.coll.out
awk '$10 > 0 {print $0}' final.inwardJ72.out  > final.inwardJ72.coll.out

sort -g -k7 -r final.inward.out | awk '{print $7,NR}' > final.inward.SFD.out

head -n $initn out.out | tail -n $initnn | sort -g -k7 -r | awk '{print $7,NR}' > init.SFD.out

#Crossover point
set tstep = `grep tstep Delivery_Drift_Coll.pl | head -n1 | awk '{print $3}' |grep -o '[0-9]\+'`
set tmax = `grep tmax Delivery_Drift_Coll.pl | head -n1 | awk '{print $3}' | grep -o '[0-9]\+'`

#echo $tstep $tmax
set tt = $tstep
set go = `echo " " | awk ' '$tt' < '$tmax' {print "1"} '$tt' >= '$tmax' {print "0"}'`

set cross = 0
set count = 0

echo -n "Crossover "

while ($tt < $tmax)
    set first = `awk '$12 == 1 && ($8) == '$tt' {print $0}' final.inward.orig.out final.inwardJ72.orig.out | wc | awk '{print $1}'`
    set nth = `awk '$12 > 1 && $8 == '$tt' {print $0}' final.inward.coll.out final.inwardJ72.coll.out  | wc  | awk '{print $1}'`
    #echo $tt $tstep $tmax $first $nth $cross $count

    if ($nth > $first) then
	if ($count == 2) then
	    echo -n $cross $amt
		break
	    endif
	if ($count == 1) then
	    set count = `echo " " | awk '{print '$count'+1}'`
	endif
	if ($count == 0) then
	    set cross = $tt
	    set amt = $nth
	    set count = `echo " " | awk '{print '$count'+1}'`
	endif
    endif

    if ($first > $nth) then
	if ($count > 0) then
	    set count = `echo " " | awk '{print 0}'`
	endif
    endif

    set tt = `echo " " | awk '{print '$tt'+'$tstep'}'`
    set go = `echo " " | awk ' '$tt' < '$tmax' {print "1"} '$tt' >= '$tmax' {print "0"}'`

end

echo  " "


set npt = `grep npt Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '.[0-9]\+'`

set rlr = `grep RadLargeRem Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '.[0-9]\+'`

set rlf = `grep RadLargeFrag Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '.[0-9]\+'`

set max = `grep MaxSize Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '.[0-9]\+'`

set sfdslope  = `grep FamSFDslope Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '[0-9].[0-9]\+'`

set ainit = `grep ainit Delivery_Drift_Coll.pl | head -n3 | awk '{print $3}' | grep -o '[0-9].[0-9]\+'`


set drift = `grep dadtbennu Delivery_Drift_Coll.pl | head -n2 | tail -n1 | awk '{print $3}' | grep -o '[a-zA-Z]\+'`

#echo "Params " $tstep $npt $rlr $rlf $max $sfdslope $ainit $drift
