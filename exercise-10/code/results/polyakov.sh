T=4
for L in 4 8 12
do
for beta in 1.00 2.00 2.20 2.26 2.28 2.30 2.32 2.34 2.40 2.60 3.00 4.00
do
	echo "$L	$beta..."
	../release/MC_SU2 $T $L $beta 5e2 1e4 1 5 > T${T}_L${L}_beta${beta}.txt
done
done
