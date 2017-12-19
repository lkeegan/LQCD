T=4
for L in 20
do
for beta in 2.295
do
	echo "$L	$beta..."
	../release/MC_SU2 $T $L $beta 1e3 1e5 5 10 > T${T}_L${L}_beta${beta}.txt
done
done
