for L in 4 6 8 10
do
	rm cL${L}.txt
	for beta in 0.00 0.05 0.10 0.11 0.12 0.13 0.14 0.1425 0.145 0.1475 0.15 0.1525 0.155 0.1575 0.16 0.17 0.18 0.19 0.20 0.25 0.30
	do
		echo "Running ${L} - ${beta}.."
		../release/MC_Ising $L ${beta} 1e3 5e3 0 1 >> cL${L}.txt
	done
done
