for L in 8 16 24 32
do
	rm L${L}.txt
	for beta in 0.300 0.325 0.350 0.375 0.400 0.425 0.430 0.435 0.437 0.439 0.441 0.443 0.445 0.447 0.450 0.475 0.500 0.525 0.550 0.575 0.600
	do
		echo "Running ${L} - ${beta}.."
		../release/MC_Ising $L ${beta} 5e3 1e5 >> L${L}.txt
	done
done
