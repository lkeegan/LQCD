for L in 40
do
	rm mL${L}.txt
	for beta in 0.300 0.325 0.350 0.375 0.400 0.425 0.435 0.445 0.455 0.465 0.475 0.500 0.525 0.550 0.575 0.600
	do
		echo "Running ${L} - ${beta}.."
		../release/MC_Ising $L ${beta} 1e4 2e6 1 0 123 >> mL${L}.txt
	done
done
