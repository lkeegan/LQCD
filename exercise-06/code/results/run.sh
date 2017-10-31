for L in 4 6 8 10
do
	rm L${L}.txt
	for beta in 0.92 0.93 0.94 0.95 0.96 0.97 0.975 0.98 0.985  0.988 0.991 0.994 0.997 1.000 1.003 1.006 1.009 1.012 1.015  1.02 1.03 1.04 1.05
	do
		echo "Running ${L} - ${beta}.."
		../release/MC_U1 $L ${beta} 5e2 1e4 1 1 >> L${L}.txt
	done
done
