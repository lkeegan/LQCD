beta=0.44068679351
for L in 8 16 32 64
do
	time ../release/MC_Ising $L ${beta} 1e3 1e4 50 10 123 > L${L}.txt
done
