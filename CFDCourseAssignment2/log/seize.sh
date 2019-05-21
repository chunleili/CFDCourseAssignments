rm seize.log
for Mx in 100 200 300
do
	for My in 100 200 300
	do
		tail -n 10 log${Mx}*${My} >>seize.log 
	done
done
sed -i '/has been/'d seize.log
sed -i 's/Converged!/**********************************************/g' seize.log
