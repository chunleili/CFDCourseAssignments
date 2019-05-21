mkdir log
for Mx in 100 200 300
do
    for My in 100 200 300
    do
        sed -i "s/Mx=..., My=.../Mx=${Mx}, My=${My}/g" main.H
        make 
        ./out | tee log/log${Mx}*${My}
        mv data/residual.dat data/residual${Mx}*${My}.dat
        mv data/finalMesh.dat data/finalMesh${Mx}*${My}.dat
    done
done
