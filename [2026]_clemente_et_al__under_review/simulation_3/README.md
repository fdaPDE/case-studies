### Simulation 3 and Simulation 4
First, if you have not generated the data for Simulation 1, run the following inside the container:
```{bash}
cd $HOME/simulation_1/
if [ ! -d "input/" ]; then
  ./compile.sh PDE_SOLUTION.cpp && ./PDE_SOLUTION.exe
fi
chown -R 1000:1000 input/
```

Outside the container, run
```{bash}
Rscript make_incidence_matrix.R
```

Finally, inside the container, run:
```{bash}
cd $HOME/simulation_3/
./compile.sh GENERATE_DATA.cpp && ./GENERATE_DATA.exe
rm GENERATE_DATA.exe
```
The `GENERATE_DATA.exe` executable generates all the data needed to run Simulation 3 and Simulation 4 avaialble in the Supplementary Material.

#### Simulation 3
To compile and run Simulation 3:
```{bash}
./compile.sh simulation_3.cpp

for sim in $(seq 0 29); do
./simulation_3.exe $sim 
done

rm simulation_3.exe
chown -R 1000:1000 simulation_3/  
```

#### Simulation 4
To compile and run Simulation 4:
```{bash}
./compile.sh simulation_3-not_unif.cpp

for sim in $(seq 0 29); do
./simulation_3-not_unif.exe $sim 
done   

rm simulation_3-not_unif.exe
chown -R 1000:1000 simulation_3-not_unif/
```

**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.

