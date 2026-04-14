### Simulation 3 and Simulation 4
Outside the container run 
```{bash}
Rscript make_incidence_matrix.R
```

Once you have runned the container:
```{bash}
cd $HOME/simulation_3/
.compile.sh GENERATE_DATA.cpp && ./GENERATE_DATA.exe
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
```

#### Simulation 4
To compile and run Simulation 4:
```{bash}
./compile.sh simulation_3-not_unif.cpp

for sim in $(seq 0 29); do
./simulation_3-not_unif.exe $sim 
done   
```

**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.

