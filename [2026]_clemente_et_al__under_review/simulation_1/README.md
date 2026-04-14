### Simulation 1, Case A, Case B 

Once you have runned the container:
```{bash}
cd $HOME/simulation_1/
./compile.sh PDE_SOLUTION.cpp && ./PDE_SOLUTION.exe
./compile.sh GENERATE_DATA.cpp && ./GENERATE_DATA.exe
rm GENERATE_DATA.exe PDE_SOLUTION.exe
chown -R 1000:1000 input/
```
The `GENERATE_DATA.exe` executable generates all the data needed to run Simulation 1 and Cases A and B.

#### Simulation 1
To compile and run Simulation 1:
```{bash}
./compile.sh simulation_1.cpp

N=250
reac=1.00 
for sim in $(seq 0 29); do
./simulation_1.exe $N $sim $reac 
done   
```

**TPS** To run the simulations related to TPS, execute run_tps.sh from a terminal outside the container. The Docker image does not include the R package mgcv, which is required.

#### Case A
To compile and run case A:
```{bash}
./compile.sh case_a.cpp

N=250
for reac in 0.20 0.40 0.60 0.80 1.00; do
for sim in $(seq 0 29); do
  ./case_a.exe $N $sim $reac
done
done  
```

#### Case B
To compile and run case B:
```{bash}
./compile.sh case_b.cpp

N=250
reac=1.00
for N in 125 250 500 1000 2000; do
for sim in $(seq 0 29); do
  ./case_b.exe $N $sim $reac
done
done  
```


**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.
