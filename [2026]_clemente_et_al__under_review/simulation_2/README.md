### Simulation 2

Once you have runned the container:
```{bash}
cd $HOME/simulation_1/
if [ ! -d "input/" ]; then
  ./compile.sh PDE_SOLUTION.cpp && ./PDE_SOLUTION.exe
fi
chown -R 1000:1000 input/

cd $HOME/simulation_2/
./compile.sh GENERATE_DATA.cpp && ./GENERATE_DATA.exe
rm GENERATE_DATA.exe
chown -R 1000:1000 input/
```
The `GENERATE_DATA.exe` executable generates all the data needed to run Simulation 2.
To compile and run Simulation 2:
```{bash}
./compile.sh simulation_2.cpp

for sim in $(seq 0 29); do
./simulation_2.exe 250 $sim 1.00 
done

rm simulation_2.exe
chown -R 1000:1000 simulation_2/  
```

**TPS** To run the simulations related to TPS, execute run_tps.sh from a terminal outside the container. The Docker image does not include the R package mgcv, which is required.


**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.

