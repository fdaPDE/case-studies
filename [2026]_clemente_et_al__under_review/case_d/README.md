### Case D

Once you have runned the container:
```{bash}
cd $HOME/case_d/
./compile.sh case_d_DATA.cpp
for N in 11, 21, 31, 41; do 
    ./case_d_DATA.exe --N $N
done
rm case_d_DATA.exe
chown -R 1000:1000 input/
```

To compile and run Case D:
```{bash}
./compile.sh case_d.cpp

for N in 11, 21, 31, 41; do
for sim in $(seq 0 29); do
  ./case_d.exe 250 $sim 1.00 $N
done
done  

rm case_d.exe
chown -R 1000:1000 case_d/  
```

**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.
