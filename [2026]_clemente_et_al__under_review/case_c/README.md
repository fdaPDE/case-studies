### Case C

Once you have runned the container:
```{bash}
cd $HOME/case_c/
./compile.sh case_c_DATA.cpp

for M in 5, 11, 21, 41; do 
    ./case_c_DATA.exe --M $M
done
rm case_c_DATA.exe
chown -R 1000:1000 input/
```

To compile and run Case C:
```{bash}
./compile.sh case_c.cpp

for M in 5 11 21 41; do
for sim in $(seq 0 29); do
  ./case_c.exe 250 $sim 1.00 $M
done
done

rm case_c.exe
chown -R 1000:1000 case_c/  
```

**Note** 
Since the simulations can take several hours, it is recommended to run them on a cluster machine.  
Once outside the container, you can collect the results and generate qualitative and quantitative outputs using the R script `results.R`.

