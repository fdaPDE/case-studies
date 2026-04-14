### Modeling Spatio-Temporal Smoothness using  Nonlinear Differential Operators

This repository contains the code necessary to reproduce the results presented in *Modeling Spatio-Temporal Smoothness Using Nonlinear Differential Operators*.  

- [`simulation_1/`, `simulation_2/`, and `simulation_3/`]  
  - Contain the C++/R files needed to reproduce the results of the simulation studies presented in the main paper.  

- [`case_c/` and `case_d/`]  
  - Contain the C++/R files needed to reproduce the results of the additional simulation studies (Cases C and D) presented in the supplementary material.  

Finally, the code to reproduce Simulation Study 4, presented in the supplementary material, is available in `simulation_3/`. The code to reproduce the results for Cases A and B of the additional simulation studies is available in `simulation_1/`.  

**Note**  
The Docker image `aldoclemente/fdapde-docker:latest` is provided to compile and run the simulation studies reported here.  
Run the container using `./run_container.sh`. Then refer to the README files in the specific subfolder of interest.
