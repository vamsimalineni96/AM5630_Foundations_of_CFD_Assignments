# Assignment 2:
### Initial Problem Definition:
  * Number of X grid points : 31
  * Number of Y grid points : 41
  * Thermal Conductivity : 380 W/m-degC
  * Thermal Diffusivity : 11.234 x 1e-5 m^2/s
  * Length : 0.3 m
  * Width : 0.4 m
 
![image](https://user-images.githubusercontent.com/98683842/165436678-86221601-8375-49ae-873d-bc8d3b08ee54.png)
  
### Problem statements :
1. `Problem statement 1:` _2D Transient Heat Conduction Equation_   using 
    * `FTCS` method
    * `ADI` method

2. `Problem statement 2:` _2D Steady State Heat Conduction Equation_   using
    * `Point Gauss Seidal` method
    * `Line Gauss Seidal` method
    * `PSOR` method
    * `LSOR` method
    * `ADI` method 

3. `Problem statement 3:` _Effect of Symmetry_ 
    * Most effective method obtained from problem 2 must be employed to a `changed boundary conditions` for the same computational domain
    * The same method has to be applied to  a `reduced computational domain of a quadrant`
   The numerical formulations :
      * `PSOR` method
      * `LSOR` method
      * `ADI` method

### General Instructions
Install the required libraries using the following command :
``` python
pip install -r requirements.txt
```
