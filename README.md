# LeeKeslerCorrelation
In this project, I use C++ along with read-only csv files. The code features the use of regex iterators to traverse through the csv files. It also features a recursive function to recursively calculate pressure and find Z values. The code approximates the pressure at a given temperature for known critical pressure and critical temperature values. 

To run the code on Windows:
  ```
  g++ LeeKesler.cpp
  ./a.exe
  ```
  
To run the code on Linux:
  ```
  g++ LeeKesler.cpp
  ./a.out
  ```
  
You will then put in every numerical value that is prompted for in this order:
1. Critical Temperature in Kelvin
2. Critical Pressure in Bar
3. Fluid Omega Value
4. System Temperature in Kelvin
5. System Volume in cm^3/mol
