# Getting started with OpenACC

## The compiler: NVIDIA HPC Software Development Kit

We will use the NVIDIA HPC SDK (PGI compiler) that can be downloaded from: [https://www.pgroup.com](https://www.pgroup.com). All information and resources can be found in: [https://www.pgroup.com/resources/accel.htm](https://www.pgroup.com/resources/accel.htm).

This compiler allows to compile both C and Fortran codes using ```pgcc``` and ```pgfortran```, respectively. For instance, to compile a serial C program we can write:

```
pgcc myprogram.c -o myprogram
```

Among the most useful compilation flags, we can find:

- ```-mp```: to compile OpenMP pragmas.
- ```-acc```: to compile OpenACC pragmas.
- ```-ta=multicore```: to run OpenACC pragmas in CPU (multi-thread).
- ```-ta=tesla:managed```: to use the managed memory.

We can also use  ```-Minfo=accel``` to see on the screen some information about the parts of the code that have been accelerated.
