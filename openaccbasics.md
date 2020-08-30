# Getting started with OpenACC

## The compiler: NVIDIA HPC Software Development Kit

We will use the NVIDIA HPC SDK (PGI compiler) that can be downloaded from: [https://www.pgroup.com](https://www.pgroup.com). All information and resources can be found in: [https://www.pgroup.com/resources/accel.htm](https://www.pgroup.com/resources/accel.htm).

This compiler allows to compile both C and Fortran codes using ```pgcc``` and ```pgfortran```, respectively. For instance, to compile a serial C program we type in the console:

```
$ pgcc myprogram.c -o myprogram
```

Among the most useful compilation flags, we can find:

- ```-mp```: to compile OpenMP pragmas.
- ```-acc```: to compile OpenACC pragmas.
- ```-ta=multicore```: to run OpenACC pragmas in CPU (multi-thread).
- ```-ta=tesla:managed```: to use the managed memory.

We can also use  ```-Minfo=accel``` to see on the screen some information about the parts of the code that have been accelerated.

## Using the unified (managed) memory

When using CUDA unified memory, both GPU and CPU memories are combined into a single pool. To tell the compiler to use the unified memory, we must use the flag ```-ta=tesla:managed```. There are three important operations in the CPU/GPU memories to be considered:

- *Array allocation*: we will use the ```enter data create()``` drective.
```
int* A=(int*)malloc(N*sizeof(int));
#pragma acc enter data create(A[0:N])
```
where N is the expected size of A.
