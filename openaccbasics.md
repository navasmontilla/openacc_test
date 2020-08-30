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

## Using CUDA Unified (managed) memory

When using CUDA Unified memory, both GPU and CPU memories are combined into a single pool. To tell the compiler to use the unified memory, we must use the flag ```-ta=tesla:managed```. There are three important operations in the CPU/GPU memories to be considered:

- *Array allocation*: we will use the ```enter data create()``` directive.
```c
int* A=(int*)malloc(N*sizeof(int));
#pragma acc enter data create(A[0:N])
```
where N is the expected size of A. Note that ```A[0:N]``` can be replaced by ```A[:N]```.

- *Array deallocation*: we will use the ```exit data delete()``` clause.
```c
#pragma acc exit data delete(A)
free(A);
```

- *Copy from host to device (CPU -> GPU)*: we will use the ```update device()``` clause.
```c
#pragma acc update device(A[0:N])
```

- *Copy from device to host (GPU -> CPU)*: we will use the ```update self()``` clause.
```c
#pragma acc update self(A[0:N])
```

When accelerating a loop, if the shared variables in the loop have already been allocated and initialized/copied in the unified memory, they will be defined inside the ```present``` clause. For example, in this loop we use matrix A that already is in the managed memory:
```c
#pragma acc parallel loop present(A[0:n])  
        for( int j = 1; j < n; j++)
        { 
          ...
        }
```
## Loop optimization

When accelerating a loop, we can use *gangs* (coarse grain), *workers* and *vectors* (fine grain). Each gang may include several workers and each worker may include several vectors. Therefore, ehen having nested loops, outer loops will be accelerated as *gangs* and inner loops as *vectors*:
```c
#pragma acc parallel loop gang
        for( int i = 1; i < n; i++)
        { 
           #pragma acc loop worker
           for( int j = 1; j < m; j++)
           {
                #pragma acc loop vector
                for( int k = 1; k < c; k++)
                 {
                        ...
                 }
           }
        }
```
Generally, when having two nested loops we can just write:
```c
#pragma acc parallel loop
        for( int i = 1; i < n; i++)
        { 
           #pragma acc loop
           for( int j = 1; j < m; j++)
           {
                ...
           }
        }
```

- *Collapse*: when having nested loops, they can be converted into a single one using:
```c
#pragma acc parallel collapse(2)
        for( int i = 1; i < n; i++)
        { 
           for( int j = 1; j < m; j++)
           {
                ...
           }
        }
```
where 2 is the number of loops.

- *Tile*: when having nested loops, they can be redefined into n x m loops to explot data locality:
```c
#pragma acc parallel tile(32,32)
        for( int i = 1; i < n; i++)
        { 
           for( int j = 1; j < m; j++)
           {
                ...
           }
        }
```

