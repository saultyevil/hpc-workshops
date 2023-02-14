# OpenMP Notes

## Parallel Regions

To use OpenMP to make a parallel part of code, we need to use a directive which the compiler will see and know to use OpenMP to parallelise the block of code. We need to use the -fopenmp flag for the gcc or gfortran compiler.

```fortran
!$OMP PARALLEL
parallel block
!$OMP END PARALLEL
```

```c
#pragma omp parallel
{
parallel block
}
```

The above syntax for C will accept any C object block, hence we do not need the curly braces if we are expecting to parallelise a single region.

To find the number of threads being used, we need to import some useful functions.

```fortran
USE OMP_LIB
INTEGER FUNCTION OMP_GET_NUM_THREADS()
```

```c
#include <omp.h>
int omp_get_num_threads(void);
```

Inside a parallel region, variables can either be *shared* or *private*.

```fortran
SHARED(list)
PRIVATE(list)
DEFAULT(SHARED|PRIVATE|NONE)
```

```c
shared(list)
private(list)
default(shared|none)
```

On entry to a parallel region, private variables are uninitialised. Variables declared inside the scope of the parallel region are automatically private. After the parallel region ends, the original variable is unaffected by any changes to private copies. **Not specifying a `DEFAULT` clause is the same as specifying `DEFAULT(SHARED)`. This is dangerous and `DEFAULT(NONE)` should always be used.**

```fortran
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(I, MYID),
!$OMP& SHARED(A, N)
myid = omp_get_thread_num() + 1
DO i = 1, n
a(i, myid) = 1.0
END DO
!$OMP END PARALLEL

```

If we declared n are being private, it would be uninitialised in the private region and hence would contain trash. If we want to initialise a private variable at the start of a parallel region, we use the `FIRSTPRVIATE` or `firstprivate` clause. This is generally not done, but below is a 'contrived' example of when it can be used.

```c
b = 23.0;
..........
#pragma omp parallel firstprivate(b), private(i,myid)
{
myid = omp_get_thread_num();
for (i=0; i<n; i++)
{
b += c[myid][i];
}
c[myid][n] = b;
}
```

## Reductions

In the below example, `n` has not been declared as being a shared variable. The dude fell into his own trap! `n` should have been declared as being a shared variable.

```fortran
b = 10
!value in original variable is saved
! each thread gets a private copy, initialised to 0
!$OMP PARALLEL REDUCTION(+:b),
!$OMP PRIVATE(I, MYID)
myid = omp_get_thread_num() + 1
! all accesses inside the parallel region are to the
! private copies
DO i = 1, n
b = b + c(i, myid)
END DO
!$OMP END PARALLEL

! At the end of the parallel region, all the private copies
! are added into the original variable
a = b
```

## Work sharing

### Parallel Loops

To call a parallel loop, use the following syntax.

```fortran
!$OMP DO [clauses]
do loop
[!$OMP END DO]
```

```c
#pragma omp for [clauses]
for loop
```

Remember that what is being calculated has to be independent for parallel loops to work.

```fortran
!$OMP PARALLEL
!$OMP DO
DO i=1,n
b(i) = (a(i)-a(i-1))*0.5
END DO
!$OMP END DO
!$OMP END PARALLEL
```

```c
#pragma omp parallel
{
#pragma omp for
for (int i=0; i<n; i++)
{
b[i] = a[i] * a[i-1] * 0.5;
}
}

```

These directives are used so commonly that a shorthand notation has been added.

```fortran
!$OMP PARALLEL DO [clauses]
do loop
[!$OMP END PARALLEL DO]
```

```c
#pragma omp parallel for [clauses]
for loop
```

`DO`/`FOR` directives can take `PRIVATE`, `FIRSTPRIVATE` and `REDUCATION` clauses which refer to the scope of the loop. Note that the parallel loop index variable is private by default. `PARALLEL` `DO`/`FOR` directives can take all clauses available to for the `PARALLEL` directive. With no additonal clauses, the `DO`/`FOR` directive will partition the iterations as equally as possible between the threads. However, this is implementation dependent, and there is still some ambiguity, e.g. 7 iterations and 3 threads could do 3+3+1 or 3+2+2.

### `SCHEDULE`

The `SCHEDULE` clause allows us to specify which loop iterations are executed by which thread. `SCHEDULE` has multiple options to split up the loop iterations.

* `STATIC` - best for load balanced loops, least overhead
* `STATIC, n` - best for mild of smooth load imbalance, but can induce overheads
* `DYNAMIC` - useful if iterations have widely varying loads, but ruins data locality
* `GUIDED` - often less expensive than `DYNAMIC`, but beware of loops where the first iterations are the most expensive
* `AUTO` - may be useful if the loop is executed many times over

### `SINGLE`

The `SINGLE` directive indicates that a block of code should be indicated only by a single thread. The first thread to reach the `SINGLE` directive will execute the block of code and the rest of the threads will skip over it, but all other blocks will wait for this thread to finish executing the `SINGLE` block.

### `MASTER`

Indicates that a block of code should be executed by the master thread (thread 0) only. There is no synchronisation at the end of this block, the other threads will skip over this block and continue.

## Synchronisation

### `BARRIER`

No thread can proceed past a barrier until all treads have arrived.

```fortran
!$OMP PARALLEL PRIVATE(I,MYID,NEIGHB)
    myid = omp_get_thread_num()
    neighb = myid - 1
    if (myid.eq.0) neighb = omp_get_num_threads()-1
    ...
    a(myid) = a(myid)*3.5
    !$OMP BARRIER
    b(myid) = a(neighb) + c
    ...
!$OMP END PARALLEL
```

The barrier forces synchronisation on `a`.

### `CRITICAL`

A critical section is a block of code which can only be executed by one thread at a time. Critical directives can be named.

```fortran
!$OMP CRITICAL [( name )]
    block
!$OMP END CRITICAL [( name )]
```

```c
#pragma omp critical [( name )]
    structured block
```

```fortran
!$OMP PARALLEL SHARED(STACK),PRIVATE(INEXT,INEW)
    ...
    !$OMP CRITICAL (STACKPROT)
        inext = getnext(stack)
    !$OMP END CRITICAL (STACKPROT)
        call work(inext,inew)
    !$OMP CRITICAL (STACKPROT)
        if (inew .gt. 0) call putnew(inew,stack)
    !$OMP END CRITICAL (STACKPROT)
    ...
!$OMP END PARALLEL
```

### `ATOMIC`

Used to protect a single update to a shared variable. Applies only to a single statement. Using `ATOMIC` is generally more efficient than using `CRITICAL`.

`CRITICAL` protects access to a block of code, whereas `ATOMIC` protects access to a block of memory.

```fortran
!$OMP ATOMIC
    ! the statement must be a binary operation
    statement
```

```c
#pragma omp parallel for
for (j=0; j<nedges; j++)
{
    #pragma omp atomic
    degree[edge[j].vertex1]++;
    #pragma omp atomic
    degree[edge[j].vertex2]++;
}
```

### `LOCK`

A lock is a special variable that may be *set* by a thread. No other thread may *set* the lock until the thread which set the lock has *unset* it. Setting a lock can either be block or non-blocking. A lock must be initialised before it is used, and may be destroyed when it is no longer required. Lock varibales should not be used for any other purpose.

```fortran
USE OMP_LIB
SUBROUTINE OMP_INIT_LOCK(OMP_LOCK_KIND var)
SUBROUTINE OMP_SET_LOCK(OMP_LOCK_KIND var)
LOGICAL FUNCTION OMP_TEST_LOCK(OMP_LOCK_KIND var)
SUBROUTINE OMP_UNSET_LOCK(OMP_LOCK_KIND var)
SUBROUTINE OMP_DESTROY_LOCK(OMP_LOCK_KIND var)
```

```c
#include <omp.h>
void omp_init_lock(omp_lock_t *lock);
void omp_set_lock(omp_lock_t *lock);
int omp_test_lock(omp_lock_t *lock);
void omp_unset_lock(omp_lock_t *lock);
void omp_destroy_lock(omp_lock_t *lock);
```

```c
for (i=0; i<nvertexes; i++)
{
    omp_init_lock(lockvar[i]);
}

#pragma omp parallel for
for (j=0; j<nedges; j++)
{
    omp_set_lock(lockvar[edge[j].vertex1]);
    degree[edge[j].vertex1]++;
    omp_unset_lock(lockvar[edge[j].vertex1]);
    omp_set_lock(lockvar[edge[j].vertex2]);
    degree[edge[j].vertex2]++;
    omp_unset_lock(lockvar[edge[j].vertex2]);
}
```

## Nested Parallelism

Nested parallelism is supported by OpenMP.

```fortran
!$OMP PARALLEL PRIVATE(myid)
myid = omp_get_thread_num()
IF (myid .eq. 0) THEN
    !$OMP PARALLEL DO
    DO i = 1,n
        x(i) = 1.0
    END DO
ELSE IF (myid .eq.1) THEN
    !$OMP PARALLEL DO
    DO j = 1,n
            y(j) = 2.0
    END DO
END IF
!$OMP END PARALLEL
```

To control the number of threads, we can change the environment variable. It takes a list of values, `export OMP_NUM_THREADS=2,4`. This will use 2 threads at the outer level and 4 threads for each of the inner teams. We can also use the `omp_set_num_threads()` or `num_threads` clause on the parallel region.

We can also control the maximum number of threads running at any one time, `export OMP_THREAD_LIMIT=64`. We can also set the maxmimum depth of nesting, `export OMP_MAX_ACTIVE_LEVELS=2` or call `omp_set_max_active_levels()`. If the thread limit is reached, the code will not stop running. Instead, OpenMP will just limit the amount of threads which the new parallel region will get.

## Nested Loops

```c
#pragma omp parallel for collapse(2)
for (int i=0; i<N; i++)
{
    for (int j=0; j<Ml k++){
        ...........
    }
}
```

The above only works for perfectly nested retagular loops. We ca parallelise multiple loops in the nest with the collapse clause. The argument is the number of loops to collapse start from the outide. It will form a single loop of length `N`x`M` and then parallelise and schedule that. Useful if `N` is O(no. of threads) so parallelising the outer loop may not have have good load balance. More efficient than using nested teams.

## Data scoping rules

When we call a subroutine from inside a parallel region,

* Variables in the argument list inherit their data scope attribute from the calling routine.
* Global varibales in C, and `COMMON` blocks or module variables in Fortran are shared, unless delcared `THREADPRIVATE`.
* `static` local variables in C and `SAVE` variables in Fortran are shared.
* All other local variables are private.

## Timing routines

OpenMP supports a portable timer. Return current wall clock time (relative to arbitrary origin) with, `DOUBLE PRECISION FUNCTION OMP_GET_WTIME()` or, `double omp_get_wtime(void);`. Return clock precision with, `DOUBLE PRECISION FUNCTION OMP_GET_WTICK()` or, `double omp_get_wtick(void);`.

## Tasks

Tasks are independent units of work. They are composed of code to execute and data to compute with. Threads are assigned to perform the work of each task.

The task construct includes a structured block of code. Inside a parallel region, a thread encuntering a task construct will package up the code block and its data for execution. Tasks can be nested: i.e. a task may itself generate tasks.

```fortan
!$OMP TASK [clauses]
    structured block
!$OMP END TASK
```

```c
#pragma omp task [clauses]
{
    structured blocked
}
```

```c
#pragma omp parallel
{
    // thread 0 packages tasks
    #pragma omp master
    {
        // tasks executed by other threads
        #pragma omp task
        fred();
        #pragma omp task
        daisy();
        #pragma omp task
        billy();
    }
}
```

Some thread in the parallel region will execute the task at some point in the order. *The tasks are run in a random order?* Very strange.

Tasks are affected by barrier regions. We can also use the `TASKWAIT` directive. This will cause a task to wait until all taasks defined in the current task have completed.

```c
#pragma omp parallel
{
    #pragma omp master
    {
        #pragma omp task
            fred();
        #pragma omp task
            daisy();
        #pragma taskwait
        #pragma omp task
            billy();
    }
}
```

In the above example, `fred()` and `daisy()` must complete before `billy()` starts.

```c
#pragma omp parallel
{
    #pragma omp master
    {
        p = listhead;
        while (p)
        {
            // makes a copy of p when the task
            // is packaged
            #pragma omp task firstprivate(p)
            {
                process(p);
            }
            p=next (p);
        }
    }
}
```

```c
#pragma omp parallel
{
    // all threads package tasks
    #pragma omp for private(p)
    for(int i =0; i <numlists; i++)
    {
        p = listheads[i];
        while(p)
        {
            #pragma omp task firstprivate(p)
            {
                process(p);
            }
            p = next(p);
        }
    }
}
```

### Data scoping with tasks

Variables can be shared, private or firstprivate with respect to the task. These concepts are a little bit different compared with threads.

* If a variable is shared on a task construct, the references to it inside the construct are to the storage with that name at the point where the task was encountered
* If a variable is private on a task construct, the references to it inside the construct are to new uninitialized storage that is created when the task is executed
* If a variable is firstprivate on a construct, the references to it inside the construct are to new storage that is created and initialized with the value of the existing storage of that name when the task is encountered

### Parallel Recursive Fibonacci

```c
int fib(int n)
{
    // x, y are local and thus private to
    // current task
    int x,y;
    if ( n < 2 ) return n;

    #pragma omp task shared(x)
        x = fib(n-1);

    #pragma omp task shared(y)
        y = fib(n-2);

    // a task cannot complete until all tasks
    // below it in the tree are complete
    // (enforced with taskwait)
    #pragma omp taskwait
        return x+y;
}

int main(void)
{
    int NN = 5000;
    #pragma omp parallel
    {
        #pragma omp master
        fib(NN);
    }
}
```

## Node Architecture

There are four principle technologies which make up HPC systems.

* Processors
* Memory
* Interconnect
* Storage

### Processors

Processors execute instructions to perform arithmetic operations. They load data from memory and store data to memory and decide which instructions to do next. The clock speed of a processor determines the rate at whch instructions are executed. Integer and floating point calculations can be done in parallel. The peak FLOP rate is just clock rate x no. of floating point operartions per clock cyce.

#### Instruction unit

* Responsible for fetching, decoding and dispatching of instructions.
* Fetches instruction from instruction caches
* Decodes instruction.
* Sends the instructions to the appropriate unit.
* May also be responsible for scheduling instructions (see later).

#### Integer unit

* Handles integer arithmetic
* Integer addition, multiplication and division.
* Logical ops (and, or, shift etc.)
* Also known as arithmetic and logic unit (ALU)

#### Floating point unit

* Handles floating point arithmetic
* Addition, multiplication, division.
* Usually the critical resource for HPC
* Machines sold by peak flop/s

#### Control unit

* Responsible for branches and jumps

#### Load/store unit

* Responsible for loading data from memory and storing it back.

#### Register file

* Local storage in the CPU
* Accessed by name (not address)
* Separate files for integers/addresses and floating point

### Pipelining

Key implementation technique for making fast CPUs. The execution of instructions is broken down into *stages*. Each stage can be executed in one CPU clock cycle. Once a stage has completed for one instuctions, it can be exeucted for the next instruction on the subsequent clock cycle. Allows one instruction to be completed per clock cycle, even though the instruction may take many clock cycles to complete. However, not all instructions are pipelined, i.e. FP square root or FP divide.

### Problems with pipelining

Any of the following can result in stopping and restarting the pipeline, and wasting cycles as a result:

* Two instructions both require the same hardware resource at the same time.
* One instruction depends on the result of another instruction further down the pipeline.
* The result of the instruction changes which instruction to execute next.

### Overcoming pipeline hazards

* Out-of-order execution: assembly code specifies an order of instructions. But the hardware chooses to reorder the instructions as they are fetched to minimise pipeline stalls. Requires some complex book keeping to ensure correctness.
* Branch preduction: the hardware tries to guess which way the next branch will go. Uses a hardware table that tracks the outcomes of recent branches in the code. Keeps the pipeline going and only stalls if the prediction is wrong.

### Instruction level parallelism

Pipeling is a form of instruction level parllelism. Multiple instructions are "in-flight" at the same time, but the maximum performance is 1 instruction per cycle.

There are two main approaches:

* Superscalar processor: parallel instructions identified at run-time, in hardware.
* SIMD (or vector) instructions: operations on multiple data items encoded in a single instruction.

### SIMD

It exists. I didn't listen that closely to what he was saying.

### Multicore chips

It's now commonplace to have multiple processors on a chip. The main difference is that processors may share caches. Typically, each core has its own Level 1 and Level 2 caches, but Level 3 cache is shared between cores. Most of a processor chip is Level 3 cache.

### SMT - simultaneous multithreading (a.k.a. Hyperthreading)

Tries to fil spare slots by mising instructions from more than one thread in the same clock cycle. Requires some replication of hardware, i.e. some hardware to store the instruction sets. Everything else is shared between the threads, such as functional units, register file, memory system (including the caches).

### Memory

Memory speed is often the limiting factor HPC applications. Keeping the CPU fed with data is the key to performance. Memory is a substatintial contributor to the cost of systems. Typically, a node will have a few Gbytes of memory per processor. It's technically possibly to have much more memory than this, but it's too expensive and power hungry.

Latency of 100s of nanoseconds and a few Gbytes/s bandwidth. Memory latencies are very long, on the order of 100s of processor ccles. Fetching data from main memory is 2 orders of magnitude slower than doing arithmetics. The solution is to use Cache memory. This memory is much faster than main storage memory, but also much smaller. Cache keeps copies of recently used memory.

### Caches

Sophisticated book keeping required. We need to know where the most up-to-date copy is kept and flush out any old data downloads to make space for new data. Caches only help performance if the application re-uses recently asscessed data. It is mostly the programmer's responsibility to order computations so as to maximise the re-use of Cache memory. Many modern memory systems also do prefetching, which makes (quite simple) guesses as to which data will be used next. This data will then be put into the chaches before the processor requests them.

#### When to cache

Always cache on reads, except in special circumstances. If a memory location is read and there isn't a copy in the cache (read miss), then cache the memory.

## Measuring performance

It's not enough to measure overall speed. You need to know where the time is going. We can use profiling tools so we can profile our code and figure out where all of the time is being spent. Many codes perform differently depending on the different input data. We have to use multiple sets of input data when measuring performance. Make sure that this data is representative of the problems where you want the code to run quickly.

### Code profiling

Code profiling is the first step for anyone interesting in performance optimisation. Profiling works by instruementing code at compile time, thus it's usually controlled by compiler flags and it can reduce performance.

### Standard Unix profilers

The standard profilers are `prof` and `gprof`. They are usually flags which you add at compile time. When the code is run, it produces an instrumentation log. We can then read the log file by using a tool to make it human readable.

### The golden rules of profiling

* Profile your code: the compiler will not do all the optimisation for you.
* Profile your code yourself: don't believe what anyone else tells you.
* Profile on the hardware you want to run on: don't profile on your laptop if you plan to run on ARCHER.
* Profile your code running the full-sized problem: the profile will almost certainly be qualitatively different for a test case.
* Keep profiling your code as you optimise
* Many proposed changes will turn out not to be useful: use some type of version control to easily revert back to the previous code.
* Create a good test framework

## OpenMP Performance

### Avoiding synchronisation overheads

Parallelise at the outermost level possible. This will minimise the frequency of barriers and may required reordering of loops and/or array indices. We can also use careful use of `NOWAIT` clause, which supresses the implicit barrier. However, doing this results it in being easy to remove barriers which are required for correctness. Atomics *may* have less overhear than citical or locks.

### Scheduling

If we create computational tasks, and realy on the runtime to assign these to the threads, then we incur some overheads. For example, non-static loop schedules will increase the overhead.

## Acknowledgements 

I would like to acknowledge financial support from the EPSRC Centre for Doctoral Training in Next Generation Computational Modelleing grant EP/L015382/1.
