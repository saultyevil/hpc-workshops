# MPI Workshop

Textbook: Using MPI *Portable Parallel Programming with the Message-Passing Interface*: William Gropp, Ewing Lusk, etc

MPI is a library of function calls. It is not a language and there is no such thing as an MPI compiler. MPI's prime goals are to provide source-code portability and to allow efficient implementation. It also offers a great deal of functionality and support for heterogenous parallel architectures.

## MPI basics

### Header files

To use MPI, we need to use the MPI library.

```c
#include <mpi.h>
```

```fortran
use mpi
```

### MPI function format

Error terms are always returned by MPI functions. However, we can handle them differently between C and FORTRAN.

```c
error = MPI_Xxxxxx(parameter, ...);
MPI_Xxxxxx(parameter, ...);
```

```fortran
CALL MPI_XXXXXXX(parameter, ..., IERROR)
```

`IERROR` is optional in the 2008 version of Fortran, but is otherwise essential.

### Initialising MPI

To use MPI, we have to call the `MPI_INIT` function before we call any other MPI functions.

```c
int MPI_Init(int *argc, char ***argv)
```

```fortran
MPI_INIT(IERROR)
INTEGER IERROR
```

For example, in an actual program, this may look like

```c
int main(int argc, char *argv[])
{
    ...
    MPI_Init(&argc, &argv);
    ...
}
```

```fortran
PROGRAM my_mpi_program integer :: ierror

    ...
    CALL MPI_INIT(IERROR)
    ...

END PROGRAM my_mpi_program
```

### Rank

We use the function `MPI_COMM_RANK` to identify the different processes. This is known as the rank and is not the physical processor number. It is a logical number which goes from 0, ..., N-1.

```c
MPI_Comm_rank(MPI_Comm comm, int *rank)
```

```c
int rank;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
printf("Hello from rank %d\n", rank);
```

### Size

To find out how many processes are contained within a communicator we use the `MPI_COMM_SIZE` function.

```c
MPI_Comm_size(MPI_Comm comm, int *size);
```

### Exiting MPI

We have to finalise MPI every time we use it. This allows MPI to clean up and send any unsent messages or recieve any messages.

```c
int MPI_Finalize();
```

```fortran
MPI_FINALIZE(IERROR)
```

We can also abort MPI, but this is not recommended.

### Machine names

We can find out the hostname of the computer which we are on. This can be useful for debugging.

```c
int namelen;
char procname[MPI_MAX_PROCESSOR_NAME];

MPI_Get_processor_name(procname, &namelen);
printf("Rank %d is on machine %s\n", rank, procname);
```

## Messages in MPI

A message contains a number of elements of some particular datatype. MPI has two datatypes, basic types and derived types. Derived types can be built up from basic types. C types are different from Fortran types.

### Point-to-point communication

This is communication between two processes. The source process sends a message to the destination process. The communication takes place in a communicator. The destination process is identified by its rank in the communicator.

### Communication modes

For the majority of cases we will want to use synchronous send.

* Synchronous send - only completes when the receive has completed.
* Buffered send - always completes (unless an error occurs), irrespective of receiver.
* Standard send - either synchronous or buffered.
* Receive - completes when a message has arrived.

### Sending a message - synchronous

```c
int MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
```

```fortran
MPI_SSEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
```

### Sending data from rank 1 to rank 3

Sending an array:

```c
int x[10];

if (rank == 1)  // to make sure only rank 1 sends to 3
{
    MPI_Ssend(x, 10, MPI_INT, 3, 0, MPI_COMM_WORLD);
}
```

Sending a scalar:

```c
int x;

if (rank == 1)
{
    MPI_Ssend(&x, 1, MPI_INT, 3, 0, MPI_COMM_WORLD);
}
```

### Receiving a message

```c
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
```

```fortran
MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

<type> BUF(*)
INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS(MPI_STATUS_SIZE), IERROR
```

### Receive data from rank 1 on rank 3

Receiving an array:

```c
int y[10];
MPI_Status status;

if (rank == 3)
{
    MPI_Recv(y, 10, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
}
```

Receiving a scalar:

```c
int y;
MPI_STATUS status;

if (rank == 3)
{
    MPI_Recv(&y, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
}
```

### Synchronous blocking message passing

Processes synchronise. Sending process specifies the synchronous mode. If blocking, both processes wait until the transfer has completed. For a communication to succeed, the sender must specify a valid destination and then receiver must specify a valid source rank. The communicator and tags must match and the receiver's buffer must be large enough.

### Wildcarding

The receiver can wildcard. To receive from any source use `MPI_ANY_SOURCE` and to receive from any tag use `MPI_ANY_TAG`. The actual source and tag will be returned in receiver's status parameter.

### Received message count

```c
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);
```

```fortran
MPI_GET_COUNT(STATUS, DATATYPE, COUNT, IERROR)
INTEGER STATUS(MPI_STATUS_SIZE), DATATYPE, COUNT, IERROR
```

## Modes, tags and communicators

### Send types

* `MPI_Ssend` - synchronous send.
* `MPI_Bsend` - buffered send. Guaranteed to be asynchronous.
* `MPI_Send` - standard send. May be implemented to be synchronous or asynchronous.

### MPI_Bsend

For `Bsend`, the user has to provide a single large block of memory. This is made available to MPI by using `MPI_Buffer_attach`. If there is another `Bsend` before a receive and there is not enough buffer space, then `Bsend` will fail. We have to allocate how much memory we required for the buffer and then also allocate the amount of memory we require for the sent metadata.

### MPI_Send

`MPI_Send` tries to solve the problems with `MPI_Ssend` and `MPI_Bsend`. Buffer space is provided by the system. `Send` will normally be asynchronous, like `Bsend`, but if the message is too big for the buffer, i.e. if it is full, then `Send` becomes synchronous, like `Ssend`. `MPI_Send` is unlikely to fail, but can cause deadlocks in your program.

### Checking for messages

MPI allows us to check if any messages have arrived. You can 'probe' for matching messages.

```c
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status, *status);
```

Status is set as if the receive took place, e.g. you can find out the size of the message and allocate space prior to receive. We have to be careful with wildcards. We can use `MPI_ANY_SOURCE` in the call to probe, but we must use specific source in receive to guarantee matching the same message.

```c
MPI_Recv(buff, count, datatype, status.MPI_SOURCE, ...)
```

### Tags

Every message can have a tag. This is a non-negative integer value. The maximum tag value can be queried using `MPI_TAG_UB` attribute. MPI guarantees to support tags of at least 32767. Many MPI programs don't use tags and set them to 0, but they have to be defined. Tags can be useful in some situations where we can choose to only receive messages from only a given tag.

### Communicators

All MPI communications take place within a communicator. A communicator is fundamentally a group of processors. There is a pre-defined communicator `MPI_COMM_WORLD` which is all of the processors. **A message can only be received from the same communicator it was sent from.** Unlike tags, wildcards can't be used with communicators. The function `MPI_Comm_dup` can be used to duplicate a communicator. Enables processes to communicate with each other safely within a piece of code. This is **essential** for people writing parallel libraries.

## Non-blocking communication

### Completion

This mode of operation determines when its operations are completed. This determines when then value will be returned???

### Non-blocking operations

Non-blocking operations returns straight away and allow the sub-program to continue to perform other work. At some later time the sub-program can test or wait for the completion of the non-blocking operation.

All non-blocking operations should have matching wait operations. Some systems cannot free resources until wait has been called. **A non-blocking operation immediately followed by a matching wait is the same as a blocking operation.**

Separate communication into three phases:

* Initiate non-blocking communication.
* Do some work (could involve other communication).
* Wait for non-blocking communication to complete.

### Non-blocking synchronous send

```c
int MPI_Issend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);

// put code in the middle here

int MPI_Wait(MPI_Request *request, MPI_Status *status);
```

### Non-blocking receive

```c
int MPI_Irecv(void *buf, int count MPI_Datatype datatype, int src, int tag, MPI_Comm comm, MPI_Request *request);

// put code in the middle

int MPI_Wait(MPI_Request *request, MPI_Status *status);
```

### Non-blocking example

```c
MPI_Request request;
MPI_Status status;

if (rank == 0)
{
    MPI_Issend(sendarray, 10, MPI_INT, 1, tag, MPI_COMM_WORLD, &request);
    // Do something else while Issend happens
    MPI_Wait(&request, &status);
}
else if (rank == 1)
{
    MPI_Irecv(recvarray, 10, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
    // Do something else while Irecv happens
    MPI_Wait(&request, &status);
}
```

### Multiple communications

* Test or wait for completion of one message.
* Test or wait for completion of all messages.
* Test or wait for completion of as many messages as possible.

### Combined send and receive

Specify all the send and receive arguments in one call. This can be used to avoid deadlocks and is useful in simple pairwise communication patterns but isn't as generally application as non-blocking.

```c
int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
```

## Collective communications

### Characteristics

All processes must communicate, thus collective communications can be slow. Synchronisation may or may not occur. Standard collective operations are blocking. There are no tags and the **receive buffers must be exactly the right size.**

### Broadcast

Broadcast will send a buffer to all of the processes. It is useful when you want to replicate data.

```c
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
```

### Scatter

We can use scatter to share data between all of the processes. The scatter is across all of the processes across the communicator.

```c
int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
```

`sendcount` and `recvcount` in almost all cases will be the same. `sendcount` is the number of "chunks" which we send to each processes.

### Gather

```c
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
```

Brings everything onto one process???

### Global reduction operations

Used to compute a result involving data distributed over a group of processes.

* global sum or product
* global maximum or minimum
* global user-defined operation

```c
int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
```

`MPI_Op` defines the type of reduction operation to be done. Below is an example of using `MPI_Reduce` to do a reduction sum.

```c
MPI_Reduce(&x, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
```

### `MPI_Allreduce`

Instead of using `MPI_Reduce` we should use `MPI_Allreduce` to return the reduction operation to all processes.

```c
int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
```

## Virtual Topologies

Virtual topologies provide convenient process naming with a naming scheme with fits the communication pattern. In principle it simplifies the writing of code and can allow MPI to optimise communications.

### Topology Types

* Cartesian - each process is "connected" to its neighbours in a virtual grid. Boundaries can be cycling or not. Processes are identified by cartesian coordinates.
* Graph topologies - general graphs??? Not covered here.

```c
int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, int reorder, MPI_Comm *comm_cart)
```

Note that this function takes in a communicator and creates a new one. If reorder is true, the ranks can be re-ordered by the hardware. In practice, it is best to use reorder as false.

MPI provides a helper function which suggests the dimensions to use for the number of processes in use.

```c
int MPI_Dims_create(int nnodes, int ndims, int *dims);
```

Be sure to set `*dims` to be zero before the call, otherwise MPI will take a non-zero value as being a restriction and will try to factor around that.

If we want to find the rank or the coordinates of a create rank in the new communicator, we can use the following functions. **We need to be sure that we provide the new Cartesian communicator.**

```c
int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);
```

```c
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords);
```

We can use `MPI_Cart_shift` to find the neighbours of a current rank.

```c
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
```

If there are no neighbouring ranks, MPI returns a `NULL` processor - `MPI_PROC_NULL`. This is a black hole as send and receive tasks will finish straight away. Send buffers vanish but the receive buffer isn't touched.

### Cartesian partitioning

Cut a grid up into "slices". A new communicator is produced for each slice. Each slice can then perform its own collective communications.

```c
int MPI_Cart_sub(MPI_Comm comm, int *remain_dims, MPI_Comm *newcomm);
```

## MPI datatypes

There are basic types and derived types - vectors, structs and others.

### Defining types

All derived types stored by MPI as a list of basic types and displacements (in bytes).

### Contiguous Data

The simplest derived datatype consists of a number of contiguous items of the same datatype.

```c
int MPI_Type_contigous(int count, MPI_Datatype oldtype, MPI_Datatype, *newtype);
```

By using contiguous, we may be making the program more readable, for example it will become more readable if we want to send equal blocks of contiguous memory to processes.

### Vector type

A vector type corresponds to a subsection of a 2D array. In C, we have to use statically allocated arrays because dynamically allocated arrays have no defined storage format.

## Acknowledgements 

I would like to acknowledge financial support from the EPSRC Centre for Doctoral Training in Next Generation Computational Modelleing grant EP/L015382/1.
