#  mpi-2dmesh-harness-instructional

This project contains the code harness for doing an MPI-based scatter-process-gather
motif.

What the code provides:  
* Loads input files, parses command line arguments, writes output  
* Performs domain decomposition using one of three different strategies: row-slab, column-slab, or tiled,   sets up tiles on every rank that contain metadata like tile size, extents, etc.  , assigns tiles to ranks  
* Executes a scattering step by iterating over all tiles and invoking a method to send
data from rank 0 to other ranks, and from non-zero ranks to receive data from rank 0.
* Executes a processing step by iterating over all tiles and sets the stage for user code to be executed on a tile  
* Executes a tile gather step, where data from all tiles not owned by rank 0 is sent back to rank 0  
* Writes output to a disk file  

# Build instructions, general

This distribution uses cmake and MPI.

After downloading, cd into the mpi-2dmesh-harness directory, then:  

% mkdir build  
% cd build  
% cmake 
% make

# Platforms

In principle, this code will build and run on both Cori@NERSC and the VM. It's
fine to use the VM for initial development and testing, but please run your
performance tests on Cori.

# Setting up your NERSC environment

% module swap PrgEnv-intel PrgEnv-gnu  
% module load openmpi  
% module load cmake

# Adding your code: 3 locations

Your coding assignment consists of adding code to implement MPI-based communication
that needs to take place during the scatter and gather phases of processing, and to
also add code that computes the Sobel filter operation on data in a distributed memory
fashion.

### 1. Adding your code to sendStridedBuffer()

This method, which is called from both scatterAllTiles() and gatherAllTiles(),
is responsible for sending data from one rank to another. 

void  
sendStridedBuffer(float \*srcBuf,   
    int srcWidth, int srcHeight,   
    int srcOffsetColumn, int srcOffsetRow,  
    int sendWidth, int sendHeight,  
    int fromRank, int toRank )   

// ADD YOUR CODE HERE  

Your code will performs sending of data using MPI\_Send(), going _fromRank_ and 
to _toRank_. The data to be sent is in _srcBuf_, which has width _srcWidth_, _srcHeight_.
Your code needs to send a subregion of _srcBuf_, where the subregion is of size
_sendWidth_ by _sendHeight_ values, and the subregion is offset from the origin of
_srcBuf_ by the values specificed by _srcOffsetColumn_, _srcOffsetRow_.



### 2. Adding your code to recvStridedBuffer()

void  
recvStridedBuffer(float \*dstBuf,  
    int dstWidth, int dstHeight,  
    int dstOffsetColumn, int dstOffsetRow,  
    int expectedWidth, int expectedHeight,  
    int fromRank, int toRank )   

// ADD YOUR CODE HERE

Your code will perform receiving of data using MPI\_Recv(), where inbound data
is coming _fromRank_ and is destined for _toRank_. The that that arrives will be of size 
_expectedWidth_ by _expectedHeight_ values.  This incoming data is to be placed into 
the subregion of _dstBuf_ that has an origin at _dstOffsetColumn, dstOffsetRow_, and 
that is _expectedWidth, expectedHeight_ in size.


### 3. Adding your sobel filtering code

The intention here is for you to transplant part of your code from HW5 into this
MPI-based code. There are two locations where you will need to add code.

First, inside the _sobelAllTiles_ routine is a doubly-nested loop that iterates over tiles.
Inside the inner loop is a conditional that checks if _t->tileRank == myrank_, and if
so, you need to add a call to your _do_sobel_filtering()_ method from your sobel_cpu.cpp
code. You will invoke that method to process a tile's worth of data, passing in the
tile's _inputBuffer_ as the input data, and the tile's _outputBuffer_ to receive the
sobel filtered results.

Note: even though your _do_sobel_filtering()_ code was parallelized using OpenMP, we are
*not* using OpenMP parallelism in this project.

Second, you need to transplant your two methods, _do_sobel_filtering()_ and _sobel_filter_pixel()_ from your sobel\_cpu.cpp code into the mpi\_2dmesh.cpp code, where they will be invoked
to do the Sobel computation.

# Information about data files

Zebra file dimensions 
* Original: 3556 2573
* 4x Augmented: 7112 5146

* zebra-gray-int8.dat - raw 8-bit grayscale pixel values from the Zebra_July_2008-1.jpg image
* zebra-gray-int8-4x.dat - raw 8-bit grayscale pixel values from the Zebra_July_2008-1.jpg image but 
augmented 2x in each direction

Source file:  Zebra_July_2008-1.jpg, obtained from Wikimedia commons, https://commons.wikimedia.org/wiki/File:Zebra_July_2008-1.jpg

# python display script

imshow.py - a python script to display the raw 8-bit pixel values in grayscale. 

Usage:  
    python imshow.py filename-of-raw-8bit-bytes int-cols-width int-rows-height

To use this script from Cori, first please enable X11-tunneling through ssh before you
log in. E.g., "ssh -Y username@cori.nersc.gov"

On Cori, before using this script, please load the Python module:  
% module load python

Some people on Mac systems have reported problems with the remote image display
not working from Cori. This is likely the result of either (1) forgetting to enable
X tunneling through ssh by forgetting to add the -Y option to ssh, or (2) a
Mac-side Quartz configuration issue. 

If this issue affects you, you might want to consider trying NX to remote desktop
to Cori. See https://docs.nersc.gov/connect/nx/ for more information.
