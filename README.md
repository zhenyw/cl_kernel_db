##cl_kernel_db

Inspired by shader-db to run standalone OCL kernel compiling to evaluate
OCL compiler.

###Run

Type *make* to build required program.

Run by

> ./run.py > new_result

More option is available from

> ./run.py -h

Currently support to specify OCL platforms for beignet or Intel OpenCL SDK.

To get report comparison by

> ./report.py old_result new_result

###Kernels

Currently kernels from luxmark are used. More will be added.

Only binary size and build time is recorded and compared for now.

