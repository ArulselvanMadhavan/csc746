## How to run
- The instructions below show how to run it on Cori. (To run it on your local, edit runsum.in and remove srun)
- O3
```bash
rm -rf build && mkdir build && cd build/ && cmake ../ && make && bash ../runsum.in > opt.log 2>&1
```bash
- O0 - The below command was used for turning off compiler optimizations
```bash
rm -rf build && mkdir build && cd build/ && cmake ../ -DCMAKE_CXX_FLAGS_RELEASE="-O0" && make && bash ../runsum.in > no_opt.log 2>&1
```
