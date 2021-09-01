## How to run
- O3
```bash
rm -rf build && mkdir build && cd build/ && cmake ../ && make && bash ../runsum.in > opt.log 2>&1
```bash
- O0 - The below command was used for turning off compiler optimizations
```bash
cmake ../ -DCMAKE_CXX_FLAGS_RELEASE="-O0"
```
