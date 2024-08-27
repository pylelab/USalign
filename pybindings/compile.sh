
c++ -O3 -Wall -shared -std=c++20 -fvisibility=hidden -I /home/russ/Scripts/git/USalign -fPIC $(python3 -m pybind11 --includes) pySOIalign.cpp -o pySOIalign$(python3-config --extension-suffix)

