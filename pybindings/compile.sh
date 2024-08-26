
c++ -O3 -Wall -shared -std=c++20 -fvisibility=hidden -I /home/russ/Scripts/git/USalign -fPIC $(python3 -m pybind11 --includes) pybind_SOIalign_main.cpp -o pybind_SOIalign_main$(python3-config --extension-suffix)

