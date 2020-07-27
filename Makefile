quick:
	g++ -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` nbody.cpp -o nbody`python3-config --extension-suffix` -D PYTHON

python:
	g++ -O3 -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` nbody.cpp -o nbody`python3-config --extension-suffix` -D PYTHON
