### Compiling the response C++ code as a dynamic library
Step 1:
```
g++ -c -fPIC response.cc -Wall -g -O2 -lfftw3 -o response.o
```
Step 2:
```
g++ -shared -Wl,-soname,libpsf.so -o libpsf.so  response.o -L/usr/lib/x86_64-linux-gnu -lfftw3
```
You will now have a `.dll` that you can import in python named `libpsf.so`

Many of the files use C++17 features and fftw3
you can compile those files like this
```
g++ -std=c++17 filename.cc -Wall -O2 -lfftw3 -o exec_name
```
