# Compile from source

There is no need to compile it from source. However, if you want to do it, you will need GCC 4.8.1+ or Clang 3.7+, and run one of the following commands, depending on the compiler you have:

`g++ -std=c++11 -lz -fopenmp -O3  main.cpp -static-libstdc++ -o profeatx`

`clang -fopenmp=libomp -O3 -o profeatx main.cpp`

# Run ProFeatX

You can run ProFeatX using the following command:

`./profeatx -i <input> -o <output> -e <encoding>`

In order to see the help:

`./profeatx --help`
