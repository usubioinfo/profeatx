# Compile from source

There is no need to compile it from source. However, if you want to do it, you will need GCC 4.8.1+ and run the following command:

`g++ -std=c++11 -lz -fopenmp -O3  main.cpp -static-libstdc++ -o profeatx`

# Run ProFeatX

You can run ProFeatX using the following command:

`./profeatx -i <input> -o <output> -e <encoding>`

In order to check the help:

`./profeatx --help`
