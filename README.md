to compile to .pd_darwin use the following command in the folder 

clang -shared -undefined dynamic_lookup -o model~.pd_darwin -I/Users/lezak/github/pure-data/src model~.c  
