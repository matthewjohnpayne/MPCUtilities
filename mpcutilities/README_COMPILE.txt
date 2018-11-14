gcc -c -Wall -Werror -fpic kepcart.c
gcc -shared -o libkepcart.so kepcart.o

