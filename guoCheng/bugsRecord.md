##bug1:
在定义print函数的时候,文件输出流没有定义成引用导致以下报错:
```shell
use of deleted function ‘std::basic_fstream<_CharT, _Traits>::basic_fstream(const std::basic_fstream<_CharT, _Traits>&) [with _CharT = char; _Traits = std::char_traits<char>]’
```
###fix:
```c++
void print(Field aField, fstream ) 改为
void print(Field aField, fstream &)
```

##bug2:
 fatal error: main.H: No such file or directory
fix:
```c++
    改
	#include<main.H>
	为
    #include"main.H"
```
##bug3:
make自动推导依赖关系是没法调试的!
下面的代码是可以编译的,但是因为没加-g, 无法调试!
```Makefile
OBJS:=$(wildcard *.o)
all:$(OBJS)
assignment4.o: assignment4.C main.H
solve.o: solve.C main.H
```

###fix:
```Makefile
OBJS:=$(wildcard *.o)
#OBJS=assignment4.o solve.o
CFLAG=-g -O2 -Wall -std=c++11 
all:$(OBJS) 
	g++ $(CFLAG) -lm $(OBJS)
assignment4.o: assignment4.C main.H
	g++ -c $(CFLAG) assignment4.C 
solve.o: solve.C main.H
	g++ -c $(CFLAG) solve.C
```