###bug1:
数组越界会报错segmenta fault! (core dumped)
写下标的时候一定要检查是否越界!

###bug2:
-nan: not a number
除以0; 开方负数; log负数等!

###bug3:
bug:
伪代码
```c++
typedef double Field[maxSpace+1][3];
typedef double Vector[3];
double WToP(Vector W)
{
    return (GAMMA-1) * ( W[2] - 0.5* W[1]*W[1]/W[0] );
}
void print( const Field  W)
{
    ofstream fout("result.dat");
    for (unsigned i = 0; i <= maxSpace; i++) 
        fout<<WToP(W[i])<<endl;
}
```

报错error: invalid conversion from ‘const double*’ to ‘double*’ [-fpermissive]
 note:   initializing argument 1 of ‘double pres(double*)’

原因: 向print函数传递的参数为 const引用, 无法再利用它传递引用

fix: 去掉void print( const Field  W)中的const

bug:

(.text+0x20): undefined reference to `main'
/tmp/ccPbhQjp.o: In function `FlowField::solve(Mesh&)':
friend.C:(.text+0x61c): undefined reference to `Mesh::getVolume(unsigned int, unsigned int)'
未定义main函数

bug:

 variable 'std:ofstream’ has initializer but incomplete type

或者是variable 'std:ifstream’ has initializer but incomplete type

其原因是因为没有包含fstream这个头文件。
