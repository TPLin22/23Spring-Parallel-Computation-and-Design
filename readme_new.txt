Author: 潘浩然 2100011091 
Update:2023.06.24


零、思路与优化

（一）
采用的并行技术：OpenMP（主要）、MPI（仅用于调用scalapack）
程序最重要的优化在下列的第（5）ii.和iii.步，
（1）到（5）完成了除scalapack外的其他任务，（6）完成了附加题2调用scalapack，

（二）
主要优化部分：在下面的步骤（5）.ii中，如果对于每一个格点(i,j,k)都要枚举所有的points来计算值的话，
时间复杂度是O(nx*ny*nz*m*m)的，这里m是points的总数。
但由于截断误差的存在，每个格点可能处在很少的points所在半径内，m*m的枚举是不必要的。
通过枚举一次m个points，就能找到真正能计算出非零值的points，符合要求的points其实很少，
时间复杂度大约能降到O(nx*ny*nz*m)。

6月24日更新：不再采用上述方式计算积分。在外层先枚举points而不是格点，
枚举两层points，仅在一对points距离小于cutoff的情况下才进行后续计算，
基础的时间复杂度仅为O(M^2)，然后枚举空间格点时，假设是两个以point为中心的立方体，
仅在立方体相交部分的i，j，k会被枚举，大大减小了遍历空间中的格点数。
经验证，优化后的程序能比原先快近百倍。


其他优化：尽可能减小了空间复杂度，开nx*ny*nz大小数组的仅有V值本身这一个数据，
插值函数在读入后就进行了预处理，储存了一阶差商、二阶差商等值，
后续求f(r)时，只需让r查找到所在的插值区间即可计算出答案，
查找区间用了二分法，这样每次插值用时仅为O(log mesh)，mesh是给定插值点的个数

一个没有实现的优化：由于V数据很大，用流读入会很慢，用fscanf可以改善速率，不过我很晚才想起来，没有实现

（三）
整个程序的流程：
（1）初始化MPI进程，下面从（1）到（5）步全部在rank 0 进程中进行：
先读入input.txt中的lx,ly,lz、diago_lib、points_path、venergy_path、distribution_path参数，
由于只做了附加题1、2，其他参数都可忽略。
（2）读入points文件，每个点的坐标都需要储存；读入V文件，储存格点的数量，动态分配内存，储存空间中每个格点的函数值；
读入分布函数文件，根据mesh个位置处的分布函数值，做三次样条插值的预处理，
将预处理数据储存在数组里，后续每次插值都只需以O(log(mesh))的复杂度找到插值位置就能简单计算出值。
（3）动态分配内存构造待求的对角矩阵H，设置好openmp线程数。

6月24日以前采用的（4）和（5）（已废弃）：
（4）用i，j，k依次枚举空间中的x、y、z方向上的格点，i的枚举用#pragma omp parallep for语句分配给各线程并行，
如此一来枚举的总量是nx*ny*nz，但每个线程都只需要处理nx/nthreads*ny*nz的任务量
（5）
i.
对于每个给定的i，j，k，这三个变量决定了此时处理的是空间中的第(i,j,k)个格点，其函数值V也是已知的，
首先根据i，j，k的值，nx，ny，nz的值，lx，ly，lz的值计算出此点真正的空间坐标(x,y,z)，
ii.
然后，先枚举所有的points，计算每个points的坐标与当前(x,y,z)的距离是否在截断半径cutoff之内。
这轮循环会统计有多少个points符合这个要求（记为numpos），
并且用一个数组记录它们分别是所有points中的第几个。
iii.
随后，用变量i2和j2枚举这numpos个points中的每一对，i2和j2表示是第i2和第j2个符合要求的point，
可以保证此时枚举的两个point都在截断半径之内，
计算V*f1(r1)*f2(r2)，f1(r1)和f2(r2)都是插值计算的，
将计算结果用omp的原子操作#pragma omp atomic累加到全局矩阵H[pos1][pos2]上。
pos1表示的是i2是所有points中的第几个，pos2表示的是j2是所有points中的第几个。


6月24日更新的（4）和（5）：
（4）用i，j枚举所有的points，对于一对points（i，j），判断两个point之间的距离是否小于等于cutoff，
仅在距离小于等于cutoff的情况下才执行（5）。
（5）假设以两个point为中心有两个长为2cutoff的正方体，这两个正方体相交的部分才有可能积分值不为0，
仅在这部分空间内枚举格点。用i1，j1，k1变量枚举空间中的格点，计算积分，累加到矩阵H[i][j]上。

（6）将矩阵H对角化。
i.
若用lapack对角化，仍然在rank0中调用lapack库，进行对角化。
ii.
若用scalapack对角化，则首先从进程0将矩阵H的维度M和矩阵H本身广播（MPI_Bcast）到其他进程中，
再在各个进程中进行如下操作，
用cblas_pinfo等函数进行scalapack上下文的预处理，用numroc获取当前进程块的行列数，
为了对应局部矩阵与全局矩阵之间的位置关系，用MPI_Allgather来收集各个进程的行、列数信息并储存在数组中，
如此就可以用求前缀和的方式来找到局部矩阵在全局的位置索引。
确认位置索引后，就可以把当前进程要处理的局部矩阵A从H中拷贝下来，
先调用一次pdsyev获取所需工作空间，再调用pdsyev将A对角化，
并将储存特征向量的局部矩阵Z根据先求所求的位置索引，对应到全局矩阵的位置，通过MPI_Reduce规约到rank 0进程，
最终将特征值和特征向量的结果由rank0储存并输出。


一、编译安装运行
（一）文件目录：
src——源代码 include——头文件 bin——可执行文件 
build——中间文件 input——输入文件 output——输出文件
（二）需要配置的环境：
对于一个新节点，进行以下步骤：
（1）用apt安装openmpi(先apt update,然后apt install openmpi-bin libopenmpi-dev)
（2）手动编译安装lapack
（即下载源码包，解压缩，编译安装，
下载的版本和命令是wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.0.tar.gz，
然后将压缩包直接tar -zxvf，进入文件夹lapack-3.11.0，将make.inc.example改为make.inc，直接运行四个命令：make blaslib,make cblaslib,make lapacklib,make lapackelib）
之后，在文件夹lapack-3.11.0下建立一个build文件夹，将lapack-3.11.0下的libcblas.a等四个后缀为.a的文件移入build文件夹中
（3）apt安装scalapack(apt install libscalapack-openmpi-dev)
（4）apt install time 
（三）编译
在cmakelist文件中，需要在第21行，根据lapack实际安装路径，修改此处的路径，
在第21行，将LAPACKE_PATH设置为编译安装的根目录；
在23行设置头文件${LAPACKE_PATH}/LAPACKE/include ${LAPACKE_PATH}/CBLAS/include ${MPI_INCLUDE_PATH}
24和25行链接库文件
其余内容不需要变动，openmp和mpi都可通过find_package自动找到。
进入build目录，先cmake.. 再make，即生成了在bin目录下的可执行文件。
（四）运行
举例：进入input目录下，运行命令：
mpirun --allow-run-as-root -np 4 --host localhost:4 ../bin/main ./INPUT.txt
最后一个参数是主输入文件，其他输入文件的路径在INPUT.txt中给出，可以在文件中修改。
输出的特征值和特征向量在output文件夹里。矩阵H也输出在文件夹里的out.log中
可能出现报错read -1 expected 20000 errno=38，
可以通过运行命令export OMPI_MCA_btl_vader_single_copy_mechanism=none解决。

读入point文件时，请以(100,100,100)这样有括号有逗号的形式读入。
在输入文件INPUT.txt输入diago_lib lapack或diago_lib scalapack决定使用哪一种库。（只要输入不是lapack，都会采用scalapack）

串行和并行执行程序的指定方式：直接在main.cpp的约77行设置nthreads为线程数量。
默认已经设定为16，如要串行执行程序，设置nthreads为1即可。

由于程序中采用文件输出，
打开输出文件时使用的是相对路径
(如main.cpp中的ofs.open("../output/scaeigenvalues.log",ios::out);)，
(涉及到文件输出的有main.cpp的line133,261,287,mathwork.cpp的line117,125，可以根据需要修改)
因此不要在根目录直接执行，在bin,input目录下执行都可以，
这样可以保证程序打开 父文件夹-output文件夹下的输出文件。


二、运行效率

在16个线程时，
V512，POINT50 运行约50s（计算部分<7s）
V1024，POINT10 运行约100s（计算部分<20s）
V1024，POINT50 运行约120s （计算部分40s）

三、数据结构设计
类Input_0读入主输入文件input.txt，
类Input_v用于读取文件v.txt并储存其值。储存空间函数值V的数组，大小最大为1024*1024*1024，会根据实际需要动态分配内存。
类Input_p用于读入point.txt，并用二维数组pos[m][4]储存了m个点的坐标，pos[i][j]，j=1,2,3时分别代表第i个点的x、y、z坐标。
类input_d用于读入分布函数在各个点的取值，
类mathwork用input_d的数据来构造，用多个一维数组储存了插值所需的一阶、二阶差商等数学数据。
