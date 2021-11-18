# cartpole_casadi_cplusplus  
小车倒立摆模型轨迹优化，模型及问题定义见https://epubs.siam.org/doi/pdf/10.1137/16M1062569  
求解工具：casadi c++接口(ipopt)  
optimal contorl problem -> nonlinear optimization problem: direct transcription: single/multiple shooting、 collocatin  
文档链接：https://zhuanlan.zhihu.com/p/391903468  
  
运行代码  
安装casadi,安装教程：https://github.com/casadi/casadi/wiki/InstallationLinux  
下载代码后进入根目录  
mkdir build && cd build && cmake ..  
make  
./mytest  
