# Frame based CWT

[English version is here](https://github.com/james34602/iir_cwt_research/readme_en.md)

#### 文件列表

1. multResByIIRInFreq_archived
   算法初形

2. multResByIIRInFreq_dirichletKernelOutterWnd
   目前所有项目中最efficient的frame based CWT

3. multResByIIRInFreq_generalizeOutterWnd
   容许自定义全局分析窗的frame based CWT版本

4. multResByIIRInFreq_FrameTheoretic
   实现帧理论变换，frame based CWT的几乎无损逆变换最小二乘解(canonical dual)

5. slidingCQT
   基于Sliding DFT的CWT，特点在于来一个sample，输出一个频谱切片，类似传统CWT

6. multResByIIRInTime
   利用窄hann窗计算频谱，在帧上面circular shift，对频域每个bins在时间轴上进行IIR滤波实现多尺度分析，一种aliased signal processing方法，算法不可能重建输入

7. minFunc_2012
   无约束梯度优化工具

8. conv_constant-q-toolbox
   前人的低效、低SNR的可逆CQT实现

9. discreteLevelWithWindowing
   类似multResByIIRInTime，但是没有在时间轴上对每个频域bins进行IIR滤波，算法不可能重建输入

10. iirLayer
    时域时变SVF滤波的导数算式和测试代码

11. lsicqs
    前人的低效又复杂的CQT实现，针对个别应用

12. ADNode.m
    由于Matlab symbolic toolbox太慢，我写了个自动(符号)微分计算器，算式复杂起来，ADNode比Matlab symbolic toolbox快过千、万倍速度，可以求解离散傅里叶、某些overlap-add的special case、某些特殊函数、SVF滤波器(Backpropagation through time)、平滑round/floor/ceil、抗缠绕相位差函数等各类的reverse mode AD导数

13. testMorletCWT.m
    代码实现Textbook CWT，里面包含了两种信号重建(逆变换)方法，用来展示只要信号完全冗余，任何Filterbank皆可完美重建，包括完全随机的Filterbank系数

##### 预计算数据

链接：[https://pan.baidu.com/s/14qUPAxsxLeBW1LjoNAea4A?pwd=uppl]() 
提取码：uppl

##### 备注

1. 我经常叫CQT(Constant Q Transform)做CWT，这个讲法应该没问题，所有连续尺度、多尺度变换都是CWT。

CQT与Morlet CWT类同，而且CQT的Q值参数相当于Morlet可调Time-bandwidth product参数

2. 花了一个月去写时域时变SVF滤波的导数算式和测试。。。
