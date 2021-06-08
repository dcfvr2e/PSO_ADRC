#include "pch.h"
#include "adrc.h"

RandomNumber r;       //随机数
int main()
{
	clock_t startTime, endTime; //定义程序开始运行时间和结束时间
	startTime = clock();  //计时开始
	pso PSO;  //定义PSO相关参数和函数

	//给定初始化粒子群速度和位置,计算粒子群适应度，初始化个体极值和全局极值
	ADRC_Init();					//ADRC init()
	PSO.Initialize_fit_extremum();

	//进入主要循环，按照公式依次迭代，直到满足精度要求
	PSO.Optimization_iteration();

	endTime = clock();//计时结束
	cout << "run time:" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}