#pragma once
// Particle_swarm.cpp : 粒子群算法实现过程。
//开发人员：陈帅    开发日期：2019年8月4日-29日   邮箱：chenshuai0614@hrbeu.edu.cn
#ifndef PCH_H
#define PCH_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>

using namespace std;
//产生随机小数或整数
class RandomNumber {
public:
	RandomNumber() {
		srand(time(0));    //构造函数，在对象创建时数据成员执行初始化操作
	}
	int integer(int begin, int end)		//随机整数
	{
		return rand() % (end - begin + 1) + begin;
	}
	double decimal(double a, double b)		//计算两者之间的随机小数
	{
		return double(rand() % 10000) / 10000 * (b - a) + a;
	}
};
//定义粒子群
class pso
{
private:
	double c1_particle = 2.0; //c1学习因子1，自我认知
	double c2_particle = 2.0; //c2学习因子2，社会认知
	double w_max = 0.9;       //最大惯性权重因子，影响着全局搜索
	double w_min = 0.6;       //最小惯性权重因子，局部深度搜索能力
	int M_particle = 40;     //最大迭代次数，200
	int D_particle = 6;       //搜索空间的维数也叫变量个数
	int N_particle = 30;      //初始化群体的个体，N很小容易陷入局部优化,N很大优化能力很好，优化的寻优精度和运行时间
public:
	vector <vector <double>>x_i;     //粒子群位置
	vector<vector <double>>v_i;      //粒子群速度
	vector<vector<double>>xp_best;   //个体最优位置
	vector<double> xg_best;          //全局最优位置
	vector<double>fp_best;           //个体最优值
	double fg_best;                  //全局最优值
	double fg_best_last;			 //上一时刻全局最优值
	double w_particle;               //权重更新
//parameters need to be tuned (beta1, beta2, b,		b1,		b2,		b3		)
	vector<double>x_low = {	1.0,	1.0,	1.0,	1.0,	200.0, 400.0	};   //优化变量下限值
	vector<double>x_high = { 50.0,	30.0,	10.0,	50.0,	500.0, 1000.0	};   //优化变量上限值
	vector<double>v_low = { -4.0,	-3.0,	-0.9,	 -5.0,	 -30.0, -100.0	};     //飞行速度下限值
	vector<double>v_high = { 4.0,	3.0,	0.9,	 5.0,	 30.0,	100.0	};     //飞行速度上限值
	double r1, r2;                    //r1、r2为增加随机搜索性
	void Initialize_fit_extremum();	//给定初始化粒子群速度和位置,计算粒子群适应度，初始化个体极值和全局极值
	void Optimization_iteration();//寻优迭代
};
double function(vector<double> x);		//目标函数
void function_1(vector<double> x);		//最终优化参数得到优化曲线
double cost_function(vector<double> y_v);	//模型函数
#endif //PCH_H