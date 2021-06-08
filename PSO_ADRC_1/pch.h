#pragma once
// Particle_swarm.cpp : ����Ⱥ�㷨ʵ�ֹ��̡�
//������Ա����˧    �������ڣ�2019��8��4��-29��   ���䣺chenshuai0614@hrbeu.edu.cn
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
//�������С��������
class RandomNumber {
public:
	RandomNumber() {
		srand(time(0));    //���캯�����ڶ��󴴽�ʱ���ݳ�Աִ�г�ʼ������
	}
	int integer(int begin, int end)		//�������
	{
		return rand() % (end - begin + 1) + begin;
	}
	double decimal(double a, double b)		//��������֮������С��
	{
		return double(rand() % 10000) / 10000 * (b - a) + a;
	}
};
//��������Ⱥ
class pso
{
private:
	double c1_particle = 2.0; //c1ѧϰ����1��������֪
	double c2_particle = 2.0; //c2ѧϰ����2�������֪
	double w_max = 0.9;       //������Ȩ�����ӣ�Ӱ����ȫ������
	double w_min = 0.6;       //��С����Ȩ�����ӣ��ֲ������������
	int M_particle = 40;     //������������200
	int D_particle = 6;       //�����ռ��ά��Ҳ�б�������
	int N_particle = 30;      //��ʼ��Ⱥ��ĸ��壬N��С��������ֲ��Ż�,N�ܴ��Ż������ܺã��Ż���Ѱ�ž��Ⱥ�����ʱ��
public:
	vector <vector <double>>x_i;     //����Ⱥλ��
	vector<vector <double>>v_i;      //����Ⱥ�ٶ�
	vector<vector<double>>xp_best;   //��������λ��
	vector<double> xg_best;          //ȫ������λ��
	vector<double>fp_best;           //��������ֵ
	double fg_best;                  //ȫ������ֵ
	double fg_best_last;			 //��һʱ��ȫ������ֵ
	double w_particle;               //Ȩ�ظ���
//parameters need to be tuned (beta1, beta2, b,		b1,		b2,		b3		)
	vector<double>x_low = {	1.0,	1.0,	1.0,	1.0,	200.0, 400.0	};   //�Ż���������ֵ
	vector<double>x_high = { 50.0,	30.0,	10.0,	50.0,	500.0, 1000.0	};   //�Ż���������ֵ
	vector<double>v_low = { -4.0,	-3.0,	-0.9,	 -5.0,	 -30.0, -100.0	};     //�����ٶ�����ֵ
	vector<double>v_high = { 4.0,	3.0,	0.9,	 5.0,	 30.0,	100.0	};     //�����ٶ�����ֵ
	double r1, r2;                    //r1��r2Ϊ�������������
	void Initialize_fit_extremum();	//������ʼ������Ⱥ�ٶȺ�λ��,��������Ⱥ��Ӧ�ȣ���ʼ�����弫ֵ��ȫ�ּ�ֵ
	void Optimization_iteration();//Ѱ�ŵ���
};
double function(vector<double> x);		//Ŀ�꺯��
void function_1(vector<double> x);		//�����Ż������õ��Ż�����
double cost_function(vector<double> y_v);	//ģ�ͺ���
#endif //PCH_H