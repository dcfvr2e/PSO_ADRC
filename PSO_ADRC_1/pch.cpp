#include "pch.h"
#include "adrc.h"


double function(vector<double> x)
{
	vector<double> y = ADRC_sim(x);	//计算该参数下模型的调节时间
	double fx = cost_function(y);			//计算该调节时间下的代价函数
	return fx;
}

void pso::Initialize_fit_extremum()
{
	//===================给定初始化粒子群速度和位置==========
	extern RandomNumber r;       //定义全局随机数
	xp_best.resize(N_particle, vector<double>(D_particle));
	xg_best.resize(D_particle);
	fp_best.resize(N_particle);
	x_i.resize(N_particle, vector<double>(D_particle));		//N_particle=20,	D_particle=2
	v_i.resize(N_particle, vector<double>(D_particle));
	for (int i = 0; i < N_particle; i++)
	{
		for (int j = 0; j < D_particle; j++)
		{
			x_i[i][j] = r.decimal(x_low[j], x_high[j]);    //随机初始化位置
			v_i[i][j] = r.decimal(v_low[j], v_high[j]);    //随机初始化速度
		}
	}
	//================计算粒子群适应度，初始化个体极值和全局极值=========
	for (int j = 0; j < N_particle; j++)
	{
		fp_best[j] = function(x_i[j]);  //先计算各个粒子的适应度
		xp_best[j] = x_i[j];
	}
	xg_best = x_i[0];
	fg_best = fp_best[0];
	fg_best_last = fg_best;
	for (int j = 1; j < N_particle; j++)
	{
		if (fp_best[j] < fg_best)
		{
			xg_best = x_i[j];   //更新最优值对应的优化变量
			fg_best = fp_best[j];
		}                   //找到全局最优变量
	}
}
//**************************************
//粒子群寻优迭代
//*************************************
void pso::Optimization_iteration()
{
	extern RandomNumber r;       //定义全局随机数
	ofstream out_beta1("PSO para beta1.csv");
	ofstream out_beta2("PSO para beta2.csv");
	ofstream out_b("PSO para b.csv");
	ofstream out_b1("PSO para b1.csv");
	ofstream out_b2("PSO para b2.csv");
	ofstream out_b3("PSO para b3.csv");
	ofstream out_res("PSO iteration result.csv");
	ofstream out_f_value("PSO f value.csv");
	double f;
	for (int k = 0; k < M_particle; k++)		//迭代次数
	{
		int flag_sum = 0;
		w_particle = w_max - k * (w_max - w_min) / M_particle;  //更新权重
		for (int j = 0; j < N_particle; j++)		//种群数量
		{
			r1 = r.decimal(0, 1.0), r2 = r.decimal(0, 1.0);    //r1、r2产生随机数
			for (int i = 0; i < D_particle; i++)
			{
				v_i[j][i] = w_particle * v_i[j][i] + c1_particle * r1*(xp_best[j][i] - x_i[j][i]) + c2_particle * r2*(xg_best[i] - x_i[j][i]);  //更新位置
			//速度越界处理，取边界值	
				if (v_i[j][i] < v_low[i])
				{
					v_i[j][i] = v_low[i];//越界更新为边界值
				}
				if (v_i[j][i] > v_high[i])
				{
					v_i[j][i] = v_high[i];//越界更新为边界值
				}
				x_i[j][i] = x_i[j][i] + v_i[j][i];   //更新位置
			//位置越界处理，取边界值
				if (x_i[j][i] < x_low[i])
				{
					x_i[j][i] = x_low[i];//越界更新为最优值
				}
				if (x_i[j][i] > x_high[i])
				{
					x_i[j][i] = x_high[i];//越界更新为边界值
				}
			}
			f = function(x_i[j]);
			if (f < fp_best[j])			//个体选优
			{
				fp_best[j] = f;
				xp_best[j] = x_i[j];
			}
			if (fp_best[j] < fg_best)           //全局最优
			{
				xg_best = xp_best[j];
				fg_best = fp_best[j];
			}
			if (abs(fp_best[j] - fg_best) < 0.03)
				flag_sum++;
			out_beta1 << setprecision(5) << x_i[j][0] << ",";//行表示种群数量20，且表示迭代数*整定参数维度20*2
			out_beta2 << setprecision(5) << x_i[j][1] << ",";
			out_b << setprecision(5) << x_i[j][2] << ",";
			out_b1 << setprecision(5) << x_i[j][3] << ",";
			out_b2 << setprecision(5) << x_i[j][4] << ",";
			out_b3 << setprecision(5) << x_i[j][5] << ",";
			out_f_value << setprecision(5) << f << ",";
		}
		out_beta1 << endl;
		out_beta2 << endl;
		out_b << endl;
		out_b1 << endl;
		out_b2 << endl;
		out_b3 << endl;
		out_f_value << endl;
		out_res << k << fixed << setw(12) << setprecision(5) << fg_best << endl;	//粒子群算法优化结果
		fg_best_last = fg_best;
		if (abs(fg_best_last - fg_best) / fg_best_last < 1e-4 && flag_sum > 0.75 * N_particle)			//设置迭代终止条件
			break;
	}
	function_1(xg_best);			//每次迭代将最优结果保存
	out_res << "最优变量:" << endl;
	for (int i = 0; i < D_particle; i++)
	{
		out_res << "x" << i << "=" << fixed << setw(12) << setprecision(5) << xg_best[i] << endl;	//输出最优变量
	}
	out_res << "最优值=" << fixed << setw(12) << setprecision(5) << fg_best << endl;
	out_res.close();
	out_beta1.close();
	out_beta2.close();
	out_b.close();
	out_b1.close();
	out_b2.close();
	out_b3.close();
	out_f_value.close();
}

double cost_function(vector<double> y_v)
{
	double Tt = 0.005;		//仿真采样周期
	double y_infty = y_v[y_v.size() - 1];
	vector<double> err;
	for (unsigned int i = 0; i < y_v.size(); i++)
		err.push_back(abs(y_v[i] - y_infty));
	reverse(err.begin(), err.end());
	double delt_err = 0.02;
	unsigned int iter;
	for (iter = 0; iter < err.size(); iter++)//逆序后首次进入稳态0.02倍的点
		if (err[iter] > y_infty * delt_err)
			break;
	double ts = (y_v.size() - iter) * Tt;		//调节时间
	double sigma;
	auto y_max_pos = max_element(y_v.begin(), y_v.end());		//峰值
	if (*y_max_pos < y_infty)
		sigma = 0;
	else
		sigma = (*y_max_pos - y_infty) / y_infty;	//超调
	double fx_v = 0.2 * (ts - 0.6) * (ts - 0.6) + 0.8 * (sigma - 0.02) * (sigma - 0.02);		//评价指标
	cout << "调节时间：" << fixed << setw(10) << setprecision(4) << ts
		<< "   超调量：" << fixed << setw(10) << setprecision(4) << sigma
		<< "   评价指标：" << fixed << setw(10) << setprecision(4) << fx_v << endl;
	return fx_v;
}
