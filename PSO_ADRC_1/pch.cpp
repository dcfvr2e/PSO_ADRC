#include "pch.h"
#include "adrc.h"


double function(vector<double> x)
{
	vector<double> y = ADRC_sim(x);	//����ò�����ģ�͵ĵ���ʱ��
	double fx = cost_function(y);			//����õ���ʱ���µĴ��ۺ���
	return fx;
}

void pso::Initialize_fit_extremum()
{
	//===================������ʼ������Ⱥ�ٶȺ�λ��==========
	extern RandomNumber r;       //����ȫ�������
	xp_best.resize(N_particle, vector<double>(D_particle));
	xg_best.resize(D_particle);
	fp_best.resize(N_particle);
	x_i.resize(N_particle, vector<double>(D_particle));		//N_particle=20,	D_particle=2
	v_i.resize(N_particle, vector<double>(D_particle));
	for (int i = 0; i < N_particle; i++)
	{
		for (int j = 0; j < D_particle; j++)
		{
			x_i[i][j] = r.decimal(x_low[j], x_high[j]);    //�����ʼ��λ��
			v_i[i][j] = r.decimal(v_low[j], v_high[j]);    //�����ʼ���ٶ�
		}
	}
	//================��������Ⱥ��Ӧ�ȣ���ʼ�����弫ֵ��ȫ�ּ�ֵ=========
	for (int j = 0; j < N_particle; j++)
	{
		fp_best[j] = function(x_i[j]);  //�ȼ���������ӵ���Ӧ��
		xp_best[j] = x_i[j];
	}
	xg_best = x_i[0];
	fg_best = fp_best[0];
	fg_best_last = fg_best;
	for (int j = 1; j < N_particle; j++)
	{
		if (fp_best[j] < fg_best)
		{
			xg_best = x_i[j];   //��������ֵ��Ӧ���Ż�����
			fg_best = fp_best[j];
		}                   //�ҵ�ȫ�����ű���
	}
}
//**************************************
//����ȺѰ�ŵ���
//*************************************
void pso::Optimization_iteration()
{
	extern RandomNumber r;       //����ȫ�������
	ofstream out_beta1("PSO para beta1.csv");
	ofstream out_beta2("PSO para beta2.csv");
	ofstream out_b("PSO para b.csv");
	ofstream out_b1("PSO para b1.csv");
	ofstream out_b2("PSO para b2.csv");
	ofstream out_b3("PSO para b3.csv");
	ofstream out_res("PSO iteration result.csv");
	ofstream out_f_value("PSO f value.csv");
	double f;
	for (int k = 0; k < M_particle; k++)		//��������
	{
		int flag_sum = 0;
		w_particle = w_max - k * (w_max - w_min) / M_particle;  //����Ȩ��
		for (int j = 0; j < N_particle; j++)		//��Ⱥ����
		{
			r1 = r.decimal(0, 1.0), r2 = r.decimal(0, 1.0);    //r1��r2���������
			for (int i = 0; i < D_particle; i++)
			{
				v_i[j][i] = w_particle * v_i[j][i] + c1_particle * r1*(xp_best[j][i] - x_i[j][i]) + c2_particle * r2*(xg_best[i] - x_i[j][i]);  //����λ��
			//�ٶ�Խ�紦��ȡ�߽�ֵ	
				if (v_i[j][i] < v_low[i])
				{
					v_i[j][i] = v_low[i];//Խ�����Ϊ�߽�ֵ
				}
				if (v_i[j][i] > v_high[i])
				{
					v_i[j][i] = v_high[i];//Խ�����Ϊ�߽�ֵ
				}
				x_i[j][i] = x_i[j][i] + v_i[j][i];   //����λ��
			//λ��Խ�紦��ȡ�߽�ֵ
				if (x_i[j][i] < x_low[i])
				{
					x_i[j][i] = x_low[i];//Խ�����Ϊ����ֵ
				}
				if (x_i[j][i] > x_high[i])
				{
					x_i[j][i] = x_high[i];//Խ�����Ϊ�߽�ֵ
				}
			}
			f = function(x_i[j]);
			if (f < fp_best[j])			//����ѡ��
			{
				fp_best[j] = f;
				xp_best[j] = x_i[j];
			}
			if (fp_best[j] < fg_best)           //ȫ������
			{
				xg_best = xp_best[j];
				fg_best = fp_best[j];
			}
			if (abs(fp_best[j] - fg_best) < 0.03)
				flag_sum++;
			out_beta1 << setprecision(5) << x_i[j][0] << ",";//�б�ʾ��Ⱥ����20���ұ�ʾ������*��������ά��20*2
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
		out_res << k << fixed << setw(12) << setprecision(5) << fg_best << endl;	//����Ⱥ�㷨�Ż����
		fg_best_last = fg_best;
		if (abs(fg_best_last - fg_best) / fg_best_last < 1e-4 && flag_sum > 0.75 * N_particle)			//���õ�����ֹ����
			break;
	}
	function_1(xg_best);			//ÿ�ε��������Ž������
	out_res << "���ű���:" << endl;
	for (int i = 0; i < D_particle; i++)
	{
		out_res << "x" << i << "=" << fixed << setw(12) << setprecision(5) << xg_best[i] << endl;	//������ű���
	}
	out_res << "����ֵ=" << fixed << setw(12) << setprecision(5) << fg_best << endl;
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
	double Tt = 0.005;		//�����������
	double y_infty = y_v[y_v.size() - 1];
	vector<double> err;
	for (unsigned int i = 0; i < y_v.size(); i++)
		err.push_back(abs(y_v[i] - y_infty));
	reverse(err.begin(), err.end());
	double delt_err = 0.02;
	unsigned int iter;
	for (iter = 0; iter < err.size(); iter++)//������״ν�����̬0.02���ĵ�
		if (err[iter] > y_infty * delt_err)
			break;
	double ts = (y_v.size() - iter) * Tt;		//����ʱ��
	double sigma;
	auto y_max_pos = max_element(y_v.begin(), y_v.end());		//��ֵ
	if (*y_max_pos < y_infty)
		sigma = 0;
	else
		sigma = (*y_max_pos - y_infty) / y_infty;	//����
	double fx_v = 0.2 * (ts - 0.6) * (ts - 0.6) + 0.8 * (sigma - 0.02) * (sigma - 0.02);		//����ָ��
	cout << "����ʱ�䣺" << fixed << setw(10) << setprecision(4) << ts
		<< "   ��������" << fixed << setw(10) << setprecision(4) << sigma
		<< "   ����ָ�꣺" << fixed << setw(10) << setprecision(4) << fx_v << endl;
	return fx_v;
}
