#include "adrc.h"
#include "pch.h"


fhanParas_TypeDef TD_fhanParas_RollRadio;
TDState_TypeDef TDState_RollRadio;
ESOParas_TypeDef ESOParas_Roll;
ESOState_TypeDef ESOState_Roll;
NLSEFState_TypeDef NLSEFState_Roll;

vector<vector<double>> data_log;
default_random_engine rand_e; //引擎
normal_distribution<double> rand_noise(0, 0.001); //均值, 方差

static void fhan_Init(void)				//fhan()函数：自抗扰技术中，离散系统最速控制综合函数的一种简化
{
	TD_fhanParas_RollRadio.h = 0.005;
	TD_fhanParas_RollRadio.r = 50;
}

static void TD_Init(void)
{
	TDState_RollRadio.h = 0.005;
	TDState_RollRadio.x1 = 0;
	TDState_RollRadio.x2 = 0;
}

static void ESO_Init(void)
{
	ESOParas_Roll.h = 0.005;
	ESOParas_Roll.b = 7;
	ESOParas_Roll.b1 = 50;
	ESOParas_Roll.b2 = 150;
	ESOParas_Roll.b3 = 800;
	ESOParas_Roll.a1 = 0.7;
	ESOParas_Roll.a2 = 0.1;
	ESOParas_Roll.d = 0.05;

	ESOState_Roll.z1 = 0;
	ESOState_Roll.z2 = 0;
	ESOState_Roll.z3 = 0;
}

static void NLSEF_Init(void)
{
	NLSEFState_Roll.b1 = 5;
	NLSEFState_Roll.b2 = 2.5;
	NLSEFState_Roll.b = 7;
	NLSEFState_Roll.a1 = 0.6;
	NLSEFState_Roll.a2 = 0.9;
	NLSEFState_Roll.d = 0.02;
	NLSEFState_Roll.u = 0;
}

void ADRC_Init(void)
{
	fhan_Init();
	TD_Init();
	ESO_Init();
	NLSEF_Init();
}

/* ------------------------- TD -------------------------------*/
static float fhan(float x1, float x2, fhanParas_TypeDef *para)
{
	float h, d, d0, y, a0, r, a;
	float res;
	r = para->r;
	h = para->h;
	d = r * h;
	d0 = h * d;

	y = x1 + h * x2;
	a0 = sqrt(d*d + 8 * r * fabs(y));

	if (fabs(y) > d0)
	{
		if (y >= 0.0f)
			a = x2 + (a0 - d) / 2;
		else
			a = x2 - (a0 - d) / 2;
	}
	else
		a = x2 + y / h;

	if (fabs(a) > d)
	{
		if (a > 0)
			res = -r;
		else
			res = r;
	}
	else
		res = -r * a / d;

	return res;
}

void TD_Atti(TDState_TypeDef *state, float v, fhanParas_TypeDef *para)
{
	float fh;
	float x1, x2;
	x1 = state->x1;
	x2 = state->x2;
	fh = fhan(x1 - v, x2, para);
	state->x1 = x1 + x2 * state->h;
	state->x2 = x2 + fh * state->h;
}

void ESO_Atti(const double y, const double u, ESOParas_TypeDef *para, ESOState_TypeDef *state)
{
	double e, fe, fe1;
	double z1, z2, z3;

	z1 = state->z1;
	z2 = state->z2;
	z3 = state->z3;

	e = state->z1 - y;
	if (abs(e) > para->d) {
		fe = pow(abs(e), para->a1) * sign(e);
		fe1 = pow(abs(e), para->a2) * sign(e);
	}
	else {
		fe = e / pow(para->d, (1 - para->a1));
		fe1 = e / pow(para->d, (1 - para->a2));
	}

	state->z1 = z1 + (z2 - para->b1 * e) * para->h;
	state->z2 = z2 + (z3 - para->b2 * fe + para->b * u) * para->h;
	state->z3 = z3 - para->b3 * fe1 * para->h;
}

void NLSEF_Atti(TDState_TypeDef *tdstate, ESOState_TypeDef *esostate, NLSEFState_TypeDef *nlsefstate)
{

	float e1, e2, u1, u2;

	e1 = tdstate->x1 - esostate->z1;
	e2 = tdstate->x2 - esostate->z2;
	if (abs(e1) > nlsefstate->d)
		u1 = pow(abs(e1), nlsefstate->a1) * sign(e1);
	else
		u1 = e1 / pow(nlsefstate->d, 1 - nlsefstate->a1);

	if (abs(e2) > nlsefstate->d)
		u2 = pow(abs(e2), nlsefstate->a2) * sign(e2);
	else
		u2 = e2 / pow(nlsefstate->d, 1 - nlsefstate->a2);

	nlsefstate->u = nlsefstate->b1 * u1 + nlsefstate->b2 * u2 - esostate->z3 / nlsefstate->b;
}
vector<double> ADRC_sim(vector<double> x_v)
{
	double Tt = 0.005;		//仿真采样时间，需要相应修改差分传递函数
	double total_t = 10;
	int N = floor(total_t / Tt);	//样本总数
	double phi_ref = 1;					//参考值
	vector<double> y(N + 1, 0);
	vector<double> u(N+1,0);
	vector<double> T = { 0,0.005 };
	//		       5					6.131e-05 z + 6.015e-05
	//G(s) = ---------------   --->  --------------------------
	//		   s^2 + 11.55s				z^2 - 1.944 z + 0.9439
	NLSEFState_Roll.b1 = x_v[0];
	NLSEFState_Roll.b2 = x_v[1];
	NLSEFState_Roll.b = x_v[2];
	ESOParas_Roll.b = x_v[2];
	for (int i = 2; i <= N; i++) {
		y[i] = 1.944 * y[i - 1] - 0.944 * y[i - 2] + 6.131e-05 * u[i - 1] + 6.015e-05 * u[i - 2] + 0.1 * rand_noise(rand_e);
		
		TD_Atti(&TDState_RollRadio, phi_ref, &TD_fhanParas_RollRadio);
		ESO_Atti(y[i-1], u[i-1], &ESOParas_Roll, &ESOState_Roll);
		NLSEF_Atti(&TDState_RollRadio, &ESOState_Roll, &NLSEFState_Roll);
		u[i] = NLSEFState_Roll.u;
		T.push_back(Tt*i);
	}
	data_log.push_back(T);
	data_log.push_back(y);
	return y;
}

void function_1(vector<double> x) {		//计算并保留最优参数下的数据
	ofstream out_y_best("PSO y_best.txt");
	ofstream out_noise("PSO noise.txt");
	double Tt = 0.005;		//仿真采样时间
	double total_t = 10;
	int N = floor(total_t / Tt);	//样本总数
	double r = 1;				//参考值
	double e = 0;
	double noise_temp = 0;
	vector<double> y(N+1, 0);
	vector<double> u(N+1, 0);
	NLSEFState_Roll.b1 = x[0];
	NLSEFState_Roll.b2 = x[1];
	NLSEFState_Roll.b = x[2];
	ESOParas_Roll.b = x[2];
	out_y_best << y[0] << endl;
	out_y_best << y[1] << endl;
	for (int i = 2; i <= N; i++) {
		noise_temp = 0.1 * rand_noise(rand_e);
		y[i] = 1.944 * y[i - 1] - 0.944 * y[i - 2] + 6.131e-05 * u[i - 1] + 6.015e-05 * u[i - 2] + noise_temp;

		TD_Atti(&TDState_RollRadio, r, &TD_fhanParas_RollRadio);
		ESO_Atti(y[i - 1], u[i - 1], &ESOParas_Roll, &ESOState_Roll);
		NLSEF_Atti(&TDState_RollRadio, &ESOState_Roll, &NLSEFState_Roll);
		u[i] = NLSEFState_Roll.u;
		out_y_best << y[i] << endl;
		out_noise << setprecision(5) << noise_temp << endl;
	}
	out_y_best.close();
	out_noise.close();
	cout << "最终参数为:" << endl;
	for (int i = 0; i < x.size(); i++)
	{
		cout << "x" << i << "=" << fixed << setw(12) << setprecision(5) << x[i] << endl;	//输出最优变量
	}
	cost_function(y);
}

void data_log_to_csv() {
	for (unsigned int j = 0; j < data_log[0].size(); j++)
		oFile << data_log[0][j] << ",";
	oFile << endl;
	for (unsigned int i = 1; i < data_log.size(); i = i+2) {
		for (unsigned int j = 0; j < data_log[i].size(); j++)
			oFile << data_log[i][j] << ",";
		oFile << endl;
	}
}
template<typename T1, typename T2>
T1 AMP_Limit(T1 value, T2 max, T2 min){
	T1 res;
	res = value > max ? max : value;
	res = res < min ? min : res;
	return res;
}
template<typename T>
T sign(T value) {
	if (value > 0)
		return 1;
	else if (value < 0)
		return -1;
	else
		return 0;
}
