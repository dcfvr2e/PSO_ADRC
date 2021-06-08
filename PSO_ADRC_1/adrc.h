#ifndef __ADRC_H_
#define __ADRC_H_
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <algorithm>
using namespace std;


/*****	 fhan	 ******
@parameters: r ���� �ٶ����ӣ�һ��Ϊ��������ı���
			 c ���� �˲�����
*/
typedef struct _fhanparas
{
	float r;
	float h;
}fhanParas_TypeDef;


/******  TD  ******
@parameters: x1 ���� ���ٵ���
			 x2 ���� ���ٵ��ٶ���
*/
typedef struct _tdstate
{
	float h;
	float x1;
	float x2;
}TDState_TypeDef;


/******  ESOPara  ******
@parameters: b1��b2��b3 ���� �۲�������
						 b ���� ��������Ӧ�ڽǼ��ٶ���������
						 c ���� �˲�����
						 a1 ���� λ�÷���fal����ָ������
						 a2 ���� �ٶȷ���fal����ָ������
*/
typedef struct _esoparas
{
	float h;
	float b;
	float b1;
	float b2;
	float b3;
	float a1;
	float a2;
	float d;
}ESOParas_TypeDef;

/******  ESOState  ******
@parameters: z1��z2��z3 ���� �۲����۲��λ�á��ٶȡ����ٶ���

*/
typedef struct _esostate
{
	float z1;
	float z2;
	float z3;
}ESOState_TypeDef;

typedef struct _nlsefpara
{
	float beta1;
	float beta2;
	float b;
	float a1;
	float a2;
	float d;
	float u;
}NLSEFState_TypeDef;

/* TD */

extern fhanParas_TypeDef TD_fhanParas_RollRadio;

extern TDState_TypeDef TDState_RollRadio;

extern ESOParas_TypeDef ESOParas_Roll;
extern ESOState_TypeDef ESOState_Roll;

/* NLSEF */
extern NLSEFState_TypeDef NLSEFState_Roll;


void ADRC_Init(void);
void TD_Atti(TDState_TypeDef *state, float v, fhanParas_TypeDef *para);
void NLSEF_Atti(TDState_TypeDef *tdstate, ESOState_TypeDef *esostate, NLSEFState_TypeDef *nlsefstate);
void ESO_Atti(const double y, const double u, ESOParas_TypeDef *para, ESOState_TypeDef *state);
vector<double> ADRC_sim(vector<double> x_v);

template<typename T1, typename T2>
T1 AMP_Limit(T1 value, T2 max, T2 min);
template<typename T>
T sign(T value);
#endif
