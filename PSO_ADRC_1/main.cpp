#include "pch.h"
#include "adrc.h"

RandomNumber r;       //�����
int main()
{
	clock_t startTime, endTime; //�������ʼ����ʱ��ͽ���ʱ��
	startTime = clock();  //��ʱ��ʼ
	pso PSO;  //����PSO��ز����ͺ���

	//������ʼ������Ⱥ�ٶȺ�λ��,��������Ⱥ��Ӧ�ȣ���ʼ�����弫ֵ��ȫ�ּ�ֵ
	ADRC_Init();					//ADRC init()
	PSO.Initialize_fit_extremum();

	//������Ҫѭ�������չ�ʽ���ε�����ֱ�����㾫��Ҫ��
	PSO.Optimization_iteration();

	endTime = clock();//��ʱ����
	cout << "run time:" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}