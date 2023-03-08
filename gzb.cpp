#include "Refractive.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string.h>
//#include <mpi.h>123

using namespace std;

string Element;

string inttoStr(int s)
{
	string c = "";

	int m = s;
	int n=0;

	while(m>0)
	{
		n++;
		m/=10;
	}

	for (int i = 0; i < n; i++)
	{
			c = (char)(s % 10 + 48) + c;
			s /= 10;
	}
	return c;

}

double lifetime(double rho)
{
	double t = 0;
	double r, p1, p2, p3, p4, p5, p6, p7;
	r = rho * 1e6;
	//cout << r <<endl;
	
	if (Element == "Si")
	{
		p1 = -0.001098;
		p2 = 0.05413;
		p3 = -1.03;
		p4 = 8.811;
		p5 = -31.2;
		
		if (rho < 1e12)
		{
			r = 1e12 * 1e6;
		}
		else if(rho > 1e20)
		{
			r = 1e20 * 1e6;
		}

		r = log10(r);
		t = p1*r*r*r*r + p2*r*r*r + p3*r*r + p4*r +p5;
		t = pow(10,t);
	}
	else if (Element == "GaP")
	{
		p1 = 0.07434;
		p2 = -7.7;
		p3 = 332;
		p4 = -7628;
		p5 = 9.851e+04;
		p6 = -6.779e+05;
		p7 = 1.942e+06;

		if (rho < 3.893418e+016)
		{
			r = 3.893418e+016 * 1e6;
		}
		else if(rho > 4.915596e+018)
		{
			r = 4.915596e+018 * 1e6;
		}

		r = log10(r);
		t = p1*r*r*r*r*r*r + p2*r*r*r*r*r + p3*r*r*r*r + p4*r*r*r + p5*r*r + p6*r +p7;
		t = pow(10,t);
	}
	else
	{
		cout << Element + "lifetime loss!!!" << endl;
		exit(0);
	}

	//cout << t <<endl;
	////exit(0);
	return t;
}

Refractive::Refractive(const char* EXP_file, ExpType file_type)
:fileType(file_type)
{
	std::ifstream input(EXP_file);
	std::string line;
	std::vector < std::vector<double>> Data;
	double code;
	
	if (input)//如果有文件
	{
		
		while (std::getline(input, line))//把文件里的内容复制到Data暂时储存
		{
			std::istringstream sstream(line);
			std::vector <double> row;
			while (sstream >> code)
				row.push_back(code);
			Data.push_back(row);
		}

		//在这考虑的是：因为文件里的时间间隔是相等的，所以用文件中的平均时间间隔作为我们的时间间隔；
		this->dt = (Data[Data.size()-1][0] - Data[0][0]) / (Data.size()-1);

		//那么我们就可以将数据按顺序填入
		switch (this->fileType)
		{
			case REFRACTIVE_DATA: //如果EXP_file文件里的第二列已经是折射率数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "折射率数据加载完毕！" << std::endl;
				break;
			case REFLECTIVE_DATA://如果EXP_file文件里的第二列已经是反射率数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "反射率数据加载完毕！" << std::endl;
				break;
			case NORMALIZED_DATA://如果EXP_file文件里第二列已经是相对光强数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "相对光强数据加载完毕！" << std::endl;
				break;
			case PHASESHIFT_DATA://如果EXP_file文件里第二列已经是相位变化数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "相位变化数据加载完毕！" << std::endl;
				break;
			default:
				break;
		}

		//初始化基本参数
		//初始X射线通量J/m^2
	}
	else
	{
		std::cout << "Data文件没找到或没数据！请再查看你的文件名！" << std::endl;
	}
}

void Refractive::measure(const char* EXP_file, ExpType EXP_type, const char* Flux_file)
{
	std::cout << "抱歉，还没有完成这方面的算法呢" << std::endl;
}

void Refractive::Progress(const char* Output_file, int jishu, double Energy,double Sample_thickness, double dx)
{
	ifstream filein;
	filein.open("ini.txt");
	if (filein.fail())
	{
	    cout << "read ini.txt failed" <<endl;
	    exit(0);
	}

	char* p = new char[24];
    double dataPos[20];
    int num1, num2;
    double Eg;
    double a[3];

    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num1 = int(dataPos[0]);
    for (int i = 0; i < num1; i++)
    {
        filein >> p;
    }
    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num2 = int(dataPos[0]);
    for (int i = 0; i < num2; i++)
    {
        filein >> p;
    }
    filein >> p;

    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->dx = double(dataPos[0]); //材料的厚度微元(m)
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->d = double(dataPos[0]); //材料厚度(m)
    unsigned int thick_size = (unsigned int)(this->d / this->dx); //计算有多少个薄片

    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->sigma = double(dataPos[0]); //每个分子的总反应截面(m^2)
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    Eg = double(dataPos[0]); //带隙(eV)
    this->alpha = Energy / Eg;
    for (int i = 0; i < 3; i++)
    {
        filein >> p; sscanf(p, "%lf", &dataPos[0]);
        a[i] = double(dataPos[0]);//E_gap参数
    }
    filein >> p;

    cout<<  a[1]<<endl;





    filein.close();

    double G_e = this->C_e / (this->C_e + this->C_l) / this->tao_e_l;
	double G_l = this->C_l / (this->C_e + this->C_l) / this->tao_e_l;
	double m_r = this->m_c * this->m_v / (this->m_c + this->m_v);

	double tao_e = this->tao_e_p / this->omega_detec;
	double tao_h = this->tao_h_p / this->omega_detec;

	this->E_pump =Energy * e_0;

	double Phy = 0;//变化的相位

	std::vector<double> n_electron_new;
	std::vector<double> T_e_new;
	std::vector<double> T_l_new;

	std::vector<double> n_new;
	std::vector<double> k_new;

	for (unsigned int i = 0; i < thick_size+1; ++i)
	{
		//初始自由电子密度设为0/m^3
		this->n_electron.push_back(0);
		n_electron_new.push_back(0);
		//初始电子温度和晶格温度设为300K
		this->T_e.push_back(300);
		T_e_new.push_back(300);
		this->T_l.push_back(300);
		T_l_new.push_back(300);

		n_new.push_back(0);
		k_new.push_back(0);
	}


	//if (this->dt >= (this->timel))
	{
		std::ofstream output((Element+string("_")+string(inttoStr(int(Energy)))+string("_")+string(inttoStr(jishu+1))+string(".dat")).c_str());
		//目前由于时间间隔dt>1ps, 所以假设T_l和T_e是同步的，
		for (unsigned int t = 0; t < EXP_Data.size(); ++t)
		{
			Phy = 0;//相位重置

			for (unsigned int i = 0; i < thick_size; ++i)
			{
				//1,2,双温模型
				//不是边界的情况
				if (i != 0 && i != (thick_size - 1))
				{
					T_l_new[i] = (this->k_l * (this->T_l[i + 1] - 2 * this->T_l[i] + this->T_l[i - 1]) / this->dx / this->dx
						+ G_l * (this->T_e[i] - this->T_l[i])) * this->dt
						+ this->T_l[i];

					T_e_new[i] = (this->k_e * (this->T_e[i + 1] - 2 * this->T_e[i] + this->T_e[i - 1]) / this->dx / this->dx
						+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
						- G_e * (this->T_e[i] - this->T_l[i])) * this->dt
						+ this->T_e[i];
				}
				//在边界上的情况(第一类边界条件)
				else
				{
					if (i == 0)
					{
						T_l_new[i] = (this->k_l * (this->T_l[i + 1] - 2 * this->T_l[i] + 300.) / this->dx / this->dx
							+ G_l * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_l[i];

						T_e_new[i] = (this->k_e * (this->T_e[i + 1] - 2 * this->T_e[i] + 300.) / this->dx / this->dx
							+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
							- G_e * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_e[i];
					}
					else
					{
						T_l_new[i] = (this->k_l * (300. - 2 * this->T_l[i] + this->T_l[i - 1]) / this->dx / this->dx
							+ G_l * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_l[i];

						T_e_new[i] = (this->k_e * (300. - 2 * this->T_e[i] + this->T_e[i - 1]) / this->dx / this->dx
							+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
							- G_e * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_e[i];
					}
				}

				//3,电子密度
				this->gamma = 1 / lifetime(n_electron[i]);
				//this->gamma = 1/1.26041e-06;

				n_electron_new[i] = (-this->gamma * this->n_electron[i]
					+ this->alpha * this->n_s * this->sigma / this->E_pump * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)) * this->dt
					+ this->n_electron[i];

				//4,带隙
				double E_gap = (a[0] - a[1] * T_l_new[i] * T_l_new[i] / (-T_l_new[i] + a[2])) * e_0;

				//5,6,Ec和Ev
				double E_c = E_gap / 2 + m_r * (hbar * this->omega_detec - E_gap) / this->m_c;
				double E_v = -E_gap / 2 - m_r * (hbar * this->omega_detec - E_gap) / this->m_v;

				//7,费米能级E_F
				double E_F = (E_c + E_v) / 2 + 3. / 4 * kB * T_l_new[i] * log(m_v / m_c);

				//8,F(E_v)-F(E_c)
				double F_F = 1. / (1 + exp((E_v - E_F) / kB / T_l_new[i])) - 1. / (1 + exp((E_c - E_F) / kB / T_l_new[i]));

				//9, <v|p|c>^2
				double VPC = m_e * m_e * E_gap / 2 / this->m_c;

				//10, k_interband
				double sqrt_term = (((hbar * this->omega_detec) > E_gap) ? sqrt(hbar * this->omega_detec - E_gap) : 0);
				double k_interband = e_0 * e_0 * pow(2 * m_r, 3. / 2) / m_e / m_e / this->n_core / hbar / hbar / hbar / this->omega_detec / this->omega_detec / eps_0
					* VPC * sqrt_term * F_F;

				//11，eps_core
				double eps_core_Re = this->n_core * this->n_core - k_interband * k_interband;
				double eps_core_Im = 2 * this->n_core * k_interband;

				//12, 电子和空穴的频率
				double omaga_pe2 = n_electron_new[i] * e_0 * e_0 / this->m_c / eps_0;
				double omaga_ph2 = n_electron_new[i] * e_0 * e_0 / this->m_v / eps_0;

				//13, eps
				double eps_Re = eps_core_Re - omaga_pe2 * this->omega_detec * this->omega_detec * tao_e * tao_e / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_e * tao_e + 1)
											- omaga_ph2 * this->omega_detec * this->omega_detec * tao_h * tao_h / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_h * tao_h + 1);
				double eps_Im = eps_core_Im + omaga_pe2 * this->omega_detec * tao_e / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_e * tao_e + 1)
											+ omaga_ph2 * this->omega_detec * tao_h / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_h * tao_h + 1);

				//14,n和k
				double theta = atan(eps_Im / eps_Re) / 2;
				double radios = sqrt(eps_Re * eps_Re + eps_Im * eps_Im);
				n_new[i] = sqrt(radios) * cos(theta);
				k_new[i] = sqrt(radios) * sin(theta);

				//更新数据
				this->n_electron[i] = n_electron_new[i];
				this->T_e[i] = T_e_new[i];
				this->T_l[i] = T_l_new[i];

				//相位
				Phy += 2 * this->omega_detec * (n_new[i]-this->n_core) * this->dx / c_0;
			}
			//反射率
			double R = ((n_new[0] - 1) * (n_new[0] - 1)) / ((n_new[0] + 1) * (n_new[0] + 1));

			output << t * dt << " "  << n_new[0] << " " << 1/this->gamma << std::endl;

			//std::cout << R << std::endl;

		}
		output.close();
	}
	n_electron_new.clear();
	T_e_new.clear();
	T_l_new.clear();

	n_new.clear();
	k_new.clear();
}

int main()
{
	Refractive* MineRefractive = new Refractive("multi.txt", PHASESHIFT_DATA);


	MineRefractive->omega_detec = 2 * Pi * c_0 / (1545.03e-9);

	ifstream filein;
	filein.open("ini.txt");//
	if (filein.fail())
	{
	    cout << "read ini.txt failed" <<endl;
	    exit(0);
	}

    char* p = new char[24];
    double dataPos[20];
    int num1, num2;
    double mu[100];
	double si[100];
	double energy[100];

    filein >> p; Element = p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num1 = int(dataPos[0]);
    for (int i = 0; i < num1; i++)
    {
        filein >> p; sscanf(p, "%lf", &dataPos[0]);
        mu[i] = double(dataPos[0]);
        si[i] = mu[i]/2.355;
    }
    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num2 = int(dataPos[0]);
    for (int i = 0; i < num2; i++)
    {
        filein >> p; sscanf(p, "%lf", &dataPos[0]);
        energy[i] = double(dataPos[0]);
    }

    filein.close();


   /* int rank=1;
    int size=1;

     MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int i = 0; i < 16;i++)
    {
        if (i % size1 == rank1)
        {
            for (unsigned int t = 0; t < MineRefractive->EXP_Data.size(); ++t)
            {

                MineRefractive->Flux_Data[t] = 0.03 * 10000 * pow(E,-(t-mu[i])*(t-mu[i])/(2*si[i]*si[i]))/(2.50662827463*si[i])/1e-12;
            }

            	MineRefractive->Progress("out.dat", i);
        }

    }

    MPI_Finalize();

    */

    for (int i = 0; i < num1;i++)
    {
        for (int j = 0; j < num2; j++)
        {
            for (unsigned int t = 0; t < MineRefractive->EXP_Data.size(); ++t)
            {

                MineRefractive->Flux_Data[t] = 0.03 * 10000 * pow(E,-(t-mu[i])*(t-mu[i])/(2*si[i]*si[i]))/(2.50662827463*si[i])/1e-12;
            }

            MineRefractive->Progress("out.dat", i, energy[j]);
        }
    }


	std::cout << "Finish!!!!!!" << std::endl;

	delete MineRefractive;
}

