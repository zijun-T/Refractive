#include "Refractive.h"

extern string Element;

string Refractive::inttoStr(int s)
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

void Refractive::measure(const char* EXP_file, ExpType EXP_type, const char* Flux_file)
{
	std::cout << "抱歉，还没有完成这方面的算法呢" << std::endl;
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
				std::cout << "Refractive index data loaded!" << std::endl;
				break;
			case REFLECTIVE_DATA://如果EXP_file文件里的第二列已经是反射率数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Reflectivity data loaded!" << std::endl;
				break;
			case NORMALIZED_DATA://如果EXP_file文件里第二列已经是相对光强数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Relative light intensity data loaded!" << std::endl;
				break;
			case PHASESHIFT_DATA://如果EXP_file文件里第二列已经是相位变化数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Phase change data loading complete!" << std::endl;
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