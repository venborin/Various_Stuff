#include<cstdio>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<algorithm>
#include<string>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<iostream>
#define KE 28706.6897
#define PI 3.14159265359
float gauss(float a,float x,float p,float s)
{
/*
a - peak intensity (Oscillator strength)
x - energy scale (eV)
p - VET, eV
s - stdev
*/
	float e=(p-x)/s;
	return KE*a*exp(-0.5*e*e)/(s*sqrt(2.0*PI));
}
std::vector<std::vector<float> > convolute(float min, float max, std::vector<double> E, std::vector<double> I,float sigma)
{
	std::vector<std::vector<float> > spect;
	std::vector<float> tmp;
	float x=min;
	float t=0;
	float delta=0.01;
	while (x <= max+delta)
	{
		tmp.push_back(x);
		for (int i=0; i<E.size();i++)
		{
			t+=gauss(I[i],x,E[i],sigma);
	
		}
		tmp.push_back(t);
		spect.push_back(tmp);
		t=0;
		tmp.clear();
		x+=delta;
	}
	return spect;
}

struct VET
{
        float E;
        float Fo;
};


std::vector<VET> read_data(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open file: " + path);
	std::exit(1);
    }

    std::vector<VET> out;
    std::string line;
    std::size_t line_no = 0;

    while (std::getline(in, line)) {
        ++line_no;

        // Trim leading whitespace for comment/blank checks
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;          // blank line
        if (line[first] == '#') continue;                  // comment line

        // Allow comma-separated columns by converting commas to spaces
        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream iss(line);
        VET r{};
        if (!(iss >> r.E >> r.Fo)) {
            throw std::runtime_error(
                "Parse error on line " + std::to_string(line_no) + ": " + line
            );
        }

        // Optional: ensure no extra non-whitespace tokens
        std::string extra;
        if (iss >> extra) {
            throw std::runtime_error(
                "Too many columns on line " + std::to_string(line_no) + ": " + line
            );
        }

        out.push_back(r);
    }

    return out;
}

int main(int argc,char** argv)
{
//	Defaults:
	std::string fname="data.dat";
	std::string units="eV";
	float fwhm = 0.15;
	float mine = 0.01;
	float maxe = 5.00;
	
	for (int i = 1; i < argc; i++)
	{
		std::string str = argv[i];
		if (str=="-h" || str=="--help")
		{
			printf("Usage: %s -d VET-file -u units -f FWHM --min MIN_E --max MAX_E\n",argv[0]);
			printf("\t-d/--data:\t a file, containing calculated VETs (VET vs Oscillator strength\n");
			printf("\t-u/--units:\t units for the broadened spectrum: eV (default), cm, or nm\n");
			printf("\t-f/--fwhm:\t the value for the FWHM. Default is 0.15 eV\n");
			printf("\t--min and --max: min and max values on the energy scale (in eV). Defaults are 0.01 and 5.00\n");
			std::exit(0);
		}
		if (str=="-d" || str=="--data")
		{
			fname = argv[++i];
		}
		else if (str=="-u" || str=="--units")
		{
			units = argv[++i];
		}
		else if (str=="-f" || str=="--fwhm")
		{
			fwhm=atof(argv[++i]);
		}
		else if (str=="--min")
		{
			mine=atof(argv[++i]);
		}
		else if (str=="--max")
		{
			maxe=atof(argv[++i]);
		}
		else
		{
			printf("Unknown argument %s!\n",str);
		}

	}
	std::transform(units.begin(), units.end(), units.begin(),[](unsigned char c){ return std::tolower(c); });
	std::vector<VET> V = read_data(fname);
	float sigma=fwhm/2.355;
	std::vector<double> En;
	std::vector<double> Fo;
	for (const VET& r : V) {
	    En.push_back(r.E);
	    Fo.push_back(r.Fo);
	}
	if (mine <=0.0)
	{
		printf("Warning! Minimum energy should be strictly greater than zero. Resetting to 0.01\n");
		mine=0.01;
	}
	std::vector<std::vector<float> > spectrum = convolute(mine,maxe,En,Fo,sigma);

	if (units=="ev")
	{
		for (int i=0;i<spectrum.size();i++)
		{
			printf("%5.2f\t%10.1f\n",spectrum[i][0],spectrum[i][1]);
		}
	}
	else if (units=="cm")
	{
		for (int i=0;i<spectrum.size();i++)
		{
			printf("%6d\t%10.1f\n",(int)spectrum[i][0]*8066,spectrum[i][1]);
		}
	}
	else if (units=="nm")
	{
		for (int i=0;i<spectrum.size();i++)
		{
			printf("%5.1f\t%10.1f\n",(1240.0/spectrum[i][0]),spectrum[i][1]);
		}
	}


}
