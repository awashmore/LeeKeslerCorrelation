#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <cmath>
using namespace std;

const double R = 0.00008314; //Gas constant cm^3*bar/mol*K

double M(double x1, double x2, double y1, double y2, double m11, double m12, double m21, double m22, double Pr, double Tr) { //double interpolation formula
	double m, m2;
	m = (m11 * (x2 - Pr)) / (x2 - x1);
	m += (m12 * (Pr - x1)) / (x2 - x1);
	m *= (y2 - Tr) / (y2 - y1);
	m2 = (m21 * (x2 - Pr)) / (x2 - x1);
	m2 += (m22 * (Pr - x1)) / (x2 - x1);
	m2 *= (Tr - y1) / (y2 - y1);
	m += m2;
	return m;
}

double SingleInterpM(double m11, double m21, double x1, double x2, double pr) { //single interpolation formula
	double m;
	m = (m11 * (x2 - pr)) / (x2 - x1);
	m += (m21 * (pr - x1)) / (x2 - x1);
	return m;
}

double Z(double Z0, double Z1, double w) { //calculates z value
	double z = Z0 + (Z1 * w);
	return z;
}

double P(double T,double Z, double V) { //calculates pressure
	double p = (R*T*Z)/V;
	return p;
}

double Table(double p, double pc, double t, double tc, double w, double v) { //recursive function that uses the Z-tables
	ifstream fZ0, fZ1;
	fZ0.open("Z0.csv");
	fZ1.open("Z1.csv");

	double pr = p / pc; //row of Z table
	double tr = t / tc; //column of Z table

	double x1, x2, y1, y2, m11, m12, m21, m22, mZ0, mZ1;
	int countX = 0, countY = 0;

	bool flagX = false, flagY = false; //flags to see if the exact values for pr and tr are in the Z tables

	string line;
	regex delim(",");
		
	getline(fZ0, line);
	auto beginX = sregex_token_iterator(line.begin(), line.end(), delim, -1);
	auto endX = sregex_token_iterator();
	auto prevX = beginX;

	if (pr < 0 || pr > 10) { //checks to see if the pr and tr are within the constraints of the Z tables
		cout << "Invalid Pr" << endl;
		return -1;
	}
	if (tr < 0 || tr > 10) {
		cout << "Invalid Tr" << endl;
		return -1;
	}

	for (sregex_token_iterator data = beginX; data != endX; ++data) { //iterates until it finds the exact pr, or it finds pr's to interpolate between
		countX++;
		if (pr == stod(*data)) {
			x1 = pr;
			x2 = pr;
			flagX = true;
			break;
		}
		else if (pr < stod(*data)) {
			x2 = stod(*data);
			x1 = stod(*prevX); 
			break;
		}
		prevX = data;
	}

	double prevY = 0.0;
	while (getline(fZ0, line)) { //iterates until it finds the exact tr, or it finds tr's to interpolate between
		countY++;
		auto tempY = sregex_token_iterator(line.begin(), line.end(), delim, -1);

		if (tr == stod(*tempY)) {
			y2 = tr;
			y1 = tr;
			flagY = true;
			break;
		}
		else if (tr < stod(*tempY)) {
			y2 = stod(*tempY);
			y1 = prevY;
			break;
		}

		prevY = stod(*tempY);
	}

	auto ym2 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
	
	//Z0 table cases
	if (!(flagY || flagX)) { //case for double interpolation
		for (int i = 0; i < countX - 2; i++)
			ym2++;

		m21 = stod(*ym2);
		m22 = stod(*(++ym2));

		fZ0.clear();
		fZ0.seekg(0, ios::beg);

		for (int i = 0; i < countY; i++)
			getline(fZ0, line);

		auto ym1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 2; i++)
			ym1++;

		m11 = stod(*ym1);
		m12 = stod(*(++ym1));
		mZ0 = M(x1, x2, y1, y2, m11, m12, m21, m22, pr, tr);
	} 
	else if (flagX && flagY){ //case for pr and tr's exact values being there
		for (int i = 0; i < countX - 1; i++)
			ym2++;

		mZ0 = stod(*ym2);	
	}
	else if (flagY) { //case for tr being exact but not pr
		for (int i = 0; i < countX - 2; i++)
			ym2++;

		m21 = stod(*ym2);
		m22 = stod(*(++ym2));
		mZ0 = SingleInterpM(m21, m22, x1, x2, pr);
	}
	else { //case for pr being exact but not tr
		for (int i = 0; i < countX - 2; i++)
			ym2++;

		m21 = stod(*(++ym2));

		fZ0.clear();
		fZ0.seekg(0, ios::beg);

		for (int i = 0; i < countY; i++)
			getline(fZ0, line);

		auto ym1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 2; i++)
			ym1++;

		m11 = stod(*(++ym1));

		mZ0 = SingleInterpM(m11, m21, y1, y2, tr);
	}

	fZ0.clear();
	fZ0.seekg(0, ios::beg);		

	getline(fZ1, line);
	
	//Z1 table cases
	if (!(flagX || flagY)) { //case for double interpolation
		for (int i = 0; i < countY - 1; i++)
			getline(fZ1, line);

		auto beginYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);

		for (int i = 0; i < countX - 2; i++)
			beginYZ1++;

		m11 = stod(*beginYZ1);
		m12 = stod(*(++beginYZ1));

		getline(fZ1, line);
		auto nextYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 2; i++)
			nextYZ1++;
		m21 = stod(*nextYZ1);
		m22 = stod(*(++nextYZ1));

		mZ1 = M(x1, x2, y1, y2, m11, m12, m21, m22, pr, tr);
	}
	else if (flagY && flagX) { //case for both pr and tr being exact values
		for (int i = 0; i < countY; i++)
			getline(fZ1, line);

		auto nextYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 1; i++)
			nextYZ1++;

		mZ1 = stod(*nextYZ1);
	} //case for tr being exact but not pr
	else if (flagY) {
		for (int i = 0; i < countY; i++)
			getline(fZ1, line);

		auto nextYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 2; i++)
			nextYZ1++;

		m21 = stod(*nextYZ1);
		m22 = stod(*(++nextYZ1));

		mZ1 = SingleInterpM(m21, m22, x1, x2, pr);
	}
	else { //case for pr being exact but not tr
		for (int i = 0; i < countY - 1; i++)
			getline(fZ1, line);

		auto beginYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);

		for (int i = 0; i < countX - 1; i++)
			beginYZ1++;

		m12 = stod(*beginYZ1);

		getline(fZ1, line);
		auto nextYZ1 = sregex_token_iterator(line.begin(), line.end(), delim, -1);
		for (int i = 0; i < countX - 1; i++)
			nextYZ1++;

		m22 = stod(*(++nextYZ1));

		mZ1 = SingleInterpM(m12, m22, y1, y2, tr);
	}

	fZ1.clear();
	fZ1.seekg(0, ios::beg);

	double z = Z(mZ0, mZ1, w); //recalulates z value

	double pNew = P(t, z, v); //calculates new pressure
	
	if (fabs(p - pNew) < 0.0001) //checks to see if there is no difference between the previous pressure and the new pressure
		return pNew;
	else
		return Table(pNew, pc, t, tc, w, v);
}

int main() {
	ifstream fZ0, fZ1; //opens both Z tables to be read from
	fZ0.open("Z0.csv");
	fZ1.open("Z1.csv");

	double tc,pc,t,v,p,w,Z0,Z1,pa;
	double z = 1;
	cout << "Please input critical temperature(K): "; //gets all user inputs
	cin >> tc; "/n";
	cout << "Please input critical pressure(bar): ";
	cin >> pc; "/n";
	cout << "Please input fluid omega value: ";
	cin >> w; " /n";
	cout << "Please input system temperature(K): ";
	cin >> t; "/n";
	cout << "Please input system volume(cm^3/mol): ";
	cin >> v; "/n";
		
	p = P(t, 1, v); //calculates pressure with an inital Z value of 1

	double tableValue = Table(p, pc, t, tc, w, v); //runs the recursive algorithm using the Z tables

	if (fabs(1 + tableValue) < 0.0001) {//checks to see if the value returned is -1, if so then either pr or tr was invalid
		cout << "Invalid values, try again with proper values" << endl;
		return -1;
	}

	std::cout << "Pressure (bar): " << tableValue << endl; //prints out the final pressure

	fZ0.close();
	fZ1.close();
	
	return 0;
}
