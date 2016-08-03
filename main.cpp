#include <iostream>
#include <fstream>
#include "const.h"
#include "func.h"
#include "adept.h"
#include "Transmittance.h"

int main() {
	std::ifstream in; 
	std::string file = "data/A457413.1.xml";
	std::string file2 = "in.txt";
	std::string file3 = "data/P20995.4.xml";
	std::string file4 = "data/A457200.1.xml";
	std::string file5 = "data/A473706.1.xml";
	std::ofstream out1("out1.txt");
	std::ofstream out2("out2.txt");
	std::ofstream out3("out3.txt");
	std::ofstream out4("out4.txt");
	Transmittance trans;
	//trans.readtrans(file);
	//trans.readtrans(file5);
	trans.readfile(file2);
	std::cout << "size " << trans.getSize() << '\n';
	trans.setThickness(97);
	//trans.readthick(file3);
	trans.setSubstrate(GLASS);
	std::cout << "thick " << trans.thickness << '\n';
	trans.print_data(std::cout);
	trans.optimization(3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 0.1, 0.1, 0.1, 10000);
	trans.print_index(std::cout);
	/*for (size_t i = 0; i < trans.getSize(); ++i) {
		std::cout << trans.wavelength[i] << ' ' << trans.n[i] << ' ' << trans.k[i] << ' ' << trans.sub[i] << ' ' \
			<< trans.calc[i] << ' ' << trans.TransTheory(trans.n[i], trans.k[i], trans.sub[i], trans.wavelength[i]) << '\n';
	}
	std::cout << trans.optimize(unknowns, 10000, out3) << '\n';
	trans.set_index(unknowns);
	for (size_t i = 0; i < trans.getSize(); ++i) {
		std::cout << trans.wavelength[i] << ' ' << trans.n[i] << ' ' << trans.k[i] << ' ' << trans.sub[i] << ' ' \
			<< trans.calc[i] << ' ' << trans.TransTheory(trans.n[i], trans.k[i], trans.sub[i], trans.wavelength[i]) << '\n';
	}
	trans.print_index(std::cout);
	trans.print_result(std::cout);
	trans.print_data(std::cout);
	//out1.close();
	out2.close();
	out3.close();
	out4.close();*/
	system("pause");
	return 0;
}