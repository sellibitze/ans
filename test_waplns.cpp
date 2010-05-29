#include <iostream>
#include "waplns.hpp"

const float k[] = {
	0.25, -0.15
};

int main()
{
	waplns ns;
	ns.set_params(0.5f, 2, k);
	std::cout << "warp_gain = " << ns.warp_gain() << '\n';
	// compute impulse response ...
	for (int i=0; i<=15; ++i) {
		float x = (i==0);
		float y = x - ns.u();
		std::cout << y << '\n';
		ns.x_was(x);
	}
}

