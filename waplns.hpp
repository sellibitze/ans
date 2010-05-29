#ifndef WAPLNS_HPP_INCLUDED
#define WAPLNS_HPP_INCLUDED

#include <cassert>

const int max_wapl_filt_order = 32;

/**
 * WAPLNS = warped all-pole lattice noise shaper
 *
 *    s : signal (to be quantized)
 *    q : quantized signal
 *    x : (unfiltered) quantization error
 *    y : filtered quantization error: q = s + y
 *    u : difference between x and y: x = y + u
 *    w : w = s - u
 *
 * Here's a signal flow graph that illustrates how noise shaping is
 * supposed to be done with this kind of filter:
 *
 *                       x
 *                       |
 *             (w)       V
 *    s -->(+)------*-->(+)---*-- q
 *          A       |         |
 *        - |       |         V
 *          |       +------->(+)
 *          |              -  |
 *      (u) |                 | (x)
 *          +--[z^-1]<--[H]<--+
 *
 * An object of the class waplns represents H combined with the unit delay.
 * You feed it with the unfiltered error (x) via x_was(). It computes and
 * remembers the next 'u' for the next quantization step. You get it
 * via u(). Noise shaping is as simple as the following pseudo loop:
 *
 *    for (i=0; i<count; ++i) {
 *       float w = s[i] - ns.u();
 *       q[i] = round( w + dither() );
 *       ns.x_was( q[i] - w );
 *    }
 *
 * To fight clipping errors you may want to do something like this:
 *
 *    for (i=0; i<count; ++i) {
 *       float w = s[i] - ns.u();
 *       q[i] = clipping_round( w + dither() );
 *       float x = q[i] - w;
 *       // restrict magnitude of x to prevent overload
 *       if (x<-thresh) x=-thresh; else if (thresh<x) x=thresh;
 *       ns.x_was(x);
 *    }
 *
 * The filter that turns x into y ("shapes x") is a frequency-warped
 * all-pole lattice filter. It is parameterized by order, k[i] (parcor
 * coefficients for 0 <= i < order) and a warping parameter lambda.
 */
class waplns
{
	// input filter parameters ...
	int order_;
	float lambda_;
	float k_[max_wapl_filt_order];

	// derived filter parameters ...
	float s1_;
	float s2_;

	// filter state
	float t_[max_wapl_filt_order];
	float next_u_;

	void precompute_derived_params();

public:
	waplns() : order_(0), lambda_(0), s1_(1), s2_(1), next_u_(0) {}

	float warp_gain() const { return 1.0 / s1_ / s2_; }
	int order() const {return order_;}
	float lambda() const {return lambda_;}
	float k(int idx) const {assert(0<=idx && idx<order_); return k_[idx];}

	void set_params(float lam, int ord, float const* newk);
	void set_params(float lam, int ord, float* newk)
	{ set_params(lam,ord,static_cast<float const*>(newk)); }
	template<class Iter>
	void set_params(float lam, int ord, Iter it);

	void reset_state();
	float u() const { return next_u_; }
	void x_was(float x);

	//void show_state(std::ostream &);
};

template<class Iter>
void waplns::set_params(float lam, int ord, Iter it)
{
	float temp[max_wapl_filt_order];
	ord = std::min(ord,max_wapl_filt_order);
	for (int i=0; i<ord; ++i) {
		temp[i] = *it;
		++it;
	}
	set_params(lam,ord,static_cast<float const*>(temp+0));
}

#if 0
inline void waplns::show_state(std::ostream & cout)
{
	cout << "k   = [";
	for (int i=0; i<order_; ++i) {
		cout << ' ' << k_[i];
	}
	cout << " ]\nlam = " << lambda_
		<< "\ns1  = " << s1_
		<< "\ns2  = " << s2_
		<< "\nt   = [";
	for (int i=0; i<order_; ++i) {
		cout << ' ' << t_[i];
	}
	cout << " ]\n";
	cout.flush();
}
#endif

#endif // WAPLNS_HPP_INCLUDED

