
/*
 * y[] --> x[] is a frequency-warped lattice FIR filter which actually
 * becomes an IIR filter due to frequency warping. Its parameters are 'lam',
 * 's1', 's2' and the parcor coefficients 'k_i'. Instead of z^-1 type deley
 * elements the frequency-warped version uses a 1st order IIR all-pass 'D':
 *
 *    y ---|>--*---------*->(+)------*->(+)-- ... -------->(+)--|>--- x
 *         s1   \         \ /         \ /                  /    s2
 *               \     -k1 X       -k2 X              -kn /
 *                \       / \         / \                /
 *                 --[D]-*->(+)--[D]-*->(+)-- ... --[D]--
 *
 *                       -lam
 *    D:  i -->(+)-----*---->(+)-- o    o[n]   = (1-lam^2) t[n] -lam i[n]
 *              A      |      A         t[n+1] = i[n] +lam t[n]
 *          lam |    [z^-1]   |
 *              |      |      |
 *              +------*------+ (t)
 *
 * We choose the scale factors s1 and s2 so that
 *
 *    x[n] = 1.0 * y[n] + u[n]  where  u[n] = sum_i h_i * t_i[n]
 *           ~~~
 *
 * which is necessary for noise shaping. Such scale factors only depend on
 * the frequency warping parameter 'lam' and the parcor coefficients 'k_i'.
 * From the above equation the value for s1*s2 follows. But I'm not yet
 * sure about how to divide this factor to minimize artefacts when changing
 * filter parameters. For now, I simply set s1=1 and s2 to satisfy the above
 * equation.
 *
 * Note that we want to use a warped IIR filter (all-pole) as noise shaper
 * and not a warped FIR filter (all-zero) since all-pole filters with their
 * resonance-type responses are expected to be better suited for
 * approximating masking curves. So, y[] is actually the filtered noise
 * (overall error) and x[] is the (unfiltered) quantization error and the
 * noise shaping filter is actually the inverse of the above structure.
 * The reverse has delay-free loops, though. The above structure doesn't
 * have delay-free loops which makes it easier to work with.
 */

/* +------+
 * | TODO |
 * -------+
 * [ ] 16*order+3 FLOPS per sample is rather high. Can we get it faster?
 *     u could be calculated in terms of t and h. This should save about
 *     6 FLOPS per sample but it requires more precomputation (h). But
 *     this probably takes O(order^2) time.
 */

#include <algorithm>
#include <cmath>
#include "waplns.hpp"
#include "tools.hpp"

namespace { // anonymous

template<class T>
inline void lattice_step(T & a, T & b,   // 4 FLOPS
	typename identity<T>::type k)
{
	T const ak = a * k;
	a -= b*k;
	b -= ak;
}

inline void apply_D_alter_t(double & io, float & t, float lambda) // 4 FLOPS
{
	float next_t = io + lambda * t;
	io = t - lambda * next_t;
	t = next_t;
}

inline void apply_D_keep_t(double & io, float t, float lambda) // 4 FLOPS
{
	float next_t = io + lambda * t;
	io = t - lambda * next_t;
}

} // anonymous namespace

void waplns::precompute_derived_params()
{
	const double negative_lam = -lambda_;
	double a = 1;
	double b = 1;
	for (int i=0; i<order_; ++i) {
		b *= negative_lam;
		lattice_step(a,b,k_[i]);
	}
	s1_ = 1.0f;
	s2_ = static_cast<float>( 1.0/a );
}

void waplns::set_params(const float lam, int ord, float const* newk)
{
	int const neword = std::min(ord,max_wapl_filt_order);
	for (int i=0; i<neword; ++i) {
		this->k_[i] = newk[i];
	}
	for (int i=order_; i<neword; ++i) {
		t_[i] = 0;
	}
	this->lambda_ = lam;
	this->order_ = neword;
	precompute_derived_params();
	double nua = 0;
	double nub = 0;
	for (int i=0; i<order_; ++i) {
		const float k = this->k_[i];
		apply_D_keep_t(nub,t_[i],lam);
		lattice_step(nua,nub,k);
	}
	next_u_ = nua * s2_;
}

void waplns::reset_state()
{
	for (int i=0; i<order_; ++i) {
		t_[i] = 0;
	}
	next_u_ = 0;
}

void waplns::x_was(float x)  // 16 * order + 3 FLOPS
{
	// y + u = x  <=>  y = x - u
	double const y = static_cast<double>(x) - next_u_;
	double a = y * s1_;
	double b = a;
	double nua = 0;
	double nub = 0;
	const float lam = lambda_;
	for (int i=0; i<order_; ++i) {
		const float k = this->k_[i];
		apply_D_alter_t( b,t_[i],lam);
		apply_D_keep_t(nub,t_[i],lam);
		lattice_step(  a,  b,k);
		lattice_step(nua,nub,k);
	}
	next_u_ = nua * s2_;
	// next_u is only a linear combination of the ts which
	// could be computed with 2*order FLOPS instead of 8*order FLOPS
	// assuming we know the weights (not yet precomputed).
}

