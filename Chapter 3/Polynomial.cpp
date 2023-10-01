#include "temp.cpp"

tcT> struct Polynomial {

    int deg;
    V<T> coeffs;

    // Constructors
    Polynomial(V<T> const& _v) {
        coeffs = _v;
        if (coeffs.empty()) coeffs.pb(0);
        deg = sz(coeffs) - 1;
    }
    Polynomial(V<pair<int, T>> const& _v) {
        each(p, _v) {
            if (sz(coeffs) <= p.f) {
                coeffs.rsz(p.f + 1);
            }
            coeffs[p.f] += p.s;
        }
        deg = sz(coeffs) - 1;
    }
    
    // Calculate the derivative
    Polynomial derivative() const {
        V<T> aux(deg);
        for (int i = 1; i < sz(coeffs); ++i) aux[i - 1] = coeffs[i] * T(i);
        return Polynomial(aux);
    }

    // Evaluate a polynomial
    T evaluate(T x) const {
        T ans = 0;
        R0F(i, sz(coeffs)) ans = ans * x + coeffs[i];
        return ans;
    }

    T evaluate_difference(V<T> const& nodes, T x) const {
        T ans = 0, val = 1;
        assert(sz(nodes) >= sz(coeffs));
        F0R(i, deg + 1) {
            ans += coeffs[i] * val;
            val *= (x - nodes[i]);
        }
        return ans;
    }

    T coef(int N) const { return N > deg ? 0 : coeffs[N]; }

};
