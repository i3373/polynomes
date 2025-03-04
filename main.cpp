#include<iostream>
#include<cstdint>
#include"polynomes.h"

long long Polynom::GF = LLONG_MAX - 1;
bool Polynom::GFEnabled = false;

vector<Polynom> h, v, gs, edf;

long long power(long long a, long long exp) {
    long long result = a;
    for(int i = 1; i < exp; i++) {
        result *= a;
    }
    return result;
}


Polynom repeated_squaring(Polynom a, long long n, Polynom f){
    Polynom result = Polynom({1});
    while(n > 0) {
        if (n % 2 == 1) {
            result *= a;
        }
        a *= a;
        n /= 2;
    }
    return result % f;
}

vector<pair<Polynom, long long>> square_free_factorization(Polynom f) {
    vector<pair<Polynom, long long>> R(0);
    Polynom c = gcd(f, f.derivative());
    Polynom w = f / c;

    int i = 1;
    while(w != 1) {
        Polynom y = gcd(w, c);
        Polynom fac = w / y;
        R.push_back({fac, i});
        w = y;
        c /= y;
        i++;
    }
    if (c != 1) {
        c = c.rootGF();
        vector<pair<Polynom, long long>> newR = square_free_factorization(c);
        for(pair<Polynom, long long> p : newR) {
            R.push_back({p.first, Polynom::getGF() * p.second});
        }
    }
    return R;
}


Polynom equal_degree_splitting(Polynom f, int d) {
    for (Polynom a = Polynom::first(d, Polynom::getGF()); !a.isLast(); a.next()) {
        if (a.degree() == 0 || a.degree() == -1) {
            continue;
        }

        if (a.wasChecked()) {
            continue;
        }

        Polynom g1 = gcd(a, f);
        if (g1 != 1) {
            return g1;
        }    

        a.markChecked();
    }    
    return Polynom({0});
}


void equal_degree_factorization(Polynom f, int d) {
    if (d == f.degree()) {
        edf.push_back(f);
        return;
    }
    Polynom g = equal_degree_splitting(f, d);

    if (g == 0 || g == f) {
        edf.push_back(f);
        return;
    }

    equal_degree_factorization(g, d);
    equal_degree_factorization(f/g, d);
    return;
}

vector<pair<Polynom, long long>> factor_polynomial_over_finite_field(Polynom f) {
    h.push_back(Polynom({0, 1}));
    v.push_back(f / f.leading_coefficient());
    int i = 0;
    vector<pair<Polynom, long long>> U;
    while(v[i] != 1) {
        i++;
        h.push_back(repeated_squaring(h[i - 1], Polynom::getGF(), f));
        Polynom test = h[i] - h[0];
        Polynom test2 = h[i];
        Polynom test0 = h[0];
        Polynom v0 = v[i-1];
        Polynom g = gcd(h[i] - h[0], v[i - 1]);
        gs.push_back(g);
        if (g == 1) {
            if(v[i - 1].isIrreducible()) {
                break;
            } else {
                v.push_back(v[i - 1]);
            }
        }
        if (g != 1) {
            equal_degree_factorization(g, i);
            vector<Polynom> factors = edf;
            edf.clear();
            v.push_back(v[i - 1]);
            for (Polynom& g_j : factors) {
                int e = 0;
                while (v[i] % g_j == Polynom({0})) {
                    v[i] /= g_j;
                    e++;
                }
                U.push_back({g_j, e});
            }
        }
    }
    return U;
}

vector<pair<Polynom, long long>> factor(Polynom f) {
    vector<pair<Polynom, long long>> result = {};
    vector<pair<Polynom, long long>> R = square_free_factorization(f);
    for(pair<Polynom, long long> p : R) {
        vector<pair<Polynom, long long>> factors = factor_polynomial_over_finite_field(p.first);
        h.clear();
        gs.clear();
        v.clear();
        for(pair<Polynom, long long> factor : factors) {
           result.push_back({factor.first, factor.second * p.second});
        
        }
        factors.clear();
    }
    return result;
}


int main(){
    Polynom::setGF(2);
    Polynom::enableGF();
    vector<long long> z(12, 0);
    z.push_back(1);
    z[9] = 1;
    z[6] = 1;
    z[3] = 1;
    z[0] = 1;
    vector<long long> f(1025, 0);
    f.push_back(1);
    f[2] = -1;
    Polynom fP(f), zp(z);
    cout << fP << endl;
    //cout << zp % Polynom({1, 1, 0, 0, 1}) << endl;
    vector<pair<Polynom, long long>> result = factor(fP);
    cout << "Divisors of f"<< ":" << endl;
    for(pair<Polynom, long long> factor : result) {
        cout << "Factor: " << factor.first << " Amount: " << factor.second << endl;
    }
    
    return 0;
}