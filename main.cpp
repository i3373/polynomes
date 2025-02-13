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

Polynom equal_degree_splitting(Polynom f, int d) {
    Polynom a = random(f.degree());
    if(a.degree() == 0) {
        return Polynom({0});
    }

    Polynom g1 = gcd(a, f);

    if(g1 != 1) {
        return g1;
    }

    Polynom b = repeated_squaring(a, (power(a.getGF(), d) - 1) / 2, f);
    Polynom g2 = gcd(b - 1, f);
    if(g2 != 1 && g2 != f) {
        return g2;
    } else {
        return Polynom({0});
    }
}


void equal_degree_factorization(Polynom f, int d) {
    if (d == f.degree()) {
        edf.push_back(f);
    }
    Polynom g({0});
    while(g == 0){
        g = equal_degree_splitting(f, d);
    }

    equal_degree_factorization(g, d);
    equal_degree_factorization(f/g, d);
    return;
}

vector<pair<Polynom, long long>> factor_polynomial_over_finite_field(Polynom f, long long q) {

    h.push_back(Polynom({0, 1}));
    v.push_back(f / f.leading_coefficient());
    int i = 0;
    vector<pair<Polynom, long long>> U;
    while(v[i] != 1) {
        i++;
        h.push_back(repeated_squaring(h[i - 1], q, f));
        Polynom test = h[i] - h[0];
        Polynom test2 = h[i];
        Polynom tesst0 = h[0];
        Polynom g = gcd(h[i] - h[0], v[i - 1]);
        if (g == 1) {
            break;
        }
        if (g != 1) {
            equal_degree_factorization(g, i);
            vector<Polynom> factors = edf;
            edf.clear();
            v.push_back(v[i - 1]);
            for (Polynom& g_j : factors) {
                int e = 0;
                while (v[i] % g_j == 0) {
                    v[i] /= g_j;
                    e++;
                }
                U.push_back({g_j, e});
            }
        }
    }
    return U;
}


int main(){
    gs.push_back(Polynom({1}));
    Polynom::setGF(31);
    Polynom::enableGF();
    Polynom a =  Polynom({1, 4, 6, 4, 1}, "a(x)");
    vector<long long> f(31, 0);
    f.push_back(1);
    Polynom fP(f);
    Polynom b =  Polynom({2, 2, 1}, "b(x)");
    cout << gcd(fP - Polynom({0, 1}), a) << endl;
    //vector<pair<Polynom, long long>> factors = factor_polynomial_over_finite_field(a, 31);
    cout << a << endl;
    cout << b << endl;
    if(a != 1) {
        cout << "Equal";
    } else {
        cout << "Not Equal";
    }
    cout << endl;
    a += 2;
    cout << Polynom(a + 2, "a(x) + b(x)");
    cout << endl;
    a -= 2;
    cout << Polynom(a - 2, "a(x) - b(x)");
    cout << endl;
    a *= 2;
    cout << Polynom(a * 2, "a(x) * b(x)");
    cout << endl;
    a /= 2;
    cout << Polynom(a / 2, "a(x) / b(x)");
    cout << endl;
    a %= 2;
    cout << Polynom(a % 2, "a(x) % b(x)") << endl;
    // for(int i =0; i < 10; i++) {
    //     cout << power(2, 3) << endl;
    // }
    
    return 0;
}