#include<iostream>
#include<cstdint>
#include"polynomes.h"

long long Polynom::GF = LLONG_MAX - 1;
bool Polynom::GFEnabled = false;

vector<Polynom> h, v, gs;

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

vector<Polynom> factor_polynomial_over_finite_field(Polynom f, long long q) {
    Polynom::setGF(q);
    return {};
}


int main(){
    gs.push_back(Polynom({1}));
    Polynom::enableGF();
    Polynom::setGF(5);
    Polynom a =  Polynom({1, 4, 6, 4, 1}, "a(x)");
    vector<long long> f(32, 0);
    f[2] = -1;
    f.push_back(1);
    Polynom fP(f);
    Polynom b =  Polynom({2, 2, 1}, "b(x)");
    cout << a << endl;
    cout << b << endl;
    if(a == b) {
        cout << "Equal";
    } else {
        cout << "Not Equal";
    }
    cout << endl;
    cout << Polynom(a + b, "a(x) + b(x)");
    cout << endl;
    cout << Polynom(a - b, "a(x) - b(x)");
    cout << endl;
    cout << Polynom(a * b, "a(x) * b(x)");
    cout << endl;
    cout << Polynom(a / b, "a(x) / b(x)");
    cout << endl;
    cout << Polynom(a % b, "a(x) % b(x)");
    return 0;
}