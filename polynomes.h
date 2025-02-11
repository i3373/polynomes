#include<iostream>
#include<vector>

using namespace std;

class Polynom {
    private:
        string name;
        vector<long long> coeffs;
        static long long GF;
        static bool GFEnabled;
    public:

        Polynom(){
            coeffs = vector<long long>();
            if(GFEnabled) normaliseGF();
        }
        Polynom(vector<long long> in, string name = "p(x)") {
            this->coeffs = in;
            this->name = name;
            if(GFEnabled) normaliseGF();
        }

        Polynom(Polynom const &other, string name = "p(x)"){
            this->coeffs = other.coeffs;
            this->name = name;
            if(GFEnabled) normaliseGF();
        }

        static void enableGF() {
            GFEnabled = true;
        }

        static void disableGF() {
            GFEnabled = false;
        }


        static void setGF(long long p) {
            GF = p;
        }

        void normaliseGF(){
            for(auto& k : coeffs){ 
                k = ((k % GF) + GF) % GF;
            }
        }

        void clearGF() {
            GF = INT32_MAX;
        }

        Polynom operator=(vector<long long> in) {
            this->coeffs = in;
            return *this;
        }

        Polynom operator=(Polynom const &other) {
            this->coeffs = other.coeffs;
            return *this;
        }


        long long& operator[](size_t i) {
            return coeffs.at(i);
        }

        bool operator==(Polynom const& other) {
            return equal(this->coeffs.begin(), this->coeffs.end(), other.coeffs.begin());
        }

        bool operator!=(Polynom const& other) {
            return !(*this == other);
        }

        Polynom operator+(Polynom const& other) {
            bool sizeComp = size(this->coeffs) < size(other.coeffs);
            Polynom result = sizeComp ? other : *this;
            Polynom lesserPolynom = sizeComp ? *this : other;
            for(long long i = 0; i <= lesserPolynom.power(); i++) {
                result[i] += lesserPolynom[i];
            }
            result.trimLeadingZeroes();
            if(GFEnabled) result.normaliseGF();
            return result; 
        }


        void operator+=(Polynom const& other) {
            *this = *this + other;    
        }

        Polynom operator-(Polynom const& other) {
            Polynom result = *this;
            if(result.coeffs.size() < other.coeffs.size()){
                result.coeffs.resize(other.coeffs.size(), 0);
            }
            for(long long i = 0; i < other.coeffs.size(); i++) {
                result[i] -= other.coeffs[i];
            }
            result.trimLeadingZeroes();
            if(GFEnabled) result.normaliseGF();
            return result; 
        }
        void operator-=(Polynom const& other) {
            *this = *this - other;    
        }

        Polynom operator* (Polynom const& other) {
            Polynom result(vector<long long>(power() + size(other.coeffs), 0));
            for(long long i = 0; i <= power(); i++) {
                for(long long j = 0; j < size(other.coeffs); j++){
                    result[i + j] = (result[i + j] + coeffs[i] * other.coeffs[j]);
                    if(GFEnabled) {
                        result[i + j] %= GF;
                    }
                }
            }
            result.trimLeadingZeroes();
            if(GFEnabled) result.normaliseGF();
            return result;
        }

        Polynom operator*(long long other) {
            return *this * Polynom({other});    
        }

        void operator*=(Polynom const& other) {
            *this = *this * other;    
        }

        Polynom operator/ (Polynom const& other) {
            return division(*this, other).first;
        }

        Polynom operator/ (long long number) {
            return *this / Polynom({number});
        }

        void operator /=(Polynom const& other) {
            *this = *this / other;    
        }

        Polynom operator% (Polynom const& other) {
            return division(*this, other).second;
        }

        void operator%=(Polynom const& other) {
            *this = *this % other;    
        }

        Polynom operator++() {
            return *this + Polynom({1});
        }

        Polynom operator--() {
            return *this - Polynom({1});
        }

        pair<Polynom, Polynom> division (const Polynom &divisable, const Polynom &divisor){
            if (divisable.coeffs.size() < divisor.coeffs.size()) {
                return {Polynom(), divisor};
            }
            Polynom div(vector<long long>(divisable.coeffs.size() - divisor.coeffs.size() + 1, 0));
            Polynom result = divisable, mod(vector<long long>(divisable.coeffs.size(), 0));
            long long dMP = divisor.coeffs.size() - 1; //divisorMaxPower
            long long inverse = calculateInverse(divisor.coeffs[dMP], GF-2, GF);
            

            for(long long i = divisable.coeffs.size() - 1; i >= 0 && dMP <= i; i--) {
                long long k1 = (result.coeffs[i] * inverse) % GF;
                div.coeffs[i - dMP] = GFEnabled ? k1 : result.coeffs[i] / divisor.coeffs[dMP];
                Polynom singlePolynom(vector<long long>(i - dMP, 0));
                singlePolynom.coeffs.push_back(GFEnabled ? 1 : result.coeffs[i] / divisor.coeffs[dMP]);
                result -=  GFEnabled ? ((singlePolynom * divisor) * k1) : singlePolynom * divisor;
                if(size(result.coeffs) == i + 1) {
                    mod.coeffs[i] += result.coeffs[i];
                    result.coeffs[i] = 0;
                }
                result.trimLeadingZeroes();
                i = result.coeffs.size();
            }
            mod += result;
            return {div, mod};
        }

        friend ostream& operator<<(ostream& os, Polynom p);

        long long calculateInverse(long long a, long long n, long long mod){
            long long result = 1;
            while(n > 0) {
                if (n % 2 == 1) {
                    result = (result * a) % mod;
                }
                a = (a * a) % mod;
                n /= 2;
            }
            return result;
        }

        void trimLeadingZeroes() {
            while (!coeffs.empty() && coeffs.back() == 0) {
                coeffs.pop_back();
            }
        }

        
        long long power() {
            return coeffs.size() - 1;
        }


        long long leading_coefficient() {
            return this->coeffs.at(power());
        }
        
        Polynom gcd(Polynom a, Polynom b) {
            while (!b.coeffs.empty()) {
                Polynom temp = a % b;
                a = b;
                b = temp;
            }
            return a;
        }

};

ostream& operator<<(ostream& os, Polynom p) {
        os << p.name << " = ";
        long long psize = size(p.coeffs) - 1;
        for(long long i = psize; i >= 0; i--){
            long long k = p.coeffs[i];
            if (k > 0) {

                i != psize ? os << "+ " : os << "";
            } else if (k < 0) {
                os << "- ";
            } else {
                continue;
            }
            k = abs(k);
            if (i == 1) {
                if(k == 1) {
                    os << "x ";
                    continue;
                }
                os << k << "x ";
                continue;
            }

            if (i == 0) {
                os << k;
                break;
            }
            k == 1 ? os << "x^" << i << " " : os << k << "x^" << i << " ";
        }
        return os;
}

