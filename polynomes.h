#include<iostream>
#include<vector>

using namespace std;

class Polynom {
    private:
        vector<long long> coeffs;
        static long long GF;
        static bool GFEnabled;
    public:
        string name;
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

        static long long getGF() {
            return GF;
        }

        void normaliseGF(){
            for(auto& k : coeffs){ 
                k = ((k % GF) + GF) % GF;
            }
        }

        void clearGF() {
            GF = INT32_MAX;
        }

        Polynom& operator=(Polynom const &other) {
            if (this != &other) { 
                this->coeffs = other.coeffs;
            }
            return *this;
        }
    
        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator=(T other) {
            this->coeffs = {static_cast<long long>(other)};
            return *this;
        }


        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator=(vector<T> other) {
            this->coeffs = other;
            return *this;
        }

        long long& operator[](size_t i) {
            return coeffs.at(i);
        }

        bool operator==(Polynom const& other) {
            return equal(this->coeffs.begin(), this->coeffs.end(), other.coeffs.begin());
        }

        template<typename T, typename = enable_if_t<is_arithmetic_v<T>>>
        bool operator==(T other) {
            if(this->coeffs.size() == 1) {
                return this->coeffs[0] == static_cast<long long>(other) ? true : false;
            }
            return false;
        }

        bool operator!=(Polynom const& other) {
            return !(*this == other);
        }

        template<typename T, typename = enable_if_t<is_arithmetic_v<T>>>
        bool operator!=(T other) {
            if(this->coeffs.size() == 1) {
                return this->coeffs[0] == static_cast<long long>(other) ? false : true;
            }
            return true;
        }

        Polynom operator+(Polynom const& other) const{
            bool sizeComp = size(this->coeffs) < size(other.coeffs);
            Polynom result = sizeComp ? other : *this;
            Polynom lesserPolynom = sizeComp ? *this : other;
            for(long long i = 0; i <= lesserPolynom.degree(); i++) {
                result[i] += lesserPolynom[i];
            }
            result.trimLeadingZeroes();
            if(GFEnabled) result.normaliseGF();
            return result; 
        }


        Polynom& operator+=(Polynom const& other) {
            *this = *this + other;
            return *this;    
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom operator+(T other) const {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            return *this + otherP;
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator+=(T other) {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            *this = *this + otherP;
            return *this;
        }

        Polynom operator-(Polynom const& other) const{
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

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom operator-(T other) const {
            const Polynom otherP = Polynom({static_cast<long long>(other)});
            return *this - otherP;
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator-=(T other) {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            *this = *this - otherP;
            return *this;
        }

        Polynom operator* (Polynom const& other) const {
            Polynom result(vector<long long>(degree() + size(other.coeffs), 0));
            for(long long i = 0; i <= degree(); i++) {
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

        Polynom& operator*=(Polynom const& other) {
            *this = *this * other;
            return *this;    
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom operator*(T other) const {
            const Polynom otherP = Polynom({static_cast<long long>(other)});
            return *this * otherP;
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator*=(T other) {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            *this = *this * otherP;
            return *this;
        }

        Polynom operator/ (Polynom const& other) const{
            return division(*this, other).first;
        }

        Polynom operator/ (long long number) const{
            return *this / Polynom({number});
        }

        void operator /=(Polynom const& other) {
            *this = *this / other;    
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom operator/(T other) const {
            const Polynom otherP = Polynom({static_cast<long long>(other)});
            return *this / otherP;
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator/=(T other) {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            *this = *this / otherP;
            return *this;
        }

        Polynom operator% (Polynom const& other) const{
            return division(*this, other).second;
        }

        void operator%=(Polynom const& other) {
            *this = *this % other;    
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom operator%(T other) const {
            const Polynom otherP = Polynom({static_cast<long long>(other)});
            return *this % otherP;
        }

        template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
        Polynom& operator%=(T other) {
            Polynom otherP = Polynom({static_cast<long long>(other)});
            *this = *this % otherP;
            return *this;
        }

        Polynom operator++() {
            return *this + Polynom({1});
        }

        Polynom operator--() {
            return *this - Polynom({1});
        }

        pair<const Polynom, const Polynom> division (const Polynom &divisable, const Polynom &divisor) const{
            if (divisable.coeffs.size() < divisor.coeffs.size()) {
                return {Polynom(), divisable};
            }
            Polynom div(vector<long long>(divisable.coeffs.size() - divisor.coeffs.size() + 1, 0));
            Polynom result = divisable, mod(vector<long long>(divisable.coeffs.size(), 0));
            long long dD = divisor.coeffs.size() - 1; //divisorDegree
            long long inverse = calculateInverse(divisor.coeffs[dD], GF-2, GF);
            

            for(long long i = divisable.coeffs.size() - 1; i >= 0 && dD <= i; i--) {
                long long k1 = (result.coeffs[i] * inverse) % GF;
                div.coeffs[i - dD] = GFEnabled ? k1 : result.coeffs[i] / divisor.coeffs[dD];
                Polynom singlePolynom(vector<long long>(i - dD, 0));
                singlePolynom.coeffs.push_back(GFEnabled ? 1 : result.coeffs[i] / divisor.coeffs[dD]);
                result -=  GFEnabled ? ((singlePolynom * divisor) * k1) : singlePolynom * divisor;
                if(size(result.coeffs) == i + 1) {
                    mod.coeffs[i] += result.coeffs[i];
                    result.coeffs[i] = 0;
                }
                result.trimLeadingZeroes();
                i = result.coeffs.size();
            }
            mod += result;
            div.trimLeadingZeroes();
            mod.trimLeadingZeroes();
            return {div, mod};
        }

        friend ostream& operator<<(ostream& os, Polynom p);

        long long calculateInverse(long long a, long long n, long long mod) const{
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

        
        long long degree() const {
            return coeffs.size() - 1;
        }


        long long leading_coefficient() {
            return this->coeffs.at(degree());
        }

        friend Polynom gcd(Polynom a, Polynom b);
        friend Polynom random(long long degree);

        static Polynom first(int degree, int q) {
            vector<long long> coeffs(degree + 1, 0);
            coeffs[1] = 1;
            return Polynom(coeffs);
        }

        bool isLast() {
            for (long long c : coeffs) {
                if (c != GF - 1) return false;
            }
            return true;
        }
    
        void next() {
            for (size_t i = 0; i < coeffs.size(); i++) {
                coeffs[i]++;
                if (coeffs[i] < GF) return;
                coeffs[i] = 0;
            }
        }

        Polynom derivative() {
            Polynom result(vector<long long>(coeffs.size() - 1, 0));
            if (coeffs.size() <= 1) {
                return Polynom({0});
            }
            for(int i = 1; i < coeffs.size(); i++) {
                result.coeffs[i - 1] = coeffs[i] * i;
            }
            result.normaliseGF();
            result.trimLeadingZeroes();
            return result;
        }

        
        Polynom rootGF(){
            int p = Polynom::getGF();
            vector<long long> new_coeffs;
            
            for (size_t i = 0; i < coeffs.size(); i += p) {
                new_coeffs.push_back(coeffs[i]);
            }

            return Polynom(new_coeffs);
        }

};

Polynom random(long long degree) {
    Polynom a;
    for(int i = rand() % degree; i > 0; i--) {
        a.coeffs.push_back(rand() % a.GF);
    }
    a.trimLeadingZeroes();
    return a;
}

Polynom gcd(Polynom a, Polynom b) {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    while (!b.coeffs.empty()) {
        Polynom temp = a % b;
        a = b;
        b = temp;
    }
    if(a.degree() == 0) {
        a.coeffs[0] = 1;
    }
    return a / a.leading_coefficient();
}

ostream& operator<<(ostream& os, Polynom p) {
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

