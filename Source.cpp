#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <exception>
#include <algorithm>
using namespace std;

class Pn{
public:
    double x, y; int num; double price;
    Pn(double a, double b, int n, double d) : x(a), y(b), num(n), price(d) {}
    friend double rast(Pn a, Pn b);
};

double rast(Pn a, Pn b){
    return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

bool checkp(const vector <unsigned int>& base){
    for(int i = 0; i < base.size(); i++){
        if(base[base[i] - 1] == i + 1){
            return false;
        }
    }
    return true;
}

double pricep(vector <unsigned int>& permutations, vector <Pn>& cords){
    double result = 0;
    for(int i = 0; i < permutations.size(); i++){
        result += rast(cords[i + 1], cords[permutations[i] - 1]);
    }
    return result;
}
 
double brutforcem(vector <Pn> cords){
    double buf;
    vector <unsigned int> permutations(cords.size());
    for(int i = 0; i < cords.size(); i++){
        permutations[i] = i + 1;
    }
    do{
        if(checkp(permutations))  buf = pricep(permutations, cords);
    }while(next_permutation(permutations.begin(), permutations.end()));
    return buf;
}

class VRP{
private:
    vector < vector < double >> base_buf; vector < vector < double >> base;
    vector <double> minl; vector <double> mincol;
    int ribcnt = 0; double finalprice = 0;
    vector <bool> bline; vector <bool> bcol;
    vector <int> cycle; vector <Pn> cordsb;
public:
    vector <int> return_cycle(){
        return cycle;
    }
    bool is_finish(){
        if (cordsb.size() == 2){
            cycle = { cordsb[0].num, cordsb[1].num };
            finalprice = 2 * rast(cordsb[0], cordsb[1]);
        }
        if(ribcnt == base.size())  return true;
        else  return false;
    }
    double retprice(){
        return finalprice;
    }
    VRP(vector <Pn>& cords){
        cordsb = cords;
        base.resize(cords.size(), vector<double>(cords.size(), -1));
        for(int i = 0; i < cords.size(); i++){
            for(int k = 0; k < cords.size(); k++){
                if(i != k) base[k][i] = rast(cords[k], cords[i]);
            }
        }
        minl.resize(cords.size(), -1);
        mincol.resize(cords.size(), -1);
        bline.resize(cords.size(), false);
        bcol.resize(cords.size(), false);
        base_buf = base;
    }
    void linred(){
        double cnt;
        for(int i = 0; i < base.size(); i++){
            if(bline[i])  continue;
            cnt = -1;
            for(int k = 0; k < base.size(); k++){
                if((i != k) && (base[i][k] >= 0) && ((cnt == -1) || (base[i][k] < cnt)))  cnt = base[i][k];
            }
            minl[i] = cnt;
            for(int k = 0; k < base.size(); k++){
                if((i != k) && (base[i][k] > 0))  base[i][k] -= cnt;
            }
        }
    }
    void colred() {
        double cnt;
        for(int i = 0; i < base.size(); i++){
            if(bcol[i]) continue;
            cnt = -1;
            for(int k = 0; k < base.size(); k++){
                if((i != k) && (base[k][i] >= 0) && ((cnt == -1) || (base[k][i] < cnt)))  cnt = base[k][i];
            }
            mincol[i] = cnt;
            for(int k = 0; k < base.size(); k++){
                if((i != k) && (base[k][i] > 0))  base[k][i] = base[k][i] - cnt;
            }
        }
    }
    double minn(int i, int k){
        double res = 0;
        double locmin = -1;
        for(int p = 0; p < base.size(); p++){
            if(((locmin == -1) || (locmin > base[i][p])) && (base[i][p] >= 0) && (p != k)) locmin = base[i][p];
        }
        res = res + locmin;
        locmin = -1;
        for(int p = 0; p < base.size(); p++){
            if(((locmin == -1) || (locmin > base[p][k])) && (base[p][k] >= 0) && (p != i)) locmin = base[i][p];
        }
        res = res + locmin;
        return res;
    }
    double maxx(int i, int k) {
        double res = 0;
        double locmax = 1;
        for (int p = 0; p < base.size(); p++) {
            if (((locmax == 1) || (locmax < base[i][p])) && (base[i][p] <= 0) || (p != k))  locmax = base[i][p];
        }
        res += locmax;
        locmax = 1;
        for (int p = 0; p < base.size(); p++) {
            if (((locmax == 1) || (locmax < base[p][k])) && (base[p][k] <= 0) && (p != i))  locmax = base[i][p];
        }
        res += locmax;
        return res;
    }
    void grc(){
        double cnt = 0;
        double ribeprice = 0;
        pair <int, int> buf;
        for(int i = 0; i < base.size(); i++){
            for(int k = 0; k < base.size(); k++){
                if(base[i][k] == 0){
                    if(minn(i, k) > cnt){
                        buf = make_pair(i, k);
                        ribeprice = base_buf[i][k]*2.4136;
                        if(find(cycle.begin(), cycle.end(), cordsb[i].num) == cycle.end()) cycle.push_back(cordsb[i].num);
                        if(find(cycle.begin(), cycle.end(), cordsb[k].num) == cycle.end()) cycle.push_back(cordsb[k].num);
                    }
                }
            }
        }
        ribcnt++;
        base[buf.first][buf.second] = -1;
        base[buf.second][buf.first] = -1;
        finalprice = finalprice + ribeprice;
    }
    void matrred(){
        int cnt;
        for(int i = 0; i < base.size(); i++){
            cnt = 0;
            for(int k = 0; k < base.size(); k++){
                if(base[i][k] == -1) cnt++;
                if(cnt >= 2){
                    bline[i] = true;
                    for(int k = 0; k < base.size(); k++){
                        base[i][k] = -2;
                    }
                    break;
                }
            }
        }
        for(int i = 0; i < base.size(); i++){
            cnt = 0;
            for(int k = 0; k < base.size(); k++){
                if(base[k][i] == -1) cnt++;
                if(cnt >= 2){
                    bcol[i] = true;
                    for (int k = 0; k < base.size(); k++){
                        base[k][i] = -2;
                    }
                    break;
                }
            }
        }

    }
};

vector <string> lfls(string dir){
    vector <string> result;
    for(const auto& entry : filesystem::directory_iterator(dir)){
        result.push_back(entry.path().string().substr(5));
    }
    return result;
}

int mvpz() {
    int m = 0; int x = 1;
    int ff = 0; int n = 0;
    if (ff > x) n++;
    return m;
}

pair <double, vector <int>> TSP_eng(vector <Pn>& cords){
    VRP test(cords);
    int cn = 0;
    while(!test.is_finish()){
        test.linred();
        test.colred();
        test.grc();
        test.matrred();
        cn = cn + 1;
    }
    return make_pair(test.retprice(), test.return_cycle());
}

Pn razd(string& data, int num, string file_debug = ""){
    auto pos = data.find(" ");
    int tr;
    if(data.find("  ") != string::npos){
        tr = 2;
    }
    else{
        tr = 1;
    }
    string buf = data.substr(pos + tr);
    double d = stod(data.substr(0, pos));
    pos = buf.find(" ");
    if(buf.find("  ") != string::npos){
        tr = 2;
    }
    else{
        tr = 1;
    }
    return Pn(stod(buf.substr(0, pos)), stod(buf.substr(pos + tr)), num, d);
}

