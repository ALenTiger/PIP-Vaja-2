#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <thread>
#include <future>
#include <chrono>
#include <cmath>

using namespace std;

void print(string text) {
    std::cout << endl << text;
}
unsigned int seed = 0;

//vector<bool> Sl;
int Si[] = { -1, +1 };

int Ck(vector<bool>* Sl, int k) {
    int sum = 0;
    for (int i = 0; i < Sl->size() - k; i++) {
        sum += Si[(*Sl)[i]] * Si[(*Sl)[i + k]];
    }
    return sum;
}

int PSL(vector<bool>* Sl) {
    int max = INT_MIN;
    for (int k = 1; k < Sl->size(); k++) {
        int temp = Ck(Sl, k);
        if (temp > max) max = temp;
    }
    return max;
}

double MF(vector<bool>* Sl) {
    int sum = 0;
    int L = Sl->size();
    for (int k = 1; k < L; k++) {
        sum += pow(Ck(Sl, k), 2);
    }
    return pow(L, 2) / (2 * sum);
}

pair<vector<bool>, int> minimizePSL(vector<bool>* Sl, bool print_neighbour_value = false) {
    int min = PSL(Sl);
    vector<bool> min_neighbour = vector<bool>(Sl->size()); //sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>

    for (int i = 0; i < Sl->size(); i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        (*Sl)[i] = !(*Sl)[i];
        int temp = PSL(Sl);
        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp <= min) {
            min = temp;
            min_neighbour = (*Sl);
        }
        (*Sl)[i] = !(*Sl)[i];
    }
    return make_pair(min_neighbour, min);
}

pair<vector<bool>, int> maximizeMF(vector<bool>* Sl, bool print_neighbour_value = false) {
    int max = MF(Sl);
    vector<bool> max_neighbour = vector<bool>(Sl->size());//sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>

    for (int i = 0; i < Sl->size(); i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        (*Sl)[i] = !(*Sl)[i];
        int temp = PSL(Sl);
        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp >= max) {
            max = temp;
            max_neighbour = (*Sl);
        }
        (*Sl)[i] = !(*Sl)[i];
    }
    return make_pair(max_neighbour, max);
}

vector<bool> evaluate(int n, int L, string type) {
    pair<vector<bool>, int> bestSl(vector<bool>(L), INT_MIN);
    pair<vector<bool>, int> cur;
    if (type == "PSL") bestSl.second = INT_MAX;
    vector<bool> Sl = vector<bool>(L);
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < L; j++) {
            Sl[j] = rand() % 2;
        }

        if (type == "PSL") {
            cur = minimizePSL(&Sl);
            if (cur.second < bestSl.second) bestSl = cur;
        }
        else if(type == "MF") {
            cur = maximizeMF(&Sl);
            if (cur.second > bestSl.second) bestSl = cur;
        }
    }
    return bestSl.first;
}

int main() {
    vector<bool> test_a = { 0, 1, 0, 0, 0 };
    vector<bool> test_d = { 1, 0, 0, 0, 1, 0, 0, 1 };
    vector<bool> test_h = { 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0 };

    unsigned int L = 250;
    string type = "PSL";
    unsigned int nfesLmt = 100;

    srand(seed);

    vector<bool> Sl = evaluate(nfesLmt, L, type);
    //Sl = test_h;


    int psl = PSL(&Sl);
    double mf = MF(&Sl);

    print("\n\nPSL: " + to_string(psl) + "   MF: " + to_string(mf));
}
