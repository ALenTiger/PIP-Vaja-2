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

int thread_n = 4;
int nfes=0;

//vector<bool> Sl;
int Si[] = { -1, +1 };

chrono::steady_clock::time_point best_begin;
chrono::steady_clock::time_point best_end;
int best_min = INT_MAX;
double best_max = INT_MIN;


int Ck(vector<bool>* Sl, int k, int start_i, int end_i) {
    int sum = 0;
    for (int i = start_i; i < end_i - k; i++) {
        sum += Si[(*Sl)[i]] * Si[(*Sl)[i + k]];
    }
    return sum;
}

int PSL(vector<bool>* Sl) {
    int max = INT_MIN;
    int L = Sl->size();
    for (int k = 1; k < L; k++) {
        int temp = Ck(Sl, k, 0, L);
        if (temp > max) max = temp;
    }
    return max;
}

int parallelPSL(vector<bool>* Sl) {
    int max = INT_MIN;
    int L = Sl->size();
    vector<future<int>> futures;

    for (int k = 1; k < thread_n; k++) {
        int strt_i = k * (L / thread_n);
        int end_i = (k + 1) * (L / thread_n);
        futures.push_back(async(launch::async, Ck, Sl, k, strt_i, end_i));
    }
    for (future<int>& f : futures) {
        int temp = f.get();
        if (temp > max) max = temp;
    }
    return max;
}

double MF(vector<bool>* Sl) {
    int sum = 0;
    int L = Sl->size();
    for (int k = 1; k < L; k++) {
        sum += pow(Ck(Sl, k, 0, L), 2);
    }
    return pow(L, 2) / (2 * sum);
}

double parallelMF(vector<bool>* Sl) {
    int sum = 0;
    int L = Sl->size();
    vector<future<int>> futures;


    for (int k = 1; k < thread_n; k++) {
        int strt_i = k * (L / thread_n);
        int end_i = (k + 1) * (L / thread_n);
        futures.push_back(async(launch::async, Ck, Sl, k, strt_i, end_i));
    }
    for (future<int>& f : futures) {
        int temp = f.get();
        sum += pow(temp, 2);
    }
    return pow(L, 2) / (2 * sum);
}

pair<vector<bool>, int> minimizePSL(vector<bool>* Sl, bool print_neighbour_value = false) {
    int min = parallelPSL(Sl);
    vector<bool> min_neighbour = vector<bool>(Sl->size()); //sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>

    for (int i = 0; i < Sl->size(); i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        (*Sl)[i] = !(*Sl)[i];
        int temp = parallelPSL(Sl);
        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp <= min) {
            min = temp;
            min_neighbour = (*Sl);
        }
        (*Sl)[i] = !(*Sl)[i];
        nfes++;
    }
    return make_pair(min_neighbour, min);
}

pair<vector<bool>, int> maximizeMF(vector<bool>* Sl, bool print_neighbour_value = false) {
    int max = parallelMF(Sl);
    vector<bool> max_neighbour = vector<bool>(Sl->size());//sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>

    for (int i = 0; i < Sl->size(); i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        (*Sl)[i] = !(*Sl)[i];
        int temp = parallelMF(Sl);
        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp >= max) {
            max = temp;
            max_neighbour = (*Sl);
        }
        (*Sl)[i] = !(*Sl)[i];
        nfes++;
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
    unsigned int nfesLmt = 20;

    srand(seed);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vector<bool> Sl = evaluate(nfesLmt, L, type);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    auto time_taken = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    //Sl = test_h;


    int psl = parallelPSL(&Sl);
    double mf = parallelMF(&Sl);

    print("\n\nPSL: " + to_string(psl) + "   MF: " + to_string(mf) + "   nfes: " + to_string(nfes));
    print("celoten cas izvajanja (" + to_string(nfes) + ")-ih ovrednotenj: " + to_string(time_taken ) + " sec");
    print("povprečen čas 1 ovrednotenja: " + to_string((time_taken / nfes)) + " sec");
}
