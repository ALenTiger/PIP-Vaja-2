#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <thread>
#include <future>
#include <chrono>
#include <cmath>
#include <omp.h>

using namespace std;

void print(string text, bool endl_ = true) {
    if(endl_) cout << endl << endl;
    cout << text;
}
unsigned int seed = 0;

int thread_n = 8;
int nfes=0;

//vector<bool> Sl;
int Si[] = { -1, +1 };

chrono::steady_clock::time_point best_begin;
chrono::steady_clock::time_point best_end;
int best_min = INT_MAX;
double best_max = INT_MIN;
bool* bestSL;
vector<int> times;
int long long total_eval_time;

int openmp_thread_n = 0;

string toSequance(bool* binary_b, int L) {
    vector<bool> binary_ = vector<bool>(L);
    for (int i = 0; i < L; i++) {
        binary_[i] = binary_b[i];
    }
    int padding = 4 - binary_.size() % 4;
    vector<bool> binary;
    for (int i = 0; i < padding; i++) binary.push_back(0);
    for (bool b : binary_) binary.push_back(b);
    int decimal = 0;
    string output = "";
    for (int i = 0; i < binary.size(); i++) {
        if (!(i % 4)) decimal = 0;
        int bit = binary[binary.size() - i - 1];
        decimal += bit * pow(2, ((i) % 4));
        if (!((i + 1) % 4)) {
            if (decimal>9) {
                output += (char)(decimal - 10 + 'A');
            }
            else {
                output = to_string(decimal) + output;
            }
        }
    }
    return "0x" + output;
}

bool** splitArr(bool* arr, int L, int thread_n_ = thread_n) {
    bool** splited = new bool*[thread_n_];
    for (int i = 0; i < thread_n_; i++) {
        splited[i] = new bool[L];
        for (int j = 0; j < L; j++) {
            splited[i][j] = arr[j];
        }
    }
    return splited;
}

bool* vectorToArr(vector<bool> vec) {
    bool* arr = new bool[vec.size()];
    for (int i = 0; i < vec.size(); i++) {
        arr[i] = vec[i];
    }
    return arr;
}

int splitedCk(bool* Sl, int k, int start_i, int end_i) {
    int sum = 0;
    for (int i = start_i; i < end_i - k; i++) {
        sum += (Sl[i] * 2 - 1) * (Sl[i + k] * 2 - 1);
    }
    return sum;
}

int Ck(bool* Sl, int k, int start_i, int end_i) {
    int sum = 0;
    for (int i = start_i; i < end_i - k; i++) {
        sum += (Sl[i] * 2 - 1)* (Sl[i + k] * 2 - 1);
    }
    return sum;
}

int PSL(bool* Sl, int L) {
    int max = INT_MIN;
    for (int k = 1; k < L; k++) {
        int temp = Ck(Sl, k, 0, L);
        if (temp > max) max = temp;
    }
    return max;
}

int splitedPSL(bool* Sl, int L, int start_i) {
    int max = INT_MIN;
    //if (start_i == 0) start_i = 1;
    for (int i = 1 + start_i; i <= L / 2; i += thread_n) {
        int temp = Ck(Sl, L - i, 0, L);
        if (temp > max) max = temp;
        if (L - i != i) {
            temp = Ck(Sl, i, 0, L);
            if (temp > max) max = temp;
        }
    }
    return max;
}

int parallelPSL(bool* Sl, int L) {
    int max = INT_MIN;
    vector<future<int>> futures;

    for (int i = 0; i < thread_n; i++) {
        futures.push_back(async(launch::async, splitedPSL, Sl, L, i));
    }
    for (future<int>& f : futures) {
        int temp = f.get();
        if (temp > max) max = temp;
    }
    return max;
}




double MF(bool* Sl, int L) {
    int sum = 0;
    for (int k = 1; k < L; k++) {
        sum += pow(Ck(Sl, k, 0, L), 2);
    }
    return pow(L, 2) / (2 * sum);
}

int splitedMF(bool* Sl, int L, int start_i, int end_i) {
    int sum = 0;
    if (start_i == 0) start_i = 1;
    for (int k = start_i; k < end_i; k++) {
        sum += pow(Ck(Sl, k, 0, L), 2);
    }
    return sum;
}

double parallelMF(bool* Sl, int L) {
    int sum = 0;
    vector<future<int>> futures;

    for (int k = 0; k < thread_n; k++) {
        int strt_i = k * (L / thread_n);
        int end_i = (k + 1) * (L / thread_n);
        if (k == thread_n - 1 && end_i != L) end_i = L;
        futures.push_back(async(launch::async, splitedMF, Sl, L, strt_i, end_i));
    }
    for (future<int>& f : futures) {
        int temp = f.get();
        sum += temp;
    }
    return pow(L, 2) / (2 * sum);
}


void minimizePSL(bool* Sl, int L, bool print_neighbour_value = false, bool print_eval_time = false) {
    int min = PSL(Sl, L);
    int best_i = 0;
    //vector<bool> min_neighbour = vector<bool>(Sl->size()); //sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>
 
    #pragma omp parallel for num_threads (thread_n)
    for (int i = 0; i < L; i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        int used_threads = omp_get_num_threads();
        if (used_threads > openmp_thread_n) openmp_thread_n = used_threads;
        Sl[i] = !Sl[i];

        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        int temp = PSL(Sl, L);
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto temp_time = chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        total_eval_time += temp_time;
        if (print_eval_time) print("1 eval time: " + to_string(temp_time));
        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp <= min) {
            min = temp;
            best_i = i;
        }
        Sl[i] = !Sl[i];
        nfes++;
    }
    Sl[best_i] = !Sl[best_i];
    if (min < best_min) {
        best_end = std::chrono::steady_clock::now();
        best_min = min;
        bestSL = Sl;
    }
}

void maximizeMF(bool* Sl, int L, bool print_neighbour_value = false) {
    int max = parallelMF(Sl, L);
    int best_i = 0;

    //vector<bool> max_neighbour = vector<bool>(Sl->size());//sprobaj če hitreje deluje če delamu tu z min_neigbour_index namesto da se kreira oz. kopira vector<bool>

    for (int i = 0; i < L; i++) { //ovrednoti PSL za vse sosede Sl in najde najmanjšega
        Sl[i] = !Sl[i];

        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        int temp = parallelMF(Sl, L);
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto temp_time = end - begin;
        total_eval_time += (int)temp_time.count();

        if (print_neighbour_value) print("neighbour at index (" + to_string(i) + ") with PSL value: " + to_string(temp));
        if (temp >= max) {
            max = temp;
            best_i = i;

        }
        Sl[i] = !Sl[i];
        nfes++;
    }
    Sl[best_i] = !Sl[best_i];
    if (max > best_max) {
        best_end = std::chrono::steady_clock::now();
        best_max = max;
        bestSL = Sl;
    }
}

bool* evaluate(int n, int L, string type, bool print_progress = true) {
    //pair<vector<bool>, int> bestSl(vector<bool>(L), INT_MIN);
    //vector<bool> cur;
    //if (type == "PSL") bestSl.second = INT_MAX;

    bool* Sl = new bool[L];

    for (int j = 0; j < L; j++) {
        Sl[j] = rand() % 2;
    }

    for (int i = 0; i < n; i++) {
        //float f = i % (n / 100);
        if (print_progress && !(i%50)) print(to_string((((float)i) / ((float)n)) * 100) + "%..    ", false);
        if (type == "PSL") {
            minimizePSL(Sl, L);
            //if (cur.second < bestSl.second) bestSl = cur;
        }
        else if(type == "MF") {
            maximizeMF(Sl, L);
            //if (cur.second > bestSl.second) bestSl = cur;
        }
    }
    return bestSL;
}

int main() {
    vector<bool> test_a = { 0, 1, 0, 0, 0 };
    vector<bool> test_d = { 1, 0, 0, 0, 1, 0, 0, 1 };
    vector<bool> test_h = { 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0 };
    unsigned int L = 200;
    string type = "PSL";
    unsigned int nfesLmt = 5000;

    srand(seed);

    bool *Sl;

    int n = thread_n;
    thread_n = 0;
    for (int i = 1; i <= n; i++) {
        best_min = INT_MAX;
        best_max = INT_MIN;
        nfes = 0; 
        total_eval_time = 0;
        print("\n\novrednotenje na: " + to_string(i) + "-ih nitih: \n\n");
        thread_n = i;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point best_begin = std::chrono::steady_clock::now();
        Sl = evaluate(nfesLmt, L, type);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto time_taken = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        auto best_time_taken = std::chrono::duration_cast<std::chrono::microseconds>(best_end - best_begin).count();
        int psl = parallelPSL(Sl, L);
        double mf = parallelMF(Sl, L);
        print("\n\nPSL: " + to_string(psl) + "   MF: " + to_string(mf) + "   nfes: " + to_string(nfes));
        auto test = total_eval_time;
        print("povprečen čas 1 ovrednotenja: " + to_string((time_taken / nfes)) + " mus");
        print("povprečen točen čas 1 ovrednotenja: " + to_string((total_eval_time / nfes)) + " mus");
        print("celoten cas izvajanja (" + to_string(nfes) + ")-ih ovrednotenj: " + to_string(time_taken) + " mus");
        print("ali: " + to_string(time_taken / (float)1000000) + " sec");
        print("Čas potreben, da smo najsli najbolslo resitev: " + to_string(best_time_taken) + " mus");
        print("ali: " + to_string(best_time_taken / (float)1000000) + " sec");
        print("peed (st ovrednotenj na sekundo)" + to_string(nfes / (total_eval_time /1000000)));
        print("sekvenca: " + toSequance(Sl, L));
        print("number of used threads with open mp: " + to_string(openmp_thread_n));
    }

    print("\n----\n\nizhod:");
    print("L:   " + to_string(L)); 
    print("nfesLmt:   " + to_string(nfesLmt));
    print("seed:   " + to_string(seed));
    /*sequence: 0x0C23529719DE19
    MF : <double>
    PSL : <unsigned int>
    //Sl = test_h;*/

    thread_n = 4;
    int psl = parallelPSL(vectorToArr(test_d), L);
    double mf = parallelMF(vectorToArr(test_d), L);
    print("\n\nPSL: " + to_string(psl) + "   MF: " + to_string(mf) + "   nfes: " + to_string(nfes));
    print("\n");
}