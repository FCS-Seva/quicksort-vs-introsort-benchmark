#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>
#include <cstdint>

std::mt19937 rng(3048927);
constexpr int INSERTION_BOUND = 16;

class ArrayGenerator {
public:
    ArrayGenerator(const size_t max_n, const int min_val, const int max_val, const uint64_t seed = 9238417)
        : max_n_(max_n), min_val_(min_val), max_val_(max_val), rng_(seed)
    {
        generate_base_arrays();
    }

    std::vector<int> get_random(const size_t n) const {
        return std::vector<int>(random_base_.begin(), random_base_.begin() + n);
    }

    std::vector<int> get_reversed(const size_t n) const {
        return std::vector<int>(reversed_base_.begin(), reversed_base_.begin() + n);
    }

    std::vector<int> get_almost_sorted(const size_t n) const {
        return std::vector<int>(almost_sorted_base_.begin(), almost_sorted_base_.begin() + n);
    }

private:
    size_t max_n_;
    int min_val_, max_val_;
    mutable std::mt19937_64 rng_;

    std::vector<int> random_base_;
    std::vector<int> reversed_base_;
    std::vector<int> almost_sorted_base_;

    void generate_base_arrays() {
        std::uniform_int_distribution<int> dist(min_val_, max_val_);

        random_base_.resize(max_n_);
        for (size_t i = 0; i < max_n_; ++i) {
            random_base_[i] = dist(rng_);
        }

        reversed_base_ = random_base_;
        std::sort(reversed_base_.begin(), reversed_base_.end());
        std::reverse(reversed_base_.begin(), reversed_base_.end());

        almost_sorted_base_ = random_base_;
        std::sort(almost_sorted_base_.begin(), almost_sorted_base_.end());

        size_t num_swaps = max_n_ / 100;
        std::uniform_int_distribution<size_t> pos_dist(0, max_n_ - 1);
        for (size_t k = 0; k < num_swaps; ++k) {
            size_t i = pos_dist(rng_);
            size_t j = pos_dist(rng_);
            std::swap(almost_sorted_base_[i], almost_sorted_base_[j]);
        }
    }
};


void InsertionSort(std::vector<int>& v, const int left, const int right) {
    for (int i = left + 1; i <= right; ++i) {
        int x = v[i];
        int j = i - 1;
        while (j >= left && v[j] > x) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = x;
    }
}

void Heapify(std::vector<int>& v, int left, int heap_size, int i) {
    while (true) {
        int largest = i;
        int l = 2 * i + 1;
        int r = 2 * i + 2;

        if (l < heap_size && v[left + l] > v[left + largest]) {
            largest = l;
        }
        if (r < heap_size && v[left + r] > v[left + largest]) {
            largest = r;
        }

        if (largest == i) {
            break;
        }

        std::swap(v[left + i], v[left + largest]);
        i = largest;
    }
}

void HeapSort(std::vector<int>& v, int left, int right) {
    const int n = right - left + 1;
    if (n <= 1) {
        return;
    }

    for (int i = n / 2 - 1; i >= 0; --i) {
        Heapify(v, left, n, i);
    }

    for (int i = n - 1; i > 0; --i) {
        std::swap(v[left], v[left + i]);
        Heapify(v, left, i, 0);
    }
}

int Partition(std::vector<int> &v, const int left, const int right) {
    std::uniform_int_distribution<int> distrib(left, right);
    const int pivot_index = distrib(rng);
    const int pivot = v[pivot_index];
    std::swap(v[pivot_index], v[right]);
    int i = left;
    for (int j = left; j < right; j++) {
        if (v[j] < pivot) {
            std::swap(v[i], v[j]);
            i++;
        }
    }
    std::swap(v[i], v[right]);
    return i;
}

void IntroSort(std::vector<int> &v, const int left, const int right, int depth) {
    if (left >= right) {
        return;
    }
    const int n = right - left + 1;

    if (n < INSERTION_BOUND) {
        InsertionSort(v, left, right);
        return;
    }

    if (depth == 0) {
        HeapSort(v, left, right);
        return;
    }

    const int pivot = Partition(v, left, right);

    IntroSort(v, left, pivot - 1, depth - 1);
    IntroSort(v, pivot + 1, right, depth - 1);
}

void QuickSortImpl(std::vector<int>& v, int left, int right) {
    if (left >= right) {
        return;
    }
    int pivot = Partition(v, left, right);
    QuickSortImpl(v, left, pivot - 1);
    QuickSortImpl(v, pivot + 1, right);
}

void QuickSort(std::vector<int>& v) {
    if (v.empty()) {
        return;
    }
    QuickSortImpl(v, 0, static_cast<int>(v.size()) - 1);
}

class SortTester {
public:
    SortTester(ArrayGenerator& gen, int repeats)
        : gen_(gen), repeats_(repeats) {}

    void run_all() {

        std::cout << "type;n;algo;threshold;time_ms\n";

        for (const std::string& type : {"random", "reversed", "almost"}) {
            for (int n = 500; n <= 100000; n += 100) {
                long long avg_ms = measure_one(type, n, "quick", 0);
                std::cout << type << ";" << n << ";quick;0;" << avg_ms << "\n";
            }
            for (int n = 500; n <= 100000; n += 100) {
                long long avg_ms = measure_one(type, n, "intro", INSERTION_BOUND);
                std::cout << type << ";" << n << ";intro;" << INSERTION_BOUND << ";" << avg_ms << "\n";
            }
        }
    }

private:
    ArrayGenerator& gen_;
    int repeats_;

    std::vector<int> make_array(const std::string& type, const int n) const {
        if (type == "random") {
            return gen_.get_random(static_cast<size_t>(n));
        }

        if (type == "reversed") {
            return gen_.get_reversed(static_cast<size_t>(n));
        }

        return gen_.get_almost_sorted(static_cast<size_t>(n));
    }

    long long measure_one(const std::string& type, const int n, const std::string& algo, int threshold) const {
        (void)threshold;
        return measure_one(type, n, algo);
    }

    long long measure_one(const std::string& type, const int n, const std::string& algo) const {
        using namespace std::chrono;
        long long total_ms = 0;

        for (int rep = 0; rep < repeats_; ++rep) {
            std::vector<int> a = make_array(type, n);

            int depth = 2 * static_cast<int>(std::log2(static_cast<double>(n)));
            auto start = high_resolution_clock::now();
            if (algo == "quick") {
                QuickSort(a);
            } else {
                IntroSort(a, 0, n - 1, depth);
            }
            auto elapsed = high_resolution_clock::now() - start;
            long long ms = duration_cast<milliseconds>(elapsed).count();
            total_ms += ms;
        }
        return total_ms / std::max(1, repeats_);
    }
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    constexpr size_t MAX_N = 100000;
    constexpr int MIN_VAL = 0;
    constexpr int MAX_VAL = 6000;

    ArrayGenerator gen(MAX_N, MIN_VAL, MAX_VAL);
    SortTester tester(gen, 5);

    tester.run_all();

    return 0;
}
