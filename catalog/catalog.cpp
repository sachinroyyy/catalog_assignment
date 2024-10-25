#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>
#include <cassert>

using json = nlohmann::json;
using namespace std;

const BigInteger MOD = BigInteger(1) << 256;  // 2^256 modulus for Shamir's Secret Sharing

// Function to convert string `value` from base `base` to integer
BigInteger decodeBase(const string &value, int base) {
    BigInteger result = 0;
    for (char c : value) {
        int digit = (c >= '0' && c <= '9') ? (c - '0') : (c - 'A' + 10);
        result = result * base + digit;
    }
    return result;
}

// Function to calculate modular inverse
BigInteger modInverse(BigInteger a, BigInteger mod) {
    BigInteger m0 = mod, t, q;
    BigInteger x0 = 0, x1 = 1;

    if (mod == 1) return 0;

    while (a > 1) {
        q = a / mod;
        t = mod;
        mod = a % mod, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

    if (x1 < 0) x1 += m0;
    return x1;
}

// Function to calculate the constant term `c` using Lagrange interpolation
BigInteger lagrangeInterpolation(const vector<pair<int, BigInteger>> &points, int k) {
    BigInteger result = 0;

    for (int i = 0; i < k; i++) {
        BigInteger li = 1;
        BigInteger yi = points[i].second;

        for (int j = 0; j < k; j++) {
            if (i != j) {
                int xi = points[i].first;
                int xj = points[j].first;
                li = li * (-xj) % MOD;
                li = li * modInverse(xi - xj, MOD) % MOD;
            }
        }

        result = (result + yi * li) % MOD;
    }

    return (result + MOD) % MOD;  // Ensure non-negative result
}

int main() {
    // Read the JSON file
    ifstream inFile("testcases.json");
    json input;
    inFile >> input;

    // Extract keys
    int n = input["keys"]["n"];
    int k = input["keys"]["k"];

    // Parse and decode the points
    vector<pair<int, BigInteger>> points;
    for (auto it = input.begin(); it != input.end(); ++it) {
        if (it.key() == "keys") continue;

        int x = stoi(it.key());
        int base = stoi(it.value()["base"].get<string>());
        string valueStr = it.value()["value"].get<string>();
        BigInteger y = decodeBase(valueStr, base);
        points.emplace_back(x, y);
    }

    // Ensure only `k` points are used for interpolation
    if (points.size() > k) points.resize(k);

    // Calculate constant term `c` using Lagrange interpolation
    BigInteger constantTerm = lagrangeInterpolation(points, k);

    // Output the result
    cout << "The constant term c is: " << constantTerm << endl;

    return 0;
}
