#include <string>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "needlemanWunsch.cpp"
#include "initData.cpp"
using namespace std;

void swapDouble(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void swapString(string* a, string* b)
{
    string t = *a;
    *a = *b;
    *b = t;
}

int partition (vector<string>& names, vector<double>& data, int low, int high)
{
    double pivot = data[high];
    int i = (low - 1);
 
    for (int j = low; j <= high- 1; j++)
    {
        if (data[j] <= pivot)
        {
            i++;   
            swapDouble(&data[i], &data[j]);
            swapString(&names[i], &names[j]);
        }
    }
    swapDouble(&data[i + 1], &data[high]);
    swapString(&names[i + 1], &names[high]);
    return (i + 1);
}

void quickSort(vector<string>& names, vector<double>& data, int low, int high)
{
    if (low < high)
    {
        int pi = partition(names, data, low, high);
        quickSort(names, data, low, pi - 1);
        quickSort(names, data, pi + 1, high);
    }
}

static int getMinSequenceLength(vector<string> names, vector<string> seqs)
{
    int currLen;
    int minLen = seqs[0].length();
    int minIndex = 0;
    for (int i = 1; i < seqs.size(); i++)
    {
        currLen = seqs[i].length();
        if (currLen < minLen)
        {
            minLen = currLen;
            minIndex = i;
        }
    }
    cout << "The shortest sequence is " << names[minIndex] << " - " << seqs[minIndex] << " and has a length of " << minLen << endl;
    return minLen;
}

static int checkFastaSameLength(vector<string> seqs)
{
    int length = seqs[0].length();
    for (string const &seq : seqs)
    {
        if (seq.length() != length)
            return -1;
    }
    return length;
}

static void removeDisallowed(string& seq, const string& allowed) 
{
    unordered_set<char> allowedSet(allowed.begin(), allowed.end());
    seq.erase
    (
        remove_if(seq.begin(), seq.end(), [&](const char c) 
        {
            return !allowedSet.count(c);
        }
    ), seq.end());
}

static vector<double> AAC(const string& seq, const string& allowed)
{
    vector<double> encoded;
    unordered_map<char, int> count;
    for (char const &c : seq)  
    {
        count[c]++;
    }
    double l = seq.length() * 1.0;
    for (char const &c : allowed)
    {
        encoded.push_back(count[c] / l);
    }
    return encoded;
}

static vector<double> EAAC(const string& seq, const string& allowed, vector<string> keys, int window, int windows)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    for (int i = 0; i < windows; i++)
    {
        for (int j = 0; j < window; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key.push_back(seq[i + j]);
            count[key]++;
        }
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / window);
    }
    return encoded;
}

static vector<double> CKSAAP(const string& seq, const string& allowed, vector<string> keys, int gap)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int seqLen = seq.length();
    for (int i = 0; i < gap + 1; i++)
    {
        for (int j = 0; j < seqLen - i; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key.push_back(seq[j]);
            key.push_back(seq[i + j + 1]);
            count[key]++;
        }
    }
    int gapCount = 0;
    int pairCount = 0;
    int allowedLen = allowed.length();
    int maxPairCount = allowedLen * allowedLen;
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (seqLen - gapCount - 1));
        pairCount++;
        if(pairCount == maxPairCount)
        {
            gapCount++;
            pairCount = 0;
        }
    }
    return encoded;
}

static vector<double> TPC(const string& seq, const string& allowed, vector<string> keys)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int i = 0; i < l - 2; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        key.push_back(seq[i + 2]);
        count[key]++;
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - 2));
    }
    return encoded;
}

static vector<double> DPC(const string& seq, const string& allowed, vector<string> keys)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key]++;
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - 1));
    }
    return encoded;
}

static vector<double> DDE(const string& seq, const string& allowed, vector<string> keys, map<string, double> tm, map<string, double> tvP)
{
    vector<double> encoded;
    map<string, double> count;
    map<string, double> tv(tvP);
    int l = seq.length();
    for (string const &key : keys)
    {
        count[key] = 0;
        tv[key] /= (l - 1);
    }

    for (int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key]++;
    }

    for (char const &c : allowed)
    {
        for (char const &d : allowed)
        {
            string key(1, c);
            key.push_back(d);
            count[key] /= (l - 1);
            count[key] = (count[key] - tm[key]) / sqrt(tv[key]);
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> GAAC(const string& seq, vector<string> keys, map<char, string> groups)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }

    int l = seq.length();
    for (char const &c : seq) 
    {
        count[groups[c]]++;
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / l);
    }
    return encoded;
}

static vector<double> EGAAC(const string& seq, vector<string> keys, map<char, string> groups, int window, int windows)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }

    int l = seq.length();
    for (int i = 0; i < windows; i++)
    {
        for (int j = 0; j < window; j++) 
        {
            string key = to_string(i);
            key.push_back('-');
            key += groups[seq[i + j]];
            count[key]++;
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / window);
    }
    return encoded;
}

static vector<double> CKSAAGP(const string& seq, vector<string> keys, int gap, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }

    int seqLen = seq.length();
    for (int i = 0; i < gap + 1; i++)
    {
        for (int j = 0; j < seqLen - i; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key += groups[seq[j]];
            key.push_back('-');
            key += groups[seq[i + j + 1]];
            count[key]++;
        }
    }
    int gapCount = 0;
    int pairCount = 0;
    int maxPairCount = 25;
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (seqLen - gapCount - 1));
        pairCount++;
        if(pairCount == maxPairCount)
        {
            gapCount++;
            pairCount = 0;
        }
    }
    return encoded;
}

static vector<double> GDPC(const string& seq, vector<string> keys, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }

    int l = seq.length();
    for (int i = 0; i < l - 1; i++) 
    {
        string key = groups[seq[i]];
        key.push_back('-');
        key += groups[seq[i + 1]];
        count[key]++;
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - 1));
    }
    return encoded;
}

static vector<double> GTPC(const string& seq, vector<string> keys, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }

    int l = seq.length();
    for (int i = 0; i < l - 2; i++) 
    {
        string key = groups[seq[i]];
        key.push_back('-');
        key += groups[seq[i + 1]];
        key.push_back('-');
        key += groups[seq[i + 2]];
        count[key]++;
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - 2));
    }
    return encoded;
}

static vector<double> binary(const string& seq, vector<string> keys)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    
    int counter = 0;
    for (char const &c : seq)
    {
        string key(1, c);
        key = to_string(counter) + key;
        count[key] = 1;
        counter++;
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> Moran(const string& seq, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<double> encoded;
    map<string, double> average;
    int seqLen = seq.length();
    int count = 0;
    for (string const &index : indexList)
    {
        average[index] = 0;
        for (char const &c : seq)
        {
            average[index] += (indices[index][order[c]] / seqLen);
        }
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            double den = 0;
            for (int j = 0; j < seqLen; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                double avg = average[index];
                if (j < seqLen - i)
                    num += (indices[index][order[c1]] - average[index]) * (indices[index][order[c2]] - average[index]);
                den += pow((indices[index][order[c1]] - average[index]), 2);
            }
            num /= (seqLen - i);
            den /= seqLen;
            encoded.push_back(num / den);
        }
    }
    return encoded;
}

static vector<double> Geary(const string& seq, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<double> encoded;
    map<string, double> average;
    int seqLen = seq.length();
    for (string const &index : indexList)
    {
        average[index] = 0;
        for (char const &c : seq)
        {
            average[index] += (indices[index][order[c]] / seqLen);
        }
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            double den = 0;
            for (int j = 0; j < seqLen; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                double avg = average[index];
                if (j < seqLen - i)
                    num += pow((indices[index][order[c1]] - indices[index][order[c2]]), 2);
                den += pow((indices[index][order[c1]] - average[index]), 2);
            }
            num /= (2*(seqLen - i));
            den /= (seqLen - 1);
            encoded.push_back(num / den);
        }
    }
    return encoded;
}

static vector<double> NMB(const string& seq, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<double> encoded;
    int seqLen = seq.length();
    for (string const &index : indexList)
    {
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            for (int j = 0; j < seqLen - i; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                num += (indices[index][order[c1]] * indices[index][order[c2]]);
            }
            encoded.push_back(num / (seqLen - i));
        }
    }
    return encoded;
}

static vector<double> CTDC(const string& seq, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    for (string const &property : properties)
    {
        for (int i = 0; i < groups.size(); i++)
        {
            for (char const &c : seq)
            {
                for (char const &d : groups[i][property])
                {
                    if (c == d)
                    {
                        string key = property;
                        key.push_back('-');
                        key += to_string(i + 1);
                        count[key]++;
                        goto endLoop;
                    }
                }
                endLoop:;
            }
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / seq.length());
    }
    return encoded;
}

static vector<double> CTDT(const string& seq, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int seqLen = seq.length();
    for (string const &property : properties)
    {
        for (int i = 0; i < groups.size(); i++)
        {
            for (int k = i + 1; k < groups.size(); k++)
            {
                for (int l = 0; l < seqLen - 1; l++)
                {
                    for (char const &c : groups[i][property])
                    {
                        for (char const &d : groups[k][property])
                        {
                            if ((seq[l] == c && seq[l + 1] == d) || (seq[l] == d && seq[l + 1] == c))
                            {
                                string key = property;
                                key.push_back('-');
                                key += to_string(i + 1);
                                key.push_back('-');
                                key += to_string(k + 1);
                                count[key]++;
                                goto endLoop;
                            }
                        }
                    }
                    endLoop:;
                }
            }
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (seqLen - 1));
    }
    return encoded;
}

static vector<double> CTDD(const string& seq, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties, array<int, 5> pcts)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int seqLen = seq.length();
    for (string const &property : properties)
    {
        for (int j = 0; j < groups.size(); j++)
        {
            int counter = 0;
            for (char const &c : seq)
                for (char const &d : groups[j][property])
                    if (c == d)
                        counter++;
            double p25 = floor(0.25 * counter);
            double p50 = floor(0.50 * counter);
            double p75 = floor(0.75 * counter);
            double p100 = counter;
            map<int, double> cutoffs = {{0, 1}, {25, p25 > 1 ? p25 : 1}, {50, p50 > 1 ? p50 : 1}, {75, p75 > 1 ? p75 : 1}, {100, p100 > 1 ? p100 : 1}};
            for (auto const &kvp : cutoffs)
            {
                string key = property;
                key.push_back('-');
                key += to_string(j + 1);
                key.push_back('-');
                key += to_string(kvp.first);
                counter = 0;
                for (int k = 0; k < seqLen; k++)
                {
                    for (char const &d : groups[j][property])
                    {
                        if (seq[k] == d)
                        {
                            counter++;
                            if (counter == kvp.second)
                            {
                                count[key] = (k + 1) * 100.0 / seqLen;
                                goto endLoop;
                            }
                        }
                    }
                }
                endLoop:;
                if (counter == 0)
                {
                    count[key] = 0;
                }
            }
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> CT(const string& seq, const string& allowed, vector<string> keys, array<string, 7> groups)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int i = 0; i < l - 2; i++) 
    {
        int g1 = 0;
        int g2 = 0;
        int g3 = 0;
        for (int j = 0; j < 7; j++)
        {
            for (char const &c : groups[j])
            {
                if (seq[i] == c)
                    g1 = j + 1;
                if (seq[i + 1] == c)
                    g2 = j + 1;
                if (seq[i + 2] == c)
                    g3 = j + 1;
            }
        }
        string key = to_string(g1);
        key.push_back('-');
        key += to_string(g2);
        key.push_back('-');
        key += to_string(g3);
        count[key]++;
    }
    double min = l;
    double max = 0;
    for (auto const &kvp : count)
    {
        min = kvp.second < min ? kvp.second : min;
        max = kvp.second > max ? kvp.second : max;
    }
    for (string const &key : keys) 
    {
        encoded.push_back((count[key] - min) / max);
    }
    return encoded;
}

static vector<double> KSCT(const string& seq, const string& allowed, vector<string> keys, array<string, 7> groups, int ks)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int currK = 0; currK <= ks; currK++)
    {
        for (int i = 0; i < l - 2; i++) 
        {
            int g1 = 0;
            int g2 = 0;
            int g3 = 0;
            for (int j = 0; j < 7; j++)
            {
                for (char const &c : groups[j])
                {
                    if (seq[i] == c)
                        g1 = j + 1;
                    if (seq[i + (1 * (currK + 1))] == c)
                        g2 = j + 1;
                    if (seq[i + (2 * (currK + 1))] == c)
                        g3 = j + 1;
                }
            }
            string key = to_string(g1);
            key.push_back('-');
            key += to_string(g2);
            key.push_back('-');
            key += to_string(g3);
            key.push_back('-');
            key += to_string(currK);
            count[key]++;
        }
        double min = l;
        double max = 0;
        for (int i = currK * 343; i < (currK + 1) * 343; i++)
        {
            min = count[keys[i]] < min ? count[keys[i]] : min;
            max = count[keys[i]] > max ? count[keys[i]] : max;
        }
        for (int i = currK * 343; i < (currK + 1) * 343; i++)
        {
            count[keys[i]] = (count[keys[i]] - min) / max;
        }
    } 
    
    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> SOCNumber(const string& seq, const string& allowed, vector<string> keys, map<char, int> order, 
    map<char, array<double, 20>> swData, map<char, array<double, 20>> gData, int lag)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int i = 1; i <= lag; i++)
    {
        double sw = 0;
        double g = 0;
        for (int j = 0; j < l - i; j++)
        {
            sw += pow(swData[seq[j]][order[seq[j + i]]], 2);
            g += pow(gData[seq[j]][order[seq[j + i]]], 2);
        }
        sw /= (l - i);
        g /= (l - i);

        string key(1, '-');
        key += to_string(i);
        count["sw" + key] = sw;
        count["g" + key] = g;
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> QSOrder(const string& seq, const string& allowed, vector<string> keys, map<char, int> order, 
    map<char, array<double, 20>> swData, map<char, array<double, 20>> gData, int lag, double weight)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    for (int i = 1; i <= lag; i++)
    {
        double sw = 0;
        double g = 0;
        for (int j = 0; j < l - i; j++)
        {
            sw += pow(swData[seq[j]][order[seq[j + i]]], 2);
            g += pow(gData[seq[j]][order[seq[j + i]]], 2);
        }

        string key(1, '-');
        key += to_string(i);
        count["sw" + key] = sw;
        count["g" + key] = g;
    }

    double sumSW = 0;
    double sumG = 0;
    for (int i = 1; i <= lag; i++)
    {
        sumSW += count["sw-" + to_string(i)];
        sumG += count["g-" + to_string(i)];
    }

    map<char, int> counts;
    for (char const &c : allowed)
    {
        counts[c] = 0;
    }
    for (char const &c : seq)
    {
        counts[c]++;
    }
    for (char const &c : allowed)
    {
        string key(1, '-');
        key.push_back(c);
        count["sw" + key] = counts[c] / (1 + weight * sumSW);
        count["g" + key] = counts[c] / (1 + weight * sumG);
    }
    for (int i = 1; i <= lag; i++)
    {
        string key(1, '-');
        key += to_string(i);
        count["sw" + key] *= weight / (1 + weight * sumSW);
        count["g" + key] *= weight / (1 + weight * sumG);
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> PAAC(const string& seq, const string& allowed, vector<string> keys, int lag, map<char, int> order,
    map<string, array<double, 20>> paac, double weight)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    int l = seq.length();
    double globalSum = 0;
    for (int i = 1; i <= lag; i++)
    {
        double res = 0;
        for (int j = 0; j < l - i; j++)
        {
            for (auto const &kvp : paac)
            {
                res += pow(kvp.second[order[seq[j]]] - kvp.second[order[seq[j + i]]], 2) / 3; // 3 for H1, H2 and SCM
            }
        }
        res /= (l - i);
        globalSum += res;
        string key = to_string(i);
        count[key] = res;
    }
    map<char, int> counts;
    for (char const &c : allowed)
    {
        counts[c] = 0;
    }
    for (char const &c : seq)
    {
        counts[c]++;
    }

    for (char const &c : allowed)
    {
        string key(1, c);
        count[key] = counts[c] / (1 + weight * globalSum);
    }
    for (int i = 1; i <= lag; i++)
    {
        string key = to_string(i);
        count[key] = (weight * count[key]) / (1 + weight * globalSum);
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> APAAC(const string& seq, const string& allowed, vector<string> keys, int lag, map<char, int> order,
    map<string, array<double, 20>> paac, string properties[], double weight)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    
    int l = seq.length();
    double globalSum = 0;
    for (int i = 1; i <= lag; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            string property = properties[j];
            double res = 0;
            for (int k = 0; k < l - i; k++)
                res += paac[property][order[seq[k]]] * paac[property][order[seq[k + i]]];
            res /= (l - i);
            globalSum += res;
            string key = property;
            key.push_back('-');
            key += to_string(i);
            count[key] = res;
        }
    }

    map<char, int> counts;  
    for (char const &c : allowed)
    {
        counts[c] = 0;
    }
    for (char const &c : seq)
    {
        counts[c]++;
    }

    for (char const &c : allowed)
    {
        string key(1, c);
        count[key] = counts[c] / (1 + weight * globalSum);
    }
    for (int i = 1; i <= lag; i++)
    {
        for (auto const &kvp : paac)
        {
            string key = kvp.first;
            key.push_back('-');
            key += to_string(i);
            count[key] = (weight * count[key]) / (1 + weight * globalSum);
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> AAINDEX(const string& seq, vector<string> keys, map<char, int> order, map<string, 
    array<double, 20>> indexData, vector<string> indexList)
{
    vector<double> encoded;
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    
    for (char const &c : seq)
    {
        for (string const &index : indexList)
        {
            encoded.push_back(indexData[index][order[c]]);
        }
    }
    return encoded;
}

static vector<double> BLOSUM62(const string& seq, vector<string> keys, map<char, array<double, 20>> blosum)
{
    vector<double> encoded;
    int l = seq.length();
    for (char const &c : seq)
    {
        for (double const &val : blosum[c])
        {
            encoded.push_back(val);
        }
    }
    return encoded;
}

static vector<double> ZSCALE(const string& seq, vector<string> keys, map<char, array<double, 5>> zscale)
{
    vector<double> encoded;
    for (char const &c : seq)
    {
        for (double const &val : zscale[c])
        {
            encoded.push_back(val);
        }
    }
    return encoded;
}

static vector<double> SSEB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file, map<char, array<int, 3>> types, const string& type)
{
    vector<double> encoded;

    string aas;
    string elements;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2

    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                elements.push_back(line[0]);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            char c = elements[i];
            if (types.find(c) != types.end())
            {
                array<int, 3> vals = types[c];
                for (int const &val : vals)
                {
                    encoded.push_back(val * 1.0);
                } 
            }
            else 
            {
                encoded.push_back(0.0000000);
                encoded.push_back(0.0000000);
                encoded.push_back(0.0000000);
            }
        }
    }
                
    return encoded;
}

static vector<double> SSEC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file, array<char, 3> types, const string& type)
{
    vector<double> encoded;

    string aas;
    string elements;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<char, int> typesCount;
    for (char const &c : types)
    {
        typesCount[c] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                elements.push_back(line[0]);
                break;
            };
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            char c = elements[i];
            if (typesCount.find(c) != typesCount.end())
                typesCount[c]++;
        }
        for (char const &c : types)
        {
            encoded.push_back(typesCount[c] / (l * 1.0));
        } 
    }
    
    return encoded;
}

static vector<double> SSPB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file, array<char, 3> types, const string& type)
{
    vector<double> encoded;

    string aas;
    map<char, vector<double>> probs;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            char currType;
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if ((counter == 3 && type.compare("ss2") == 0) || (counter == 7 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['H'].push_back(val);
            }
            else if ((counter == 4 && type.compare("ss2") == 0) || (counter == 5 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['E'].push_back(val);
            }
            else if ((counter == 5 && type.compare("ss2") == 0) || (counter == 6 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['C'].push_back(val);
                break;
            }
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (char const &c : types)
        {
            for (char const &d : types)
            {
                string key(1, c);
                key.push_back(d);
                for (int i = found; i < found + l - n; i++)
                {
                    count[key] += probs[c][i] * probs[d][i + n];
                }
            }
        }

        for (string const &key : keys)
        {
            encoded.push_back(count[key]);
        } 
    }
    
    return encoded;
}

static vector<double> SSPAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file, array<char, 3> types, const string& type)
{
    vector<double> encoded;

    string aas;
    map<char, vector<double>> probs;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<string, double> count;
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            char currType;
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if ((counter == 3 && type.compare("ss2") == 0) || (counter == 7 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['H'].push_back(val);
            }
            else if ((counter == 4 && type.compare("ss2") == 0) || (counter == 5 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['E'].push_back(val);
            }
            else if ((counter == 5 && type.compare("ss2") == 0) || (counter == 6 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['C'].push_back(val);
                break;
            }
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = 1; i <= n; i++)
        {
            string iString = to_string(i);
            for (char const &c : types)
            {
                string key = iString;
                key.push_back(c);
                for (int j = found; j < found + l - i; j++)
                { 
                    count[key] += probs[c][j] * probs[c][j + i];
                }
            }
        }

        for (string const &key : keys)
        {
            encoded.push_back(count[key] / l);
        } 
    }
    
    return encoded;
}

static vector<double> Disorder(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read disorder file
    string line;
    string aas;
    vector<double> scores;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                double val = stof(line);
                scores.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            encoded.push_back(scores[i]);
        }
    }

    return encoded;
}

static vector<double> DisorderC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read disorder file
    string line;
    string aas;
    string types;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                char val = line[0];
                types.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        double order = 0;
        double disorder = 0;
        for (int i = found; i < found + l; i++)
        {
            if (types[i] == 'D')
                disorder++;
            else
                order++;
        }
        encoded.push_back(disorder / l);
        encoded.push_back(order / l);
    }

    return encoded;
}

static vector<double> DisorderB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read disorder file
    string line;
    string aas;
    string types;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                char val = line[0];
                types.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            if (types[i] == 'D')
            {
                encoded.push_back(0.0000000);
                encoded.push_back(1.0000000);
            }
            else
            {
                encoded.push_back(1.0000000);
                encoded.push_back(0.0000000);
            }
        }
    }

    return encoded;
}

static vector<double> TA(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read disorder file
    string line;
    string aas;
    vector<double> phis;
    vector<double> psis;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line);
                phis.push_back(val);
            }
            else if (counter == 4)
            {
                double val = stof(line);
                psis.push_back(val);
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            encoded.push_back(phis[i]);
            encoded.push_back(psis[i]);
        }
    }

    return encoded;
}

static vector<double> TAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        double phiSinVal = 0;
        double phiCosVal = 0;
        double psiSinVal = 0;
        double psiCosVal = 0;
        for (int i = found; i < found + l; i++)
        {
            phiSinVal += phiSin[i];
            phiCosVal += phiCos[i];
            psiSinVal += psiSin[i];
            psiCosVal += psiCos[i];
        }
        encoded.push_back(phiSinVal / l);
        encoded.push_back(phiCosVal / l);
        encoded.push_back(psiSinVal / l);
        encoded.push_back(psiCosVal / l);
    }

    return encoded;
}

static vector<double> TAB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file)
{
    vector<double> encoded;

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        map<string, double> values;
        int l = seq.length();
        for (int i = found; i < found + l - n; i++)
        {
            values["phiSin-phiSin"] = phiSin[i] * phiSin[i + n];
            values["phiSin-phiCos"] = phiSin[i] * phiCos[i + n];
            values["phiSin-psiSin"] = phiSin[i] * psiSin[i + n];
            values["phiSin-psiCos"] = phiSin[i] * psiCos[i + n];
            values["phiCos-phiCos"] = phiCos[i] * phiCos[i + n];
            values["phiCos-psiSin"] = phiCos[i] * psiSin[i + n];
            values["phiCos-psiCos"] = phiCos[i] * psiCos[i + n];
            values["psiSin-psiSin"] = psiSin[i] * psiSin[i + n];
            values["psiSin-psiCos"] = psiSin[i] * psiCos[i + n];
            values["psiCos-psiCos"] = psiCos[i] * psiCos[i + n];
        }

        for (string const &key : keys)
        {
            encoded.push_back(values[key] / l);
        }
    }

    return encoded;
}

static vector<double> TAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file)
{
    vector<double> encoded;

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        map<string, double> values;
        int l = seq.length();
        for (int i = 1; i <= n; i++)
        {
            string iString = to_string(i);
            for (int j = found; j < found + l - i; j++)
            {
                values[iString + "-phiSin"] = phiSin[i] * phiSin[i + n];
                values[iString + "-phiCos"] = phiCos[i] * phiCos[i + n];
                values[iString + "-psiSin"] = psiSin[i] * psiSin[i + n];
                values[iString + "-psiCos"] = psiCos[i] * psiCos[i + n];
            }
        }

        for (string const &key : keys)
        {
            encoded.push_back(values[key] / l);
        }
    }

    return encoded;
}

static vector<double> ASA(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read disorder file
    string line;
    string aas;
    vector<double> asas;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 10)
            {
                double val = stof(line);
                asas.push_back(val);
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            encoded.push_back(asas[i] / l);
        }
    }

    return encoded;
}

static vector<double> KNNpeptide(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 20>> blosum, map<char, int> order,
    vector<string> trainNames, vector<string> trainSeqs, map<string, string> labelsData, set<string> labelsP, vector<int> kNums)
{
    vector<double> encoded;
    
    vector<string> labels;
    vector<double> distances;
    for (int i = 0; i < trainSeqs.size(); i++)
    {
        if (seqName != trainNames[i])
        {
            labels.push_back(labelsData[trainNames[i]]);
            double minValue = -4.0;
            double maxValue = 11.0;
            double similarity = 0;
            for (int j = 0; j < trainSeqs[i].length(); j++)
            {
                similarity += (blosum[trainSeqs[i][j]][order[seq[j]]] - minValue) / (maxValue - minValue);
            }
            similarity = 1 - similarity / trainSeqs[i].length();
            distances.push_back(similarity);
        }
    }

    quickSort(labels, distances, 0, distances.size() - 1);
    for (int const &num : kNums)
    {
        map<string, double> contents;
        for (string const &label : labelsP)
        {
            contents[label] = 0;
        }
        for (int i = 0; i < num; i++)
        {
            contents[labels[i]]++;
        }
        for (string const &label : labelsP)
        {
            encoded.push_back(contents[label] /= num);
        }
    }
    
    return encoded;
}

static vector<double> KNNprotein(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 20>> blosum, map<char, int> order,
    vector<string> trainNames, vector<string> trainSeqs, map<string, string> labelsData, set<string> labelsP, vector<int> kNums)
{
    vector<double> encoded;
    
    vector<string> labels;
    vector<double> distances;
    for (int i = 0; i < trainSeqs.size(); i++)
    {
        if (seqName != trainNames[i])
        {
            labels.push_back(labelsData[trainNames[i]]);
            double similarity = 0;
            array<string, 2> nw = calculateNeedlemanWunsch(trainSeqs[i], seq, blosum, order, -10, -1);
            for (int j = 0; j < nw[0].length(); j++)
            {
                if (nw[0][j] == nw[1][j])
                    similarity++;
            }
            similarity *= 2 / (trainSeqs[i].length() + seq.length());
            distances.push_back(similarity);
        }
    }

    quickSort(labels, distances, 0, distances.size() - 1);
    for (int const &num : kNums)
    {
        map<string, double> contents;
        for (string const &label : labelsP)
        {
            contents[label] = 0;
        }
        for (int i = 0; i < num; i++)
        {
            contents[labels[i]]++;
        }
        for (string const &label : labelsP)
        {
            encoded.push_back(contents[label] / num);
        }
    }

    return encoded;
}

static vector<double> AAPAS(const string& seq, const string& allowed, vector<string> keys, map<string, double> freqs)
{
    vector<double> encoded;

    map<string, double> count;
    
    int l = seq.length();
    for (int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key] += freqs[key];
    }
    for (string const &key : keys) 
    {
        encoded.push_back(count[key]);
    }
    return encoded;
}

static vector<double> TVD(const string& seq, const string& allowed, vector<string> keys, map<char, array<double, 10>> tvd)
{
    vector<double> encoded;
    
    int l = seq.length();
    for (char const &c : seq) 
    {
        for (double const &val : tvd[c])
        {
            encoded.push_back(val);
        }
    }
    return encoded;
}

static vector<double> CMV(const string& seq, const string& allowed, vector<string> keys)
{
    vector<double> encoded;

    map<char, double> count;
    for (char const &c : allowed)
    {
        count[c] = 0;
    }

    int l = seq.length();
    for (int i = 0; i < l; i++)
    {
        count[seq[i]] += i + 1;
    }

    for (int i = 0; i < l; i++) 
    {
        for (char const &c : allowed)
        {
            encoded.push_back(count[c] / (l * (l - 1)));
        }
    }
    return encoded;
}

static vector<double> EBGW(const string& seq, const string& allowed, vector<string> keys, map<char, array<int, 3>> groups, int n)
{
    vector<double> encoded;

    array<vector<double>, 3> count;

    int l = seq.length();

    for (int k = 1; k <= n; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            count[i].push_back(0);
        }
        int subSize = floor(k * l / (n * 1.0));
        if (subSize == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                count[i][k - 1] = 0;
            }
        }
        else
        {
            double subSizeD = subSize * 1.0;
            string sub = seq.substr(0, subSize);
            for (char const &c : sub)
            {
                for (int i = 0; i < 3; i++)
                {
                    count[i][k - 1] += groups[c][i];
                }
            }
            for (int i = 0; i < 3; i++)
            {
                count[i][k - 1] /= subSizeD;
            }
        }
    }

    for (vector<double> const &v : count)
    {
        for (double const &d : v)
        {
            encoded.push_back(d);
        }
    }

    return encoded;
}

static vector<double> PSSM(const string& seq, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<double> encoded;

    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char const &c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                encoded.push_back(stof(line));
            };
            counter++;
        }
        endLoop:;
    }
    return encoded;
}

static vector<double> PSSMAAC(const string& seq, vector<string> keys, const string& orderString, ifstream &file)
{
    vector<double> encoded;

    map<char, double> count;
    int l = seq.length();
    for (char const &c : orderString)
    {
        count[c] = 0;
    }
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char const &c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                count[orderString[counter - 2]] += val;
            };
            counter++;
        }
        endLoop:;
    }
    for (char const &c : orderString)
    {
        encoded.push_back(count[c] / (l * 1.0));
    }
    return encoded;
}

static vector<double> BiPSSM(const string& seq, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<double> encoded;

    map<string, double> count;
    int l = seq.length();
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    // Read PSSM file
    string line;
    int lines = 0;
    queue<array<int, 20>> prevLines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        array<int, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char const &c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                if (lines >= n)
                {
                    char c1 = orderString[counter - 2];
                    for (int i = 0; i < 20; i++)
                    {
                        char c2 = orderString[i];
                        int val2 = prevLines.front()[i];
                        string key(1, c1);
                        key.push_back(c2);
                        count[key] += (val * val2);
                    }
                }
            };
            counter++;
            if (counter == 22)
            {
                prevLines.push(newLine);
                if(lines >= n)
                    prevLines.pop();
                lines++;
            }
        }
        endLoop:;
    }

    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - (n * 1.0)));
    }
    return encoded;
}

static vector<double> PSSMAC(const string& seq, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<double> encoded;

    map<char, double> avg;
    map<string, double> count;
    int l = seq.length();
    for (char const &c : orderString)
    {
        avg[c] = 0;
    }
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    // Read PSSM file
    string line;
    vector<array<int, 20>> lines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        array<int, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char const &c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                avg[orderString[counter - 2]] += val;
            };
            counter++;
            if (counter == 22)
            {
                lines.push_back(newLine);
            }
        }
        endLoop:;
    }

    for (char const &c : orderString)
    {
        avg[c] /= l;
    }

    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < lines.size() - i; j++)
        {
            for (int k = 0; k < 20; k++)
            {
                char c = orderString[k];
                string key = to_string(i);
                key.push_back(c);
                count[key] += ((lines[j][k] - avg[c]) * (lines[j + i][k] - avg[c]));
            }
        }
    }
    
    for (string const &key : keys) 
    {
        encoded.push_back(count[key] / (l - (n * 1.0)));
    }
    return encoded;
}

static vector<double> PPSSM(const string& seq, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<double> encoded;
    
    vector<double> avg;
    map<char, double> avgChar;
    map<string, double> count;
    int l = seq.length();
    for (char const &c : seq)
    {
        avg.push_back(0);
    }
    for (char const &c : orderString)
    {
        avgChar[c] = 0;
    }
    for (string const &key : keys)
    {
        count[key] = 0;
    }
    // Read PSSM file
    string line;
    vector<array<double, 20>> lines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    int lineCount = 0;
    while (getline(file, line))
    {
        array<double, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for (string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char const &c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                avg[lineCount] += val;
            };
            counter++;
            if (counter == 22)
            {
                lineCount++;
                lines.push_back(newLine);
            }
        }
        endLoop:;
    }

    for (int i = 0; i < l; i++)
    {
        avg[i] /= 20;
        double den = 0;
        for (double const &val : lines[i])
        {
            den += pow(val - avg[i], 2);
        }
        den = sqrt(den / 20);
        for (int j = 0; j < 20; j++)
        {
            lines[i][j] -= avg[i];
            lines[i][j] /= den;
            string key(1, orderString[j]);
            count["M" + key] += lines[i][j];
            if (i >= n)
                count["G" + key] += pow(lines[i - n][j] - lines[i][j], 2);
        }
    }

    int counter = 0;
    for (string const &key : keys) 
    {
        if (counter < 20)
            encoded.push_back(count[key] / l);
        else
            encoded.push_back(count[key] / (l - n));
        counter++;
    }
    
    return encoded;
}

static vector<double> PseKRAAC(const string& seq, const string& allowed, vector<string> keys, map<char, int> aaMap, const string& ktuple, 
    const string& subtype, int gapLambda)
{
    vector<double> encoded;

    int increment = subtype.compare("g-gap") == 0 ? gapLambda + 1 : 1;
    int rangeCheck = subtype.compare("g-gap") == 0 ? 1 : gapLambda;

    map<string, double> counts;

    if (ktuple == "1")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            counts[to_string(aaMap[seq[i]] + 1)]++;
        }
    }
    else if (ktuple == "2")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            if (i + rangeCheck < seq.length())
            {
                counts[to_string(aaMap[seq[i]] + 1) + "_" + to_string(aaMap[seq[i + rangeCheck]] + 1)]++;
            }
        }
    }
    else if (ktuple == "3")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            if (i + rangeCheck < seq.length() && i + (2 * rangeCheck) < seq.length())
            {
                counts[to_string(aaMap[seq[i]] + 1) + "_" + to_string(aaMap[seq[i + rangeCheck]] + 1) + 
                    to_string(aaMap[seq[i + (2 * rangeCheck)]] + 1)]++;
            }
        }
    }

    for (string const &key : keys) 
    {
        encoded.push_back(counts[key]);
    }
    return encoded;
}