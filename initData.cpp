#include <map>
#include <array>

#include <unordered_set>

using namespace std;

static map<string, char> getCodonTable()
{
	map<string, char> codonTable;

	codonTable["ATA"] = 'I';
	codonTable["ATC"] = 'I';
	codonTable["ATT"] = 'I';
	codonTable["ATG"] = 'M';
    codonTable["ACA"] = 'T';
	codonTable["ACC"] = 'T';
	codonTable["ACG"] = 'T';
	codonTable["ACT"] = 'T';
    codonTable["AAC"] = 'N';
	codonTable["AAT"] = 'N';
	codonTable["AAA"] = 'K';
	codonTable["AAG"] = 'K';
    codonTable["AGC"] = 'S';
	codonTable["AGT"] = 'S';
	codonTable["AGA"] = 'R';
	codonTable["AGG"] = 'R';          
    codonTable["CTA"] = 'L';
	codonTable["CTC"] = 'L';
	codonTable["CTG"] = 'L';
	codonTable["CTT"] = 'L';
    codonTable["CCA"] = 'P';
	codonTable["CCC"] = 'P';
	codonTable["CCG"] = 'P';
	codonTable["CCT"] = 'P';
    codonTable["CAC"] = 'H';
	codonTable["CAT"] = 'H';
	codonTable["CAA"] = 'Q';
	codonTable["CAG"] = 'Q';
    codonTable["CGA"] = 'R';
	codonTable["CGC"] = 'R';
	codonTable["CGG"] = 'R';
	codonTable["CGT"] = 'R';
    codonTable["GTA"] = 'V';
	codonTable["GTC"] = 'V';
	codonTable["GTG"] = 'V';
	codonTable["GTT"] = 'V';
    codonTable["GCA"] = 'A';
	codonTable["GCC"] = 'A';
	codonTable["GCG"] = 'A';
	codonTable["GCT"] = 'A';
    codonTable["GAC"] = 'D';
	codonTable["GAT"] = 'D';
	codonTable["GAA"] = 'E';
	codonTable["GAG"] = 'E';
    codonTable["GGA"] = 'G';
	codonTable["GGC"] = 'G';
	codonTable["GGG"] = 'G';
	codonTable["GGT"] = 'G';
    codonTable["TCA"] = 'S';
	codonTable["TCC"] = 'S';
	codonTable["TCG"] = 'S';
	codonTable["TCT"] = 'S';
    codonTable["TTC"] = 'F';
	codonTable["TTT"] = 'F';
	codonTable["TTA"] = 'L';
	codonTable["TTG"] = 'L';
    codonTable["TAC"] = 'Y';
	codonTable["TAT"] = 'Y';
	// codonTable["TAA"] = '_';
	// codonTable["TAG"] = '_';
    codonTable["TGC"] = 'C';
	codonTable["TGT"] = 'C';
	// codonTable["TGA"] = '_';
	codonTable["TGG"] = 'W';

	return codonTable;
}

static map<char, array<double, 20>> getBlosum()
{
    map<char, array<double, 20>> blosum;

    blosum['A'] = { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0}; 
    blosum['R'] = {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3}; 
    blosum['N'] = {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3}; 
    blosum['D'] = {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3}; 
    blosum['C'] = { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1}; 
    blosum['Q'] = {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2}; 
    blosum['E'] = {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2}; 
    blosum['G'] = { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3}; 
    blosum['H'] = {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3}; 
    blosum['I'] = {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3}; 
    blosum['L'] = {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1}; 
    blosum['K'] = {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2}; 
    blosum['M'] = {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1}; 
    blosum['F'] = {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1}; 
    blosum['P'] = {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2}; 
    blosum['S'] = { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2}; 
    blosum['T'] = { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0}; 
    blosum['W'] = {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3}; 
    blosum['Y'] = {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1}; 
    blosum['V'] = { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}; 

    return blosum;
}

static string helpPseKRAAC()
{
    string help = "The 'raactype' value for each subtype descriptor could be chosen from:\n"
        "type 1    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 2	  [2, 3, 4, 5, 6,    8,                        15,                 20]\n"
        "type 3A   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 3B   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 4	  [         5,       8, 9,     11,     13,                         20]\n"
        "type 5    [   3, 4,          8,    10,                 15,                 20]\n"
        "type 6A   [      4, 5,                                                     20]\n"
        "type 6B   [         5,                                                       ]\n"
        "type 6C   [         5,                                                       ]\n"
        "type 7    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 8    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 9    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 10   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 11   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 12   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20]\n"
        "type 13   [      4,                        12,                 17,         20]\n"
        "type 14   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]\n"
        "type 15   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]\n"
        "type 16   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]";
    return help;
}

static map<string, vector<int>> usagePseKRAAC()
{
    map<string, vector<int>> usage;

    usage["1"]  = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["2"]  = {2, 3, 4, 5, 6,    8,                        15,                 20};
	usage["3A"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["3B"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["4"]  = {         5,       8, 9,     11,     13,                         20};
	usage["5"]  = {   3, 4,          8,    10,                 15,                 20};
	usage["6A"] = {      4, 5,                                                     20};
	usage["6B"] = {         5,                                                       };
	usage["6C"] = {         5,                                                       };
	usage["7"]  = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["8"]  = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["9"]  = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["10"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["11"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["12"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20};
	usage["13"] = {      4,                        12,                 17,         20};
	usage["14"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	usage["15"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20};
	usage["16"] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20};
 
 	return usage;
}

static map<int, vector<string>> type1()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] = {"CMFILVWY", "AGTSNQDEHRKP"};
    aaGroups[3] = {"CMFILVWY", "AGTSP", "NQDEHRK"};
    aaGroups[4] = {"CMFWY", "ILV", "AGTS", "NQDEHRKP"};
    aaGroups[5] = {"WFYH", "MILV", "CATSP", "G", "NQDERK"};
    aaGroups[6] = {"WFYH", "MILV", "CATS", "P", "G", "NQDERK"};
    aaGroups[7] = {"WFYH", "MILV", "CATS", "P", "G", "NQDE", "RK"};
    aaGroups[8] = {"WFYH", "MILV", "CA", "NTS", "P", "G", "DE", "QRK"};
    aaGroups[9] = {"WFYH", "MI", "LV", "CA", "NTS", "P", "G", "DE", "QRK"};
    aaGroups[10] = {"WFY", "ML", "IV", "CA", "TS", "NH", "P", "G", "DE", "QRK"};
    aaGroups[11] = {"WFY", "ML", "IV", "CA", "TS", "NH", "P", "G", "D", "QE", "RK"};
    aaGroups[12] = {"WFY", "ML", "IV", "C", "A", "TS", "NH", "P", "G", "D", "QE", "RK"};
    aaGroups[13] = {"WFY", "ML", "IV", "C", "A", "T", "S", "NH", "P", "G", "D", "QE", "RK"};
    aaGroups[14] = {"WFY", "ML", "IV", "C", "A", "T", "S", "NH", "P", "G", "D", "QE", "R", "K"};
    aaGroups[15] = {"WFY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "QE", "R", "K"};
    aaGroups[16] = {"W", "FY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "QE", "R", "K"};
    aaGroups[17] = {"W", "FY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
    aaGroups[18] = {"W", "FY", "M", "L", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
    aaGroups[19] = {"W", "F", "Y", "M", "L", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
    aaGroups[20] = {"W", "F", "Y", "M", "L", "I", "V", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};

    return aaGroups;
}

static map<int, vector<string>> type2()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] =  {"LVIMCAGSTPFYW", "EDNQKRH"};
	aaGroups[3] =  {"LVIMCAGSTP", "FYW", "EDNQKRH"};
	aaGroups[4] =  {"LVIMC", "AGSTP", "FYW", "EDNQKRH"};
	aaGroups[5] =  {"LVIMC", "AGSTP", "FYW", "EDNQ", "KRH"};
	aaGroups[6] =  {"LVIM", "AGST", "PHC", "FYW", "EDNQ", "KR"};
	aaGroups[8] =  {"LVIMC", "AG", "ST", "P", "FYW", "EDNQ", "KR", "H"};
	aaGroups[15] = {"LVIM", "C", "A", "G", "S", "T", "P", "FY", "W", "E", "D", "N", "Q", "KR", "H"};
	aaGroups[20] = {"L", "V", "I", "M", "C", "A", "G", "S", "T", "P", "F", "Y", "W", "E", "D", "N", "Q", "K", "R", "H"};

    return aaGroups;
}

static map<int, vector<string>> type3A()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2]  = {"AGSPDEQNHTKRMILFYVC", "W"};
	aaGroups[3]  = {"AGSPDEQNHTKRMILFYV", "W", "C"};
	aaGroups[4]  = {"AGSPDEQNHTKRMIV", "W", "YFL", "C"};
	aaGroups[5]  = {"AGSPDEQNHTKR", "W", "YF", "MIVL", "C"};
	aaGroups[6]  = {"AGSP", "DEQNHTKR", "W", "YF", "MIL", "VC"};
	aaGroups[7]  = {"AGP", "DEQNH", "TKRMIV", "W", "YF", "L", "CS"};
	aaGroups[8]  = {"AG", "DEQN", "TKRMIV", "HY", "W", "L", "FP", "CS"};
	aaGroups[9]  = {"AG", "P", "DEQN", "TKRMI", "HY", "W", "F", "L", "VCS"};
	aaGroups[10]  = {"AG", "P", "DEQN", "TKRM", "HY", "W", "F", "I", "L", "VCS"};
	aaGroups[11]  = {"AG", "P", "DEQN", "TK", "RI", "H", "Y", "W", "F", "ML", "VCS"};
	aaGroups[12]  = {"FAS", "P", "G", "DEQ", "NL", "TK", "R", "H", "W", "Y", "IM", "VC"};
	aaGroups[13]  = {"FAS", "P", "G", "DEQ", "NL", "T", "K", "R", "H", "W", "Y", "IM", "VC"};
	aaGroups[14]  = {"FA", "P", "G", "T", "DE", "QM", "NL", "K", "R", "H", "W", "Y", "IV", "CS"};
	aaGroups[15]  = {"FAS", "P", "G", "T", "DE", "Q", "NL", "K", "R", "H", "W", "Y", "M", "I", "VC"};
	aaGroups[16]  = {"FA", "P", "G", "ST", "DE", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "VC"};
	aaGroups[17]  = {"FA", "P", "G", "S", "T", "DE", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "VC"};
	aaGroups[18]  = {"FA", "P", "G", "S", "T", "DE", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "V", "C"};
	aaGroups[19]  = {"FA", "P", "G", "S", "T", "D", "E", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "V", "C"};
	aaGroups[20]  = {"F", "A", "P", "G", "S", "T", "D", "E", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "V", "C"};

    return aaGroups;
}

static map<int, vector<string>> type3B()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] = {"HRKQNEDSTGPACVIM", "LFYW"};
	aaGroups[3] = {"HRKQNEDSTGPACVIM", "LFY", "W"};
	aaGroups[4] = {"HRKQNEDSTGPA", "CIV", "MLFY", "W"};
	aaGroups[5] = {"HRKQNEDSTGPA", "CV", "IML", "FY", "W"};
	aaGroups[6] = {"HRKQNEDSTPA", "G", "CV", "IML", "FY", "W"};
	aaGroups[7] = {"HRKQNEDSTA", "G", "P", "CV", "IML", "FY", "W"};
	aaGroups[8] = {"HRKQSTA", "NED", "G", "P", "CV", "IML", "FY", "W"};
	aaGroups[9] = {"HRKQ", "NED", "ASTG", "P", "C", "IV", "MLF", "Y", "W"};
	aaGroups[10] = {"RKHSA", "Q", "NED", "G", "P", "C", "TIV", "MLF", "Y", "W"};
	aaGroups[11] = {"RKQ", "NG", "ED", "AST", "P", "C", "IV", "HML", "F", "Y", "W"};
	aaGroups[12] = {"RKQ", "ED", "NAST", "G", "P", "C", "IV", "H", "ML", "F", "Y", "W"};
	aaGroups[13] = {"RK", "QE", "D", "NG", "HA", "ST", "P", "C", "IV", "ML", "F", "Y", "W"};
	aaGroups[14] = {"R", "K", "QE", "D", "NG", "HA", "ST", "P", "C", "IV", "ML", "F", "Y", "W"};
	aaGroups[15] = {"R", "K", "QE", "D", "NG", "HA", "ST", "P", "C", "IV", "M", "L", "F", "Y", "W"};
	aaGroups[16] = {"R", "K", "Q", "E", "D", "NG", "HA", "ST", "P", "C", "IV", "M", "L", "F", "Y", "W"};
	aaGroups[17] = {"R", "K", "Q", "E", "D", "NG", "HA", "S", "T", "P", "C", "IV", "M", "L", "F", "Y", "W"};
	aaGroups[18] = {"R", "K", "Q", "E", "D", "NG", "HA", "S", "T", "P", "C", "I", "V", "M", "L", "F", "Y", "W"};
	aaGroups[19] = {"R", "K", "Q", "E", "D", "NG", "H", "A", "S", "T", "P", "C", "I", "V", "M", "L", "F", "Y", "W"};
	aaGroups[20] = {"R", "K", "Q", "E", "D", "N", "G", "H", "A", "S", "T", "P", "C", "I", "V", "M", "L", "F", "Y", "W"};

    return aaGroups;
}

static map<int, vector<string>> type4()
{
    map<int, vector<string>> aaGroups;

	aaGroups[5] = {"G", "IVFYW", "ALMEQRK", "P", "NDHSTC"};
	aaGroups[8] = {"G", "IV", "FYW", "ALM", "EQRK", "P", "ND", "HSTC"};
	aaGroups[9] = {"G", "IV", "FYW", "ALM", "EQRK", "P", "ND", "HS", "TC"};
	aaGroups[11] = {"G", "IV", "FYW", "A", "LM", "EQRK", "P", "ND", "HS", "T", "C"};
	aaGroups[13] = {"G", "IV", "FYW", "A", "L", "M", "E", "QRK", "P", "ND", "HS", "T", "C"};
	aaGroups[20] = {"G", "I", "V", "F", "Y", "W", "A", "L", "M", "E", "Q", "R", "K", "P", "N", "D", "H", "S", "T", "C"}; 

    return aaGroups;
}

static map<int, vector<string>> type5()
{
    map<int, vector<string>> aaGroups;

	aaGroups[3] = {"FWYCILMVAGSTPHNQ", "DE", "KR"};
	aaGroups[4] = {"FWY", "CILMV", "AGSTP", "EQNDHKR"};
	aaGroups[8] = {"FWY", "CILMV", "GA", "ST", "P", "EQND", "H", "KR"};
	aaGroups[10] = {"G", "FYW", "A", "ILMV", "RK", "P", "EQND", "H", "ST", "C"};
	aaGroups[15] = {"G", "FY", "W", "A", "ILMV", "E", "Q", "RK", "P", "N", "D", "H", "S", "T", "C"};
	aaGroups[20] = {"G", "I", "V", "F", "Y", "W", "A", "L", "M", "E", "Q", "R", "K", "P", "N", "D", "H", "S", "T", "C"};

    return aaGroups;
}

static map<int, vector<string>> type6A()
{
    map<int, vector<string>> aaGroups;

    aaGroups[4] = {"AGPST", "CILMV", "DEHKNQR", "FYW"};
	aaGroups[5] = {"AHT", "CFILMVWY", "DE", "GP", "KNQRS"};
	aaGroups[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

    return aaGroups;
}

static map<int, vector<string>> type6B()
{
    map<int, vector<string>> aaGroups;
    
	aaGroups[5] = {"AEHKQRST", "CFILMVWY", "DN", "G", "P"};

    return aaGroups;
}

static map<int, vector<string>> type6C()
{
    map<int, vector<string>> aaGroups;
    
	aaGroups[5] = {"AG", "C", "DEKNPQRST", "FILMVWY", "H"};

    return aaGroups;
}

static map<int, vector<string>> type7()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] = {"C", "MFILVWYAGTSNQDEHRKP"};
	aaGroups[3] = {"C", "MFILVWYAKR", "GTSNQDEHP"};
	aaGroups[4] = {"C", "KR", "MFILVWYA", "GTSNQDEHP"};
	aaGroups[5] = {"C", "KR", "MFILVWYA", "DE", "GTSNQHP"};
	aaGroups[6] = {"C", "KR", "WYA", "MFILV", "DE", "GTSNQHP"};
	aaGroups[7] = {"C", "KR", "WYA", "MFILV", "DE", "QH", "GTSNP"};
	aaGroups[8] = {"C", "KR", "WYA", "MFILV", "D", "E", "QH", "GTSNP"};
	aaGroups[9] = {"C", "KR", "WYA", "MFILV", "D", "E", "QH", "TP", "GSN"};
	aaGroups[10] = {"C", "KR", "WY", "A", "MFILV", "D", "E", "QH", "TP", "GSN"};
	aaGroups[11] = {"C", "K", "R", "WY", "A", "MFILV", "D", "E", "QH", "TP", "GSN"};
	aaGroups[12] = {"C", "K", "R", "WY", "A", "MFILV", "D", "E", "QH", "TP", "GS", "N"};
	aaGroups[13] = {"C", "K", "R", "W", "Y", "A", "MFILV", "D", "E", "QH", "TP", "GS", "N"};
	aaGroups[14] = {"C", "K", "R", "W", "Y", "A", "FILV", "M", "D", "E", "QH", "TP", "GS", "N"};
	aaGroups[15] = {"C", "K", "R", "W", "Y", "A", "FILV", "M", "D", "E", "Q", "H", "TP", "GS", "N"};
	aaGroups[16] = {"C", "K", "R", "W", "Y", "A", "FILV", "M", "D", "E", "Q", "H", "TP", "G", "S", "N"};
	aaGroups[17] = {"C", "K", "R", "W", "Y", "A", "FI", "LV", "M", "D", "E", "Q", "H", "TP", "G", "S", "N"};
	aaGroups[18] = {"C", "K", "R", "W", "Y", "A", "FI", "LV", "M", "D", "E", "Q", "H", "T", "P", "G", "S", "N"};
	aaGroups[19] = {"C", "K", "R", "W", "Y", "A", "F", "I", "LV", "M", "D", "E", "Q", "H", "T", "P", "G", "S", "N"};
	aaGroups[20] = {"C", "K", "R", "W", "Y", "A", "F", "I", "L", "V", "M", "D", "E", "Q", "H", "T", "P", "G", "S", "N"};

    return aaGroups;
}

static map<int, vector<string>> type8()
{
    map<int, vector<string>> aaGroups;
    
    aaGroups[2] = {"ADEGKNPQRST", "CFHILMVWY"};
	aaGroups[3] = {"ADEGNPST", "CHKQRW", "FILMVY"};
	aaGroups[4] = {"AGNPST", "CHWY", "DEKQR", "FILMV"};
	aaGroups[5] = {"AGPST", "CFWY", "DEN", "HKQR", "ILMV"};
	aaGroups[6] = {"APST", "CW", "DEGN", "FHY", "ILMV", "KQR"};
	aaGroups[7] = {"AGST", "CW", "DEN", "FY", "HP", "ILMV", "KQR"};
	aaGroups[8] = {"AST", "CG", "DEN", "FY", "HP", "ILV", "KQR", "MW"};
	aaGroups[9] = {"AST", "CW", "DE", "FY", "GN", "HQ", "ILV", "KR", "MP"};
	aaGroups[10] = {"AST", "CW", "DE", "FY", "GN", "HQ", "IV", "KR", "LM", "P"};
	aaGroups[11] = {"AST", "C", "DE", "FY", "GN", "HQ", "IV", "KR", "LM", "P", "W"};
	aaGroups[12] = {"AST", "C", "DE", "FY", "G", "HQ", "IV", "KR", "LM", "N", "P", "W"};
	aaGroups[13] = {"AST", "C", "DE", "FY", "G", "H", "IV", "KR", "LM", "N", "P", "Q", "W"};
	aaGroups[14] = {"AST", "C", "DE", "FL", "G", "H", "IV", "KR", "M", "N", "P", "Q", "W", "Y"};
	aaGroups[15] = {"AST", "C", "DE", "F", "G", "H", "IV", "KR", "L", "M", "N", "P", "Q", "W", "Y"};
	aaGroups[16] = {"AT", "C", "DE", "F", "G", "H", "IV", "KR", "L", "M", "N", "P", "Q", "S", "W", "Y"};
	aaGroups[17] = {"AT", "C", "DE", "F", "G", "H", "IV", "K", "L", "M", "N", "P", "Q", "R", "S", "W", "Y"};
	aaGroups[18] = {"A", "C", "DE", "F", "G", "H", "IV", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y"};
	aaGroups[19] = {"A", "C", "D", "E", "F", "G", "H", "IV", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y"};
	aaGroups[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "V", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y"};

    return aaGroups;
}

static map<int, vector<string>> type9()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] = {"ACDEFGHILMNPQRSTVWY", "K"};
	aaGroups[3] = {"ACDFGMPQRSTW", "EHILNVY", "K"};
	aaGroups[4] = {"AGPT", "CDFMQRSW", "EHILNVY", "K"};
	aaGroups[5] = {"AGPT", "CDQ", "EHILNVY", "FMRSW", "K"};
	aaGroups[6] = {"AG", "CDQ", "EHILNVY", "FMRSW", "K", "PT"};
	aaGroups[7] = {"AG", "CDQ", "EHNY", "FMRSW", "ILV", "K", "PT"};
	aaGroups[8] = {"AG", "C", "DQ", "EHNY", "FMRSW", "ILV", "K", "PT"};
	aaGroups[9] = {"AG", "C", "DQ", "EHNY", "FMW", "ILV", "K", "PT", "RS"};
	aaGroups[10] = {"A", "C", "DQ", "EHNY", "FMW", "G", "ILV", "K", "PT", "RS"};
	aaGroups[11] = {"A", "C", "DQ", "EHNY", "FM", "G", "ILV", "K", "PT", "RS", "W"};
	aaGroups[12] = {"A", "C", "DQ", "EHNY", "FM", "G", "IL", "K", "PT", "RS", "V", "W"};
	aaGroups[13] = {"A", "C", "DQ", "E", "FM", "G", "HNY", "IL", "K", "PT", "RS", "V", "W"};
	aaGroups[14] = {"A", "C", "D", "E", "FM", "G", "HNY", "IL", "K", "PT", "Q", "RS", "V", "W"};
	aaGroups[15] = {"A", "C", "D", "E", "FM", "G", "HNY", "IL", "K", "PT", "Q", "R", "S", "V", "W"};
	aaGroups[16] = {"A", "C", "D", "E", "F", "G", "HNY", "IL", "K", "M", "PT", "Q", "R", "S", "V", "W"};
	aaGroups[17] = {"A", "C", "D", "E", "F", "G", "HNY", "IL", "K", "M", "P", "Q", "R", "S", "T", "V", "W"};
	aaGroups[18] = {"A", "C", "D", "E", "F", "G", "HNY", "I", "K", "L", "M", "P", "Q", "R", "S", "T", "V", "W"};
	aaGroups[19] = {"A", "C", "D", "E", "F", "G", "HN", "I", "K", "L", "M", "P", "Q", "R", "S", "T", "V", "W", "Y"};
	aaGroups[20] = {"A", "C", "D", "E", "F", "G", "H", "N", "I", "K", "L", "M", "P", "Q", "R", "S", "T", "V", "W", "Y"};

    return aaGroups;
}

static map<int, vector<string>> type10()
{
    map<int, vector<string>> aaGroups;

	aaGroups[2] = {"CMFILVWY", "AGTSNQDEHRKP"};
	aaGroups[3] = {"CMFILVWY", "AGTSP", "NQDEHRK"};
	aaGroups[4] = {"CMFWY", "ILV", "AGTS", "NQDEHRKP"};
	aaGroups[5] = {"FWYH", "MILV", "CATSP", "G", "NQDERK"};
	aaGroups[6] = {"FWYH", "MILV", "CATS", "P", "G", "NQDERK"};
	aaGroups[7] = {"FWYH", "MILV", "CATS", "P", "G", "NQDE", "RK"};
	aaGroups[8] = {"FWYH", "MILV", "CA", "NTS", "P", "G", "DE", "QRK"};
	aaGroups[9] = {"FWYH", "ML", "IV", "CA", "NTS", "P", "G", "DE", "QRK"};
	aaGroups[10] = {"FWY", "ML", "IV", "CA", "TS", "NH", "P", "G", "DE", "QRK"};
	aaGroups[11] = {"FWY", "ML", "IV", "CA", "TS", "NH", "P", "G", "D", "QE", "RK"};
	aaGroups[12] = {"FWY", "ML", "IV", "C", "A", "TS", "NH", "P", "G", "D", "QE", "RK"};
	aaGroups[13] = {"FWY", "ML", "IV", "C", "A", "T", "S", "NH", "P", "G", "D", "QE", "RK"};
	aaGroups[14] = {"FWY", "ML", "IV", "C", "A", "T", "S", "NH", "P", "G", "D", "QE", "R", "K"};
	aaGroups[15] = {"FWY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "QE", "R", "K"};
	aaGroups[16] = {"W", "FY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "QE", "R", "K"};
	aaGroups[17] = {"W", "FY", "ML", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
	aaGroups[18] = {"W", "FY", "M", "L", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
	aaGroups[19] = {"W", "F", "Y", "M", "L", "IV", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};
	aaGroups[20] = {"W", "F", "Y", "M", "L", "I", "V", "C", "A", "T", "S", "N", "H", "P", "G", "D", "Q", "E", "R", "K"};

    return aaGroups;
}

static map<int, vector<string>> type11()
{
    map<int, vector<string>> aaGroups;

	aaGroups[2] = {"CFYWMLIV", "GPATSNHQEDRK"};
	aaGroups[3] = {"CFYWMLIV", "GPATS", "NHQEDRK"};
	aaGroups[4] = {"CFYW", "MLIV", "GPATS", "NHQEDRK"};
	aaGroups[5] = {"CFYW", "MLIV", "G", "PATS", "NHQEDRK"};
	aaGroups[6] = {"CFYW", "MLIV", "G", "P", "ATS", "NHQEDRK"};
	aaGroups[7] = {"CFYW", "MLIV", "G", "P", "ATS", "NHQED", "RK"};
	aaGroups[8] = {"CFYW", "MLIV", "G", "P", "ATS", "NH", "QED", "RK"};
	aaGroups[9] = {"CFYW", "ML", "IV", "G", "P", "ATS", "NH", "QED", "RK"};
	aaGroups[10] = {"C", "FYW", "ML", "IV", "G", "P", "ATS", "NH", "QED", "RK"};
	aaGroups[11] = {"C", "FYW", "ML", "IV", "G", "P", "A", "TS", "NH", "QED", "RK"};
	aaGroups[12] = {"C", "FYW", "ML", "IV", "G", "P", "A", "TS", "NH", "QE", "D", "RK"};
	aaGroups[13] = {"C", "FYW", "ML", "IV", "G", "P", "A", "T", "S", "NH", "QE", "D", "RK"};
	aaGroups[14] = {"C", "FYW", "ML", "IV", "G", "P", "A", "T", "S", "N", "H", "QE", "D", "RK"};
	aaGroups[15] = {"C", "FYW", "ML", "IV", "G", "P", "A", "T", "S", "N", "H", "QE", "D", "R", "K"};
	aaGroups[16] = {"C", "FY", "W", "ML", "IV", "G", "P", "A", "T", "S", "N", "H", "QE", "D", "R", "K"};
	aaGroups[17] = {"C", "FY", "W", "ML", "IV", "G", "P", "A", "T", "S", "N", "H", "Q", "E", "D", "R", "K"};
	aaGroups[18] = {"C", "FY", "W", "M", "L", "IV", "G", "P", "A", "T", "S", "N", "H", "Q", "E", "D", "R", "K"};
	aaGroups[19] = {"C", "F", "Y", "W", "M", "L", "IV", "G", "P", "A", "T", "S", "N", "H", "Q", "E", "D", "R", "K"};
	aaGroups[20] = {"C", "F", "Y", "W", "M", "L", "I", "V", "G", "P", "A", "T", "S", "N", "H", "Q", "E", "D", "R", "K"};

    return aaGroups;
}

static map<int, vector<string>> type12()
{
    map<int, vector<string>> aaGroups;

    aaGroups[2] = {"IVMLFWYC", "ARNDQEGHKPST"};
	aaGroups[3] = {"IVLMFWC", "YA", "RNDQEGHKPST"};
	aaGroups[4] = {"IVLMFW", "C", "YA", "RNDQEGHKPST"};
	aaGroups[5] = {"IVLMFW", "C", "YA", "G", "RNDQEHKPST"};
	aaGroups[6] = {"IVLMF", "WY", "C", "AH", "G", "RNDQEKPST"};
	aaGroups[7] = {"IVLMF", "WY", "C", "AH", "GP", "R", "NDQEKST"};
	aaGroups[8] = {"IVLMF", "WY", "C", "A", "G", "R", "Q", "NDEHKPST"};
	aaGroups[9] = {"IVLMF", "WY", "C", "A", "G", "P", "H", "K", "RNDQEST"};
	aaGroups[10] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "RN", "DQEKPST"};
	aaGroups[11] = {"IVLMF", "W", "Y", "C", "A", "H", "G", "R", "N", "Q", "DEKPST"};
	aaGroups[12] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "T", "RDEKPS"};
	aaGroups[13] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "DEKST"};
	aaGroups[14] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "DEST"};
	aaGroups[15] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "D", "EST"};
	aaGroups[16] = {"IVLM", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "DE"};
	aaGroups[17] = {"IVL", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "DE"};
	aaGroups[18] = {"IVL", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "D", "E"};
	aaGroups[20] = {"I", "V", "L", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "D", "E"};

    return aaGroups;
}

static map<int, vector<string>> type13()
{
    map<int, vector<string>> aaGroups;

    aaGroups[4] = {"ADKERNTSQ", "YFLIVMCWH", "G", "P"};
	aaGroups[12] = {"A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"};
	aaGroups[17] = {"A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"};
	aaGroups[20] = {"A", "D", "K", "E", "R", "N", "T", "S", "Q", "Y", "F", "L", "I", "V", "M", "C", "W", "H", "G", "P"};

    return aaGroups;
}

static map<int, vector<string>> type14()
{
    map<int, vector<string>> aaGroups;

	aaGroups[2] = {"ARNDCQEGHKPST", "ILMFWYV"};
	aaGroups[3] = {"ARNDQEGHKPST", "C", "ILMFWYV"};
	aaGroups[4] = {"ARNDQEGHKPST", "C", "ILMFYV", "W"};
	aaGroups[5] = {"AGPST", "RNDQEHK", "C", "ILMFYV", "W"};
	aaGroups[6] = {"AGPST", "RNDQEK", "C", "H", "ILMFYV", "W"};
	aaGroups[7] = {"ANDGST", "RQEK", "C", "H", "ILMFYV", "P", "W"};
	aaGroups[8] = {"ANDGST", "RQEK", "C", "H", "ILMV", "FY", "P", "W"};
	aaGroups[9] = {"AGST", "RQEK", "ND", "C", "H", "ILMV", "FY", "P", "W"};
	aaGroups[10] = {"AGST", "RK", "ND", "C", "QE", "H", "ILMV", "FY", "P", "W"};
	aaGroups[11] = {"AST", "RK", "ND", "C", "QE", "G", "H", "ILMV", "FY", "P", "W"};
	aaGroups[12] = {"AST", "RK", "ND", "C", "QE", "G", "H", "IV", "LM", "FY", "P", "W"};
	aaGroups[13] = {"AST", "RK", "N", "D", "C", "QE", "G", "H", "IV", "LM", "FY", "P", "W"};
	aaGroups[14] = {"AST", "RK", "N", "D", "C", "Q", "E", "G", "H", "IV", "LM", "FY", "P", "W"};
	aaGroups[15] = {"A", "RK", "N", "D", "C", "Q", "E", "G", "H", "IV", "LM", "FY", "P", "ST", "W"};
	aaGroups[16] = {"A", "RK", "N", "D", "C", "Q", "E", "G", "H", "IV", "LM", "F", "P", "ST", "W", "Y"};
	aaGroups[17] = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "IV", "LM", "K", "F", "P", "ST", "W", "Y"};
	aaGroups[18] = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "IV", "LM", "K", "F", "P", "S", "T", "W", "Y"};
	aaGroups[19] = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "IV", "L", "K", "M", "F", "P", "S", "T", "W", "Y"};
	aaGroups[20] = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "V", "L", "K", "M", "F", "P", "S", "T", "W", "Y"};

    return aaGroups;
}

static map<int, vector<string>> type15()
{
    map<int, vector<string>> aaGroups;

	aaGroups[2] = {"MFILVAW", "CYQHPGTSNRKDE"};
	aaGroups[3] = {"MFILVAW", "CYQHPGTSNRK", "DE"};
	aaGroups[4] = {"MFILV", "ACW", "YQHPGTSNRK", "DE"};
	aaGroups[5] = {"MFILV", "ACW", "YQHPGTSN", "RK", "DE"};
	aaGroups[6] = {"MFILV", "A", "C", "WYQHPGTSN", "RK", "DE"};
	aaGroups[7] = {"MFILV", "A", "C", "WYQHP", "GTSN", "RK", "DE"};
	aaGroups[8] = {"MFILV", "A", "C", "WYQHP", "G", "TSN", "RK", "DE"};
	aaGroups[9] = {"MF", "ILV", "A", "C", "WYQHP", "G", "TSN", "RK", "DE"};
	aaGroups[10] = {"MF", "ILV", "A", "C", "WYQHP", "G", "TSN", "RK", "D", "E"};
	aaGroups[11] = {"MF", "IL", "V", "A", "C", "WYQHP", "G", "TSN", "RK", "D", "E"};
	aaGroups[12] = {"MF", "IL", "V", "A", "C", "WYQHP", "G", "TS", "N", "RK", "D", "E"};
	aaGroups[13] = {"MF", "IL", "V", "A", "C", "WYQHP", "G", "T", "S", "N", "RK", "D", "E"};
	aaGroups[14] = {"MF", "I", "L", "V", "A", "C", "WYQHP", "G", "T", "S", "N", "RK", "D", "E"};
	aaGroups[15] = {"MF", "IL", "V", "A", "C", "WYQ", "H", "P", "G", "T", "S", "N", "RK", "D", "E"};
	aaGroups[16] = {"MF", "I", "L", "V", "A", "C", "WYQ", "H", "P", "G", "T", "S", "N", "RK", "D", "E"};
	aaGroups[20] = {"M", "F", "I", "L", "V", "A", "C", "W", "Y", "Q", "H", "P", "G", "T", "S", "N", "R", "K", "D", "E"};

    return aaGroups;
}

static map<int, vector<string>> type16()
{
    map<int, vector<string>> aaGroups;

	aaGroups[2] = {"IMVLFWY", "GPCASTNHQEDRK"};
	aaGroups[3] = {"IMVLFWY", "GPCAST", "NHQEDRK"};
	aaGroups[4] = {"IMVLFWY", "G", "PCAST", "NHQEDRK"};
	aaGroups[5] = {"IMVL", "FWY", "G", "PCAST", "NHQEDRK"};
	aaGroups[6] = {"IMVL", "FWY", "G", "P", "CAST", "NHQEDRK"};
	aaGroups[7] = {"IMVL", "FWY", "G", "P", "CAST", "NHQED", "RK"};
	aaGroups[8] = {"IMV", "L", "FWY", "G", "P", "CAST", "NHQED", "RK"};
	aaGroups[9] = {"IMV", "L", "FWY", "G", "P", "C", "AST", "NHQED", "RK"};
	aaGroups[10] = {"IMV", "L", "FWY", "G", "P", "C", "A", "STNH", "RKQE", "D"};
	aaGroups[11] = {"IMV", "L", "FWY", "G", "P", "C", "A", "STNH", "RKQ", "E", "D"};
	aaGroups[12] = {"IMV", "L", "FWY", "G", "P", "C", "A", "ST", "N", "HRKQ", "E", "D"};
	aaGroups[13] = {"IMV", "L", "F", "WY", "G", "P", "C", "A", "ST", "N", "HRKQ", "E", "D"};
	aaGroups[14] = {"IMV", "L", "F", "WY", "G", "P", "C", "A", "S", "T", "N", "HRKQ", "E", "D"};
	aaGroups[15] = {"IMV", "L", "F", "WY", "G", "P", "C", "A", "S", "T", "N", "H", "RKQ", "E", "D"};
	aaGroups[16] = {"IMV", "L", "F", "W", "Y", "G", "P", "C", "A", "S", "T", "N", "H", "RKQ", "E", "D"};
	aaGroups[20] = {"I", "M", "V", "L", "F", "W", "Y", "G", "P", "C", "A", "S", "T", "N", "H", "R", "K", "Q", "E", "D"};

    return aaGroups;
}