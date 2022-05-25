#include <limits>
#include <vector>
#include <array>
#include <map>
#include <algorithm>

using namespace std;

// Struct for matrix cells
typedef struct Cell 
{
    int score;
    int up;
    int left;
} Cell;

static int getGap(int gapUp, int gapLeft, int length)
{
    if (length > 0)
        return gapUp + (length - 1) * gapLeft;
    else   
        return 0;
}

static array<string, 2> calculateNeedlemanWunsch(string seq1, string seq2, map<char, array<double, 20>> blosum, map<char, int> order, int gapUp, int gapLeft)
{
    const int negInf = numeric_limits<int>::max() * -1;
    vector<vector<Cell>> matrix;
    for (int i = 0; i < seq1.length() + 1; i++)
    {
        vector<Cell> row;
        for (int j = 0; j < seq2.length() + 1; j++)
        {
            Cell cell = { 0, negInf, negInf };
            row.push_back(cell);
        }
        matrix.push_back(row);
    }
    
    // Fill first row
    for (int i = 1; i < matrix.size(); i++)
    {
        int gap = getGap(gapUp, gapLeft, i);
        Cell cell = { gap, negInf, gap };
        matrix[i][0] = cell;
    }

    // Fill first column
    for (int i = 1; i < matrix[0].size(); i++)
    {
        int gap = getGap(gapUp, gapLeft, i);
        Cell cell = { gap, gap, negInf };
        matrix[0][i] = cell;
    }

    // Fill the matrix
    for (int i = 1; i < matrix.size(); i++)
    {
        for (int j = 1; j < matrix[i].size(); j++)
        {
            int diag = matrix[i - 1][j - 1].score + blosum[seq1[i - 1]][order[seq2[j - 1]]];
            int up = max(matrix[i][j - 1].score + gapUp, matrix[i][j - 1].up + gapLeft);
            int left = max(matrix[i - 1][j].score + gapUp, matrix[i - 1][j].left + gapLeft);
            int score = max(diag, up);
            score = max(score, left);
            Cell cell = { score, up, left };
            matrix[i][j] = cell;
        }
    }

    // Get aligned strings of sequences
    string res1 = "";
    string res2 = "";
    int i = seq1.length();
    int j = seq2.length();
    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0 && matrix[i][j].score == matrix[i - 1][j - 1].score + blosum[seq1[i - 1]][order[seq2[j - 1]]])
        {
            i--;
            j--;
            res1.insert(0, 1, seq1[i]);
            res2.insert(0, 1, seq2[j]);
        }
        else if (i > 0 && matrix[i][j].score == matrix[i][j].left)
        {
            i--;
            res1.insert(0, 1, seq1[i]);
            res2.insert(0, 1, '-');
        }
        else if (j > 0 && matrix[i][j].score == matrix[i][j].up)
        {
            j--;
            res1.insert(0, 1, '-');
            res2.insert(0, 1, seq2[j]);
        }
    }

    array<string, 2> nw = { res1, res2 };
    return nw;
}