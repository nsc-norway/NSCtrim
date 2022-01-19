using namespace std;


class DegenerateMatch {
    public:
    static bool match(char templ, char test) {
        switch (templ) {
            case 'A': return test == 'A';
            case 'C': return test == 'C';
            case 'G': return test == 'G';
            case 'T': return test == 'T';
            case 'U': return test == 'U';
            case 'W': return test == 'A' || test == 'T';
            case 'S': return test == 'C' || test == 'G';
            case 'M': return test == 'A' || test == 'C';
            case 'K': return test == 'G' || test == 'T';
            case 'R': return test == 'A' || test == 'G';
            case 'Y': return test == 'C' || test == 'T';
            case 'B': return test == 'C' || test == 'G' || test == 'T';
            case 'D': return test == 'A' || test == 'C' || test == 'T';
            case 'V': return test == 'A' || test == 'C' || test == 'G';
            case 'N': return true;
            default:
                return false;
        }
    }
};

class ExactMatch {
    public:
    static bool match(char templ, char test) {
        return templ == test;
    }
};



/** 
 * mismatch_and_alignment
 *
 * This function performs two operations:
 * - Computes the Levenshtein distance between s1 and s2, up to a maximum bound
 * on edit operations

 * If there are more characters left in s2 after all characters in s1 are used,
 * there is no penalty for these extra characters.
 *
 * This algorithm is faster than a naive approach, as it applies the max_val
 * bounds at every step, excluding a lot of computations which are guaranteed
 * to yield a distance above or equal to max_val. It aborts as soon as there
 * is no chance of getting a value below max_val.
 *
 * Parameters:
 *  max_val: The max value to return for Levenshtein distance of 
 *                   string (all higher values are clipped). Such a
 *                   threshold is relevant only for performance reasons.
 *  s1:              The first string, representing the known primer sequence.
 *  s2:              The data to examine for alignment and mismatches. The full
 *                   length read can be passed without performance impact.
 *                   
 *
 * Returns:
 *  first:  The edit distance between s2 and s1 using the best matching substring
            of s2 (starting at 0). If the distance exceeds max_val, that value is
            returned instead.
 *  second: The length of s2 when aligned to the full sequence of s1. The length
 *          of the match with fewest edit operations is returned. If several
 *          matches have equal distance, the shortest match is returned (please
 *          see loop before return statement, this could change). If the edit
 *          distance of the primer is greater or equal to max_val, the
 *          value of second is undefined.
 *
 */

template <typename NucleotideMatch>
pair<int,int> mismatch_and_alignment(int max_val, const string& s1, const string& s2)
{
    size_t s1len = s1.length();
    size_t s2len = min(s2.length(), s1len + (size_t)max_val);

    // Buffer for two columns of mismatch value results
    int buf1[s2len+1], buf2[s2len+1];

    int* column = buf1, * prevcol = buf2;

    // Columns include a dummy cell at index 0. Every time an element in
    // the column is indexed, one has to add 1 to the index i.

    // Initialise "prev" column before the actual data. This encodes the 
    // penalties for insertions at the beginning of the aligned string.
    // As a special case, prevcol[0] is used in case there is a match at
    // the first position of s1 and s2.
    for (int j=0; j<s2len+1; j++) {
        prevcol[j] = j;
    }
    
    // - first points to the first row to analyse. It is greater than 0 when
    // rows are excluded by the max_val bound. first is the actual position
    // in s2, not the index into the column array (which is first+1).
    // - last is the index of the last row to analyse (position in s2). Points
    // to the last element in column with value below max_val. The initial value
    // represents max_val-1 insertions before s2. Each column opens up another
    // row for last.
    int first = 0, last = min(max_val-1, (int)s2len-1);
    int primer_mismatches = max_val;
    for (int j = 0; j < s1len; j++) { // Column!

        int i; // i is an index of the row / position in s2, does not include the
               // offset of 1, as in the column arrays.

        int minval = max_val;
        
        if (first > 0) {
            // If first is nonzero, some rows have already been deemed excluded
            // Row first here correspond to position first - 1 (due to one-based
            // offset in the array)
            column[first] = max_val;
        }
        else {
            // The first element of column represents deletions at the start of the
            // aligned string. Not an actual position in s2.
            column[0] = j+1;
        }

        int last_below = -1, first_tmp = first;
        first = s2len; // unless found, set to position past end of list

        // Second loop: find new values in range; keep track of last row below
        // threshold, so we can possibly constrain the "last"
        for (i=first_tmp; i <= last && i < s2len; i++) {
            auto possibilities = {
                prevcol[1 + i] + 1,
                column[i] + 1,
                prevcol[i] + (NucleotideMatch::match(s1[j], s2[i]) ? 0 : 1)
            };
            int best_possibility = std::min(possibilities);
            column[i+1] = best_possibility;
            if (best_possibility < max_val) {
                if (first == s2len) first = i;
                last_below = i;
                minval = min(best_possibility, minval);
            }
        }
        // Possibly, we allow to go down one diagonally in case we're still within
        // the threshold
        if (i == last + 1 && i < s2len) {
            auto possibilities = {
                column[i] + 1,
                prevcol[i] + (NucleotideMatch::match(s1[j], s2[i]) ? 0 : 1)
            };
            int best_possibility = std::min(possibilities);
            column[i + 1] = best_possibility;
            if (best_possibility < max_val) {
                if (first == s2len) first = i;
                last_below = i;
                minval = min(best_possibility, minval);
            }
        }

        // Next round 
        last = last_below;
        
        if (minval >= max_val) {
            return make_pair(max_val, 0);
        }
        
        swap(column, prevcol);
        if (first > last) {
            return make_pair(primer_mismatches, 0);
        }
    }
    int best_len = 0, minval = max_val;
    for (int i=first; i<=last; ++i) { // Prefer long match when equal score
    //for (int i=last; i>=first; --i) { // Prefer short match when equal score
        if (prevcol[i+1] < minval) {
            minval = prevcol[i+1];
            best_len = i+1;
        }
    }
    return make_pair(minval, best_len);
}

