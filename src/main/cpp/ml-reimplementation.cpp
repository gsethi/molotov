/*
**    Copyright (C) 2003-2011 Institute for Systems Biology
**                            Seattle, Washington, USA.
**
**    This library is free software; you can redistribute it and/or
**    modify it under the terms of the GNU Lesser General Public
**    License as published by the Free Software Foundation; either
**    version 2.1 of the License, or (at your option) any later version.
**
**    This library is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**    Lesser General Public License for more details.
**
**    You should have received a copy of the GNU Lesser General Public
**    License along with this library; if not, write to the Free Software
**    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
*/

/**
 * Molotov (Motif Locator with a Vengeance)
 *
 *  Reimplementation of a reverse-engineered R implementation of the motif
 *  location scoring algorithm described in G. Thijs 2003,
 *  "Probabilistic methods to search for regulatory elements in sets of
 *   coregulated genes." 
 *
 *  Full thesis here:
 *   ftp://ftp.esat.kuleuven.ac.be/sista/thijs/reports/phd.pdf
 *
 *  This implementation requires the Armadillo linear algebra library to
 *  run. It can be downloaded here:
 *   http://arma.sourceforge.net
 *
 * Initial Revision: June 15th, 2011
 * Author: Thomas Robinson, with code borrowed from Vesteinn Thorsson's
 *  implementation in R
 *
 * Usage: simply run the tool, supplying your input sequence, precomputed
 *  background model, and thresholds file similar to the original Motif
 *  Locator. Three additional switches have been provided for your use that
 *  modify the output of the tool:
 *
 * + -filter:   explicitly filter results using the supplied threshold
 *               file. Otherwise, ALL results will be output to stdout.
 * 
 * + -quiet:   disables all unnecessary print statements in favor of plain,
 *              elegant, parseable output. Errors will still be output on stderr,
 *              of course.
 * 
 * + -reverse: instead of the given input, compute the reverse strand. This
 *              will take the given input and calculate its reverse complement,
 *              which is extremely useful when this data is otherwise
 *              unavailable. 
 *
 * Conventions:
 * + m0: zeroth order model, as vector
 * + m1: first order model, as vector, following Thijs' convention
 * + m2: second order matrix of conditional probabilities, following Thijs'
 *        convention
 * + ind: index vector of dinucleotide
 * + PWM: position weight matrix
 *
 * Notes:
 * + R uses one-based collections, whereas C++ uses a zero-based system. We
 *    follow C++'s conventions here, updating indices where appropriate.
 *
 * Todos:
 * + Support an extended background model for FASTA multiples other than "N". This
 *    tool takes a naive averaging approach to the HMM background calculations of
 *    multiples, which lacks robustness against certain (pathological) types of
 *    sequence data.
 * + Short-circuit when a stream containing all variants of 'N' or 'X' is
 *    specified. Currently, we suboptimize to exit quickly, but not immediately.
 */

#include <assert.h>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>
#include <armadillo>

using namespace std;
using namespace arma;

#define u8 unsigned char
#define u32 unsigned int
#define __LINE_MAX 2048
#define __READ_MAX 268435456
#define FASTA_BIT_A 0x1
#define FASTA_BIT_C 0x2
#define FASTA_BIT_G 0x4
#define FASTA_BIT_T 0x8

/**
 * State globals, for keeping our data tidy and organized. Trivially
 * replaceable via pointer, but we treat these as singletons anyway. 
 */
//unordered_map<string, pair<double, double>, hash<string>, strequiv> motifScores;
map<string, pair<mat, double> > pwmMap;
bool quiet = false;
bool filter = false;
bool reverseComplement = false;
bool warnedAboutUracil = false;
bool prevWasMeaningless = false;
char strand = '+';
 
/**
 * Populate the positional weight matrix map with all of the expected
 * threshold and matrix files.
 */
void populatePwmData(ifstream* ifs, string path) {
    assert(ifs != NULL);

    char buf[__LINE_MAX];
    char seps[] = "\t";

    
    while (ifs->good()) {
        ifs->getline(buf, __LINE_MAX);
        stringstream filePath;
        stringstream matStream;
        
        if (buf[0] != 0) { // Check if the string is greater than length 0
            char* token = strtok(buf, seps);
            if (token == NULL) {
                cerr << "Row length was insufficient for the supplied threshold file. Aborting." << endl;
            }

            string fileName(token);
            filePath << path << "/" << fileName;
            token = strtok(NULL, seps);
            if (token == NULL) {
                cerr << "Row length was insufficient for the supplied threshold file. Aborting." << endl;
            }
            string threshold(token);

            // Skip any lines beginning with a comment, then load the remainder
            // as a matrix.
            ifstream pwmIfs(filePath.str().c_str());
            if (pwmIfs == NULL) {
                filePath << ".dat";
                pwmIfs.open(filePath.str().c_str());
                if (pwmIfs == NULL) {
                    cerr << "File could not be loaded: " << path << "/"
                         << fileName << endl; 
                }
            }

            char buf[__LINE_MAX];
            while (pwmIfs.good()) {
                buf[0] = 0x00;
                pwmIfs.getline(buf, __LINE_MAX);
                if (buf[0] != 0x00 && buf[0] != '#') {
                    matStream << buf << endl;
                }
            }
            pwmIfs.close();

            if (!quiet) { cout << "Processing positional weight matrix " << filePath.str() << endl; }
            mat pwm; pwm.load(matStream);

            pwmMap[fileName] = pair<mat,double>(pwm, atof(threshold.c_str()));
        }
    }
}

/**
 * Scan a FASTA file, reading the contents into memory until we hit our maximum
 * read length. Beware: we assume that each file represents exactly one sequence,
 * ignoring comments and headers.
 */
string fasta2Seq(ifstream* ifs) {
    assert(ifs != NULL);

    char buf[__LINE_MAX] = {0};
    stringstream ret;

    while (!ifs->eof() && !ifs->fail()) {
        buf[0] = 0x00;
        ifs->getline(buf, __LINE_MAX);
        if (buf[0] != 0x00 && buf[0] != '>' && buf[0] != '#' && buf[0] != ' ') {
            ret << buf;
        } else if (buf[0] == '>') {
            // Inane: check if a sequence is of negative strand. This check is
            //  non-robust, can easily lead to false positives, and is the only
            //  method given to determine the strand of a given sequence.
            //
            // We further assume that the strandedness is the same for the
            //  entirety of a sequence. If it is not, the sequence should be
            //  further subdivided until this is true.
            for (u32 i = 1; i < __LINE_MAX; ++i) {
                if (buf[i] == '-') {
                    strand = '-';
                    break;
                }
            }
        }
    }
    return ret.str();
}

/**
 * Get the bitfield for a given FASTA character, accounting for the reverse
 * complement.
 */
u32 getBitfieldFromFastaChar(char c) {
    u32 ret = 0;
    switch(c) {
    case 'A':       // Bit 0x1
    case 'a': 
        ret = FASTA_BIT_A;
        break;
    case 'C':       // Bit 0x2
    case 'c':
        ret = FASTA_BIT_C;
        break;
    case 'M':
    case 'm':
        ret = FASTA_BIT_A | FASTA_BIT_C;
        break;
    case 'G':       // Bit 0x4
    case 'g':
        ret = FASTA_BIT_G;
        break;
    case 'R':
    case 'r':
        ret = FASTA_BIT_A | FASTA_BIT_G;
        break;
    case 'S':
    case 's':
        ret = FASTA_BIT_C | FASTA_BIT_G;
        break;
    case 'V':
    case 'v':
        ret = FASTA_BIT_A | FASTA_BIT_C | FASTA_BIT_G;
        break;
    case 'T':       // Bit 0x8
    case 't':
        ret = FASTA_BIT_T;
        break;
    case 'W':
    case 'w':
        ret = FASTA_BIT_A | FASTA_BIT_T;
        break;
    case 'Y':
    case 'y':
        ret = FASTA_BIT_C | FASTA_BIT_T;
        break;
    case 'H':
    case 'h':
        ret = FASTA_BIT_A | FASTA_BIT_C | FASTA_BIT_T;
        break;
    case 'K':
    case 'k':
        ret = FASTA_BIT_G | FASTA_BIT_T;
        break;
    case 'D':
    case 'd':
        ret = FASTA_BIT_A | FASTA_BIT_G | FASTA_BIT_T;
        break;
    case 'B':
    case 'b':
        ret = FASTA_BIT_C | FASTA_BIT_G | FASTA_BIT_T;
        break;
    case 'N':
    case 'n':
        ret = FASTA_BIT_A | FASTA_BIT_C | FASTA_BIT_G | FASTA_BIT_T;
        break;
    case 'U':
    case 'u':
        if (!warnedAboutUracil) {
            cerr << "WARNING: Uracil characters will be demoted to 'X' for "
                 << "sequence processing. This tool is NOT optimized for RNA sequences. "
                 << "You may experience spurious results!" << endl;
            warnedAboutUracil = true;
        }
        ret = 0x0;
        break;
    case 'X':
    case 'x':
        ret = 0x0;
        break;
    case '-':
        cerr << "The presence of gaps of indeterminate length (-) was detected "
             << "in the sequence. Please split the input sequence along these gaps to get "
             << "usable results from this tool. Aborting." << endl;
        exit(-1);
            
    default:
        cerr << "Received unexpected char " << c << " from input. Aborting." << endl;
        exit(-1);
    }

    // If the given character is lowercase, set bit 5 signaling such
    if (c >= 'a' && c <= 'z') {
        ret |= 0x10;
    }

    // If we're taking the reverse complement, swap A with T, G with C.
    if (reverseComplement) {
        ret = ((ret & 0x8) >> 3) + ((ret & 0x4) >> 1) + ((ret & 0x2) << 1) + ((ret & 0x1) << 3);
    }
    return ret;
}

/**
 * Inversion of the above, only performing the character conversion. The
 * reverse complement is not reprocessed here.
 */
char getFastaCharFromBitfield(u32 u) {
    char ret = '?';
    ret = "XACMGRSVTWYHKDBN"[u & 0xf]; // Translate back to the character representation

    // If a lowercase character was originally specified, "promote" the return
    // character into it.
    if (u & 0x10) {
        u += ('a' - 'A');
    }

    return ret;
}

/**
 * Check if a subsequence only contains meaningless data
 */
inline bool hasMeaningfulSignal(const char* given, size_t length, bool prevWasMeaningless) {
    assert(given != NULL);
    
    // Ignore sequences containing all "N" or "X", as no signal can
    // meaningfully be derived from them.

    // Optimization: if we know the previous previous char-length step of the
    // sequence was meaningless, we only need to check the last character.
    if (prevWasMeaningless) {
        char c = given[length - 1];
        if (c != 'N' && c != 'n' && c != 'X' && c != 'x') {
            return false;
        }
    }
    
    size_t i;
    for (i = 0; i < length; ++i) {
        char c = given[i];
        if (c != 'N' && c != 'n' && c != 'X' && c != 'x') {
            break;
        }
    }
    
    if (i == length) {
        return false;
    } return true;
}

/**
 * Perform a conversion of nucleotide chars into a bitfield representation of
 * each character. This will be decomposed into each component within the set
 * {0,1,2,3} later.
 */
uvec seq2numvec(string *seq) {
    assert(seq != NULL);

    uvec ret(seq->size());
    const char* inseq = seq->c_str();

    for (size_t i = 0; i < seq->size(); ++i) {
        ret[i] = getBitfieldFromFastaChar(inseq[i]);
    }

    // Additionally, reverse complements should have opposite facing.
    if (reverseComplement) {
        ret = flipud(ret);
    }
    return ret;
}

/**
 * Return the average probability of a given sequence
 */
double avgProb(u32 nuc, rowvec* given) {
    assert(given != NULL);
    
    double prob = 0.0;
    u32 denom = 0;
    if (nuc & 0x1) {
        prob += (*given)[0]; //prob *= log2(given[0]);
        ++denom;
    }
    if (nuc & 0x2) {
        prob += (*given)[1]; //prob *= log2(given[1]);
        ++denom;
    }
    if (nuc & 0x4) {
        prob += (*given)[2]; //prob *= log2(given[2]);
        ++denom;
    }
    if (nuc & 0x8) {
        prob += (*given)[3]; //prob *= log2(given[3]);
        ++denom;
    }

    if (denom > 0) {
        //prob = pow(2,prob);
        prob /= (double)denom;
    }
    return prob;
}

/**
 * Calculate the PWM score of query sequences in an index vector
 */
double matscore(uvec* query, mat* matrix) {
    assert(query != NULL);
    assert(matrix != NULL);
    assert(query->n_elem >= matrix->n_rows);
    
    double outscore = 1.0;
    for (u32 j = 0; j < matrix->n_rows; ++j) {
        rowvec row = matrix->row(j);
        double prob = avgProb((*query)[j], &row);
        outscore *= prob;
    }
    return outscore;
}

/**
 * Calculate the average probability of query sequences
 */
double matprob(uvec* query, mat* matrix) {
    assert(query != NULL);
    assert(matrix != NULL);
    assert(query->n_elem >= matrix->n_rows);
    
    double outprob = 0.0;
    u32 count = 0;
    for (u32 j = 0; j < matrix->n_rows; ++j) {
        rowvec row = matrix->row(j);
        outprob += avgProb((*query)[j], &row);
        //outprob *= log2(avgProb((*query)[j], &row));
        ++count;
    }
    if (count > 0) {
        //outprob = pow(2,outprob);
        outprob /= (double)count;
    }
    return outprob;
}

/**
 * Returns the index to address within m1 and m2
 */
uvec dinucindex(uvec* ind) {
    assert(ind != NULL);
    assert(ind->n_elem >= 2);

    u32 left = (((*ind)[0] & 0x8) >> 3) + (((*ind)[0] & 0x4) >> 2) +  (((*ind)[0] & 0x2) >> 1) + ((*ind)[0] & 0x1); 
    u32 right = (((*ind)[1] & 0x8) >> 3) + (((*ind)[1] & 0x4) >> 2) +  (((*ind)[1] & 0x2) >> 1) + ((*ind)[1] & 0x1); 
    uvec ret(left * right);

    u32 k = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if ((*ind)[0] & (1 << i) && (*ind)[1] & (1 << j)) {
                ret[k] = (i * 4) + j;
                ++k;
            }
        }
    }
    
    return ret;
}

/**
 * Calculates background score per position
 */
vec bgvec(string* seq, rowvec* m0, colvec* m1, mat* m2) {
    assert(seq != NULL);
    assert(m0 != NULL);
    assert(m1 != NULL);
    assert(m2 != NULL);
    assert(m0->n_elem == 4);
    assert(m1->n_elem == 16);
    assert(m2->n_rows == 16);
    assert(m2->n_cols == 4);

    uvec uSeq = seq2numvec(seq);
    u32 n = uSeq.n_elem;

    vec outvec(n);
    outvec[0] = avgProb(uSeq[0], m0);

    uvec dnt(2);
    dnt[0] = uSeq[0];
    dnt[1] = uSeq[1];

    uvec inds = dinucindex(&dnt);
    double prob = 0.0;
    u32 count = 0;
    for(u32 i = 0; i < inds.n_elem; ++i) {
        prob += (*m1)[inds[i]]; // prob *= log2((*m1)[inds[i]]);
        ++count;
    }
    if (count > 0) {
        //prob = pow(2,prob);
        prob /= (double)count;
    }
    outvec[1] = prob / avgProb(uSeq[1], m0);

    for(u32 i = 2; i < n; ++i) {
        dnt[0] = uSeq[i-2];
        dnt[1] = uSeq[i-1];

        inds = dinucindex(&dnt);
        prob = 0.0;
        count = 0;
        for(u32 j = 0; j < inds.n_elem; ++j) {
            if (uSeq[i] & 0x1) {
                prob += (*m2)(inds[j],0); //prob *= log2(m2(inds[j],0));
                ++count;
            }
            if (uSeq[i] & 0x2) {
                prob += (*m2)(inds[j],1); //prob *= log2((*m2)(inds[j],1));
                ++count;
            }
            if (uSeq[i] & 0x4) {
                prob += (*m2)(inds[j],2); //prob *= log2((*m2)(inds[j],2));
                ++count;
            }
            if (uSeq[i] & 0x8) {
                prob += (*m2)(inds[j],3); //prob *= log2((*m2)(inds[j],3));
                ++count;
            }
        }
        if (count > 0) {
            //prob = pow(2,prob);
            prob /= (double)count;
        }
        outvec[i] = prob;
    }

    return outvec;
}

/**
 * Main entry point
 */
int main(int argc, char** argv) {

    // Assert our numeric boundaries are valid
    assert(numeric_limits<double>::has_infinity);
    assert(numeric_limits<double>::has_quiet_NaN);
    assert(numeric_limits<double>::has_signaling_NaN);
    
    if (argc < 5) {
        cout << "Usage: " << argv[0]
             << " <sequence file> <background file> <threshold file> <matrix lib> [-quiet] [-filter] [-reverse]" << endl
             << "For example: " << argv[0]
             << " E2F_Q2_concensus.fa bModel.fa q999.recomputed.tsv matrices"
             << endl;
        exit(0xe6a5u);
    }

    if (argc > 5) {
        for (int i = 5; i < argc; ++i) {
            quiet |= strcmp(argv[i],"-quiet") == 0;
            filter |= strcmp(argv[i],"-filter") == 0;
            reverseComplement |= strcmp(argv[i],"-reverse") == 0;
        }
    }

    // Open the sequence file
    ifstream seqFile(argv[1]);
    if (seqFile == NULL) {
        cerr << "Sequence file failed to load. Aborting." << endl;
        exit(-1);
    }
    
    // Read in the sequence from stdin.
    string seq = fasta2Seq(&seqFile);
    size_t n = seq.size();
    if (n == 0) {
        cerr << "Input sequence was empty: " << argv[1] << ". Aborting." << endl;
        exit(-1);
    }

    // Open the background file
    ifstream backFile(argv[2]);
    if (backFile == NULL) {
        cerr << "Background file failed to load. Aborting." << endl;
        exit(-1);
    }

    // Open the threshold file
    ifstream threshFile(argv[3]);
    if (threshFile == NULL) {
        cerr << "Threshold file failed to load. Aborting." << endl;
        exit(-1);
    }

    // Read in the path to our positional weight matrices as a string
    string pwmPath(argv[4]);

    // Populate the positional weight matrices
    populatePwmData(&threshFile, pwmPath);

    // Retrieve data from our background file
    stringstream m0ss, m1ss, m2ss;
    char buf[__LINE_MAX];
    u32 index = 0;

    while (!backFile.eof() && !backFile.fail()) {
        
        buf[0] = 0x00;
        backFile.getline(buf, __LINE_MAX);
        if (buf[0] != 0x00 && buf[0] != '#' && buf[0] != ' ') {
            if (index == 0) {
                m0ss << buf << endl;
            } else if (index > 0 && index < 17) {
                m1ss << buf << endl;
            } else {
                m2ss << buf << endl;
            }
            ++index;
        }
    }

    // Close the input file handles, since we no longer need to access these files
    seqFile.close();
    backFile.close();
    threshFile.close();
    
    // Read m0. Input file has values only.
    rowvec m0; m0.load(m0ss);
    if (m0.n_elem != 4) {
        cerr << "Vector m0 was of incorrect length: " << m0.n_elem
             << " (should be 4). Aborting." << endl;
        exit(-1);
    }

    // Read m1. Input file has values only.
    colvec m1; m1.load(m1ss);
    if (m1.n_elem != 16) {
        cerr << "Vector m1 was of incorrect length: " << m1.n_elem
             << " (should be 16). Aborting." << endl;
        exit(-1);
    }

    // Read m2. Input file has values only.
    mat m2; m2.load(m2ss);
    if (m2.n_rows != 16) {
        cerr << "Matrix m2 was of incorrect row length: " << m2.n_rows
             << " (should be 16). Aborting." << endl;
        exit(-1);
    }
    if (m2.n_cols != 4) {
        cerr << "Matrix m2 was of incorrect column length: " << m2.n_cols
             << " (should be 4). Aborting." << endl;
        exit(-1);
    }

    if (!quiet) { cout << "Now processing background vector. Please wait..." << endl; }
    vec bgv = bgvec(&seq,&m0,&m1,&m2);
    
    if (!quiet) { cout << "Processing motifs..." << endl; }
    for (map<string, pair<mat, double> >::const_iterator it = pwmMap.begin(); it != pwmMap.end(); ++it) {
        string name = it->first;
        pair<mat, double> comp = it->second;
        
        mat pwm(comp.first);
        double threshold = comp.second;
        
        u32 m = pwm.n_rows;
        
        if (!quiet) { cout << "Processing: " << name << endl; }
        if (m > n) {
            cerr << "Candidate motifs were larger than the sequence! " <<
                "Please check both to verify if this is correct. Skipping."
                 << endl;
            continue;
        }
        if (m == 0) {
            cerr << "Input matrix was empty: " << name << ". Skipping." << endl;
            continue;
        }

        // Scan a sequence
        double max = 0.0; //(- HUGE_VAL);
        double min = 0.0; //HUGE_VAL // Assume a minimum of 0, allowing us to
                                     //  transitively prune the output.

        // Calculate the minimum and maximum scores for the entire sequence,
        //  such that we can compute the threshold value thereafter
        if (!quiet) { cout << "Now processing the local threshold values. Please wait..." << endl; }
        for (u32 i = 0; i <= n-m; ++i) {
            string subseq = seq.substr(i,m);                
            const char* cSubseq = subseq.c_str();

            if (filter) {
                if (!hasMeaningfulSignal(cSubseq,m,prevWasMeaningless)) {
                    prevWasMeaningless = true;
                    continue;
                } else {
                    prevWasMeaningless = false;
                }
            }
            
            uvec indvec = seq2numvec(&subseq);
            double prob = matprob(&indvec,&pwm);
            double numerator = matscore(&indvec,&pwm);
            double denominator = 1.0;
            for (u32 j = i; j < i+m; ++j ) {
                denominator *= bgv[j];
            }
                
            double pscore = numerator / denominator;                
            double logpscore = log(pscore);
            if (logpscore == numeric_limits<double>::infinity()
                || logpscore == -numeric_limits<double>::infinity()
                || logpscore == numeric_limits<double>::quiet_NaN()
                || logpscore == numeric_limits<double>::signaling_NaN()) { // Floor nonsensical values and div-by-zero

                logpscore = 0.0;
            }
            if (max < logpscore) {
                max = logpscore;
            }

            // Filter the output. Note that scores < 1.0 are squashed by this
            //  calculation, which we believe is reasonable as these offer
            //  extremely low confidence. The raw logarithmic scores are
            //  available when the filter switch is turned off.
            double norm = filter ? (logpscore - min) / (max - min) : logpscore;
            if (norm == numeric_limits<double>::infinity()
                || norm == -numeric_limits<double>::infinity()
                || norm == numeric_limits<double>::quiet_NaN()
                || norm == numeric_limits<double>::signaling_NaN()) { // Floor nonsensical values and div-by-zero

                norm = 0.0;
            }
            if (!filter || norm >= threshold) {
                stringstream retSubseq;
                char lStrand = strand;
                
                // If the reverse complement is desired
                if (reverseComplement) {    
                    // Reverse the strand orientation 
                    lStrand = (strand == '+') ? '-' : '+';
                
                    // Invert the character set
                    size_t j = m - 1;
                    do {
                        retSubseq << getFastaCharFromBitfield(
                            getBitfieldFromFastaChar(cSubseq[j]));
                    } while (j-- > 0); // Defined in this way to prevent
                                       //  reliance on an underflow of
                                       //  potentially unsigned size_t. 
                    
                } else {
                    retSubseq << subseq;
                }
                cout << argv[1] << "\t"
                     << i << "\t"
                     << (i + m - 1) << "\t"
                     << lStrand << "\t"
                     << prob << "\t"
                     << norm << "\t"
                     << pscore << "\t"
                     << name << "\t"
                     << retSubseq.str() << endl;
            }
        }
    }
}
